using System.Diagnostics;
using System.Net;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace DarkUniverse;

/// <summary>An immutable run directory and explicit terminal outcome, with success-only latest publication.</summary>
public sealed class RunSession : IDisposable
{
    static readonly JsonSerializerOptions JsonOptions = new()
    {
        WriteIndented = true,
        NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals
    };
    readonly SafetyPaths paths;
    readonly Stopwatch timer = Stopwatch.StartNew();
    readonly DateTimeOffset started = DateTimeOffset.UtcNow;
    bool terminal;

    public string RunId { get; }
    public string DirectoryPath { get; }
    public string ReportPath => Path.Combine(DirectoryPath, "run_report.json");
    public string Stage { get; private set; } = "initializing";
    public string Command { get; }
    public Dictionary<string, object?> Metadata { get; } = [];

    RunSession(SafetyPaths paths, string command)
    {
        this.paths = paths;
        Command = command;
        RunId = started.ToString("yyyyMMdd'T'HHmmssfffffff'Z'") + "-" + Guid.NewGuid().ToString("N");
        DirectoryPath = Path.Combine(paths.OutputRoot, "runs", RunId);
        var analysisAssembly = typeof(RunSession).Assembly;
        var entryAssembly = Assembly.GetEntryAssembly();
        Metadata["assembly_sha256"] = AssemblyHash(analysisAssembly);
        Metadata["assembly_version"] = analysisAssembly.GetName().Version?.ToString();
        Metadata["entry_assembly_sha256"] = AssemblyHash(entryAssembly);
        Metadata["runtime"] = Environment.Version.ToString();
        Metadata["framework"] = RuntimeInformation.FrameworkDescription;
        Metadata["os_platform"] = RuntimeInformation.OSDescription;
        Metadata["data_directory"] = paths.DataRoot;
        Metadata["requested_output_directory"] = paths.OutputRoot;
    }

    /// <summary>Call only after path validation. Manifest validation belongs inside the caller's try/catch so its failure is recorded.</summary>
    public static RunSession Start(SafetyPaths paths, string command)
    {
        // Revalidate rather than trusting a publicly constructible paths record.
        SafetyPaths checkedPaths = ExecutionSafety.ValidateRoots(paths.DataRoot, paths.OutputRoot);
        var session = new RunSession(checkedPaths, command);
        ExecutionSafety.AssertSafeWritePath(checkedPaths.OutputRoot, session.ReportPath);
        if (Directory.Exists(session.DirectoryPath) || File.Exists(session.DirectoryPath))
            throw new IOException("A run identifier already exists; no output will be reused.");
        Directory.CreateDirectory(session.DirectoryPath);
        try
        {
            session.WriteReport(session.BaseReport("running", false));
            session.WriteAttempt("running", false);
            return session;
        }
        catch (Exception error)
        {
            // Start has not returned to the caller yet, so it owns the failure
            // transition if publishing the initial attempt record fails.
            session.Fail(error);
            throw;
        }
    }

    public void SetStage(string name)
    {
        if (terminal) throw new InvalidOperationException("A terminal run cannot enter another stage.");
        if (string.IsNullOrWhiteSpace(name)) throw new ArgumentException("A stage name is required.");
        Stage = name;
        WriteReport(BaseReport("running", false));
        WriteAttempt("running", false);
    }

    /// <summary>Seal a successful report, then atomically replace latest pointers; previous runs and galleries are never overwritten.</summary>
    public void Complete(IDictionary<string, object?> report)
    {
        if (terminal) throw new InvalidOperationException("This run already has a terminal outcome.");
        foreach (var (key, value) in report) Metadata[key] = value;
        Stage = "publish_latest";
        ExecutionSafety.AssertNoReparsePoints(DirectoryPath);
        string gallery = Path.Combine(DirectoryPath, "index.html");
        if (!File.Exists(gallery))
        {
            string label = WebUtility.HtmlEncode(Command);
            AtomicWrite(gallery, "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Completed run</title>" +
                $"<h1>Completed {label} run</h1><p>Run {RunId}</p><p><a href=\"run_report.json\">Run report and results</a></p></html>");
        }
        Stage = "completed";
        var completed = BaseReport("completed", true);
        AddEnd(completed);
        WriteReport(completed);

        string relativeGallery = "runs/" + RunId + "/index.html";
        string html = "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Latest successful analysis</title>" +
            $"<meta http-equiv=\"refresh\" content=\"0;url={relativeGallery}\"><h1>Latest successful run</h1>" +
            $"<p><a href=\"{relativeGallery}\">{WebUtility.HtmlEncode(Command)} — {RunId}</a></p>" +
            "<p><a href=\"latest-attempt.json\">Most recent attempt, including failures</a></p></html>";
        var pointer = Pointer("completed", true);
        pointer["completed_at_utc"] = completed["ended_at_utc"];
        // Replacement is atomic per file. If any promotion fails, restore the
        // old pointers before propagating the error to the caller's Fail path.
        var replacements = new (string Path, string Text)[]
        {
            (Path.Combine(paths.OutputRoot, "index.html"), html),
            (Path.Combine(paths.OutputRoot, "latest.json"), Serialize(pointer)),
            (Path.Combine(paths.OutputRoot, "latest-attempt.json"), Serialize(pointer))
        };
        foreach (var item in replacements) ExecutionSafety.AssertSafeWritePath(paths.OutputRoot, item.Path);
        var previous = replacements.Select(item => (item.Path, Bytes: File.Exists(item.Path) ? File.ReadAllBytes(item.Path) : null)).ToArray();
        int promoted = 0;
        try
        {
            foreach (var item in replacements) { AtomicWrite(item.Path, item.Text); promoted++; }
        }
        catch
        {
            for (int i = promoted - 1; i >= 0; i--)
            {
                try
                {
                    if (previous[i].Bytes is null)
                    {
                        ExecutionSafety.AssertSafeWritePath(paths.OutputRoot, previous[i].Path);
                        File.Delete(previous[i].Path);
                    }
                    else AtomicWriteBytes(previous[i].Path, previous[i].Bytes!);
                }
                catch (Exception rollbackError) { Console.Error.WriteLine("Could not restore a latest-run pointer: " + rollbackError.Message); }
            }
            Stage = "publish_latest";
            throw;
        }
        terminal = true;
        timer.Stop();
    }

    /// <summary>Record a failure inside its own run; never modify the latest successful report or gallery.</summary>
    public void Fail(Exception exception, IDictionary<string, object?>? report = null)
    {
        if (terminal) return;
        if (report is not null) foreach (var (key, value) in report) Metadata[key] = value;
        var failure = BaseReport("failed", false);
        AddEnd(failure);
        failure["failed_stage"] = Stage;
        failure["error"] = new { type = exception.GetType().FullName, message = exception.Message, details = exception.ToString() };
        terminal = true;
        timer.Stop();
        try { WriteReport(failure); }
        catch (Exception writeError) { Console.Error.WriteLine("Could not persist the failed-run report: " + writeError.Message); }
        try { WriteAttempt("failed", false, failure["ended_at_utc"]); }
        catch (Exception pointerError) { Console.Error.WriteLine("Could not update the latest-attempt pointer: " + pointerError.Message); }
    }

    Dictionary<string, object?> BaseReport(string status, bool passed)
    {
        var value = new Dictionary<string, object?>(Metadata)
        {
            ["run_id"] = RunId,
            ["command"] = Command,
            ["status"] = status,
            ["passed"] = passed,
            ["stage"] = Stage,
            ["started_at_utc"] = started.ToString("O"),
            ["output_directory"] = DirectoryPath
        };
        return value;
    }

    void AddEnd(Dictionary<string, object?> report)
    {
        report["ended_at_utc"] = DateTimeOffset.UtcNow.ToString("O");
        report["duration_seconds"] = timer.Elapsed.TotalSeconds;
    }

    Dictionary<string, object?> Pointer(string status, bool passed) => new()
    {
        ["run_id"] = RunId, ["command"] = Command, ["status"] = status, ["passed"] = passed,
        ["stage"] = Stage, ["started_at_utc"] = started.ToString("O"),
        ["run_report"] = "runs/" + RunId + "/run_report.json",
        ["gallery"] = "runs/" + RunId + "/index.html"
    };

    void WriteAttempt(string status, bool passed, object? end = null)
    {
        var pointer = Pointer(status, passed);
        if (end is not null) pointer["ended_at_utc"] = end;
        AtomicWrite(Path.Combine(paths.OutputRoot, "latest-attempt.json"), Serialize(pointer));
    }

    void WriteReport(object value) => AtomicWrite(ReportPath, Serialize(value));
    static string Serialize(object value) => JsonSerializer.Serialize(value, JsonOptions) + "\n";
    void AtomicWrite(string path, string text) => AtomicWriteBytes(path, new UTF8Encoding(false).GetBytes(text));

    void AtomicWriteBytes(string path, byte[] bytes)
    {
        ExecutionSafety.AssertSafeWritePath(paths.OutputRoot, path);
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);
        string temporary = path + "." + Guid.NewGuid().ToString("N") + ".tmp";
        ExecutionSafety.AssertSafeWritePath(paths.OutputRoot, temporary);
        try
        {
            using (var stream = new FileStream(temporary, FileMode.CreateNew, FileAccess.Write, FileShare.None))
            {
                stream.Write(bytes);
                stream.Flush(true);
            }
            ExecutionSafety.AssertSafeWritePath(paths.OutputRoot, path);
            File.Move(temporary, path, true);
        }
        finally
        {
            if (File.Exists(temporary))
            {
                ExecutionSafety.AssertSafeWritePath(paths.OutputRoot, temporary);
                File.Delete(temporary);
            }
        }
    }

    static string? AssemblyHash(Assembly? assembly)
    {
        if (assembly is null || string.IsNullOrEmpty(assembly.Location)) return null;
        using var stream = File.OpenRead(assembly.Location);
        return Convert.ToHexStringLower(SHA256.HashData(stream));
    }

    public void Dispose()
    {
        if (!terminal) Fail(new InvalidOperationException("The run ended without an explicit success or failure outcome."));
    }
}
