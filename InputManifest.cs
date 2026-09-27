using System.Collections.ObjectModel;
using System.Security.Cryptography;
using System.Text.Json;

namespace DarkUniverse;

/// <summary>A complete, verified input snapshot. Files excludes the manifest itself.</summary>
public sealed class InputSnapshot
{
    public string RootPath { get; }
    public string ManifestSha256 { get; }
    public IReadOnlyDictionary<string, string> Files { get; }

    internal InputSnapshot(string root, string manifestSha256, Dictionary<string, string> files)
    {
        RootPath = root;
        ManifestSha256 = manifestSha256;
        Files = new ReadOnlyDictionary<string, string>(files);
    }

    public void VerifyUnchanged()
    {
        InputSnapshot current = InputManifest.Verify(RootPath);
        if (current.ManifestSha256 != ManifestSha256 || current.Files.Count != Files.Count ||
            Files.Any(pair => !current.Files.TryGetValue(pair.Key, out string? value) || value != pair.Value))
            throw new InvalidDataException("Inputs or their manifest changed during the run.");
    }
}

public static class InputManifest
{
    public const string FileName = "input_manifest.json";

    /// <summary>Require a nonempty complete manifest, ordinary files, portable unique paths and exact SHA-256 matches.</summary>
    public static InputSnapshot Verify(string dataRoot)
    {
        string root = Path.GetFullPath(dataRoot);
        var ordinaryFiles = ExecutionSafety.EnumerateRegularFiles(root);
        var actual = new Dictionary<string, string>(StringComparer.Ordinal);
        var insensitiveNames = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (string path in ordinaryFiles)
        {
            string name = Path.GetRelativePath(root, path).Replace(Path.DirectorySeparatorChar, '/');
            ValidateName(name);
            if (!insensitiveNames.Add(name)) throw new InvalidDataException("Case-aliased input paths: " + name);
            ExecutionSafety.AssertUnaliasedInputFile(path);
            actual.Add(name, path);
        }
        if (!actual.Remove(FileName, out string? manifestPath))
            throw new InvalidDataException("A complete input_manifest.json is required for every analysis.");
        byte[] bytes = File.ReadAllBytes(manifestPath);
        string manifestHash = Convert.ToHexStringLower(SHA256.HashData(bytes));
        using var document = JsonDocument.Parse(bytes);
        if (document.RootElement.ValueKind != JsonValueKind.Object)
            throw new InvalidDataException("The input manifest must be a path-to-SHA-256 object.");
        var expected = new Dictionary<string, string>(StringComparer.Ordinal);
        var manifestNames = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (var property in document.RootElement.EnumerateObject())
        {
            ValidateName(property.Name);
            if (property.Name.Equals(FileName, StringComparison.OrdinalIgnoreCase))
                throw new InvalidDataException("The input manifest must not list itself.");
            if (!manifestNames.Add(property.Name))
                throw new InvalidDataException("Duplicate or case-aliased manifest path: " + property.Name);
            if (property.Value.ValueKind != JsonValueKind.String)
                throw new InvalidDataException("A manifest digest must be a hexadecimal SHA-256 string: " + property.Name);
            string hash = property.Value.GetString()!;
            if (hash.Length != 64 || hash.Any(c => !Uri.IsHexDigit(c)))
                throw new InvalidDataException("Invalid SHA-256 digest: " + property.Name);
            expected.Add(property.Name, hash.ToLowerInvariant());
        }
        if (expected.Count == 0) throw new InvalidDataException("The input manifest must not be empty.");
        var missing = expected.Keys.Except(actual.Keys, StringComparer.Ordinal).Order().ToArray();
        var extra = actual.Keys.Except(expected.Keys, StringComparer.Ordinal).Order().ToArray();
        if (missing.Length > 0 || extra.Length > 0)
            throw new InvalidDataException("Input manifest coverage differs. Missing: " + string.Join(", ", missing) + "; unlisted: " + string.Join(", ", extra));
        foreach (var (name, hash) in expected)
        {
            using var stream = File.OpenRead(actual[name]);
            string observed = Convert.ToHexStringLower(SHA256.HashData(stream));
            if (observed != hash) throw new InvalidDataException("Input manifest mismatch: " + name);
        }
        return new InputSnapshot(root, manifestHash, expected);
    }

    static void ValidateName(string name)
    {
        if (string.IsNullOrEmpty(name) || Path.IsPathRooted(name) || name.Contains('\\') || name.Contains(':'))
            throw new InvalidDataException("Manifest paths must be portable relative paths with forward slashes: " + name);
        foreach (string component in name.Split('/'))
        {
            if (component.Length == 0 || component is "." or ".." || component.EndsWith('.') || component.EndsWith(' ') ||
                component.Any(c => c < 32 || "<>:\"|?*".Contains(c)))
                throw new InvalidDataException("Unsafe or ambiguous input path: " + name);
            string device = component.Split('.')[0].ToUpperInvariant();
            if (device is "CON" or "PRN" or "AUX" or "NUL" ||
                device.Length == 4 && (device.StartsWith("COM", StringComparison.Ordinal) || device.StartsWith("LPT", StringComparison.Ordinal)) && device[3] is >= '1' and <= '9')
                throw new InvalidDataException("Reserved device name in input path: " + name);
        }
    }
}
