using System.Globalization;
using System.Security;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace DarkUniverse;

/// <summary>
/// Descriptive replay of explicitly prescribed formation experiments. This does
/// not fit observations or assign a probability measure to real galaxies.
/// </summary>
public static class FormationStudy
{
    public const string Scope = "conditional-design";
    public const string ScopeLabel = "Conditional design experiments; no observed galaxy population inference";
    static readonly JsonSerializerOptions Options = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
        UnmappedMemberHandling = JsonUnmappedMemberHandling.Disallow
    };

    public static FormationSummary Run(string experimentPath, string outputRoot)
    {
        string input = ExecutionSafety.CanonicalizePath(experimentPath);
        ExecutionSafety.AssertUnaliasedInputFile(input);
        if (ExecutionSafety.ContainsPath(outputRoot, input))
            throw new InvalidDataException("The formation bundle must be outside the output run.");
        string inputHash = Data.FileSha256(input);
        byte[] bytes = File.ReadAllBytes(input);
        using (var raw = JsonDocument.Parse(bytes)) RejectDuplicateProperties(raw.RootElement);
        FormationBundle bundle = JsonSerializer.Deserialize<FormationBundle>(bytes, Options)
            ?? throw new InvalidDataException("Empty formation bundle.");
        Validate(bundle);
        var sources = LoadSources(bundle, Path.GetDirectoryName(input)!, outputRoot);
        foreach (var study in bundle.Studies)
            foreach (var item in study.Cases)
                VerifyValue(item.OutcomeSource, item.Outcome, sources, study.Id + "/" + item.Id);

        var summaries = bundle.Studies.SelectMany(study => study.Cases.GroupBy(item => item.Group)
            .Select(group => Summarize(study.Id, group.Key, group.ToArray()))).ToArray();
        var report = new FormationSummary(1, Scope, false, bundle.Title, inputHash,
            sources.ToDictionary(pair => pair.Key, pair => pair.Value.Hash), summaries);
        string directory = Path.Combine(outputRoot, "formation");
        ExecutionSafety.AssertSafeWritePath(outputRoot, Path.Combine(directory, "summary.json"));
        Directory.CreateDirectory(directory);
        Data.SaveJson(Path.Combine(directory, "summary.json"), report);
        Csv.Write(Path.Combine(directory, "case_statistics.csv"), bundle.Studies.SelectMany(study => study.Cases.Select(item =>
            new Dictionary<string, object?>
            {
                ["study"] = study.Id, ["case"] = item.Id, ["group"] = item.Group,
                ["control"] = item.Control, ["outcome"] = item.Outcome, ["weight"] = item.Weight,
                ["outcome_units"] = study.OutcomeUnits, ["scope"] = Scope
            })));
        foreach (var study in bundle.Studies) WritePlot(study, directory);
        WriteIndex(bundle, summaries, directory);
        File.WriteAllText(Path.Combine(outputRoot, "index.html"),
            "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta http-equiv=\"refresh\" content=\"0;url=formation/index.html\"><title>Conditional formation experiments</title><a href=\"formation/index.html\">Open the conditional formation experiment report</a></html>");
        if (Data.FileSha256(input) != inputHash) throw new InvalidDataException("Formation bundle changed during execution.");
        foreach (var source in sources.Values)
            if (Data.FileSha256(source.Path) != source.Hash)
                throw new InvalidDataException("Formation source changed during execution: " + source.Path);
        return report;
    }

    public static void Validate(FormationBundle bundle)
    {
        if (bundle.SchemaVersion != 1 || bundle.Scope != Scope || bundle.ObservationalValidation)
            throw new InvalidDataException("Formation schema 1 requires scope=conditional-design and observational_validation=false.");
        RequireText(bundle.Title, "title");
        if (bundle.Provenance is null || bundle.Provenance.Length == 0 || bundle.Studies is null || bundle.Studies.Length == 0)
            throw new InvalidDataException("At least one provenance source and study are required.");
        Unique(bundle.Provenance.Select(item => item.Id), "source");
        foreach (var source in bundle.Provenance)
        {
            RequireText(source.Path, "source path");
            if (Path.IsPathRooted(source.Path) || source.Path.Contains('\\') || source.Path.Contains(':'))
                throw new InvalidDataException("Formation source paths must be relative and use forward slashes.");
            if (source.Sha256 is null || source.Sha256.Length != 64 || source.Sha256.Any(c => !Uri.IsHexDigit(c)))
                throw new InvalidDataException("Formation sources require SHA-256 digests.");
        }
        Unique(bundle.Studies.Select(item => item.Id), "study");
        foreach (var study in bundle.Studies)
        {
            foreach (var value in new[] { study.Title, study.Method, study.ControlLabel, study.OutcomeLabel, study.OutcomeUnits }) RequireText(value, "study label");
            if (study.Cases is null || study.Cases.Length == 0) throw new InvalidDataException("A formation study must contain cases.");
            if (study.OutcomeMinimum is { } minimum) Finite(minimum, "outcome minimum");
            if (study.OutcomeMaximum is { } maximum) Finite(maximum, "outcome maximum");
            if (study.OutcomeMinimum > study.OutcomeMaximum) throw new InvalidDataException("Outcome bounds are reversed.");
            Unique(study.Cases.Select(item => item.Id), "case");
            foreach (var item in study.Cases)
            {
                RequireText(item.Group, "case group");
                Finite(item.Control, "control"); Finite(item.Outcome, "outcome"); Finite(item.Weight, "weight");
                if (item.Weight <= 0) throw new InvalidDataException("Design weights must be positive.");
                if (item.Outcome < study.OutcomeMinimum || item.Outcome > study.OutcomeMaximum)
                    throw new InvalidDataException("Outcome is outside its declared range: " + item.Id);
                if (item.OutcomeSource is null || !bundle.Provenance.Any(source => source.Id == item.OutcomeSource.SourceId))
                    throw new InvalidDataException("Every outcome must identify a declared provenance source.");
                if (item.OutcomeSource.JsonPointer is null || !item.OutcomeSource.JsonPointer.StartsWith('/'))
                    throw new InvalidDataException("Outcome source must contain an absolute JSON pointer.");
                if (item.Diagnostics is null) throw new InvalidDataException("Diagnostics must be an object, possibly empty.");
                foreach (var (name, value) in item.Diagnostics) { RequireText(name, "diagnostic key"); Finite(value, "diagnostic"); }
            }
        }
    }

    public static FormationGroupSummary Summarize(string study, string group, FormationCase[] cases)
    {
        if (cases.Length == 0) throw new InvalidDataException("Cannot summarize an empty design.");
        foreach (var item in cases)
        {
            Finite(item.Outcome, "outcome"); Finite(item.Control, "control"); Finite(item.Weight, "weight");
            if (item.Weight <= 0) throw new InvalidDataException("Design weights must be positive.");
        }
        // Normalize by the largest weight to avoid overflow when the measure is
        // rescaled. Weighted Welford is stable for closely spaced outcomes.
        double maxWeight = cases.Max(item => item.Weight), total = 0, mean = 0, sumSquares = 0;
        foreach (var item in cases)
        {
            double weight = item.Weight / maxWeight;
            if (weight == 0) continue; // Negligible underflowed relative weight.
            double next = total + weight;
            double delta = item.Outcome - mean;
            double nextMean = mean + delta * (weight / next);
            sumSquares += weight * delta * (item.Outcome - nextMean);
            total = next; mean = nextMean;
        }
        double variance = Math.Max(0, sumSquares / total);
        Finite(mean, "design mean"); Finite(variance, "design variance");
        var diagnostics = cases.SelectMany(item => item.Diagnostics).GroupBy(pair => pair.Key)
            .ToDictionary(values => values.Key, values => values.Max(pair => Math.Abs(pair.Value)));
        return new FormationGroupSummary(study, group, cases.Length, mean, Math.Sqrt(variance),
            cases.Min(item => item.Outcome), cases.Max(item => item.Outcome),
            cases.Min(item => item.Control), cases.Max(item => item.Control), diagnostics);
    }

    sealed record SourceData(string Path, string Hash, JsonElement Root);

    static Dictionary<string, SourceData> LoadSources(FormationBundle bundle, string baseDirectory, string outputRoot)
    {
        var sources = new Dictionary<string, SourceData>(StringComparer.Ordinal);
        foreach (var source in bundle.Provenance)
        {
            string path = ExecutionSafety.CanonicalizePath(Path.Combine(baseDirectory, source.Path.Replace('/', Path.DirectorySeparatorChar)));
            if (ExecutionSafety.ContainsPath(outputRoot, path)) throw new InvalidDataException("Formation sources must be outside the output run.");
            ExecutionSafety.AssertUnaliasedInputFile(path);
            string hash = Data.FileSha256(path);
            if (!hash.Equals(source.Sha256, StringComparison.OrdinalIgnoreCase)) throw new InvalidDataException("Formation provenance mismatch: " + source.Id);
            using var raw = JsonDocument.Parse(File.ReadAllBytes(path));
            RejectDuplicateProperties(raw.RootElement);
            sources.Add(source.Id, new SourceData(path, hash, raw.RootElement.Clone()));
        }
        return sources;
    }

    static void VerifyValue(FormationValueSource reference, double expected, Dictionary<string, SourceData> sources, string context)
    {
        JsonElement value = sources[reference.SourceId].Root;
        foreach (string encoded in reference.JsonPointer[1..].Split('/'))
        {
            if (encoded.Contains('~') && encoded.Replace("~0", "").Replace("~1", "").Contains('~'))
                throw new InvalidDataException("Invalid JSON pointer escape: " + context);
            string token = encoded.Replace("~1", "/").Replace("~0", "~");
            if (value.ValueKind == JsonValueKind.Object && value.TryGetProperty(token, out var property)) value = property;
            else if (value.ValueKind == JsonValueKind.Array && int.TryParse(token, NumberStyles.None, CultureInfo.InvariantCulture, out int index) && index >= 0 && index < value.GetArrayLength()) value = value[index];
            else throw new InvalidDataException("Unresolved outcome JSON pointer: " + context);
        }
        if (value.ValueKind != JsonValueKind.Number || !value.TryGetDouble(out double actual) || !double.IsFinite(actual) || actual != expected)
            throw new InvalidDataException("Outcome differs from its provenance source: " + context);
    }

    static void RejectDuplicateProperties(JsonElement element)
    {
        if (element.ValueKind == JsonValueKind.Object)
        {
            var names = new HashSet<string>(StringComparer.Ordinal);
            foreach (var property in element.EnumerateObject())
            {
                if (!names.Add(property.Name)) throw new InvalidDataException("Duplicate JSON property: " + property.Name);
                RejectDuplicateProperties(property.Value);
            }
        }
        else if (element.ValueKind == JsonValueKind.Array)
            foreach (var child in element.EnumerateArray()) RejectDuplicateProperties(child);
    }

    static void WritePlot(FormationExperiment study, string directory)
    {
        var document = new PlotDocument(study.Title, 2, 1) { Scope = ScopeLabel };
        var response = new PlotPanel("Prescribed formation response", study.ControlLabel, study.OutcomeLabel) { XMin = study.Cases.Min(item => item.Control) == 0 ? 0 : null, YMin = study.OutcomeMinimum, YMax = study.OutcomeMaximum };
        var distribution = new PlotPanel("Distribution over the declared design", study.OutcomeLabel, "Cumulative design weight") { XMin = study.OutcomeMinimum, XMax = study.OutcomeMaximum, YMin = 0, YMax = 1 };
        string[] colors = ["#235a81", "#ad5822", "#41845f", "#8064a2", "#a84762", "#575757"];
        int index = 0;
        foreach (var group in study.Cases.GroupBy(item => item.Group))
        {
            var ordered = group.OrderBy(item => item.Control).ThenBy(item => item.Id).ToArray();
            string color = colors[index++ % colors.Length];
            response.Add(group.Key, ordered.Select(item => item.Control).ToArray(), ordered.Select(item => item.Outcome).ToArray(), color, true);
            var outcomes = group.OrderBy(item => item.Outcome).ToArray();
            double max = outcomes.Max(item => item.Weight), total = outcomes.Sum(item => item.Weight / max), running = 0;
            var x = new List<double> { study.OutcomeMinimum ?? outcomes[0].Outcome }; var y = new List<double> { 0 };
            foreach (var item in outcomes) { running += item.Weight / max; x.Add(item.Outcome); y.Add(running / total); }
            distribution.Add(group.Key, x.ToArray(), y.ToArray(), color).Step = true;
        }
        document.Panels.Add(response); document.Panels.Add(distribution);
        document.Save(Path.Combine(directory, study.Id + ".svg"));
    }

    static void WriteIndex(FormationBundle bundle, FormationGroupSummary[] summaries, string directory)
    {
        static string E(string text) => SecurityElement.Escape(text) ?? "";
        static string N(double value) => value.ToString("G8", CultureInfo.InvariantCulture);
        var html = new StringBuilder("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Conditional formation experiments</title><style>body{font:16px system-ui;max-width:1250px;margin:2rem auto;padding:0 1rem;line-height:1.5}table{border-collapse:collapse}td,th{padding:.4rem .7rem;border:1px solid #ccc;text-align:right}th:first-child,td:first-child{text-align:left}img{width:100%}.scope{padding:1rem;background:#edf3f6}</style>");
        html.Append("<h1>").Append(E(bundle.Title)).Append("</h1><p class=\"scope\">").Append(ScopeLabel)
            .Append(". Every control, initial state and design weight is prescribed. Means and population standard deviations summarize only those design weights; they are not uncertainty intervals, p-values or inferred galaxy frequencies. Existing publication fits and statistical tests are unchanged.</p><p><a href=\"summary.json\">Machine-readable statistics</a> | <a href=\"case_statistics.csv\">Cases and design weights</a></p>");
        foreach (var study in bundle.Studies)
        {
            html.Append("<h2>").Append(E(study.Title)).Append("</h2><p>").Append(E(study.Method)).Append("</p><p>Outcome: ")
                .Append(E(study.OutcomeLabel)).Append("; units: ").Append(E(study.OutcomeUnits)).Append(".</p><img alt=\"")
                .Append(E(study.Title + ": conditional design response and cumulative design weight")).Append("\" src=\"").Append(study.Id).Append(".svg\"><table><tr><th>Group</th><th>Cases</th><th>Design mean</th><th>Design SD</th><th>Minimum</th><th>Maximum</th></tr>");
            foreach (var row in summaries.Where(row => row.Study == study.Id))
                html.Append("<tr><td>").Append(E(row.Group)).Append("</td><td>").Append(row.Count).Append("</td><td>")
                    .Append(N(row.WeightedMean)).Append("</td><td>").Append(N(row.WeightedPopulationStandardDeviation)).Append("</td><td>")
                    .Append(N(row.Minimum)).Append("</td><td>").Append(N(row.Maximum)).Append("</td></tr>");
            html.Append("</table><p><a href=\"").Append(study.Id).Append(".series.csv\">Plotted coordinates</a></p>");
        }
        html.Append("</html>"); File.WriteAllText(Path.Combine(directory, "index.html"), html.ToString());
    }

    static void Finite(double value, string name) { if (!double.IsFinite(value)) throw new InvalidDataException("Nonfinite formation " + name + "."); }
    static void RequireText(string? value, string name) { if (string.IsNullOrWhiteSpace(value)) throw new InvalidDataException("Missing formation " + name + "."); }
    static void Unique(IEnumerable<string> values, string name)
    {
        var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (string value in values)
        {
            if (string.IsNullOrWhiteSpace(value) || value.Any(c => !(char.IsAsciiLetterOrDigit(c) || c is '_' or '-' or '.')) || value is "." or ".." || value.EndsWith('.') || !seen.Add(value))
                throw new InvalidDataException("Formation " + name + " IDs must be unique portable identifiers.");
        }
    }
}

public sealed record FormationBundle
{
    public required int SchemaVersion { get; init; }
    public required string Scope { get; init; }
    public required bool ObservationalValidation { get; init; }
    public required string Title { get; init; }
    public required FormationProvenance[] Provenance { get; init; }
    public required FormationExperiment[] Studies { get; init; }
}
public sealed record FormationProvenance
{
    public required string Id { get; init; }
    public required string Path { get; init; }
    public required string Sha256 { get; init; }
}
public sealed record FormationExperiment
{
    public required string Id { get; init; }
    public required string Title { get; init; }
    public required string Method { get; init; }
    public required string ControlLabel { get; init; }
    public required string OutcomeLabel { get; init; }
    public required string OutcomeUnits { get; init; }
    public double? OutcomeMinimum { get; init; }
    public double? OutcomeMaximum { get; init; }
    public required FormationCase[] Cases { get; init; }
}
public sealed record FormationCase
{
    public required string Id { get; init; }
    public required string Group { get; init; }
    public required double Control { get; init; }
    public required double Outcome { get; init; }
    public required double Weight { get; init; }
    public required FormationValueSource OutcomeSource { get; init; }
    public required Dictionary<string, double> Diagnostics { get; init; }
}
public sealed record FormationValueSource
{
    public required string SourceId { get; init; }
    public required string JsonPointer { get; init; }
}
public sealed record FormationGroupSummary(string Study, string Group, int Count, double WeightedMean,
    double WeightedPopulationStandardDeviation, double Minimum, double Maximum, double ControlMinimum,
    double ControlMaximum, Dictionary<string, double> MaximumAbsoluteDiagnostics);
public sealed record FormationSummary(int SchemaVersion, string Scope, bool ObservationalValidation, string Title,
    string ExperimentSha256, Dictionary<string, string> SourceSha256, FormationGroupSummary[] Groups);
