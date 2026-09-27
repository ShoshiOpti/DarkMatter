using System.IO.Compression;
using System.Text;
using System.Text.Json;
using System.Xml;
using System.Xml.Linq;

namespace DarkUniverse;

/// <summary>
/// Replays all current manuscript figures from audited frozen vector scenes.
/// The scenes retain exact glyphs, paths, axes, colors and artist coordinates.
/// This is deliberately separate from plots made from newly computed statistics:
/// changed scientific inputs require scene regeneration with the retained source
/// renderers, not silently displaying an old figure as a fresh calculation.
/// </summary>
public static class PublicationPlots
{
    public const string Scope = "Frozen publication figure replay from retained vector primitives and scientific artist coordinates; no new fitting, inference or field evolution.";
    static readonly XNamespace Svg = "http://www.w3.org/2000/svg";
    static readonly JsonSerializerOptions JsonOptions = new() { WriteIndented = true };
    static readonly HashSet<string> ForbiddenElements = ["script", "foreignObject", "iframe", "object", "embed"];

    public sealed record Report(int FigureCount, int AxisCount, int SeriesCount, long CoordinateRows,
        string ReportPath, string[] Files);

    /// <summary>
    /// Read dataRoot/publication/figures and write outputRoot/publication/figures.
    /// SVGs preserve the original renderer's vector scene rather than approximating
    /// logarithmic transforms, error bars, colorbars or mathematical typography.
    /// Numerical sidecars expose scientific coordinates independently of display paths.
    /// </summary>
    public static Report Run(string dataRoot, string outputRoot)
    {
        PublicationCompatibility.Verify(dataRoot, outputRoot);
        string input = Path.Combine(dataRoot, "publication", "figures");
        string output = Path.Combine(outputRoot, "publication", "figures");
        if (Path.GetFullPath(input).Equals(Path.GetFullPath(output), StringComparison.OrdinalIgnoreCase))
            throw new InvalidOperationException("Publication figure output must differ from immutable input data.");
        using var manifestDocument = JsonDocument.Parse(File.ReadAllText(Path.Combine(input, "manifest.json")));
        var manifest = manifestDocument.RootElement;
        Require(manifest.GetProperty("schema_version").GetInt32() == 1, "Unsupported publication scene schema");
        var figures = manifest.GetProperty("figures").EnumerateArray().ToArray();
        Require(figures.Length == 25 && manifest.GetProperty("figure_count").GetInt32() == 25,
            "The consolidated main paper and supplement require exactly 25 distinct figures");
        Require(figures.Select(f => f.GetProperty("name").GetString()).Distinct().Count() == 25,
            "Publication figure names must be unique");
        // Validate every frozen file before creating any output.
        foreach (var entry in figures)
        {
            string file = entry.GetProperty("file").GetString()!;
            Require(Path.GetFileName(file) == file, "A scene filename must not leave the input directory");
            string path = Path.Combine(input, file);
            Require(Data.FileSha256(path).Equals(entry.GetProperty("sha256").GetString(), StringComparison.OrdinalIgnoreCase),
                $"Publication scene checksum mismatch: {file}");
        }
        Directory.CreateDirectory(output);
        var files = new List<string>();
        var records = new List<object>();
        int axes = 0, series = 0;
        long coordinates = 0;
        foreach (var entry in figures)
        {
            string name = entry.GetProperty("name").GetString()!;
            Require(Path.GetFileName(name) == name, "A figure name must not leave the output directory");
            string source = Path.Combine(input, entry.GetProperty("file").GetString()!);
            using var sourceStream = File.OpenRead(source);
            using var gzip = new GZipStream(sourceStream, CompressionMode.Decompress);
            using var document = JsonDocument.Parse(gzip);
            var scene = document.RootElement;
            Require(scene.GetProperty("schema_version").GetInt32() == 1 && scene.GetProperty("name").GetString() == name,
                $"Publication scene identity mismatch: {name}");
            var axisRecords = scene.GetProperty("axes").EnumerateArray().ToArray();
            var seriesRecords = scene.GetProperty("series").EnumerateArray().ToArray();
            long rowCount = seriesRecords.Sum(s => (long)s.GetProperty("points").GetArrayLength());
            Require(axisRecords.Length == entry.GetProperty("axes").GetInt32() &&
                seriesRecords.Length == entry.GetProperty("series").GetInt32() &&
                rowCount == entry.GetProperty("coordinate_rows").GetInt64(), $"Publication coordinate inventory mismatch: {name}");
            var root = RenderNode(scene.GetProperty("svg"));
            Require(root.Name == Svg + "svg", $"Publication scene must have an SVG root: {name}");
            int sceneElements = root.DescendantsAndSelf().Count();
            root.SetAttributeValue("role", "img");
            root.AddFirst(new XElement(Svg + "desc", name + ". " + Scope));
            string svgPath = Path.Combine(output, name + ".svg");
            var settings = new XmlWriterSettings { Encoding = new UTF8Encoding(false), Indent = false, OmitXmlDeclaration = false };
            using (var writer = XmlWriter.Create(svgPath, settings))
                new XDocument(new XDeclaration("1.0", "utf-8", null), root).Save(writer);
            string csvPath = Path.Combine(output, name + ".series.csv");
            Csv.Write(csvPath, CoordinateRows(axisRecords, seriesRecords));
            string axesPath = Path.Combine(output, name + ".axes.json");
            WriteJson(axesPath, new
            {
                name, scope = Scope, axes = axisRecords,
                series = seriesRecords.Select((s, i) => new
                {
                    id = i,
                    metadata = s.EnumerateObject().Where(p => p.Name != "points").ToDictionary(p => p.Name, p => p.Value),
                    coordinate_rows = s.GetProperty("points").GetArrayLength()
                }),
                source_script = scene.GetProperty("source_script"),
                source_command = scene.GetProperty("source_command"),
                inputs_sha256 = scene.GetProperty("inputs_sha256"),
                archived_pdf_sha256 = scene.GetProperty("archived_pdf_sha256"),
                extraction_versions = scene.GetProperty("extraction_versions")
            });
            // Check serialization preserved the full typed scene and all internal
            // references. No browser or Python process is involved in replay.
            var saved = XDocument.Load(svgPath);
            Require(saved.Root!.DescendantsAndSelf().Count() == sceneElements + 1,
                $"SVG element count changed during replay: {name}");
            var ids = saved.Descendants().Attributes("id").Select(a => a.Value).ToHashSet();
            foreach (var href in saved.Descendants().Attributes().Where(a => a.Name.LocalName == "href"))
                Require(href.Value.StartsWith('#') ? ids.Contains(href.Value[1..]) : IsEmbeddedPng(href.Value), $"Unresolved SVG reference in {name}");
            files.AddRange([svgPath, csvPath, axesPath]);
            records.Add(new
            {
                name, source_file = Path.GetFileName(source), source_sha256 = Data.FileSha256(source),
                archived_pdf_sha256 = scene.GetProperty("archived_pdf_sha256").GetString(),
                source_script = scene.GetProperty("source_script").GetString(),
                axes = axisRecords.Length, series = seriesRecords.Length, coordinate_rows = rowCount,
                svg_elements = sceneElements, source_inputs = scene.GetProperty("inputs_sha256").EnumerateObject().Count(),
                outputs_sha256 = new[] { svgPath, csvPath, axesPath }.ToDictionary(path => Path.GetFileName(path)!, Data.FileSha256)
            });
            axes += axisRecords.Length; series += seriesRecords.Length; coordinates += rowCount;
        }
        string reportPath = Path.Combine(outputRoot, "publication", "figure_replay_report.json");
        WriteJson(reportPath, new
        {
            schema_version = 1, status = "pass", scope = Scope, frozen_scene = true,
            dynamic_science_regeneration = false, figure_count = figures.Length,
            result_compatibility_report = "result_compatibility.json",
            axis_count = axes, series_count = series, coordinate_rows = coordinates,
            input_manifest_sha256 = Data.FileSha256(Path.Combine(input, "manifest.json")),
            manuscript_sha256 = manifest.GetProperty("manuscript_sha256"),
            maintenance_recipe = "tools/extract_publication_scenes.py --reference <retained-reference-root> --output <data-root>/publication/figures",
            source_regeneration = "Run the retained Python renderers and re-extract scenes after changes to fitted states, predictions or statistics. The frozen scene replay does not propagate such changes.",
            figures = records
        });
        files.Add(reportPath);
        return new Report(figures.Length, axes, series, coordinates, reportPath, files.ToArray());
    }

    static IEnumerable<Dictionary<string, object?>> CoordinateRows(JsonElement[] axes, JsonElement[] series)
    {
        for (int id = 0; id < series.Length; id++)
        {
            var record = series[id];
            int axis = record.GetProperty("axis").GetInt32();
            Require(axis >= 0 && axis < axes.Length, "Coordinate refers to missing axis");
            var colors = record.TryGetProperty("color_values", out var c) ? c.EnumerateArray().ToArray() : [];
            int index = 0;
            foreach (var point in record.GetProperty("points").EnumerateArray())
            {
                Require(point.GetArrayLength() == 2, "Publication coordinate must have x and y");
                yield return new Dictionary<string, object?>
                {
                    ["axis"] = axis, ["axis_title"] = axes[axis].GetProperty("title").GetString(),
                    ["series"] = id, ["kind"] = record.GetProperty("kind").GetString(),
                    ["label"] = record.GetProperty("label").GetString(),
                    ["coordinate_system"] = record.GetProperty("coordinate_system").GetString(),
                    ["point"] = index, ["x"] = NumberOrNull(point[0]), ["y"] = NumberOrNull(point[1]),
                    ["color_value"] = index < colors.Length ? NumberOrNull(colors[index]) : null
                };
                index++;
            }
        }
    }

    static double? NumberOrNull(JsonElement value) => value.ValueKind == JsonValueKind.Null ? null : value.GetDouble();

    static XElement RenderNode(JsonElement node)
    {
        XName name = XName.Get(node.GetProperty("tag").GetString()!);
        Require(!ForbiddenElements.Contains(name.LocalName), "An immutable publication scene cannot contain active content");
        var element = new XElement(name);
        foreach (var attribute in node.GetProperty("attributes").EnumerateObject())
        {
            XName key = XName.Get(attribute.Name);
            string value = attribute.Value.GetString()!;
            Require(!key.LocalName.StartsWith("on", StringComparison.OrdinalIgnoreCase), "SVG event attributes are not permitted");
            if (key.LocalName == "href") Require(value.StartsWith('#') || IsEmbeddedPng(value), "External SVG resources are not permitted");
            element.SetAttributeValue(key, value);
        }
        if (node.GetProperty("text").ValueKind == JsonValueKind.String)
            element.Add(new XText(node.GetProperty("text").GetString()!));
        foreach (var child in node.GetProperty("children").EnumerateArray())
        {
            element.Add(RenderNode(child));
            if (child.GetProperty("tail").ValueKind == JsonValueKind.String)
                element.Add(new XText(child.GetProperty("tail").GetString()!));
        }
        return element;
    }

    // Matplotlib retains one self-contained raster colorbar; it is part of the
    // original figure scene and never causes a network or filesystem lookup.
    static bool IsEmbeddedPng(string value)
    {
        const string prefix = "data:image/png;base64,";
        if (!value.StartsWith(prefix, StringComparison.Ordinal)) return false;
        try
        {
            byte[] bytes = Convert.FromBase64String(value[prefix.Length..]);
            return bytes.AsSpan().StartsWith(new byte[] { 137, 80, 78, 71, 13, 10, 26, 10 });
        }
        catch (FormatException) { return false; }
    }

    static void WriteJson(string path, object value) => File.WriteAllText(path, JsonSerializer.Serialize(value, JsonOptions) + "\n", new UTF8Encoding(false));
    static void Require(bool condition, string message)
    {
        if (!condition) throw new InvalidDataException(message);
    }
}
