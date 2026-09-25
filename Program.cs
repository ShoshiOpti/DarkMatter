using System.Globalization;
using System.Text;
using System.Text.Json;
using DarkUniverse;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
CultureInfo.CurrentUICulture = CultureInfo.InvariantCulture;

try
{
    if (args.Contains("--help") || args.FirstOrDefault() == "help")
    {
        Console.WriteLine(
            "DarkUniverse.Console [all|import|statistics|pointwise|plots|verify] [--data PATH] [--out PATH]\n" +
            "Import only: --catalogue FILE.mrt|.zip --mass-models FILE.mrt|.zip\n" +
            "Default: all. Native C#; saved-fit reconstruction, no population refitting.");
        return 0;
    }

    var options = ParseOptions(args);
    string command = options.Command;
    string dataRoot = options.DataRoot;
    string outputRoot = options.OutputRoot;
    Directory.CreateDirectory(outputRoot);

    string reportPath = Path.Combine(outputRoot, "run_report.json");
    Data.SaveJson(reportPath, new
    {
        command,
        passed = false,
        status = "running"
    });
    var inputHashes = CheckInputManifest(dataRoot);
    var report = new Dictionary<string, object?>
    {
        ["command"] = command,
        ["runtime"] = Environment.Version.ToString(),
        ["mathnet_numerics"] = typeof(MathNet.Numerics.SpecialFunctions).Assembly.GetName().Version!.ToString(),
        ["data_directory"] = dataRoot,
        ["output_directory"] = outputRoot,
        ["scope"] = "Native import, saved-fit statistical reconstruction and plotting; no new population fit or PDE solve"
    };

    if (command is "all" or "import" or "verify")
    {
        Console.WriteLine("Importing SPARC tables...");
        var data = Sparc.Run(dataRoot, outputRoot, options.CataloguePath, options.MassModelsPath);
        Console.WriteLine($"Selected {data.SelectedNames.Length} galaxies / {data.SelectedPoints.Count()} radii.");
        if (options.CataloguePath == null && options.MassModelsPath == null)
            report["sparc"] = SparcChecks.Run(data, dataRoot, outputRoot);
    }

    string currentSummaryHtml = "";
    if (command is "all" or "statistics" or "pointwise" or "plots" or "verify")
    {
        Console.WriteLine("Reconstructing the current 25 September pointwise comparison...");
        report["current_analysis"] = "pointwise_comparison_2026_09_25";
        report["pointwise_inputs"] = PointwiseInputs.Verify(dataRoot, outputRoot);
        report["pointwise_statistics"] = PointwiseStatistics.Run(dataRoot, outputRoot);
        currentSummaryHtml = PointwiseReport.Write(outputRoot);
    }

    if (command is "all" or "statistics" or "verify")
    {
        Console.WriteLine("Reconstructing historical sign tests and auditing their saved predictions...");
        Statistics.Run(dataRoot, outputRoot);
        report["baryon_bootstrap"] = BaryonBootstrap.Run(dataRoot, outputRoot);
        report["statistics"] = "passed";
    }

    if (command is "all" or "plots" or "pointwise" or "verify")
    {
        Console.WriteLine("Generating native SVG plots and numerical sidecars...");
        bool currentOnly = command == "pointwise";
        var figures = RenderFigures(dataRoot, outputRoot, currentOnly);
        report["plots"] = figures.Select(path => RelativePath(outputRoot, path)).ToArray();
        WriteGallery(outputRoot, figures, currentSummaryHtml, currentOnly);
        Console.WriteLine($"Generated {figures.Count} figures.");
    }

    foreach (var (file, hash) in inputHashes)
    {
        if (Data.FileSha256(Path.Combine(dataRoot, file)) != hash)
            throw new InvalidDataException("Input changed during run: " + file);
    }

    report["consumed_csv_inputs"] = Csv.Inputs
        .Where(entry => entry.Key.StartsWith(dataRoot + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase))
        .ToDictionary(entry => RelativePath(dataRoot, entry.Key), entry => entry.Value);
    report["input_files_unchanged"] = inputHashes.Count;
    report["passed"] = true;
    Data.SaveJson(reportPath, report);
    Console.WriteLine("Completed. Results: " + outputRoot);
    return 0;
}
catch (Exception exception)
{
    Console.Error.WriteLine(exception.ToString());
    return 1;
}

static Options ParseOptions(string[] arguments)
{
    bool hasCommand = arguments.Length > 0 && !arguments[0].StartsWith('-');
    string command = hasCommand ? arguments[0] : "all";
    if (command is not ("all" or "import" or "statistics" or "pointwise" or "plots" or "verify"))
        throw new ArgumentException("Unknown command: " + command);

    var values = new Dictionary<string, string>();
    for (int i = hasCommand ? 1 : 0; i < arguments.Length; i += 2)
    {
        string name = arguments[i];
        if (name is not ("--data" or "--out" or "--catalogue" or "--mass-models"))
            throw new ArgumentException("Unknown option: " + name);
        if (i + 1 >= arguments.Length || arguments[i + 1].StartsWith("--"))
            throw new ArgumentException("Missing value for " + name);
        if (command != "import" && name is "--catalogue" or "--mass-models")
            throw new ArgumentException("Custom raw tables are supported by the import command.");
        values.TryAdd(name, arguments[i + 1]);
    }

    string project = FindProject();
    string dataRoot = Path.GetFullPath(values.GetValueOrDefault("--data") ?? Path.Combine(project, "data"));
    string outputRoot = Path.GetFullPath(values.GetValueOrDefault("--out") ?? Path.Combine(project, "output"));
    if (!Directory.Exists(dataRoot))
        throw new DirectoryNotFoundException(dataRoot);
    if (outputRoot.Equals(dataRoot, StringComparison.OrdinalIgnoreCase) ||
        outputRoot.StartsWith(dataRoot + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase))
        throw new ArgumentException("Output must be outside the input data directory.");

    return new Options(command, dataRoot, outputRoot,
        values.GetValueOrDefault("--catalogue"), values.GetValueOrDefault("--mass-models"));
}

static string FindProject()
{
    const string projectFile = "DarkUniverse.Console.csproj";
    foreach (string start in new[] { AppContext.BaseDirectory, Environment.CurrentDirectory })
    {
        for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
        {
            if (File.Exists(Path.Combine(directory.FullName, projectFile)))
                return directory.FullName;

            string nested = Path.Combine(directory.FullName, "DarkUniverse.Console");
            if (File.Exists(Path.Combine(nested, projectFile)))
                return nested;
        }
    }
    return Environment.CurrentDirectory;
}

static string RelativePath(string root, string path) => Path.GetRelativePath(root, path).Replace('\\', '/');

static Dictionary<string, string> CheckInputManifest(string dataRoot)
{
    var hashes = Directory.EnumerateFiles(dataRoot, "*", SearchOption.AllDirectories)
        .ToDictionary(path => RelativePath(dataRoot, path), Data.FileSha256);
    string manifestPath = Path.Combine(dataRoot, "input_manifest.json");
    if (!File.Exists(manifestPath))
        return hashes;

    var expected = JsonSerializer.Deserialize<Dictionary<string, string>>(File.ReadAllText(manifestPath))!;
    foreach (var (file, hash) in expected)
    {
        if (!hashes.TryGetValue(file, out var actual) || actual != hash)
            throw new InvalidDataException("Input manifest mismatch: " + file);
    }
    return hashes;
}

static List<string> RenderFigures(string dataRoot, string outputRoot, bool currentOnly)
{
    var figures = new List<string>(PointwisePlots.Run(dataRoot, outputRoot));
    if (!currentOnly)
    {
        string figureRoot = Path.Combine(outputRoot, "figures");
        figures.AddRange(ObservationPlots.Run(dataRoot, outputRoot));
        figures.AddRange(ControlledPlots.Run(dataRoot, figureRoot));
        figures.AddRange(GalaxyClassPlots.Run(dataRoot, outputRoot));
        figures.Add(LegacyPilot.Run(dataRoot, figureRoot));
    }

    int expectedCount = currentOnly ? 3 : 47;
    if (figures.Count != expectedCount || figures.Distinct().Count() != expectedCount)
        throw new InvalidDataException($"Expected {expectedCount} distinct plot variants.");
    return figures;
}

static void WriteGallery(string outputRoot, IEnumerable<string> figures, string currentSummaryHtml, bool currentOnly)
{
    var html = new StringBuilder(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Dark Universe — C# reconstruction</title>" +
        "<style>body{font:16px system-ui;max-width:1250px;margin:40px auto;padding:0 24px;color:#202b35}" +
        "section{border-top:1px solid #ddd;padding-top:20px;margin-top:32px}img{width:100%;height:auto}" +
        "a{color:#24658c}p{line-height:1.6}</style><h1>Dark Universe — C# reconstruction</h1>" +
        "<p>Native C# rendering of the supplied scientific results. Statistical tests are recomputed from saved fits " +
        "and per-radius predictions. Figures retain their observational, controlled, or historical scope. " +
        "Each CSV contains the plotted numerical series.</p>");
    html.Append(currentSummaryHtml);
    if (!currentOnly)
        html.Append("<p>The three current pointwise figures appear first. The remaining 44 figures preserve the earlier manuscript workflows; their sign tests answer a different question from the current mean-loss tests.</p>");

    foreach (string file in figures)
    {
        string relative = RelativePath(outputRoot, file);
        string name = Path.GetFileNameWithoutExtension(file);
        string title = System.Net.WebUtility.HtmlEncode(name);
        string series = Path.ChangeExtension(relative, "series.csv");
        html.Append($"<section><h2>{title}</h2><a href=\"{relative}\">SVG</a> · " +
            $"<a href=\"{series}\">Numerical series</a><img loading=\"lazy\" src=\"{relative}\" alt=\"{name}\"></section>");
    }

    html.Append("</html>");
    File.WriteAllText(Path.Combine(outputRoot, "index.html"), html.ToString());
}

sealed record Options(string Command, string DataRoot, string OutputRoot, string? CataloguePath, string? MassModelsPath);
