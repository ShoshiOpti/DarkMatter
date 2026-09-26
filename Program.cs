using System.CommandLine;
using System.Globalization;
using System.Text;
using System.Text.Json;
using DarkUniverse;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
CultureInfo.CurrentUICulture = CultureInfo.InvariantCulture;

return CreateCommandLine().Parse(args).Invoke();

static RootCommand CreateCommandLine()
{
    var dataOption = new Option<DirectoryInfo?>("--data")
    {
        Description = "Directory containing the input data. Defaults to the project's data directory.",
        Recursive = true
    };
    var outputOption = new Option<DirectoryInfo?>("--out")
    {
        Description = "Directory for generated reports and figures. Defaults to the project's output directory.",
        Recursive = true
    };

    var root = new RootCommand(
        "Reconstruct the Dark Universe analyses from saved fits and render their figures.")
    {
        dataOption,
        outputOption
    };
    root.SetAction(result => Run(() =>
        CreateOptions("all", result.GetValue(dataOption), result.GetValue(outputOption))));

    (string Name, string Description)[] commands =
    [
        ("all", "Run the complete reconstruction."),
        ("statistics", "Reconstruct the current and historical statistical analyses."),
        ("pointwise", "Reconstruct the current pointwise analysis and its three figures."),
        ("plots", "Reconstruct the current analysis and render all figures."),
        ("verify", "Run the complete reconstruction with all verification checks.")
    ];
    foreach (var (name, description) in commands)
    {
        var command = new Command(name, description);
        command.SetAction(result => Run(() =>
            CreateOptions(name, result.GetValue(dataOption), result.GetValue(outputOption))));
        root.Subcommands.Add(command);
    }

    var catalogueOption = new Option<FileInfo?>("--catalogue")
    {
        Description = "SPARC catalogue in MRT or ZIP format."
    };
    var massModelsOption = new Option<FileInfo?>("--mass-models")
    {
        Description = "SPARC mass-model table in MRT or ZIP format."
    };
    var importCommand = new Command("import", "Import the SPARC tables and verify the selected sample.")
    {
        catalogueOption,
        massModelsOption
    };
    importCommand.SetAction(result => Run(() => CreateOptions(
        "import",
        result.GetValue(dataOption),
        result.GetValue(outputOption),
        result.GetValue(catalogueOption),
        result.GetValue(massModelsOption))));
    root.Subcommands.Add(importCommand);

    return root;
}

static int Run(Func<Options> createOptions)
{
    try
    {
        Options options = createOptions();
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
    catch (ArgumentException exception)
    {
        Console.Error.WriteLine(exception.Message);
        return 1;
    }
    catch (DirectoryNotFoundException exception)
    {
        Console.Error.WriteLine($"Data directory not found: {exception.Message}");
        return 1;
    }
    catch (Exception exception)
    {
        Console.Error.WriteLine(exception.ToString());
        return 1;
    }
}

static Options CreateOptions(
    string command,
    DirectoryInfo? dataDirectory,
    DirectoryInfo? outputDirectory,
    FileInfo? catalogue = null,
    FileInfo? massModels = null)
{
    string project = FindProject();
    string dataRoot = dataDirectory?.FullName ?? Path.Combine(project, "data");
    string outputRoot = outputDirectory?.FullName ?? Path.Combine(project, "output");
    if (!Directory.Exists(dataRoot))
        throw new DirectoryNotFoundException(dataRoot);
    if (outputRoot.Equals(dataRoot, StringComparison.OrdinalIgnoreCase) ||
        outputRoot.StartsWith(dataRoot + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase))
        throw new ArgumentException("Output must be outside the input data directory.");

    return new Options(command, dataRoot, outputRoot, catalogue?.FullName, massModels?.FullName);
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
