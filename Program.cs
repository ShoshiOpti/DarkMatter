using System.Globalization;
using System.Text;
using System.Text.Json;
using DarkUniverse;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
CultureInfo.CurrentUICulture = CultureInfo.InvariantCulture;

return Cli.Invoke(args, Run);

static int Run(CommandOptions options)
{
    var paths = ExecutionSafety.ValidateRoots(options.DataRoot.FullName, options.OutputRoot.FullName);
    using var session = RunSession.Start(paths, options.Command);
    Csv.Inputs.Clear();
    try
    {
        session.SetStage("verify_inputs");
        var inputs = InputManifest.Verify(paths.DataRoot);
        var report = new Dictionary<string, object?>
        {
            ["command"] = options.Command,
            ["runtime"] = Environment.Version.ToString(),
            ["mathnet_numerics"] = typeof(MathNet.Numerics.SpecialFunctions).Assembly.GetName().Version!.ToString(),
            ["data_directory"] = paths.DataRoot,
            ["output_directory"] = session.DirectoryPath,
            ["input_manifest_sha256"] = inputs.ManifestSha256,
            ["scope"] = "Frozen publication replay or explicitly versioned research; no independent observational validation is implied."
        };
        if (options.ProtocolPath is not null)
            report["study_protocol_sha256"] = Data.FileSha256(options.ProtocolPath);
        foreach (var (key, value) in report) session.Metadata[key] = value;
        if (options.Command is "calibrate" or "resolution" || options.Command.StartsWith("research-"))
        {
            session.SetStage(options.Command);
            object result = options.Command switch
            {
                "calibrate" => CalibrationStudy.Run(paths.DataRoot, session.DirectoryPath, ReadStudyProtocol(options.ProtocolPath)),
                "resolution" => MechanismResolution.Run(paths.DataRoot, session.DirectoryPath, ReadStudyProtocol(options.ProtocolPath)),
                "research-create" => ResearchWorkflow.CreateProtocol(paths.DataRoot, session.DirectoryPath),
                "research-lock" => ResearchWorkflow.LockProtocol(options.ProtocolPath!),
                "research-validate" => ResearchWorkflow.ValidateProtocol(options.ProtocolPath!),
                "research-evaluate" => ResearchWorkflow.Evaluate(options.ProtocolPath!, options.PredictionsPath!, options.FitArtifactPath!, session.DirectoryPath),
                "research-pilot" => ResearchWorkflow.RunPilot(options.ProtocolPath!, options.ReferencePath!, session.DirectoryPath, options.PythonExecutable!, options.PackagesPath),
                _ => throw new ArgumentException("Unknown research command.")
            };
            report["research"] = result;
            Data.SaveJson(Path.Combine(session.DirectoryPath, "research_command_report.json"), result);
            WriteResearchIndex(session.DirectoryPath, options.Command);
        }
        else RunPublication(options, session, report);
        session.SetStage("verify_preservation");
        inputs.VerifyUnchanged();
        if (options.ProtocolPath is not null && Data.FileSha256(options.ProtocolPath) != (string)report["study_protocol_sha256"]!)
            throw new InvalidDataException("Study protocol changed during execution.");
        report["consumed_csv_inputs"] = Csv.Inputs.ToDictionary(entry => entry.Key, entry => entry.Value);
        report["input_files_unchanged"] = inputs.Files.Count + 1;
        session.Complete(report);
        Console.WriteLine("Completed. Results: " + session.DirectoryPath);
        return 0;
    }
    catch (Exception exception)
    {
        session.Fail(exception);
        throw;
    }
}

static AnalysisStudyProtocol? ReadStudyProtocol(string? path) => path is null ? null :
    JsonSerializer.Deserialize<AnalysisStudyProtocol>(File.ReadAllText(path), new JsonSerializerOptions
    {
        PropertyNameCaseInsensitive = true,
        UnmappedMemberHandling = System.Text.Json.Serialization.JsonUnmappedMemberHandling.Disallow
    }) ?? throw new InvalidDataException("Empty study protocol.");

static void RunPublication(CommandOptions options, RunSession session, Dictionary<string, object?> report)
{
    string command = options.Command, dataRoot = options.DataRoot.FullName, outputRoot = session.DirectoryPath;
    if (command is "all" or "import" or "verify")
    {
        session.SetStage("sparc_import");
        Console.WriteLine("Importing SPARC tables...");
        var data = Sparc.Run(dataRoot, outputRoot, options.CataloguePath, options.MassModelsPath);
        Console.WriteLine($"Selected {data.SelectedNames.Length} galaxies / {data.SelectedPoints.Count()} radii.");
        if (options.CataloguePath == null && options.MassModelsPath == null)
            report["sparc"] = SparcChecks.Run(data, dataRoot, outputRoot);
    }
    string currentSummaryHtml = "";
    if (command is "publication" or "all" or "statistics" or "pointwise" or "plots" or "verify")
    {
        session.SetStage("empirical_statistics");
        Console.WriteLine("Reconstructing the original empirical pointwise comparison...");
        report["current_analysis"] = "pointwise_comparison_2026_09_25";
        report["pointwise_inputs"] = PointwiseInputs.Verify(dataRoot, outputRoot);
        report["pointwise_statistics"] = PointwiseStatistics.Run(dataRoot, outputRoot);
        currentSummaryHtml = PointwiseReport.Write(outputRoot);
    }
    if (command is "publication" or "all" or "statistics" or "plots" or "verify")
    {
        session.SetStage("publication_statistics");
        Console.WriteLine("Recomputing consolidated envelope, pressure and charged comparison families...");
        report["current_analysis"] = "consolidated_scalar_gr_2026_09_26";
        report["publication_statistics"] = PublicationStatistics.Run(dataRoot, outputRoot);
    }
    if (command is "publication" or "all" or "diagnostics" or "plots" or "verify")
    {
        session.SetStage("publication_diagnostics");
        Console.WriteLine("Recomputing frozen-target scalar and parent compatibility diagnostics...");
        report["publication_diagnostics"] = PublicationDiagnostics.Run(dataRoot, outputRoot);
    }
    if (command is "all" or "statistics" or "verify")
    {
        session.SetStage("historical_statistics");
        Statistics.Run(dataRoot, outputRoot);
        report["baryon_bootstrap"] = BaryonBootstrap.Run(dataRoot, outputRoot);
        report["statistics"] = "passed";
    }
    if (command is "all" or "plots" or "pointwise" or "verify")
    {
        session.SetStage("historical_plots");
        bool currentOnly = command == "pointwise";
        var figures = RenderFigures(dataRoot, outputRoot, currentOnly);
        report["plots"] = figures.Select(path => RelativePath(outputRoot, path)).ToArray();
        WriteGallery(outputRoot, figures, currentSummaryHtml, currentOnly);
    }
    if (command is "publication" or "all" or "plots" or "verify")
    {
        session.SetStage("current_result_plots");
        report["current_result_plots"] = CurrentResultsPlots.Run(outputRoot);
        session.SetStage("publication_figure_compatibility");
        var publicationFigures = PublicationPlots.Run(dataRoot, outputRoot);
        report["publication_figures"] = publicationFigures;
        WritePublicationGallery(outputRoot, publicationFigures.Files, currentSummaryHtml);
    }
}

static string RelativePath(string root, string path) => Path.GetRelativePath(root, path).Replace('\\', '/');

static void WriteResearchIndex(string outputRoot, string command)
{
    var html = new StringBuilder("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Research workflow</title><style>body{font:16px system-ui;max-width:950px;margin:40px auto;padding:0 20px;line-height:1.6}</style>");
    html.Append("<h1>" + System.Net.WebUtility.HtmlEncode(command) + "</h1><p>Versioned research output, separate from the frozen publication baseline. Synthetic calibration and perturbation pilots do not constitute independent observational validation. Read each report's scope and failure diagnostics.</p><ul>");
    foreach (string file in Directory.EnumerateFiles(outputRoot, "*", SearchOption.AllDirectories).Where(p => Path.GetExtension(p) is ".json" or ".csv" or ".md"))
    {
        string relative = RelativePath(outputRoot, file);
        html.Append($"<li><a href=\"{System.Net.WebUtility.HtmlEncode(relative)}\">{System.Net.WebUtility.HtmlEncode(relative)}</a></li>");
    }
    html.Append("</ul></html>");
    File.WriteAllText(Path.Combine(outputRoot, "index.html"), html.ToString());
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
    File.WriteAllText(Path.Combine(outputRoot, currentOnly ? "index.html" : "historical_figures.html"), html.ToString());
}

static void WritePublicationGallery(string outputRoot, IEnumerable<string> figures, string empiricalSummaryHtml)
{
    var html = new StringBuilder("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Scalar reduction: publication reproduction</title>" +
        "<style>body{font:16px system-ui;max-width:1250px;margin:40px auto;padding:0 24px;color:#202b35}p{line-height:1.6}section{border-top:1px solid #ddd;margin-top:32px;padding-top:20px}img{width:100%;height:auto}a{color:#24658c}</style>" +
        "<h1>Scalar reduction and internal stress in two-field models</h1>" +
        "<p>All three declared statistical families are recomputed natively from frozen predictions. Inverse algebra uses retained density derivatives. " +
        "The 25 manuscript figures replay audited frozen vector scenes with exact scientific-coordinate sidecars. " +
        "A changed fit requires regenerated scenes; these displays do not represent new physical evolution.</p>" +
        "<p><a href=\"run_report.json\">Run report</a> · <a href=\"publication_diagnostics/report.json\">Inverse and persistence audit</a> · " +
        "<a href=\"publication/figure_replay_report.json\">Figure provenance</a> · " +
        "<a href=\"publication_statistics/validation.json\">All three statistical families</a> · " +
        "<a href=\"publication_statistics/README.md\">Statistical interpretation</a></p>" +
        "<p><a href=\"current_results/index.html\">Current calculated comparison panels and tables</a> · " +
        "<a href=\"publication/result_compatibility.json\">Figure/result compatibility</a></p>" +
        "<h2>Original empirical-envelope comparison</h2>");
    html.Append(empiricalSummaryHtml);
    foreach (var path in figures.Where(path => path.EndsWith(".svg", StringComparison.OrdinalIgnoreCase)))
    {
        string relative = RelativePath(outputRoot, path);
        string label = System.Net.WebUtility.HtmlEncode(Path.GetFileNameWithoutExtension(path));
        html.Append($"<section><h2>{label}</h2><a href=\"{relative}\">SVG</a> · <a href=\"{Path.ChangeExtension(relative, "series.csv")}\">Scientific coordinates</a> · " +
            $"<a href=\"{Path.ChangeExtension(relative, "axes.json")}\">Axes and source metadata</a><img loading=\"lazy\" src=\"{relative}\" alt=\"{label}\"></section>");
    }
    html.Append("</html>");
    File.WriteAllText(Path.Combine(outputRoot, "index.html"), html.ToString());
}
