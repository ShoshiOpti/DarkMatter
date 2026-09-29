using System.CommandLine;
using DarkUniverse;

internal static class Cli
{
    /// <summary>Parses command-line arguments and runs the selected analysis.</summary>
    public static int Invoke(string[] args, Func<CommandOptions, int> action) =>
        CreateRootCommand(action).Parse(args).Invoke();

    /// <summary>Defines the available analyses and the options accepted by each one.</summary>
    private static RootCommand CreateRootCommand(Func<CommandOptions, int> action)
    {
        var dataOption = new Option<DirectoryInfo?>("--data")
        {
            Description = "Directory containing the input data. Defaults to the project's data directory.",
            Recursive = true
        };
        var outputOption = new Option<DirectoryInfo?>("--out")
        {
            Description = "Directory for generated reports and figures. Defaults to the project's output directory; each run gets a new subdirectory.",
            Recursive = true
        };

        var root = new RootCommand(
            "Reproduce the consolidated scalar-reduction paper from retained inputs.")
        {
            dataOption,
            outputOption
        };

        root.SetAction(result => Execute(action, () =>
            ResolveOptions("publication", result.GetValue(dataOption), result.GetValue(outputOption))));

        (string Name, string Description)[] commands =
        [
            ("sparc-full", "Refit all 175 SPARC galaxies natively; reproduce full-profile figures, outer inference and green-fit significance bounds."),
            ("publication", "Recompute three current statistical families and inverse diagnostics; replay all 25 publication figures."),
            ("diagnostics", "Recompute frozen-target inverse algebra and summarize retained persistence audits."),
            ("reconstruction", "Check parent/relative field equations, hidden stress and saved-curve density bookkeeping."),
            ("all", "Run publication and historical reconstruction, plus the latest full-SPARC native refits and statistics."),
            ("statistics", "Reconstruct the current and historical statistical analyses."),
            ("pointwise", "Reconstruct the original empirical pointwise analysis and its three figures."),
            ("plots", "Recompute current inference and reproduce current and historical figures."),
            ("verify", "Run complete historical/publication checks and the latest full-SPARC refit, plot and statistical verification.")
        ];
        foreach (var (name, description) in commands)
        {
            var command = new Command(name, description);
            command.SetAction(result => Execute(action, () =>
                ResolveOptions(name, result.GetValue(dataOption), result.GetValue(outputOption))));
            root.Subcommands.Add(command);
        }

        var experimentOption = new Option<FileInfo?>("--experiment")
        {
            Description = "Versioned conditional formation experiment bundle with hashed source JSON files.",
            Required = true
        };
        var formationCommand = new Command("formation", "Recompute descriptive statistics and plots for prescribed formation experiments; no galaxy population inference.") { experimentOption };
        formationCommand.SetAction(result => Execute(action, () =>
            ResolveOptions("formation", result.GetValue(dataOption), result.GetValue(outputOption)) with
            {
                ExperimentPath = result.GetValue(experimentOption)!.FullName
            }));
        root.Subcommands.Add(formationCommand);

        foreach (string name in new[] { "calibrate", "resolution", "research-create", "research-lock", "research-validate", "research-evaluate", "research-pilot" })
        {
            var protocolOption = new Option<FileInfo?>("--protocol") { Description = "Versioned study or research protocol JSON file.", Required = name.StartsWith("research-") && name != "research-create" };
            var predictionsOption = new Option<FileInfo?>("--predictions") { Description = "Fresh prediction CSV for the locked research protocol.", Required = name == "research-evaluate" };
            var artifactOption = new Option<FileInfo?>("--fit-artifact") { Description = "Fitter provenance and training contract JSON.", Required = name == "research-evaluate" };
            var referenceOption = new Option<DirectoryInfo?>("--reference") { Description = "Retained scientific reference tree (required for the Python pilot)." };
            var pythonOption = new Option<string?>("--python") { Description = "Python executable for the optional empirical-refit pilot." };
            var packagesOption = new Option<DirectoryInfo?>("--packages") { Description = "Optional installed Python scientific-package directory." };
            var descriptions = new Dictionary<string, string>
            {
                ["calibrate"] = "Run versioned synthetic null/coverage calibration; does not alter published statistics.",
                ["resolution"] = "Measure charged/free differences against explicitly configured numerical/practical budgets.",
                ["research-create"] = "Create an editable development protocol and separated training/evaluation inputs.",
                ["research-lock"] = "Validate and hash-lock a research protocol and its input/split files.",
                ["research-validate"] = "Validate the locked research contract without fitting.",
                ["research-evaluate"] = "Evaluate fresh predictions against the locked protocol and fitter artifact.",
                ["research-pilot"] = "Run the original Python empirical fit pipeline with paired uncertainty perturbations."
            };
            var command = new Command(name, descriptions[name]) { protocolOption };
            if (name == "research-evaluate") { command.Options.Add(predictionsOption); command.Options.Add(artifactOption); }
            if (name == "research-pilot") { command.Options.Add(referenceOption); command.Options.Add(pythonOption); command.Options.Add(packagesOption); }
            command.SetAction(result => Execute(action, () => ResolveOptions(name, result.GetValue(dataOption), result.GetValue(outputOption)) with
            {
                ProtocolPath = result.GetValue(protocolOption)?.FullName,
                PredictionsPath = name == "research-evaluate" ? result.GetValue(predictionsOption)?.FullName : null,
                FitArtifactPath = name == "research-evaluate" ? result.GetValue(artifactOption)?.FullName : null,
                ReferencePath = name == "research-pilot" ? result.GetValue(referenceOption)?.FullName ?? Path.GetFullPath(Path.Combine(FindProject(), "..", "reference")) : null,
                PythonExecutable = name == "research-pilot" ? result.GetValue(pythonOption) ?? "python" : null,
                PackagesPath = name == "research-pilot" ? result.GetValue(packagesOption)?.FullName : null
            }));
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
        importCommand.SetAction(result => Execute(action, () => ResolveOptions(
            "import",
            result.GetValue(dataOption),
            result.GetValue(outputOption),
            result.GetValue(catalogueOption)?.FullName,
            result.GetValue(massModelsOption)?.FullName)));
        root.Subcommands.Add(importCommand);

        return root;
    }

    /// <summary>Runs an analysis and converts validation or runtime failures into CLI exit codes.</summary>
    private static int Execute(Func<CommandOptions, int> action, Func<CommandOptions> createOptions)
    {
        try
        {
            return action(createOptions());
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

    /// <summary>Resolves CLI paths and rejects missing input or output paths that could overwrite it.</summary>
    private static CommandOptions ResolveOptions(
        string command,
        DirectoryInfo? dataDirectory,
        DirectoryInfo? outputDirectory,
        string? cataloguePath = null,
        string? massModelsPath = null)
    {
        string project = FindProject();
        var dataRoot = dataDirectory ?? new DirectoryInfo(Path.Combine(project, "data"));
        var outputRoot = outputDirectory ?? new DirectoryInfo(Path.Combine(project, "output"));

        if (!dataRoot.Exists)
            throw new DirectoryNotFoundException(dataRoot.FullName);
        var safe = ExecutionSafety.ValidateRoots(dataRoot.FullName, outputRoot.FullName);
        dataRoot = new DirectoryInfo(safe.DataRoot);
        outputRoot = new DirectoryInfo(safe.OutputRoot);

        return new CommandOptions(command, dataRoot, outputRoot, cataloguePath, massModelsPath);
    }

    /// <summary>Locates the project directory from either the executable or current working directory.</summary>
    private static string FindProject()
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Environment.CurrentDirectory })
        {
            for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
            {
                if (File.Exists(Path.Combine(directory.FullName, "DarkUniverse.Console.csproj")) ||
                    File.Exists(Path.Combine(directory.FullName, "data", "input_manifest.json")))
                    return directory.FullName;
            }
        }

        return Environment.CurrentDirectory;
    }
}

internal sealed record CommandOptions
{
    public string Command { get; }
    public DirectoryInfo DataRoot { get; }
    public DirectoryInfo OutputRoot { get; }
    public string? CataloguePath { get; }
    public string? MassModelsPath { get; }
    public string? ProtocolPath { get; init; }
    public string? PredictionsPath { get; init; }
    public string? FitArtifactPath { get; init; }
    public string? ReferencePath { get; init; }
    public string? PythonExecutable { get; init; }
    public string? PackagesPath { get; init; }
    public string? ExperimentPath { get; init; }

    /// <summary>Captures the validated paths and inputs for one requested analysis.</summary>
    public CommandOptions(string command, DirectoryInfo dataRoot, DirectoryInfo outputRoot, string? cataloguePath = null, string? massModelsPath = null)
    {
        Command = command;
        DataRoot = dataRoot;
        OutputRoot = outputRoot;
        CataloguePath = cataloguePath;
        MassModelsPath = massModelsPath;
    }

}
