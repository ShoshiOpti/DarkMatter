using System.CommandLine;

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
            Description = "Directory for generated reports and figures. Defaults to the project's output directory.",
            Recursive = true
        };

        var root = new RootCommand(
            "Reconstruct the Dark Universe analyses from saved fits and render their figures.")
        {
            dataOption,
            outputOption
        };

        root.SetAction(result => Execute(action, () =>
            ResolveOptions("all", result.GetValue(dataOption), result.GetValue(outputOption))));

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
            command.SetAction(result => Execute(action, () =>
                ResolveOptions(name, result.GetValue(dataOption), result.GetValue(outputOption))));
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
        if (ContainsPath(dataRoot, outputRoot))
            throw new ArgumentException("Output must be outside the input data directory.");

        return new CommandOptions(command, dataRoot, outputRoot, cataloguePath, massModelsPath);
    }

    /// <summary>Determines whether a candidate path is the given directory or one of its descendants.</summary>
    private static bool ContainsPath(DirectoryInfo directory, DirectoryInfo candidate)
    {
        string relativePath = Path.GetRelativePath(directory.FullName, candidate.FullName);
        return relativePath == "." ||
            (!Path.IsPathRooted(relativePath) &&
             relativePath != ".." &&
             !relativePath.StartsWith(".." + Path.DirectorySeparatorChar));
    }

    /// <summary>Locates the project directory from either the executable or current working directory.</summary>
    private static string FindProject()
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Environment.CurrentDirectory })
        {
            for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
            {
                if (File.Exists(Path.Combine(directory.FullName, "DarkUniverse.Console.csproj")))
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
