using System.ComponentModel;
using System.Runtime.InteropServices;
using Microsoft.Win32.SafeHandles;

namespace DarkUniverse;

public sealed record SafetyPaths(string DataRoot, string OutputRoot);

/// <summary>Filesystem boundaries for immutable inputs and isolated run outputs.</summary>
public static class ExecutionSafety
{
    static StringComparison PathComparison => OperatingSystem.IsWindows() ? StringComparison.OrdinalIgnoreCase : StringComparison.Ordinal;

    /// <summary>Resolve existing path components, reject overlapping roots and inspect existing output entries without following links.</summary>
    public static SafetyPaths ValidateRoots(string dataRoot, string outputRoot)
    {
        string originalData = Path.GetFullPath(dataRoot);
        if (!Directory.Exists(originalData)) throw new DirectoryNotFoundException(originalData);
        RejectReparseEntry(originalData);
        string input = CanonicalizePath(originalData), output = CanonicalizePath(outputRoot);
        if (ContainsPath(input, output) || ContainsPath(output, input))
            throw new ArgumentException("Input and output directories must not overlap, including resolved links and ancestor paths.");
        if (File.Exists(output)) throw new ArgumentException("The output directory is an existing file.");
        if (Directory.Exists(output)) AssertNoReparsePoints(output);
        return new SafetyPaths(input, output);
    }

    /// <summary>Resolve junctions/symbolic links in every existing component, including a nonexistent destination's ancestors.</summary>
    public static string CanonicalizePath(string path) => CanonicalizePath(path, 0);

    static string CanonicalizePath(string path, int depth)
    {
        if (depth > 40) throw new IOException("Too many filesystem links while resolving a path.");
        if (string.IsNullOrWhiteSpace(path)) throw new ArgumentException("A filesystem path is required.");
        string full = Path.GetFullPath(path);
        // Device namespaces give the same volume another root spelling. Reject
        // them before lexical containment checks rather than accept an alias
        // that could place a generated run inside the input directory.
        if (OperatingSystem.IsWindows() && (full.StartsWith(@"\\?\", StringComparison.Ordinal) || full.StartsWith(@"\\.\", StringComparison.Ordinal)))
            throw new ArgumentException("Windows device-namespace paths are not supported; use an ordinary drive or UNC path.");
        string root = Path.GetPathRoot(full)!;
        string current = root;
        string[] parts = full[root.Length..].Split([Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar], StringSplitOptions.RemoveEmptyEntries);
        for (int i = 0; i < parts.Length; i++)
        {
            current = Path.Combine(current, parts[i]);
            var attributes = ExistingAttributes(current);
            if (attributes is null) continue;
            if ((attributes.Value & FileAttributes.ReparsePoint) != 0)
            {
                FileSystemInfo info = (attributes.Value & FileAttributes.Directory) != 0 ? new DirectoryInfo(current) : new FileInfo(current);
                var target = info.ResolveLinkTarget(true) ?? throw new IOException("Unsupported filesystem reparse point: " + current);
                current = CanonicalizePath(target.FullName, depth + 1);
                if (!Directory.Exists(current) && !File.Exists(current)) throw new IOException("Dangling filesystem link: " + info.FullName);
            }
            if (i < parts.Length - 1 && File.Exists(current)) throw new IOException("A path component is a file: " + current);
        }
        return Path.TrimEndingDirectorySeparator(current);
    }

    public static bool ContainsPath(string directory, string candidate)
    {
        string parent = Path.TrimEndingDirectorySeparator(Path.GetFullPath(directory));
        string child = Path.TrimEndingDirectorySeparator(Path.GetFullPath(candidate));
        string prefix = Path.EndsInDirectorySeparator(parent) ? parent : parent + Path.DirectorySeparatorChar;
        return child.Equals(parent, PathComparison) || child.StartsWith(prefix, PathComparison);
    }

    /// <summary>Enumerate ordinary files only; never recurse through a link.</summary>
    public static IReadOnlyList<string> EnumerateRegularFiles(string root)
    {
        root = Path.GetFullPath(root);
        RejectReparseEntry(root);
        if (!Directory.Exists(root)) throw new DirectoryNotFoundException(root);
        var files = new List<string>();
        var pending = new Stack<string>();
        pending.Push(root);
        while (pending.Count > 0)
        {
            foreach (string entry in Directory.EnumerateFileSystemEntries(pending.Pop()))
            {
                FileAttributes attributes = File.GetAttributes(entry);
                if ((attributes & FileAttributes.ReparsePoint) != 0)
                    throw new IOException("Filesystem links are not permitted inside immutable inputs or output trees: " + entry);
                if ((attributes & FileAttributes.Directory) != 0) pending.Push(entry);
                else files.Add(entry);
            }
        }
        return files;
    }

    public static void AssertNoReparsePoints(string root) => _ = EnumerateRegularFiles(root);

    /// <summary>Ensure a prospective output remains below its run root and that existing components are ordinary filesystem entries.</summary>
    public static void AssertSafeWritePath(string outputRoot, string path)
    {
        string root = Path.TrimEndingDirectorySeparator(Path.GetFullPath(outputRoot));
        string full = Path.GetFullPath(path);
        if (!ContainsPath(root, full) || full.Equals(root, PathComparison))
            throw new IOException("An output file must stay inside the selected output directory.");
        string canonicalRoot = CanonicalizePath(root);
        if (!canonicalRoot.Equals(root, PathComparison))
            throw new IOException("The output directory changed to a filesystem link: " + root);
        string current = root;
        RejectReparseEntry(current);
        foreach (string part in Path.GetRelativePath(root, full).Split(Path.DirectorySeparatorChar))
        {
            current = Path.Combine(current, part);
            RejectReparseEntry(current);
        }
    }

    public static void RejectReparseEntry(string path)
    {
        var attributes = ExistingAttributes(path);
        if (attributes is not null && (attributes.Value & FileAttributes.ReparsePoint) != 0)
            throw new IOException("Filesystem links are not permitted at this path: " + path);
    }

    static FileAttributes? ExistingAttributes(string path)
    {
        try { return File.GetAttributes(path); }
        catch (FileNotFoundException) { return null; }
        catch (DirectoryNotFoundException) { return null; }
    }

    /// <summary>NTFS hardlinks can alias files outside the declared input root. Reject them instead of treating path spelling as isolation.</summary>
    public static void AssertUnaliasedInputFile(string path)
    {
        RejectReparseEntry(path);
        if (!OperatingSystem.IsWindows()) return;
        using SafeFileHandle handle = File.OpenHandle(path, FileMode.Open, FileAccess.Read, FileShare.Read);
        if (!GetFileInformationByHandle(handle, out var information))
            throw new IOException("Cannot establish input file identity: " + path, new Win32Exception(Marshal.GetLastWin32Error()));
        if (information.NumberOfLinks != 1)
            throw new IOException("Hardlinked input files are not permitted: " + path);
    }

    [StructLayout(LayoutKind.Sequential)]
    struct FileInformation
    {
        public uint Attributes;
        public uint CreationTimeLow, CreationTimeHigh, AccessTimeLow, AccessTimeHigh, WriteTimeLow, WriteTimeHigh;
        public uint VolumeSerial, SizeHigh, SizeLow, NumberOfLinks, IndexHigh, IndexLow;
    }

    [DllImport("kernel32.dll", SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    static extern bool GetFileInformationByHandle(SafeFileHandle file, out FileInformation information);
}
