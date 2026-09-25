using System.Globalization;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.Json;

namespace DarkUniverse;

public sealed record SparcGalaxy(string Name, int Type, double DistanceMpc, double DistanceErrorMpc,
    int DistanceMethod, double InclinationDeg, double InclinationErrorDeg, double Luminosity36,
    double LuminosityError36, double EffectiveRadiusKpc, double EffectiveSurfaceBrightness,
    double DiskScaleKpc, double DiskSurfaceBrightness, double HiMass, double HiRadiusKpc,
    double FlatSpeedKms, double FlatSpeedErrorKms, int Quality, string References);

public sealed record SparcPoint(string Galaxy, double DistanceMpc, double RadiusKpc,
    double ObservedSpeedKms, double ErrorSpeedKms, double GasSpeedKms, double DiskSpeedKms,
    double BulgeSpeedKms, double DiskSurfaceBrightness, double BulgeSurfaceBrightness)
{
    // SPARC signed components represent inward or outward radial acceleration: V * |V|.
    public double BaryonsSquared(double diskMassToLight = .5, double bulgeMassToLight = .7) =>
        GasSpeedKms * Math.Abs(GasSpeedKms) + diskMassToLight * DiskSpeedKms * Math.Abs(DiskSpeedKms)
        + bulgeMassToLight * BulgeSpeedKms * Math.Abs(BulgeSpeedKms);
}

public sealed record SparcDataset(IReadOnlyDictionary<string, SparcGalaxy> Catalogue,
    IReadOnlyDictionary<string, SparcPoint[]> Points, string[] SelectedNames,
    string CatalogueSha256, string MassModelsSha256)
{
    public IEnumerable<SparcPoint> SelectedPoints => SelectedNames.SelectMany(n => Points[n]);
    public string[] DisplayNames
    {
        get
        {
            var displayIndices = Enumerable.Range(0, 6)
                .Select(index => (int)Math.Floor((index + .5) * SelectedNames.Length / 6)).ToHashSet();
            return SelectedNames.OrderBy(name => Sparc.Median(Points[name].Select(point => point.ObservedSpeedKms)))
                .ThenBy(name => name, StringComparer.Ordinal)
                .Where((_, index) => displayIndices.Contains(index)).ToArray();
        }
    }
}

public static class Sparc
{
    static readonly CultureInfo Invariant = CultureInfo.InvariantCulture;
    static readonly JsonSerializerOptions JsonOptions = new() { WriteIndented = true };
    public const string CatalogueFile = "SPARC_Lelli2016c.mrt";
    public const string MassModelsFile = "MassModels_Lelli2016c.mrt";

    public static SparcDataset Load(string dataRoot, string? cataloguePath = null, string? massModelsPath = null)
    {
        var catalogueBytes = ReadTable(cataloguePath ?? FindFile(dataRoot, "numerics/observations/data/" + CatalogueFile), CatalogueFile);
        var massBytes = ReadTable(massModelsPath ?? FindFile(dataRoot, "numerics/observations/data/" + MassModelsFile), MassModelsFile);
        var catalogue = new Dictionary<string, SparcGalaxy>(StringComparer.Ordinal);
        bool started = false;
        foreach (var line in Lines(catalogueBytes))
        {
            if (string.IsNullOrWhiteSpace(line))
                continue;
            string[] p = Fields(line);
            if (!IsCatalogueRow(p))
            {
                var fixedFields = FixedFields(line, [11, 2, 6, 5, 2, 4, 4, 7, 7, 5, 8, 5, 8, 7, 5, 5, 5, 3, 14]);
                if (!IsCatalogueRow(fixedFields))
                {
                    if (started)
                        throw new InvalidDataException("Malformed catalogue row: " + line);
                    continue;
                }
                p = fixedFields;
            }
            started = true;
            var row = new SparcGalaxy(p[0], I(p[1]), D(p[2]), D(p[3]), I(p[4]), D(p[5]), D(p[6]), D(p[7]),
                D(p[8]), D(p[9]), D(p[10]), D(p[11]), D(p[12]), D(p[13]), D(p[14]), D(p[15]), D(p[16]), I(p[17]), p[18]);
            if (!catalogue.TryAdd(row.Name, row))
                throw new InvalidDataException("Duplicate catalogue galaxy: " + row.Name);
        }
        if (catalogue.Count == 0)
            throw new InvalidDataException("No SPARC catalogue rows found.");
        var points = catalogue.Keys.ToDictionary(n => n, _ => new List<SparcPoint>(), StringComparer.Ordinal);
        started = false;
        foreach (var line in Lines(massBytes))
        {
            if (string.IsNullOrWhiteSpace(line))
                continue;
            string[] p = Fields(line);
            if (p.Length != 10 || !catalogue.ContainsKey(p[0]) || !p.Skip(1).All(IsNumber))
            {
                var fixedFields = MassFields(line);
                if (fixedFields.Length != 10 || !catalogue.ContainsKey(fixedFields[0]) || !fixedFields.Skip(1).All(IsNumber))
                {
                    if (started)
                        throw new InvalidDataException("Malformed mass-model row or unknown galaxy: " + line);
                    continue;
                }
                p = fixedFields;
            }
            started = true;
            var row = new SparcPoint(p[0], D(p[1]), D(p[2]), D(p[3]), D(p[4]), D(p[5]), D(p[6]), D(p[7]), D(p[8]), D(p[9]));
            if (row.RadiusKpc <= 0 || row.ErrorSpeedKms <= 0 || row.DistanceMpc <= 0)
                throw new InvalidDataException("Non-positive radius, error, or distance for " + row.Galaxy);
            if (Math.Abs(row.DistanceMpc - catalogue[row.Galaxy].DistanceMpc) > 1e-8)
                throw new InvalidDataException("Catalogue and radial distances disagree for " + row.Galaxy);
            var list = points[row.Galaxy];
            if (list.Count != 0 && list[^1].RadiusKpc >= row.RadiusKpc)
                throw new InvalidDataException("Radial rows must be strictly increasing for " + row.Galaxy);
            list.Add(row);
        }
        if (points.Any(p => p.Value.Count == 0))
            throw new InvalidDataException("Missing mass-model rows for a catalogue galaxy.");
        var selected = catalogue.Values.Where(m => m.Quality < 3 && m.InclinationDeg >= 30)
            .Select(m => m.Name).Order(StringComparer.Ordinal).ToArray();
        return new SparcDataset(catalogue, points.ToDictionary(p => p.Key, p => p.Value.ToArray(), StringComparer.Ordinal), selected,
            Convert.ToHexStringLower(SHA256.HashData(catalogueBytes)), Convert.ToHexStringLower(SHA256.HashData(massBytes)));
    }

    public static SparcDataset Run(string dataRoot, string outputRoot, string? cataloguePath = null, string? massModelsPath = null)
    {
        var data = Load(dataRoot, cataloguePath, massModelsPath);
        string folder = Path.Combine(outputRoot, "sparc");
        Directory.CreateDirectory(folder);
        SaveJson(Path.Combine(folder, "catalogue.json"), data.Catalogue.Values);
        SaveJson(Path.Combine(folder, "selected_metadata.json"), data.SelectedNames.Select(n =>
        {
            var m = data.Catalogue[n];
            return new
            {
                name = m.Name,
                type = m.Type,
                D = m.DistanceMpc,
                e_D = m.DistanceErrorMpc,
                inc = m.InclinationDeg,
                e_inc = m.InclinationErrorDeg,
                L36 = m.Luminosity36,
                Q = m.Quality,
                refs = m.References
            };
        }));
        using (var writer = new StreamWriter(Path.Combine(folder, "selected_points.csv")))
        {
            writer.WriteLine("galaxy,D_Mpc,r_kpc,Vobs_kms,e_Vobs_kms,Vgas_kms,Vdisk_kms,Vbul_kms,SBdisk,SBbul,Vbar2_kms2,Vbar_kms,partition");
            foreach (string name in data.SelectedNames)
            {
                var rows = data.Points[name];
                int holdout = HoldoutCount(rows.Length);
                for (int i = 0; i < rows.Length; i++)
                {
                    var p = rows[i];
                    double b = p.BaryonsSquared();
                    writer.WriteLine(string.Join(',', p.Galaxy, F(p.DistanceMpc), F(p.RadiusKpc), F(p.ObservedSpeedKms), F(p.ErrorSpeedKms),
                        F(p.GasSpeedKms), F(p.DiskSpeedKms), F(p.BulgeSpeedKms), F(p.DiskSurfaceBrightness), F(p.BulgeSurfaceBrightness),
                        F(b), b >= 0 ? F(Math.Sqrt(b)) : "", holdout == 0 ? "not_eligible" : i < rows.Length - holdout ? "train" : "test"));
                }
            }
        }
        using (var writer = new StreamWriter(Path.Combine(folder, "holdout_split.csv")))
        {
            writer.WriteLine("galaxy,n,n_train,n_test,eligible");
            foreach (string name in data.SelectedNames)
            {
                int n = data.Points[name].Length, test = HoldoutCount(n);
                writer.WriteLine($"{name},{n},{n - test},{test},{(test > 0 ? "True" : "False")}");
            }
        }
        SaveJson(Path.Combine(folder, "summary.json"), new
        {
            catalogue_galaxies = data.Catalogue.Count,
            catalogue_rows = data.Points.Values.Sum(p => p.Length),
            selected_galaxies = data.SelectedNames.Length,
            selected_rows = data.SelectedPoints.Count(),
            holdout_galaxies = data.SelectedNames.Count(n => HoldoutCount(data.Points[n].Length) > 0),
            holdout_rows = data.SelectedNames.Sum(n => HoldoutCount(data.Points[n].Length)),
            negative_gas_rows = data.Points.Values.Sum(p => p.Count(r => r.GasSpeedKms < 0)),
            negative_selected_baryon_rows = data.SelectedPoints.Count(p => p.BaryonsSquared() < 0),
            selection = "Q < 3 and inclination >= 30 degrees; no velocity-error cut",
            holdout = "N >= 8; outer max(2, ceil(0.2*N)) radii held out",
            displayed = data.DisplayNames,
            catalogue_sha256 = data.CatalogueSha256,
            mass_models_sha256 = data.MassModelsSha256
        });
        return data;
    }

    // Freeze the original outer-radius split before any prediction or error comparison.
    public static int HoldoutCount(int rowCount) => rowCount >= 8 ? Math.Max(2, (int)Math.Ceiling(.2 * rowCount)) : 0;
    public static double Median(IEnumerable<double> values)
    {
        var a = values.Order().ToArray();
        if (a.Length == 0)
            throw new ArgumentException("Median requires a nonempty sequence.");
        return a.Length % 2 == 1 ? a[a.Length / 2] : (a[a.Length / 2 - 1] + a[a.Length / 2]) / 2;
    }
    internal static string FindFile(string root, string relative)
    {
        string path = Path.Combine(root, relative.Replace('/', Path.DirectorySeparatorChar));
        if (File.Exists(path))
            return path;
        string direct = Path.Combine(root, Path.GetFileName(relative));
        if (File.Exists(direct))
            return direct;
        throw new FileNotFoundException("Required SPARC input was not found: " + path, path);
    }
    internal static void SaveJson(string path, object value) => File.WriteAllText(path, JsonSerializer.Serialize(value, JsonOptions) + Environment.NewLine);
    static string F(double value) => value.ToString("G17", Invariant);
    static double D(string value) => double.Parse(value, NumberStyles.Float, Invariant);
    static int I(string value) => int.Parse(value, NumberStyles.Integer, Invariant);
    static bool IsNumber(string value) => double.TryParse(value, NumberStyles.Float, Invariant, out var x) && double.IsFinite(x);
    static bool IsCatalogueRow(string[] p) => p.Length == 19 && p[0].Length > 0 && p.Skip(1).Take(17).All(IsNumber)
        && int.TryParse(p[1], out _) && int.TryParse(p[4], out _) && int.TryParse(p[17], out _);
    static string[] Fields(string line) => line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
    static IEnumerable<string> Lines(byte[] bytes)
    {
        using var reader = new StreamReader(new MemoryStream(bytes));
        while (reader.ReadLine() is { } line)
            yield return line;
    }
    static byte[] ReadTable(string path, string tableName)
    {
        if (!Path.GetExtension(path).Equals(".zip", StringComparison.OrdinalIgnoreCase))
            return File.ReadAllBytes(path);
        using var archive = ZipFile.OpenRead(path);
        var entries = archive.Entries.Where(e => e.Name.Equals(tableName, StringComparison.OrdinalIgnoreCase)).ToArray();
        if (entries.Length != 1)
            throw new InvalidDataException($"ZIP must contain exactly one {tableName} table.");
        using var source = entries[0].Open();
        using var target = new MemoryStream();
        source.CopyTo(target);
        return target.ToArray();
    }
    static string[] FixedFields(string line, int[] widths)
    {
        if (line.Length < widths.Sum())
            return [];
        int offset = 0;
        return widths.Select(width =>
        {
            string field = line.Substring(offset, width).Trim();
            offset += width;
            return field;
        }).ToArray();
    }
    static string[] MassFields(string line)
    {
        int[] starts = [0, 12, 19, 26, 33, 39, 46, 53, 60, 68], lengths = [11, 6, 6, 6, 5, 6, 6, 6, 7, 8];
        if (line.Length < 76)
            return [];
        return starts.Select((start, i) => line.Substring(start, lengths[i]).Trim()).ToArray();
    }
}
