using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using Microsoft.VisualBasic.FileIO;

namespace DarkUniverse;

public static class Csv
{
    public static readonly Dictionary<string, string> Inputs = new(StringComparer.OrdinalIgnoreCase);

    public static List<Dictionary<string, string>> Read(string path)
    {
        Inputs[Path.GetFullPath(path)] = Data.FileSha256(path);
        using var parser = new TextFieldParser(path, Encoding.UTF8)
        {
            TextFieldType = FieldType.Delimited,
            HasFieldsEnclosedInQuotes = true
        };
        parser.SetDelimiters(",");
        var headers = parser.ReadFields() ?? throw new InvalidDataException($"Empty CSV: {path}");
        if (headers.Distinct().Count() != headers.Length)
            throw new InvalidDataException($"Duplicate CSV headers: {path}");

        var rows = new List<Dictionary<string, string>>();
        while (!parser.EndOfData)
        {
            var fields = parser.ReadFields()!;
            if (fields.Length != headers.Length)
                throw new InvalidDataException($"Wrong column count: {path}:{parser.LineNumber}");
            rows.Add(headers.Zip(fields).ToDictionary(field => field.First, field => field.Second));
        }
        return rows;
    }

    public static string Text(IReadOnlyDictionary<string, string> row, string key) =>
        row.TryGetValue(key, out var value) ? value : throw new InvalidDataException($"Missing column: {key}");

    public static double Number(IReadOnlyDictionary<string, string> row, string key)
    {
        string value = Text(row, key).Trim();
        return value switch
        {
            "" or "None" or "null" or "nan" or "NaN" => double.NaN,
            "inf" or "+inf" or "Infinity" => double.PositiveInfinity,
            "-inf" or "-Infinity" => double.NegativeInfinity,
            _ => double.Parse(value, CultureInfo.InvariantCulture)
        };
    }

    public static bool Flag(IReadOnlyDictionary<string, string> row, string key) =>
        Text(row, key).Trim().ToLowerInvariant() is "true" or "1";

    public static void Write(string path, IEnumerable<Dictionary<string, object?>> rows)
    {
        Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(path))!);
        using var writer = new StreamWriter(path, false, new UTF8Encoding(false));
        string[]? columns = null;
        foreach (var row in rows)
        {
            if (columns is null)
            {
                columns = row.Keys.ToArray();
                writer.WriteLine(string.Join(',', columns.Select(Escape)));
            }
            writer.WriteLine(string.Join(',', columns.Select(column => Escape(Format(row.GetValueOrDefault(column))))));
        }
    }

    static string Format(object? value) => value is IFormattable number
        ? number.ToString(null, CultureInfo.InvariantCulture)
        : value?.ToString() ?? "";

    static string Escape(string value) => value.IndexOfAny([',', '"', '\r', '\n']) >= 0
        ? "\"" + value.Replace("\"", "\"\"") + "\""
        : value;
}

public static class Data
{
    static readonly JsonSerializerOptions JsonOptions = new()
    {
        WriteIndented = true,
        NumberHandling = System.Text.Json.Serialization.JsonNumberHandling.AllowNamedFloatingPointLiterals
    };

    public static string FileSha256(string path) => Convert.ToHexStringLower(SHA256.HashData(File.ReadAllBytes(path)));

    public static double Median(IEnumerable<double> values) => Quantile(values, .5);

    public static double Quantile(IEnumerable<double> values, double probability)
    {
        var sorted = values.Where(value => !double.IsNaN(value)).Order().ToArray();
        if (sorted.Length == 0)
            return double.NaN;

        // Linear interpolation between adjacent order statistics.
        double position = probability * (sorted.Length - 1);
        int index = (int)position;
        if (index == sorted.Length - 1 || position == index)
            return sorted[index];
        return sorted[index] + (sorted[index + 1] - sorted[index]) * (position - index);
    }

    public static double[] Linspace(double min, double max, int count) =>
        Enumerable.Range(0, count).Select(index => min + (max - min) * index / (count - 1)).ToArray();

    public static double[] Column(this IEnumerable<Dictionary<string, string>> rows, string key) =>
        rows.Select(row => Csv.Number(row, key)).ToArray();

    public static JsonElement Json(string path)
    {
        using var document = JsonDocument.Parse(File.ReadAllText(path));
        return document.RootElement.Clone();
    }

    public static void SaveJson(string path, object data)
    {
        Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(path))!);
        File.WriteAllText(path, JsonSerializer.Serialize(data, JsonOptions));
    }
}
