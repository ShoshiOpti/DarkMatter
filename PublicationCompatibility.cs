using System.Globalization;
using System.Text.Json;
using Row = System.Collections.Generic.Dictionary<string, string>;

namespace DarkUniverse;

/// <summary>Refuses to display the frozen atlas beside incompatible newly computed results.</summary>
public static class PublicationCompatibility
{
    public static Dictionary<string, object?> Verify(string dataRoot, string outputRoot)
    {
        string contractPath = Path.Combine(dataRoot, "publication", "figures", "result_contract.json");
        var contract = Data.Json(contractPath);
        Require(contract.GetProperty("schema_version").GetInt32() == 1, "Unsupported figure/result contract.");
        foreach (var input in contract.GetProperty("input_sha256").EnumerateObject())
            Require(Data.FileSha256(Local(dataRoot, input.Name)) == input.Value.GetString(), "Figure source baseline changed: " + input.Name);
        var records = new List<object>();
        long checkedValues = 0;
        foreach (var table in contract.GetProperty("tables").EnumerateArray())
        {
            string referenceName = table.GetProperty("reference").GetString()!;
            string resultName = table.GetProperty("result").GetString()!;
            string expectedPath = Local(dataRoot, referenceName), actualPath = Local(outputRoot, resultName);
            string[] keys = table.GetProperty("keys").EnumerateArray().Select(x => x.GetString()!).ToArray();
            string[] fields = table.GetProperty("fields").EnumerateArray().Select(x => x.GetString()!).ToArray();
            var expected = Csv.Read(expectedPath);
            var actual = Csv.Read(actualPath);
            string Key(Row row) => string.Join("\u001f", keys.Select(k => Normalize(row[k])));
            Require(expected.Count == actual.Count, "Figure/result row count changed: " + resultName);
            var indexed = actual.ToDictionary(Key);
            foreach (var row in expected)
            {
                Require(indexed.TryGetValue(Key(row), out var got), "Figure/result row missing: " + resultName + "/" + Key(row));
                foreach (string field in fields)
                {
                    Require(got!.ContainsKey(field), "Figure/result field missing: " + field);
                    Require(Equivalent(row[field], got[field], field), $"Frozen figure disagrees with current result: {resultName}/{Key(row)}/{field}");
                    checkedValues++;
                }
            }
            records.Add(new { result = resultName, reference = referenceName, rows = actual.Count,
                result_sha256 = Data.FileSha256(actualPath), reference_sha256 = Data.FileSha256(expectedPath) });
        }
        var report = new Dictionary<string, object?>
        {
            ["status"] = "pass", ["baseline"] = contract.GetProperty("baseline").GetString(),
            ["contract_sha256"] = Data.FileSha256(contractPath), ["checked_values"] = checkedValues,
            ["scope"] = "Frozen scene/input identity and numerical agreement of newly generated statistical and closure tables with the scene baseline. Does not turn frozen physical plots into new simulations.",
            ["tolerance"] = "Counts/flags exact. Probabilities relative 2e-10. Numeric values absolute 2e-10 plus relative 2e-10; masses absolute 1e-35; six endpoint-sensitive phase fields absolute 3.1e-8 (archived CSV parser rounding).",
            ["tables"] = records
        };
        Data.SaveJson(Path.Combine(outputRoot, "publication", "result_compatibility.json"), report);
        return report;
    }

    public static bool Equivalent(string expected, string actual, string field)
    {
        if (expected == actual) return true;
        bool Missing(string x) => x.Length == 0 || x.Equals("nan", StringComparison.OrdinalIgnoreCase);
        if (Missing(expected) || Missing(actual)) return Missing(expected) && Missing(actual);
        if (bool.TryParse(expected, out var flag)) return bool.TryParse(actual, out var other) && flag == other;
        if (!double.TryParse(expected, NumberStyles.Float, CultureInfo.InvariantCulture, out double a) ||
            !double.TryParse(actual, NumberStyles.Float, CultureInfo.InvariantCulture, out double b)) return false;
        if (a == b) return true;
        if (!double.IsFinite(a) || !double.IsFinite(b)) return false;
        bool integer = field is "index" or "wins" or "losses" or "ties" or "joint_failures" or "informative_n" or "informative_signs"
            or "wild_draws" or "wild_exceedances" or "n_outer" or "wrong_sign_rows" or "observed_interval_node_crossings"
            || field.StartsWith("n_") || field.EndsWith("_galaxies") || field.EndsWith("_points") || field.EndsWith("_reversals");
        if (integer) return false;
        bool phase = field is "optimized_component1_over_sqrt_density" or "optimized_component2_over_sqrt_density"
            or "optimized_component2_fraction" or "quadrature_excess_contact_normalized" or "optimized_excess_contact_normalized"
            or "fractional_total_contact_reduction";
        double absolute = phase ? 3.1e-8 : field.Contains("mass_eV") ? 1e-35 :
            field.StartsWith("p_") || field.EndsWith("_p") ? 1e-45 : 2e-10;
        return Math.Abs(a - b) <= absolute + 2e-10 * Math.Abs(a);
    }

    static string Normalize(string value) => double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out double d)
        ? d.ToString("G17", CultureInfo.InvariantCulture) : value;
    static string Local(string root, string relative)
    {
        string path = Path.GetFullPath(Path.Combine(root, relative));
        Require(!Path.IsPathRooted(relative) && path.StartsWith(Path.GetFullPath(root) + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase), "Nonlocal figure contract path.");
        return path;
    }
    static void Require(bool condition, string message) { if (!condition) throw new InvalidDataException(message); }
}
