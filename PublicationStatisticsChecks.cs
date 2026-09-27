using System.Text.Json;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public static partial class PublicationStatistics
{
    static bool BlockFailed(Prediction[] predictions) => predictions.Any(p => p.Failed) ||
        predictions.Any(p => p.Observation.Outer && (!double.IsFinite(p.Value) || p.Value < 0));

    // Small adversarial checks protect the scientific failure and multiplicity conventions.
    // They do not repeat the frozen results used by the independent table audits.
    static Result SelfCheck(JsonElement rng)
    {
        var outer = new Observation("fixture", 1, 2, 3, 1, true);
        var inner = new Observation("fixture", 0, 1, 2, 1, false);
        Require(!BlockFailed([new(outer, "m", 3, false), new(inner, "m", 2, false)]), "Finite block failed.");
        Require(BlockFailed([new(outer, "m", 3, false), new(inner, "m", 2, true)]), "Failed training fit was hidden by finite outer predictions.");
        Require(BlockFailed([new(outer, "m", -1, false)]), "Negative outer speed was accepted.");
        Require(BlockFailed([new(outer, "m", double.NaN, false)]), "Missing outer speed was accepted.");
        Require(BlockFailed([new(outer, "m", double.PositiveInfinity, false)]), "Infinite outer speed was accepted.");
        var sign = Sign([1, double.PositiveInfinity, double.PositiveInfinity, 2, 2, 2],
            [double.PositiveInfinity, 1, double.PositiveInfinity, 2.00001, 3, 1], 1e-4);
        Require((int)sign["wins"]! == 2 && (int)sign["losses"]! == 2 && (int)sign["ties"]! == 1 &&
            (int)sign["joint_failures"]! == 1 && (int)sign["informative_signs"]! == 4, "Failure-aware sign accounting changed.");
        Require((double)Sign([1, 2], [1, 2], 1e-4)["p_sign_raw"]! == 1, "All-tie sign probability changed.");
        Require(!EqualObservation(outer, outer with { Sigma = 2 }), "Changed measurement error was accepted.");
        Require(!EqualObservation(outer, outer with { Outer = false }), "Changed experimental split was accepted.");
        var adjusted = Statistics.Holm([.01, 1, .03, .02, .04, .05]);
        Require(Near(adjusted[0], .06) && adjusted[1] == 1, "Undefined-contrast Holm placeholder removed a family slot.");
        bool refused = false;
        try { _ = Bootstrap([[1, 1]], [1, 1], 0, null, rng); }
        catch (InvalidDataException ex) when (ex.Message.Contains("constant")) { refused = true; }
        Require(refused, "Nonzero constant effects must fail studentized inference.");
        return new() { ["status"] = "pass", ["checks"] = 11,
            ["coverage"] = "training failure propagation; negative/nonfinite outer speeds; finite-vs-failed and joint failures; ties; metadata/split integrity; six-slot Holm; undefined studentization" };
    }
}
