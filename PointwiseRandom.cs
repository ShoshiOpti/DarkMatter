using System.Globalization;
using System.Security.Cryptography;
using System.Text.Json;

namespace DarkUniverse;

// Frozen protocol adapter over the shared PCG64 core; NumPy states and buffer semantics are unchanged.
internal sealed class PointwiseRandom
{
    readonly Pcg64Random random;

    public PointwiseRandom(ulong seed)
    {
        var initial = seed switch
        {
            2026092501 => ("213948995852192707599246433140438501816", "321904506844260870564164000364087464915"),
            2026092502 => ("9976519398804675029817934456828511402", "112449609693538579548892286257500921437"),
            2026092503 => ("171364352295095532094269204547027915199", "231027942869066986000708987154309789453"),
            _ => throw new ArgumentOutOfRangeException(nameof(seed), "Use one of the three frozen PROTOCOL.json seeds.")
        };
        random = new Pcg64Random(UInt128.Parse(initial.Item1, CultureInfo.InvariantCulture), UInt128.Parse(initial.Item2, CultureInfo.InvariantCulture));
    }
    public ulong Next64() => random.Next64();
    public int NextInt(int upper) => random.NextInt(upper);
    public void StartInt8Batch() => random.StartInt8Batch();
    public int NextInt8(int upper) => random.NextInt8(upper);

    public static Dictionary<string, object?> SelfCheck(string referencePath)
    {
        using var reference = JsonDocument.Parse(File.ReadAllText(referencePath));
        long regenerated = 0;
        foreach (var stream in reference.RootElement.GetProperty("streams").EnumerateArray())
        {
            ulong seed = stream.GetProperty("seed").GetUInt64();
            var raw = new PointwiseRandom(seed);
            foreach (var value in stream.GetProperty("raw_uint64").EnumerateArray())
                if (raw.Next64() != ulong.Parse(value.GetString()!, CultureInfo.InvariantCulture))
                    throw new InvalidDataException($"PCG64 raw sequence differs for seed {seed}.");
            var rng = new PointwiseRandom(seed);
            int upper = stream.GetProperty("upper").GetInt32(), galaxies = stream.GetProperty("galaxies").GetInt32();
            int draws = stream.GetProperty("draws").GetInt32(), batch = stream.GetProperty("batch_draws").GetInt32();
            bool int8 = stream.GetProperty("dtype").GetString() == "int8";
            using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
            var values = new byte[batch * galaxies];
            for (int start = 0; start < draws; start += batch)
            {
                int length = Math.Min(batch, draws - start) * galaxies;
                rng.StartInt8Batch();
                for (int i = 0; i < length; i++)
                    values[i] = (byte)(int8 ? rng.NextInt8(upper) : rng.NextInt(upper));
                hash.AppendData(values.AsSpan(0, length));
                regenerated += length;
            }
            string actual = Convert.ToHexStringLower(hash.GetHashAndReset());
            if (actual != stream.GetProperty("full_stream_sha256_uint8").GetString())
                throw new InvalidDataException($"PCG64 full bootstrap stream differs for seed {seed}.");
        }
        foreach (var check in reference.RootElement.GetProperty("rejection_and_buffer_checks").EnumerateArray())
        {
            var rng = new PointwiseRandom(check.GetProperty("seed").GetUInt64());
            int upper = check.GetProperty("upper").GetInt32();
            bool int8 = check.GetProperty("dtype").GetString() == "int8";
            foreach (var batch in check.GetProperty("batches").EnumerateArray())
            {
                rng.StartInt8Batch();
                foreach (var value in batch.EnumerateArray())
                    if ((int8 ? rng.NextInt8(upper) : rng.NextInt(upper)) != value.GetInt32())
                        throw new InvalidDataException("PCG64 rejection or byte buffer check differs from NumPy.");
            }
        }
        return new()
        {
            ["status"] = "pass",
            ["streams"] = 3,
            ["regenerated_values"] = regenerated,
            ["numpy_version"] = reference.RootElement.GetProperty("numpy_version").GetString(),
            ["method"] = "Native PCG64; full draw-stream SHA256 plus raw, rejection and batch-buffer tests"
        };
    }
}

