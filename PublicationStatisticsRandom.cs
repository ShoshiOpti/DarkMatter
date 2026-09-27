using System.Globalization;
using System.Security.Cryptography;
using System.Text.Json;

namespace DarkUniverse;

// Frozen NumPy SeedSequence adapter; complete draw hashes verify the shared PCG64 core.
internal sealed class PublicationStatisticsRandom
{
    readonly Pcg64Random random;

    internal PublicationStatisticsRandom(JsonElement reference, ulong seed)
    {
        var stream = reference.GetProperty("streams").EnumerateArray().Single(s => s.GetProperty("seed").GetUInt64() == seed);
        random = new Pcg64Random(UInt128.Parse(stream.GetProperty("state").GetString()!, CultureInfo.InvariantCulture),
            UInt128.Parse(stream.GetProperty("increment").GetString()!, CultureInfo.InvariantCulture));
    }
    internal ulong Next64() => random.Next64();
    internal int NextInt(int upper) => random.NextInt(upper);
    internal void StartInt8Batch() => random.StartInt8Batch();
    internal int NextSignBit() => random.NextInt8(2);

    internal static Dictionary<string, object?> SelfCheck(JsonElement reference)
    {
        long regenerated = 0;
        foreach (var stream in reference.GetProperty("streams").EnumerateArray())
        {
            ulong seed = stream.GetProperty("seed").GetUInt64();
            var raw = new PublicationStatisticsRandom(reference, seed);
            foreach (var value in stream.GetProperty("raw_uint64").EnumerateArray())
                if (raw.Next64() != ulong.Parse(value.GetString()!, CultureInfo.InvariantCulture))
                    throw new InvalidDataException($"Publication PCG64 raw stream differs for seed {seed}.");
            var rng = new PublicationStatisticsRandom(reference, seed);
            int batch = stream.GetProperty("batch_draws").GetInt32(), galaxies = stream.GetProperty("galaxies").GetInt32();
            int draws = stream.GetProperty("draws").GetInt32(), upper = stream.GetProperty("upper").GetInt32();
            bool int8 = stream.GetProperty("dtype").GetString() == "int8";
            using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
            byte[] buffer = new byte[batch * galaxies];
            for (int start = 0; start < draws; start += batch)
            {
                int length = Math.Min(batch, draws - start) * galaxies;
                rng.StartInt8Batch();
                for (int i = 0; i < length; i++) buffer[i] = (byte)(int8 ? rng.NextSignBit() : rng.NextInt(upper));
                hash.AppendData(buffer.AsSpan(0, length));
                regenerated += length;
            }
            if (Convert.ToHexStringLower(hash.GetHashAndReset()) != stream.GetProperty("full_stream_sha256_uint8").GetString())
                throw new InvalidDataException($"Publication PCG64 full bootstrap stream differs for seed {seed}.");
        }
        return new() { ["status"] = "pass", ["streams"] = reference.GetProperty("streams").GetArrayLength(),
            ["regenerated_values"] = regenerated, ["numpy_version"] = reference.GetProperty("numpy_version").GetString(),
            ["method"] = "Native PCG64 raw uint64 values and complete bootstrap draw stream SHA-256 checks" };
    }
}

