using System.Globalization;
using System.Numerics;
using System.Security.Cryptography;
using System.Text.Json;

namespace DarkUniverse;

internal sealed class PointwiseRandom
{
    static readonly UInt128 Multiplier = ((UInt128)2549297995355413924UL << 64) | 4865540595714422341UL;
    UInt128 state;
    readonly UInt128 increment;
    uint cached32, byteBuffer;
    bool hasCached32;
    int bytesRemaining;

    public PointwiseRandom(ulong seed)
    {
        // Frozen NumPy SeedSequence states for the three protocol seeds.
        var initial = seed switch
        {
            2026092501 => ("213948995852192707599246433140438501816", "321904506844260870564164000364087464915"),
            2026092502 => ("9976519398804675029817934456828511402", "112449609693538579548892286257500921437"),
            2026092503 => ("171364352295095532094269204547027915199", "231027942869066986000708987154309789453"),
            _ => throw new ArgumentOutOfRangeException(nameof(seed), "Use one of the three frozen PROTOCOL.json seeds.")
        };
        state = UInt128.Parse(initial.Item1, CultureInfo.InvariantCulture);
        increment = UInt128.Parse(initial.Item2, CultureInfo.InvariantCulture);
    }

    public ulong Next64()
    {
        state = unchecked(state * Multiplier + increment);
        return BitOperations.RotateRight((ulong)(state >> 64) ^ (ulong)state, (int)(state >> 122));
    }

    uint Next32()
    {
        if (hasCached32)
        {
            hasCached32 = false;
            return cached32;
        }
        ulong value = Next64();
        cached32 = (uint)(value >> 32);
        hasCached32 = true;
        return (uint)value;
    }

    public int NextInt(int upper)
    {
        if (upper <= 0)
            throw new ArgumentOutOfRangeException(nameof(upper));
        if (upper == 1)
            return 0;
        // Lemire bounded integers, with rejection to avoid modulo bias.
        uint bound = (uint)upper, threshold = unchecked(0u - bound) % bound;
        ulong product;
        do
        {
            product = (ulong)Next32() * bound;
        } while ((uint)product < threshold);
        return (int)(product >> 32);
    }

    // NumPy discards unused bytes between separate integers(dtype=int8) calls.
    public void StartInt8Batch() => bytesRemaining = 0;

    int NextByte()
    {
        if (bytesRemaining == 0)
        {
            byteBuffer = Next32();
            bytesRemaining = 4;
        }
        int value = (int)(byteBuffer & 255);
        byteBuffer >>= 8;
        bytesRemaining--;
        return value;
    }

    public int NextInt8(int upper)
    {
        if (upper <= 0 || upper > 128)
            throw new ArgumentOutOfRangeException(nameof(upper));
        if (upper == 1)
            return 0;
        int threshold = 256 % upper, product;
        do
        {
            product = NextByte() * upper;
        } while ((product & 255) < threshold);
        return product >> 8;
    }

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
