using System.Globalization;

namespace DarkUniverse;

public interface IFullSparcProfile
{
    double Mass(double radius);
    double Depth { get; }
}

/// <summary>Positive cubic density interpolation with analytic volume integrals and an infinite tail.</summary>
public sealed class FullSparcDensity : IFullSparcProfile
{
    readonly double[] x, prefix;
    readonly double[][] mass;
    readonly double total, tail, decay;
    const double Join = 30;
    public double Depth { get; }

    public FullSparcDensity(double[] radii, double[] radialDensity, double k)
    {
        if (radii.Length < 3 || radii.Length != radialDensity.Length || k <= 0 ||
            radii.Where((r, i) => !double.IsFinite(r) || r <= 0 || (i > 0 && r <= radii[i - 1])).Any() ||
            radialDensity.Any(v => !double.IsFinite(v) || v < 0) || radii[^1] < Join)
            throw new ArgumentException("Require ordered positive radii, positive density, and a finite tail rate.");
        decay = k;
        var rho = radii.Select((r, i) => radialDensity[i] / (r * r)).ToArray();
        int count = radii.Count(r => r < Join);
        double central = Math.Max(0, (rho[0] * radii[1] * radii[1] - rho[1] * radii[0] * radii[0]) /
            (radii[1] * radii[1] - radii[0] * radii[0]));
        x = [0, .. radii.Take(count), Join];
        var full = Cubics(radii, rho);
        int segment = Math.Clamp(Array.BinarySearch(radii, Join) is var idx && idx >= 0 ? idx : ~idx - 1, 0, radii.Length - 2);
        double last = Math.Max(0, Evaluate(full[segment], Join - radii[segment]));
        double[] y = [central, .. rho.Take(count), last];
        var pieces = Cubics(x, y);
        mass = new double[pieces.Length][]; prefix = new double[x.Length];
        double depth = 0;
        for (int i = 0; i < pieces.Length; i++)
        {
            mass[i] = Integral(Product(pieces[i], [x[i] * x[i], 2 * x[i], 1]));
            prefix[i + 1] = prefix[i] + Evaluate(mass[i], x[i + 1] - x[i]);
            depth += Evaluate(Integral(Product(pieces[i], [x[i], 1])), x[i + 1] - x[i]);
        }
        tail = last * Join * Join / (2 * k);
        total = prefix[^1] + tail;
        // Integral exp(-u)/(2*k*Join+u) du = exp(a) E1(a); no overflowing exp(a).
        double a = 2 * k * Join, sum = 1 / a;
        const int intervals = 1000; const double h = .1;
        for (int j = 1; j < intervals; j++) sum += (j % 2 == 0 ? 2 : 4) * Math.Exp(-j * h) / (a + j * h);
        sum += Math.Exp(-intervals * h) / (a + intervals * h);
        Depth = (depth + last * Join * Join * sum * h / 3) / total;
        if (!(total > 0 && double.IsFinite(Depth) && Depth > 0)) throw new InvalidDataException("Invalid density normalization.");
    }

    public double Mass(double radius)
    {
        if (double.IsNaN(radius) || radius < 0) throw new ArgumentOutOfRangeException(nameof(radius));
        if (radius >= Join) return (prefix[^1] + tail * -double.ExpM1(-2 * decay * (radius - Join))) / total;
        int i = Array.BinarySearch(x, radius); i = Math.Clamp(i >= 0 ? i : ~i - 1, 0, mass.Length - 1);
        return (prefix[i] + Evaluate(mass[i], radius - x[i])) / total;
    }

    static double[][] Cubics(double[] x, double[] y)
    {
        int n = x.Length; var h = new double[n - 1]; var slope = new double[n - 1]; var d = new double[n];
        for (int i = 0; i < n - 1; i++) { h[i] = x[i + 1] - x[i]; slope[i] = (y[i + 1] - y[i]) / h[i]; }
        for (int i = 1; i < n - 1; i++)
            if (slope[i - 1] * slope[i] > 0)
            {
                double w1 = 2 * h[i] + h[i - 1], w2 = h[i] + 2 * h[i - 1];
                d[i] = (w1 + w2) / (w1 / slope[i - 1] + w2 / slope[i]);
            }
        static double End(double h0, double h1, double m0, double m1)
        {
            double value = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
            if (Math.Sign(value) != Math.Sign(m0)) return 0;
            return Math.Sign(m0) != Math.Sign(m1) && Math.Abs(value) > 3 * Math.Abs(m0) ? 3 * m0 : value;
        }
        d[0] = End(h[0], h[1], slope[0], slope[1]);
        d[^1] = End(h[^1], h[^2], slope[^1], slope[^2]);
        return Enumerable.Range(0, n - 1).Select(i => new[] { y[i], d[i],
            (3 * slope[i] - 2 * d[i] - d[i + 1]) / h[i], (d[i] + d[i + 1] - 2 * slope[i]) / (h[i] * h[i]) }).ToArray();
    }
    static double[] Product(double[] a, double[] b)
    {
        var c = new double[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++) for (int j = 0; j < b.Length; j++) c[i + j] += a[i] * b[j];
        return c;
    }
    static double[] Integral(double[] a) => [0, .. a.Select((v, i) => v / (i + 1))];
    static double Evaluate(double[] a, double x)
    {
        double value = 0; for (int j = a.Length - 1; j >= 0; j--) value = value * x + a[j]; return value;
    }
}

public sealed class FullSparcEnvelope(IFullSparcProfile core, double dilation, double fraction, double a, double b, bool exact = false) : IFullSparcProfile
{
    readonly double scale = (a + b) / 2, e = (b - a) / (b + a);
    static readonly double[] Nodes = [-.9602898564975363, -.7966664774136267, -.525532409916329, -.1834346424956498,
        .1834346424956498, .525532409916329, .7966664774136267, .9602898564975363];
    static readonly double[] Weights = [.1012285362903763, .2223810344533745, .3137066458778873, .362683783378362,
        .362683783378362, .3137066458778873, .2223810344533745, .1012285362903763];
    public double Depth => (1 - fraction) * core.Depth / dilation + fraction * 2 / (Math.PI * scale) *
        (exact && e != 0 ? Math.Atanh(e) / e : 1 + e * e / 3);
    public static double F0(double z)
    {
        if (z < 1e-3)
        {
            double sum = 0; for (int j = 1; j <= 6; j++) sum += (j % 2 == 1 ? 1 : -1) * 2.0 * j / (2 * j + 1) * Math.Pow(z, 2 * j + 1);
            return 2 / Math.PI * sum;
        }
        return 2 / Math.PI * (Math.Atan(z) - z / (1 + z * z));
    }
    public static double F2(double z, double q) => F0(z) + 8 * q / (3 * Math.PI) * Math.Pow(z / (1 + z * z), 3);
    public double Mass(double radius)
    {
        double z = radius / scale;
        double env = exact ? Nodes.Select((u, i) => Weights[i] * F0(z / (1 + e * u))).Sum() / 2 : F2(z, e * e);
        return (1 - fraction) * core.Mass(radius / dilation) + fraction * env;
    }
}

public sealed record FullSparcGalaxy(string Name, double[] R, double[] Y, double[] Error, double[] BaryonSquared, bool Quality, string Group)
{
    public int InnerCount => R.Length - Math.Max(2, (int)Math.Ceiling(.2 * R.Length));
}

public sealed record FullSparcFit(double Length, double TotalMass, double Amplitude, double Chi2, double[] Predicted,
    bool LengthBound, bool PotentialBound, double ExteriorMassFraction, double CentralPotential);

public static class FullSparcFitter
{
    public const double G = 4.300917270e-6, C2 = 299792.458 * 299792.458;
    public static FullSparcFit Fit(FullSparcGalaxy d, IFullSparcProfile p, bool inner, int gridCount = 65)
    {
        int n = inner ? d.InnerCount : d.R.Length;
        if (n < 1 || gridCount < 3) throw new ArgumentException("Invalid fit mask/grid.");
        double lo = Math.Log(d.R[0] / 100), hi = Math.Log(d.R[^1] * 100), upper = .001 * C2 / p.Depth;
        (double cost, double amp, double length, double[] H) Objective(double logL)
        {
            double L = Math.Exp(logL); var H = d.R.Select(r => p.Mass(r / L) / (r / L)).ToArray();
            double lower = Math.Max(0, Enumerable.Range(0, H.Length).Max(i => -d.BaryonSquared[i] / H[i]));
            lower = lower * (1 + 1e-12) + (lower > 0 ? 1e-10 : 0);
            if (lower >= upper) return (double.PositiveInfinity, lower, L, H);
            double Derivative(double a)
            {
                double value = 0;
                for (int i = 0; i < n; i++) value += H[i] / (d.Error[i] * d.Error[i]) *
                    (1 - d.Y[i] / Math.Sqrt(Math.Max(d.BaryonSquared[i] + a * H[i], 1e-30)));
                return value;
            }
            double amplitude;
            if (Derivative(lower) >= 0) amplitude = lower;
            else if (Derivative(upper) <= 0) amplitude = upper;
            else
            {
                double l = lower, u = upper;
                for (int iter = 0; iter < 110; iter++)
                {
                    double mid = (l + u) / 2;
                    if (Derivative(mid) > 0) u = mid; else l = mid;
                    if (u - l < 1e-11 + 1e-14 * Math.Abs(mid)) break;
                }
                amplitude = (l + u) / 2;
            }
            double cost = 0;
            for (int i = 0; i < n; i++) cost += Math.Pow((Math.Sqrt(Math.Max(d.BaryonSquared[i] + amplitude * H[i], 0)) - d.Y[i]) / d.Error[i], 2);
            return (cost, amplitude, L, H);
        }
        var grid = Data.Linspace(lo, hi, gridCount); var costs = grid.Select(s => Objective(s).cost).ToArray();
        int best = costs[0] <= costs[^1] ? 0 : grid.Length - 1; double bestS = grid[best], bestCost = costs[best];
        for (int i = 1; i < grid.Length - 1; i++)
        {
            if (costs[i] > costs[i - 1] || costs[i] > costs[i + 1]) continue;
            double l = grid[i - 1], u = grid[i + 1], ratio = (Math.Sqrt(5) - 1) / 2;
            double s1 = u - ratio * (u - l), s2 = l + ratio * (u - l), c1 = Objective(s1).cost, c2 = Objective(s2).cost;
            for (int iter = 0; iter < 90 && u - l > 5e-10; iter++)
            {
                if (c1 <= c2) { u = s2; s2 = s1; c2 = c1; s1 = u - ratio * (u - l); c1 = Objective(s1).cost; }
                else { l = s1; s1 = s2; c1 = c2; s2 = l + ratio * (u - l); c2 = Objective(s2).cost; }
            }
            double s = c1 <= c2 ? s1 : s2, cost = Math.Min(c1, c2);
            if (cost < bestCost) { bestCost = cost; bestS = s; }
        }
        var solution = Objective(bestS);
        if (!double.IsFinite(solution.cost)) throw new InvalidDataException("No feasible fit: " + d.Name);
        var pred = solution.H.Select((h, i) => Math.Sqrt(Math.Max(d.BaryonSquared[i] + solution.amp * h, 0))).ToArray();
        return new(solution.length, solution.amp * solution.length / G, solution.amp, solution.cost, pred,
            Math.Abs(bestS - lo) < 1e-6 || Math.Abs(bestS - hi) < 1e-6, solution.amp >= upper * (1 - 1e-8),
            1 - p.Mass(d.R[^1] / solution.length), solution.amp * p.Depth / C2);
    }
}
