using System.Numerics;
using MathNet.Numerics.IntegralTransforms;

namespace DarkUniverse;

/// <summary>
/// Independent instantaneous checks of the retained original two-field GPP action.
/// Dimensionless m = hbar = 1; no galaxy fit, Poisson solve or time integration.
/// </summary>
public static class FieldReconstructionChecks
{
    public const double EquationTolerance = 2e-10;
    // Differentiating the complex RHS twice more amplifies double-precision
    // Fourier roundoff. This separately declared bound is not a fit tolerance.
    public const double AccelerationTolerance = 2e-9;
    const double L1 = 1, L2 = .64, Lx = .48, Lp = 1.12;
    const string Scope = "Dimensionless periodic 1D original two-field retained quartic GPP action; instantaneous identities only; no halo solution, time integration or observational inference.";

    public sealed record GridCheck(int GridPoints, bool Passed,
        IReadOnlyDictionary<string, double> MaximumAbsoluteErrors);
    public sealed record HiddenStressCheck(int GridPoints, double PhaseAmplitude, bool Passed,
        IReadOnlyDictionary<string, double> MaximumAbsoluteErrors,
        double MeanInternalGradientEnergy, double MaximumDensityAccelerationChange);
    public sealed record Report(bool Passed, string Scope, string Units, double[] QuarticCoefficients,
        double EquationTolerance, double AccelerationTolerance,
        GridCheck[] Grids, HiddenStressCheck[] HiddenStress, string[] Artifacts);

    sealed record HiddenFixture(double[] X, double[] Density, double[] Current,
        double[] Acceleration, double[] DeltaAcceleration, double[] InternalEnergy,
        double[] DeltaFlux, Dictionary<string, double> Errors);

    /// <summary>Checks one full-band Fourier grid without writing files.</summary>
    public static GridCheck CheckGrid(int gridPoints) => EvaluateGrid(gridPoints, null);

    /// <summary>Checks the three fixed-density phase fixtures without writing files.</summary>
    public static HiddenStressCheck[] CheckHiddenStress(int gridPoints) =>
        EvaluateHidden(gridPoints, null).Checks;

    /// <summary>Writes report, coordinate CSVs and SVGs directly in outputFolder.</summary>
    public static Report Run(string outputFolder)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(outputFolder);
        var coordinates = new List<Dictionary<string, object?>>();
        var hiddenCoordinates = new List<Dictionary<string, object?>>();
        GridCheck[] grids = [EvaluateGrid(64, coordinates), EvaluateGrid(128, coordinates)];
        var hidden64 = EvaluateHidden(64, hiddenCoordinates);
        var hidden128 = EvaluateHidden(128, hiddenCoordinates);
        HiddenStressCheck[] hidden = [.. hidden64.Checks, .. hidden128.Checks];
        Directory.CreateDirectory(outputFolder);
        Csv.Write(Path.Combine(outputFolder, "reconstruction_coordinates.csv"), coordinates);
        Csv.Write(Path.Combine(outputFolder, "hidden_stress_coordinates.csv"), hiddenCoordinates);
        SavePlots(outputFolder, grids, hidden128.Checks, hidden128.Fixtures);
        string[] artifacts = ["report.json", "reconstruction_coordinates.csv", "hidden_stress_coordinates.csv",
            "reconstruction_checks.svg", "reconstruction_checks.series.csv",
            "hidden_stress_response.svg", "hidden_stress_response.series.csv"];
        var report = new Report(true, Scope, "m = hbar = 1; x in [0, 2*pi); full-band unfiltered Fourier derivatives",
            [L1, L2, Lx, Lp], EquationTolerance, AccelerationTolerance, grids, hidden, artifacts);
        Data.SaveJson(Path.Combine(outputFolder, "report.json"), report);
        return report;
    }

    static GridCheck EvaluateGrid(int count, List<Dictionary<string, object?>>? rows)
    {
        ValidateGrid(count);
        var x = Grid(count);
        var n = Map(x, z => 1 + .1 * Math.Cos(z));
        var f = Map(x, z => .35 + .07 * Math.Sin(2 * z));
        var theta = Map(x, z => .2 * Math.Sin(z));
        var delta = Map(x, z => 1.1 + .15 * Math.Cos(3 * z));
        var phi = Map(x, z => .1 * Math.Cos(2 * z));
        var p1 = new Complex[count]; var p2 = new Complex[count];
        for (int j = 0; j < count; j++)
        {
            p1[j] = Complex.FromPolarCoordinates(Math.Sqrt(n[j] * (1 - f[j])), theta[j]);
            p2[j] = Complex.FromPolarCoordinates(Math.Sqrt(n[j] * f[j]), theta[j] + delta[j]);
        }
        var (t1, t2) = ParentRhs(p1, p2, phi);
        var p1x = Derivative(p1); var p2x = Derivative(p2);
        var fx = Derivative(f); var tx = Derivative(theta); var dx = Derivative(delta);
        var sqrtN = Map(n, Math.Sqrt); var sqrtNx = Derivative(sqrtN); var sqrtNxx = Derivative(sqrtN, 2);
        var vel = new double[count]; var ff = new double[count]; var fp = new double[count];
        var contact = new double[count]; var un = new double[count]; var uf = new double[count]; var ud = new double[count];
        var nDot = new double[count]; var fDot = new double[count]; var deltaDot = new double[count]; var thetaDot = new double[count];
        var current = new double[count]; var energy = new double[count]; var berry = new double[count]; var polarBerry = new double[count];
        var parentContact = new double[count];
        for (int j = 0; j < count; j++)
        {
            double nj = n[j], fj = f[j];
            ff[j] = fj * (1 - fj); fp[j] = 1 - 2 * fj;
            vel[j] = tx[j] + fj * dx[j];
            double angular = 4 * Lx + 2 * Lp * Math.Cos(2 * delta[j]);
            double a = L1 * (1 - fj) * (1 - fj) + L2 * fj * fj + angular * ff[j];
            contact[j] = nj * nj * a / 16;
            un[j] = nj * a / 8;
            uf[j] = nj * nj * (-2 * L1 * (1 - fj) + 2 * L2 * fj + angular * fp[j]) / 16;
            ud[j] = -Lp * nj * nj * ff[j] * Math.Sin(2 * delta[j]) / 4;
            nDot[j] = 2 * (Complex.Conjugate(p1[j]) * t1[j] + Complex.Conjugate(p2[j]) * t2[j]).Real;
            fDot[j] = (2 * (Complex.Conjugate(p2[j]) * t2[j]).Real - fj * nDot[j]) / nj;
            deltaDot[j] = (t2[j] / p2[j] - t1[j] / p1[j]).Imaginary;
            thetaDot[j] = (t1[j] / p1[j]).Imaginary;
            current[j] = (Complex.Conjugate(p1[j]) * p1x[j] + Complex.Conjugate(p2[j]) * p2x[j]).Imaginary;
            energy[j] = (NormSquared(p1x[j]) + NormSquared(p2x[j])) / 2;
            parentContact[j] = Contact(p1[j], p2[j]);
            berry[j] = -(Complex.Conjugate(p1[j]) * t1[j] + Complex.Conjugate(p2[j]) * t2[j]).Imaginary;
            polarBerry[j] = -nj * (thetaDot[j] + fj * deltaDot[j]);
        }
        var nv = Index(count, j => n[j] * vel[j]);
        var relFluxX = Derivative(Index(count, j => n[j] * ff[j] * dx[j]));
        var compFluxX = Derivative(Index(count, j => n[j] * fx[j] / (4 * ff[j])));
        var nNew = Map(Derivative(nv), z => -z);
        var fNew = new double[count]; var dNew = new double[count]; var tNew = new double[count]; var energyNew = new double[count];
        for (int j = 0; j < count; j++)
        {
            fNew[j] = -vel[j] * fx[j] - relFluxX[j] / n[j] + ud[j] / n[j];
            dNew[j] = -vel[j] * dx[j] + compFluxX[j] / n[j] + fp[j] * fx[j] * fx[j] / (8 * ff[j] * ff[j])
                - fp[j] * dx[j] * dx[j] / 2 - uf[j] / n[j];
            double geometry = fx[j] * fx[j] / (4 * ff[j]) + ff[j] * dx[j] * dx[j];
            tNew[j] = -f[j] * dNew[j] - vel[j] * vel[j] / 2 - phi[j] - un[j] - geometry / 2 + sqrtNxx[j] / (2 * sqrtN[j]);
            energyNew[j] = (sqrtNx[j] * sqrtNx[j] + n[j] * vel[j] * vel[j] + n[j] * geometry) / 2;
        }
        Dictionary<string, double> errors = new()
        {
            ["continuity"] = Error(nDot, nNew), ["composition"] = Error(fDot, fNew),
            ["relative_phase"] = Error(deltaDot, dNew), ["common_phase"] = Error(thetaDot, tNew),
            ["gradient_energy"] = Error(energy, energyNew), ["contact_energy"] = Error(parentContact, contact),
            ["berry_term"] = Error(berry, polarBerry), ["common_current"] = Error(current, nv)
        };
        ValidateErrors(errors, EquationTolerance, "reconstruction N=" + count);
        if (rows is not null)
            for (int j = 0; j < count; j++)
                rows.Add(new()
                {
                    ["grid_points"] = count, ["x"] = x[j], ["density"] = n[j], ["composition"] = f[j],
                    ["theta"] = theta[j], ["delta"] = delta[j], ["external_potential"] = phi[j],
                    ["psi1_real"] = p1[j].Real, ["psi1_imag"] = p1[j].Imaginary,
                    ["psi2_real"] = p2[j].Real, ["psi2_imag"] = p2[j].Imaginary,
                    ["n_dot_parent"] = nDot[j], ["n_dot_polar"] = nNew[j],
                    ["f_dot_parent"] = fDot[j], ["f_dot_polar"] = fNew[j],
                    ["delta_dot_parent"] = deltaDot[j], ["delta_dot_polar"] = dNew[j],
                    ["theta_dot_parent"] = thetaDot[j], ["theta_dot_polar"] = tNew[j],
                    ["gradient_energy_parent"] = energy[j], ["gradient_energy_polar"] = energyNew[j],
                    ["contact_energy_parent"] = parentContact[j], ["contact_energy_polar"] = contact[j],
                    ["berry_parent"] = berry[j], ["berry_polar"] = polarBerry[j],
                    ["current_parent"] = current[j], ["current_polar"] = nv[j]
                });
        return new(count, true, errors);
    }

    static (HiddenStressCheck[] Checks, HiddenFixture[] Fixtures) EvaluateHidden(int count,
        List<Dictionary<string, object?>>? rows)
    {
        ValidateGrid(count);
        double[] amplitudes = [0, .25, .5];
        var fixtures = amplitudes.Select(a => Hidden(count, a)).ToArray();
        var checks = new HiddenStressCheck[amplitudes.Length];
        for (int aIndex = 0; aIndex < amplitudes.Length; aIndex++)
        {
            var item = fixtures[aIndex]; var baseline = fixtures[0]; double a = amplitudes[aIndex];
            var actualChange = Index(count, j => item.Acceleration[j] - baseline.Acceleration[j]);
            item.Errors["hidden_stress_acceleration"] = Error(actualChange, item.DeltaAcceleration);
            double expectedMeanEnergy = .35 * .65 * a * a / 4;
            item.Errors["mean_internal_energy"] = Math.Abs(item.InternalEnergy.Average() - expectedMeanEnergy);
            foreach (var (name, value) in item.Errors)
                ValidateError(value, name is "density_acceleration" or "hidden_stress_acceleration"
                    ? AccelerationTolerance : EquationTolerance, $"hidden N={count}, A={a}: {name}");
            double change = actualChange.Max(Math.Abs);
            if (a > 0 && !(change > .01)) throw new InvalidDataException("Hidden-stress fixture did not change instantaneous density acceleration.");
            checks[aIndex] = new(count, a, true, item.Errors, item.InternalEnergy.Average(), change);
            if (rows is not null)
                for (int j = 0; j < count; j++)
                    rows.Add(new()
                    {
                        ["grid_points"] = count, ["phase_amplitude"] = a, ["x"] = item.X[j],
                        ["density"] = item.Density[j], ["composition"] = .35,
                        ["theta"] = -.35 * a * Math.Sin(item.X[j]),
                        ["delta"] = Math.PI / 2 + a * Math.Sin(item.X[j]),
                        ["current"] = item.Current[j], ["internal_gradient_energy"] = item.InternalEnergy[j],
                        ["delta_momentum_flux"] = item.DeltaFlux[j], ["density_acceleration_parent"] = item.Acceleration[j],
                        ["density_acceleration_change_parent"] = actualChange[j],
                        ["density_acceleration_change_stress"] = item.DeltaAcceleration[j]
                    });
        }
        return (checks, fixtures);
    }

    static HiddenFixture Hidden(int count, double amplitude)
    {
        double f = .35, ff = f * (1 - f);
        var x = Grid(count); var n = Map(x, z => 1 + .1 * Math.Cos(z));
        var p1 = new Complex[count]; var p2 = new Complex[count];
        for (int j = 0; j < count; j++)
        {
            double angle = amplitude * Math.Sin(x[j]);
            p1[j] = Complex.FromPolarCoordinates(Math.Sqrt(n[j] * (1 - f)), -f * angle);
            p2[j] = Complex.FromPolarCoordinates(Math.Sqrt(n[j] * f), Math.PI / 2 + (1 - f) * angle);
        }
        var (t1, t2) = ParentRhs(p1, p2, new double[count]);
        var p1x = Derivative(p1); var p2x = Derivative(p2);
        var t1x = Derivative(t1); var t2x = Derivative(t2); var nxx = Derivative(n, 2);
        var sqrtNx = Derivative(Map(n, Math.Sqrt));
        var current = new double[count]; var jt = new double[count]; var flux = new double[count];
        var deltaFlux = new double[count]; var internalEnergy = new double[count]; var reconstructedDensity = new double[count];
        var complexInternalEnergy = new double[count]; var densityVelocity = new double[count];
        for (int j = 0; j < count; j++)
        {
            reconstructedDensity[j] = NormSquared(p1[j]) + NormSquared(p2[j]);
            densityVelocity[j] = 2 * (Complex.Conjugate(p1[j]) * t1[j] + Complex.Conjugate(p2[j]) * t2[j]).Real;
            complexInternalEnergy[j] = (NormSquared(p1x[j]) + NormSquared(p2x[j]) - sqrtNx[j] * sqrtNx[j]) / 2;
            current[j] = (Complex.Conjugate(p1[j]) * p1x[j] + Complex.Conjugate(p2[j]) * p2x[j]).Imaginary;
            jt[j] = (Complex.Conjugate(t1[j]) * p1x[j] + Complex.Conjugate(p1[j]) * t1x[j]
                + Complex.Conjugate(t2[j]) * p2x[j] + Complex.Conjugate(p2[j]) * t2x[j]).Imaginary;
            double u = Contact(p1[j], p2[j]);
            flux[j] = NormSquared(p1x[j]) + NormSquared(p2x[j]) - nxx[j] / 4 + u;
            double u0 = n[j] * n[j] * (L1 * (1 - f) * (1 - f) + L2 * f * f + (4 * Lx - 2 * Lp) * ff) / 16;
            double relativeGradient = n[j] * ff * amplitude * amplitude * Math.Pow(Math.Cos(x[j]), 2);
            internalEnergy[j] = relativeGradient / 2;
            deltaFlux[j] = relativeGradient + u - u0;
        }
        var acceleration = Map(Derivative(jt), z => -z);
        Dictionary<string, double> errors = new()
        {
            ["density"] = Error(reconstructedDensity, n), ["common_current"] = Error(current, new double[count]),
            ["density_velocity"] = Error(densityVelocity, new double[count]),
            ["internal_gradient_energy"] = Error(complexInternalEnergy, internalEnergy),
            ["momentum"] = Error(jt, Map(Derivative(flux), z => -z)),
            ["density_acceleration"] = Error(acceleration, Derivative(flux, 2))
        };
        return new(x, n, current, acceleration, Derivative(deltaFlux, 2), internalEnergy, deltaFlux, errors);
    }

    static (Complex[] First, Complex[] Second) ParentRhs(Complex[] p1, Complex[] p2, double[] phi)
    {
        var p1xx = Derivative(p1, 2); var p2xx = Derivative(p2, 2);
        var t1 = new Complex[p1.Length]; var t2 = new Complex[p2.Length];
        for (int j = 0; j < p1.Length; j++)
        {
            double n1 = NormSquared(p1[j]), n2 = NormSquared(p2[j]);
            var u1 = ((L1 * n1 + 2 * Lx * n2) * p1[j] + Lp * Complex.Conjugate(p1[j]) * p2[j] * p2[j]) / 8;
            var u2 = ((L2 * n2 + 2 * Lx * n1) * p2[j] + Lp * Complex.Conjugate(p2[j]) * p1[j] * p1[j]) / 8;
            t1[j] = -Complex.ImaginaryOne * (-p1xx[j] / 2 + phi[j] * p1[j] + u1);
            t2[j] = -Complex.ImaginaryOne * (-p2xx[j] / 2 + phi[j] * p2[j] + u2);
        }
        return (t1, t2);
    }

    static double Contact(Complex a, Complex b)
    {
        double n1 = NormSquared(a), n2 = NormSquared(b);
        return (L1 * n1 * n1 + L2 * n2 * n2 + 4 * Lx * n1 * n2
            + 2 * Lp * (Complex.Conjugate(a) * Complex.Conjugate(a) * b * b).Real) / 16;
    }

    static Complex[] Derivative(Complex[] values, int order = 1)
    {
        var modes = (Complex[])values.Clone();
        Fourier.Forward(modes, FourierOptions.Matlab);
        for (int j = 0; j < modes.Length; j++)
        {
            int k = j < modes.Length / 2 ? j : j - modes.Length;
            modes[j] *= order == 1 ? new Complex(0, k) : new Complex(-k * k, 0);
        }
        Fourier.Inverse(modes, FourierOptions.Matlab);
        return modes;
    }
    static double[] Derivative(double[] values, int order = 1) =>
        Derivative(values.Select(v => new Complex(v, 0)).ToArray(), order).Select(v => v.Real).ToArray();
    static double[] Grid(int count) => Index(count, j => 2 * Math.PI * j / count);
    static double[] Map(double[] input, Func<double, double> transform) => input.Select(transform).ToArray();
    static double[] Index(int count, Func<int, double> transform) => Enumerable.Range(0, count).Select(transform).ToArray();
    static double NormSquared(Complex value) => value.Real * value.Real + value.Imaginary * value.Imaginary;
    static double Error(double[] first, double[] second)
    {
        double maximum = 0;
        for (int j = 0; j < first.Length; j++)
        {
            if (!double.IsFinite(first[j]) || !double.IsFinite(second[j]))
                throw new InvalidDataException("Nonfinite reconstruction value.");
            maximum = Math.Max(maximum, Math.Abs(first[j] - second[j]));
        }
        return maximum;
    }
    static void ValidateGrid(int count)
    {
        if (count < 16 || (count & (count - 1)) != 0)
            throw new ArgumentOutOfRangeException(nameof(count), "Use a power-of-two periodic grid with at least 16 points.");
    }
    static void ValidateErrors(IReadOnlyDictionary<string, double> errors, double tolerance, string scope)
    {
        foreach (var (name, error) in errors) ValidateError(error, tolerance, scope + ": " + name);
    }
    static void ValidateError(double error, double tolerance, string label)
    {
        if (!double.IsFinite(error) || error > tolerance)
            throw new InvalidDataException($"Field reconstruction failed: {label}: {error:R} exceeds {tolerance:R}.");
    }

    static void SavePlots(string folder, GridCheck[] grids, HiddenStressCheck[] hidden, HiddenFixture[] fixtures)
    {
        var checksPlot = new PlotDocument("Original two-field reconstruction: independent instantaneous checks", 2, 1)
        { Scope = "Synthetic dimensionless fixtures; no halo fit, Poisson solve or time evolution." };
        var errorsPanel = new PlotPanel("Parent complex equations vs polar equations", "Checked identity", "Maximum absolute difference")
        { YLog = true, XMin = -.3, XMax = 7.3, YMin = 1e-17, YMax = 1e-9,
            XTicks = new() { [0] = "n", [1] = "f", [2] = "delta", [3] = "theta", [4] = "grad", [5] = "U", [6] = "Berry", [7] = "j" },
            Note = "Display floor 1e-17 only; report retains raw errors. Bound 2e-10." };
        string[] keys = ["continuity", "composition", "relative_phase", "common_phase", "gradient_energy", "contact_energy", "berry_term", "common_current"];
        string[] colors = ["#24658c", "#bd5e39", "#38916b"];
        for (int i = 0; i < grids.Length; i++)
            errorsPanel.Add($"N = {grids[i].GridPoints}", Index(keys.Length, j => j),
                keys.Select(k => Math.Max(1e-17, grids[i].MaximumAbsoluteErrors[k])).ToArray(), colors[i], true);
        errorsPanel.Add("Declared tolerance", [0, 7], [EquationTolerance, EquationTolerance], "#888888").Dash = "5,4";
        var energyPanel = new PlotPanel("Same density and zero common current", "Relative-phase amplitude A", "Mean internal gradient energy")
        { Note = "n = 1 + 0.1 cos(x), f = 0.35; energy = f(1-f) A^2 /4." };
        energyPanel.Add("Exact phase-family energy", hidden.Select(h => h.PhaseAmplitude).ToArray(),
            hidden.Select(h => h.MeanInternalGradientEnergy).ToArray(), "#38916b", true);
        checksPlot.Panels.Add(errorsPanel); checksPlot.Panels.Add(energyPanel);
        checksPlot.Save(Path.Combine(folder, "reconstruction_checks.svg"));

        var response = new PlotDocument("Fixed density does not determine its later response", 2, 1)
        { Scope = "Instantaneous second derivative only; no integration, galaxy model or observational claim." };
        var densityPanel = new PlotPanel("Identical leading density / gravity source", "Periodic coordinate x", "n(x)")
        { Note = "All A have this density and j = 0 initially; same gravity for same boundaries." };
        densityPanel.Add("A = 0, 0.25, 0.5 (coincident)", fixtures[0].X, fixtures[0].Density);
        var accelerationPanel = new PlotPanel("Different instantaneous density acceleration", "Periodic coordinate x", "n_tt(A) - n_tt(0)")
        { Note = "Phi =0; parent complex RHS. Stress identity checked separately." };
        for (int i = 0; i < hidden.Length; i++)
            accelerationPanel.Add($"A = {hidden[i].PhaseAmplitude:G}", fixtures[i].X,
                Index(fixtures[i].X.Length, j => fixtures[i].Acceleration[j] - fixtures[0].Acceleration[j]), colors[i]);
        response.Panels.Add(densityPanel); response.Panels.Add(accelerationPanel);
        response.Save(Path.Combine(folder, "hidden_stress_response.svg"));
    }
}
