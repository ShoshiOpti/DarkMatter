using System.Globalization;
using System.Security;
using System.Text;

namespace DarkUniverse;

public sealed class PlotSeries(string label, double[] x, double[] y, string color, bool points)
{
    public string Label = label, Color = color, Dash = "";
    public double[] X = x, Y = y;
    public double[]? Lower, Upper, XLower, XUpper;
    public bool Points = points, Hollow, Step, Bars, Fill;
}

public sealed class PlotPanel(string title, string xLabel, string yLabel)
{
    public string Title = title, XLabel = xLabel, YLabel = yLabel, Note = "";
    public bool XLog, YLog, XSymLog, YSymLog, ReverseY, TextOnly;
    public double SymLogThreshold = .1, XSymLogThreshold = 1, XSymLogScale = 1;
    public double? XMin, XMax, YMin, YMax;
    public List<PlotSeries> Series = [];
    public Dictionary<double, string>? XTicks, YTicks;
    public List<(double Start, double End, string Color)> Shades = [];

    public PlotSeries Add(string label, double[] x, double[] y, string color = "#24658c", bool points = false)
    {
        if (x.Length != y.Length)
            throw new InvalidDataException($"Series length mismatch: {Title}/{label}");
        var series = new PlotSeries(label, x, y, color, points);
        Series.Add(series);
        return series;
    }

    public PlotSeries ErrorBars(string label, double[] x, double[] y, double[] error, string color = "#222222")
    {
        var series = Add(label, x, y, color, true);
        series.Lower = y.Zip(error, (a, b) => a - b).ToArray();
        series.Upper = y.Zip(error, (a, b) => a + b).ToArray();
        return series;
    }

    public PlotSeries Band(string label, double[] x, double[] lower, double[] upper, string color = "#24658c")
    {
        var series = Add(label, x, lower.Zip(upper, (a, b) => (a + b) / 2).ToArray(), color);
        series.Lower = lower;
        series.Upper = upper;
        series.Fill = true;
        return series;
    }

    public void Ecdf(string label, IEnumerable<double> values, string color, bool origin = false)
    {
        var all = values.ToArray();
        var x = all.Where(double.IsFinite).Order().ToArray();
        if (x.Length == 0)
            return;

        // Failed predictions remain in the denominator as mass at infinity.
        var y = Enumerable.Range(1, x.Length).Select(i => (double)i / all.Length).ToArray();
        if (origin)
        {
            x = [x[0], .. x];
            y = [0, .. y];
        }
        Add(label, x, y, color).Step = true;
    }
}

public sealed class PlotDocument(string title, int columns, int rows)
{
    public List<PlotPanel> Panels = [];
    public string Scope = "Saved-result reconstruction; no refitting";

    public void Save(string path)
    {
        if (Panels.Count > columns * rows)
            throw new InvalidDataException("Too many panels");
        Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(path))!);
        int rowHeight = 390 + Panels.Select(p => LegendHeight(p) + (p.Note.Length > 0 ? 20 : 0)).DefaultIfEmpty(0).Max();
        int width = columns * 570;
        int height = rows * rowHeight + 110;
        var svg = new StringBuilder($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/><style>text{{font-family:Arial,sans-serif;fill:#222}}.tick{{font-size:11px}}</style>");
        var titleLines = Wrap(title, Math.Max(35, width / 10));
        for (int i = 0; i < titleLines.Length; i++)
            Text(svg, width / 2, 25 + i * 21, titleLines[i], 18, "middle");
        for (int i = 0; i < Panels.Count; i++)
            Draw(svg, Panels[i], i, (i % columns) * 570, (i / columns) * rowHeight + 65);
        Text(svg, width / 2, height - 12, Scope, 11, "middle");
        svg.Append("</svg>");
        File.WriteAllText(path, svg.ToString());
        Csv.Write(Path.ChangeExtension(path, "series.csv"), SeriesRows());
    }

    IEnumerable<Dictionary<string, object?>> SeriesRows()
    {
        for (int panelIndex = 0; panelIndex < Panels.Count; panelIndex++)
        {
            var panel = Panels[panelIndex];
            for (int seriesIndex = 0; seriesIndex < panel.Series.Count; seriesIndex++)
            {
                var series = panel.Series[seriesIndex];
                for (int i = 0; i < series.X.Length; i++)
                {
                    yield return new Dictionary<string, object?>
                    {
                        ["panel"] = panelIndex,
                        ["title"] = panel.Title,
                        ["series"] = seriesIndex,
                        ["label"] = series.Label,
                        ["index"] = i,
                        ["x"] = series.X[i],
                        ["y"] = series.Y[i],
                        ["y_lower"] = series.Lower?[i],
                        ["y_upper"] = series.Upper?[i],
                        ["x_lower"] = series.XLower?[i],
                        ["x_upper"] = series.XUpper?[i]
                    };
                }
            }
        }
    }

    static string[] Wrap(string text, int width)
    {
        var lines = new List<string>();
        string line = "";
        foreach (string word in text.Split(' ', StringSplitOptions.RemoveEmptyEntries))
        {
            if (line.Length > 0 && line.Length + word.Length + 1 > width)
            {
                lines.Add(line);
                line = "";
            }
            line += (line.Length > 0 ? " " : "") + word;
        }
        if (line.Length > 0)
            lines.Add(line);
        return lines.ToArray();
    }

    static PlotSeries[] LegendEntries(PlotPanel panel) => panel.Series
        .Where(series => series.Label.Length > 0)
        .GroupBy(series => series.Label)
        .Select(group => group.First())
        .ToArray();

    static int LegendHeight(PlotPanel panel)
    {
        var labels = LegendEntries(panel);
        int height = 0;
        for (int i = 0; i < labels.Length; i += 2)
        {
            int leftLines = Wrap(labels[i].Label, 34).Length;
            int rightLines = i + 1 < labels.Length ? Wrap(labels[i + 1].Label, 34).Length : 1;
            height += Math.Max(leftLines, rightLines) * 12 + 5;
        }
        return height;
    }

    static string F(double x) => x.ToString("0.###", CultureInfo.InvariantCulture);
    static string E(string text) => SecurityElement.Escape(text) ?? "";
    static void Text(StringBuilder svg, double x, double y, string text, int size = 12, string anchor = "start") =>
        svg.Append($"<text x=\"{F(x)}\" y=\"{F(y)}\" font-size=\"{size}\" text-anchor=\"{anchor}\">{E(text)}</text>");
    static void Line(StringBuilder svg, double x1, double y1, double x2, double y2, string color, string dash = "") =>
        svg.Append($"<line x1=\"{F(x1)}\" y1=\"{F(y1)}\" x2=\"{F(x2)}\" y2=\"{F(y2)}\" stroke=\"{color}\" stroke-width=\"1\" stroke-dasharray=\"{dash}\"/>");

    static void Draw(StringBuilder svg, PlotPanel panel, int id, double ox, double oy)
    {
        double left = ox + 88, top = oy + 45, width = 450, height = 255;
        if (panel.TextOnly)
        {
            Text(svg, left, oy + 12, panel.Title, 14);
            var lines = panel.Note.Split('\n');
            for (int i = 0; i < lines.Length; i++)
                Text(svg, left, top + 25 + i * 25, lines[i], 13);
            return;
        }

        double TransformX(double value)
        {
            if (panel.XLog)
                return value > 0 ? Math.Log10(value) : double.NaN;
            if (!panel.XSymLog)
                return value;
            return Math.Abs(value) <= panel.XSymLogThreshold
                ? value / panel.XSymLogThreshold * panel.XSymLogScale * 10 / 9
                : Math.Sign(value) * (panel.XSymLogScale * 10 / 9 + Math.Log10(Math.Abs(value) / panel.XSymLogThreshold));
        }
        double TransformY(double value)
        {
            if (panel.YLog)
                return value > 0 ? Math.Log10(value) : double.NaN;
            if (!panel.YSymLog)
                return value;
            return Math.Abs(value) <= panel.SymLogThreshold
                ? value / panel.SymLogThreshold * 10 / 9
                : Math.Sign(value) * (10.0 / 9 + Math.Log10(Math.Abs(value) / panel.SymLogThreshold));
        }
        double InverseX(double value)
        {
            if (panel.XLog)
                return Math.Pow(10, value);
            if (!panel.XSymLog)
                return value;
            return Math.Abs(value) <= panel.XSymLogScale * 10 / 9
                ? value * panel.XSymLogThreshold * 9 / (10 * panel.XSymLogScale)
                : Math.Sign(value) * panel.XSymLogThreshold * Math.Pow(10, Math.Abs(value) - panel.XSymLogScale * 10 / 9);
        }
        double InverseY(double value)
        {
            if (panel.YLog)
                return Math.Pow(10, value);
            if (!panel.YSymLog)
                return value;
            return Math.Abs(value) <= 10.0 / 9
                ? value * panel.SymLogThreshold * 9 / 10
                : Math.Sign(value) * panel.SymLogThreshold * Math.Pow(10, Math.Abs(value) - 10.0 / 9);
        }

        var xs = panel.Series.SelectMany(s => s.X.Concat(s.XLower ?? []).Concat(s.XUpper ?? [])).Select(TransformX).Where(double.IsFinite).ToArray();
        var ys = panel.Series.SelectMany(s => s.Y.Concat(s.Lower ?? []).Concat(s.Upper ?? [])).Select(TransformY).Where(double.IsFinite).ToArray();
        double xmin = panel.XMin is { } xMinimum ? TransformX(xMinimum) : xs.DefaultIfEmpty(0).Min();
        double xmax = panel.XMax is { } xMaximum ? TransformX(xMaximum) : xs.DefaultIfEmpty(1).Max();
        double ymin = panel.YMin is { } yMinimum ? TransformY(yMinimum) : ys.DefaultIfEmpty(0).Min();
        double ymax = panel.YMax is { } yMaximum ? TransformY(yMaximum) : ys.DefaultIfEmpty(1).Max();
        if (xmin == xmax)
        {
            xmin -= .5;
            xmax += .5;
        }
        if (ymin == ymax)
        {
            ymin -= .5;
            ymax += .5;
        }
        double xspan = xmax - xmin;
        if (panel.XMin is null)
            xmin -= xspan * .04;
        if (panel.XMax is null)
            xmax += xspan * .04;
        if (panel.YMin is null)
            ymin -= (ymax - ymin) * .04;
        if (panel.YMax is null)
            ymax += (ymax - ymin) * .06;
        double ProjectX(double x) => left + (TransformX(x) - xmin) / (xmax - xmin) * width;
        double ProjectY(double y) => top + (panel.ReverseY ? (TransformY(y) - ymin) / (ymax - ymin) : 1 - (TransformY(y) - ymin) / (ymax - ymin)) * height;

        var titleLines = Wrap(panel.Title, 58);
        for (int i = 0; i < titleLines.Length; i++)
            Text(svg, left, oy + 12 + i * 16, titleLines[i], 14);
        svg.Append($"<defs><clipPath id=\"c{id}\"><rect x=\"{F(left)}\" y=\"{F(top)}\" width=\"{width}\" height=\"{height}\"/></clipPath></defs>");
        svg.Append($"<g clip-path=\"url(#c{id})\">");
        foreach (var shade in panel.Shades)
            svg.Append($"<rect x=\"{F(ProjectX(shade.Start))}\" y=\"{F(top)}\" width=\"{F(ProjectX(shade.End) - ProjectX(shade.Start))}\" height=\"{height}\" fill=\"{shade.Color}\"/>");
        svg.Append("</g>");

        var xticks = panel.XTicks ?? Enumerable.Range(0, 6).ToDictionary(
            i => InverseX(xmin + (xmax - xmin) * i / 5),
            i => InverseX(xmin + (xmax - xmin) * i / 5).ToString("G3", CultureInfo.InvariantCulture));
        var yticks = panel.YTicks ?? Enumerable.Range(0, 6).ToDictionary(
            i => InverseY(ymin + (ymax - ymin) * i / 5),
            i => InverseY(ymin + (ymax - ymin) * i / 5).ToString("G3", CultureInfo.InvariantCulture));
        foreach (var (x, label) in xticks)
        {
            if (!(ProjectX(x) >= left - 1 && ProjectX(x) <= left + width + 1))
                continue;
            Line(svg, ProjectX(x), top, ProjectX(x), top + height, "#ededed");
            Text(svg, ProjectX(x), top + height + 19, label, 10, "middle");
        }
        foreach (var (y, label) in yticks)
        {
            if (!(ProjectY(y) >= top - 1 && ProjectY(y) <= top + height + 1))
                continue;
            Line(svg, left, ProjectY(y), left + width, ProjectY(y), "#ededed");
            Text(svg, left - 8, ProjectY(y) + 4, label, 10, "end");
        }

        svg.Append($"<g clip-path=\"url(#c{id})\">");
        foreach (var series in panel.Series)
            DrawSeries(svg, series, ProjectX, ProjectY);
        svg.Append("</g>");
        Line(svg, left, top + height, left + width, top + height, "#444");
        Line(svg, left, top, left, top + height, "#444");
        Text(svg, left + width / 2, top + height + 41, panel.XLabel, 12, "middle");
        svg.Append($"<text transform=\"translate({F(ox + 18)},{F(top + height / 2)}) rotate(-90)\" font-size=\"12\" text-anchor=\"middle\">{E(panel.YLabel)}</text>");
        DrawLegend(svg, panel, left, top + height + 62);
    }

    static void DrawSeries(StringBuilder svg, PlotSeries series, Func<double, double> projectX, Func<double, double> projectY)
    {
        bool Valid(int i) => double.IsFinite(projectX(series.X[i])) && double.IsFinite(projectY(series.Y[i]));
        if (series.Fill && series.Lower is not null && series.Upper is not null)
        {
            var indices = Enumerable.Range(0, series.X.Length).Where(i => Valid(i) && double.IsFinite(projectY(series.Lower[i])) && double.IsFinite(projectY(series.Upper[i]))).ToArray();
            var points = indices.Select(i => $"{F(projectX(series.X[i]))},{F(projectY(series.Lower[i]))}")
                .Concat(indices.Reverse().Select(i => $"{F(projectX(series.X[i]))},{F(projectY(series.Upper[i]))}"));
            svg.Append($"<polygon points=\"{string.Join(' ', points)}\" fill=\"{series.Color}\" opacity=\".17\"/>");
        }

        var path = new StringBuilder();
        bool previous = false;
        double priorY = 0;
        for (int i = 0; i < series.X.Length; i++)
        {
            if (!Valid(i))
            {
                previous = false;
                continue;
            }
            double x = projectX(series.X[i]);
            double y = projectY(series.Y[i]);
            if (series.Points)
                svg.Append($"<circle cx=\"{F(x)}\" cy=\"{F(y)}\" r=\"2.4\" fill=\"{(series.Hollow ? "white" : series.Color)}\" stroke=\"{series.Color}\"/>");
            else if (series.Bars)
            {
                double zero = projectY(0);
                svg.Append($"<rect x=\"{F(x - 9)}\" y=\"{F(Math.Min(y, zero))}\" width=\"18\" height=\"{F(Math.Abs(zero - y))}\" fill=\"{series.Color}\"/>");
            }
            else if (!previous)
                path.Append($"M{F(x)},{F(y)}");
            else
            {
                if (series.Step)
                    path.Append($"L{F(x)},{F(priorY)}");
                path.Append($"L{F(x)},{F(y)}");
            }
            if (!series.Fill && series.Lower is not null && series.Upper is not null && double.IsFinite(projectY(series.Lower[i])) && double.IsFinite(projectY(series.Upper[i])))
            {
                Line(svg, x, projectY(series.Lower[i]), x, projectY(series.Upper[i]), series.Color);
                Line(svg, x - 3, projectY(series.Lower[i]), x + 3, projectY(series.Lower[i]), series.Color);
                Line(svg, x - 3, projectY(series.Upper[i]), x + 3, projectY(series.Upper[i]), series.Color);
            }
            if (series.XLower is not null && series.XUpper is not null && double.IsFinite(projectX(series.XLower[i])) && double.IsFinite(projectX(series.XUpper[i])))
                Line(svg, projectX(series.XLower[i]), y, projectX(series.XUpper[i]), y, series.Color);
            previous = true;
            priorY = y;
        }
        if (path.Length > 0)
            svg.Append($"<path d=\"{path}\" fill=\"none\" stroke=\"{series.Color}\" stroke-width=\"1.6\" stroke-dasharray=\"{series.Dash}\"/>");
    }

    static void DrawLegend(StringBuilder svg, PlotPanel panel, double left, double legendY)
    {
        var labels = LegendEntries(panel);
        for (int i = 0; i < labels.Length; i += 2)
        {
            int lines = 1;
            for (int j = i; j < Math.Min(i + 2, labels.Length); j++)
            {
                double x = left + (j % 2) * 225;
                var wrapped = Wrap(labels[j].Label, 34);
                lines = Math.Max(lines, wrapped.Length);
                Line(svg, x, legendY - 4, x + 15, legendY - 4, labels[j].Color, labels[j].Dash);
                for (int k = 0; k < wrapped.Length; k++)
                    Text(svg, x + 21, legendY + k * 12, wrapped[k], 10);
            }
            legendY += lines * 12 + 5;
        }
        if (panel.Note.Length > 0)
            Text(svg, left, legendY + 2, panel.Note, 10);
    }
}
