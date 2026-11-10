using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using Cidep.Core;
using Microsoft.Win32;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Wpf;

namespace Cidep.App;

public sealed class MainWindow : Window
{
    private readonly Dictionary<string, FrameworkElement> configInputs = [];
    private readonly Dictionary<string, FrameworkElement> shapeInputs = [];
    private readonly TextBox[] gInputs = [new(), new()];
    private readonly DataGrid[] groupGrids = [new(), new()];
    private readonly ComboBox mode = new() { ItemsSource = new[] { "Sticks", "Relaxation", "Slow exchange" }, SelectedIndex = 0 };
    private readonly ComboBox view = new() { MinWidth = 160, Margin = new(8) };
    private readonly TextBlock status = new() { Margin = new(10, 5, 10, 5), VerticalAlignment = System.Windows.VerticalAlignment.Center };
    private readonly TextBlock readout = new() { Margin = new(12), FontFamily = new("Consolas"), TextWrapping = TextWrapping.Wrap };
    private readonly PlotView plot = new() { Focusable = true, Margin = new(8) };
    private readonly StackPanel inputs = new() { Margin = new(12) };
    private readonly Menu menu = new();
    private readonly Button calculate = IconButton("\uE768", "Calculate");
    private SpectrumDocument? current;
    private Dictionary<string, double[]>? data;
    private int cursorIndex;
    private LineAnnotation? cursor;
    private bool busy;

    public MainWindow()
    {
        Title = "CIDEP - Electron Polarization";
        Width = 1200;
        Height = 820;
        MinWidth = 760;
        MinHeight = 600;
        Background = Brushes.White;
        var root = new DockPanel();
        Content = root;
        DockPanel.SetDock(menu, Dock.Top);
        root.Children.Add(menu);
        var file = new MenuItem { Header = "_File" };
        file.Items.Add(MenuAction("_Open JSON...", (_, _) => Open(false)));
        file.Items.Add(MenuAction("_Save JSON...", (_, _) => Save()));
        file.Items.Add(MenuAction("_Import Legacy Binary...", (_, _) => Open(true)));
        file.Items.Add(new Separator());
        file.Items.Add(MenuAction("E_xit", (_, _) => Close()));
        menu.Items.Add(file);
        var simulation = new MenuItem { Header = "_Simulation" };
        simulation.Items.Add(MenuAction("_Calculate", async (_, _) => await Calculate()));
        simulation.Items.Add(MenuAction("_Restore Defaults", (_, _) => SetInputs(new(), null)));
        simulation.Items.Add(MenuAction("C_lear", (_, _) => Clear()));
        menu.Items.Add(simulation);
        var toolbar = new DockPanel { Background = new SolidColorBrush(Color.FromRgb(241, 245, 242)), LastChildFill = true };
        DockPanel.SetDock(toolbar, Dock.Top);
        root.Children.Add(toolbar);
        var title = new TextBlock { Text = "CIDEP", FontFamily = new("Georgia"), FontSize = 26,
            Foreground = new SolidColorBrush(Color.FromRgb(18, 105, 93)), Margin = new(16, 8, 24, 8) };
        toolbar.Children.Add(title);
        calculate.Click += async (_, _) => await Calculate();
        toolbar.Children.Add(calculate);
        var reset = IconButton("\uE777", "Restore defaults", false);
        reset.Click += (_, _) => SetInputs(new(), null);
        toolbar.Children.Add(reset);
        toolbar.Children.Add(status);
        var footer = new Border { BorderBrush = Brushes.LightGray, BorderThickness = new(0, 1, 0, 0), Child = readout };
        DockPanel.SetDock(footer, Dock.Bottom);
        root.Children.Add(footer);
        var body = new Grid();
        body.ColumnDefinitions.Add(new() { Width = new(350) });
        body.ColumnDefinitions.Add(new() { Width = new(1, GridUnitType.Star) });
        root.Children.Add(body);
        body.Children.Add(new ScrollViewer { Content = inputs, VerticalScrollBarVisibility = ScrollBarVisibility.Auto,
            HorizontalScrollBarVisibility = ScrollBarVisibility.Disabled, Background = new SolidColorBrush(Color.FromRgb(248, 249, 248)) });
        var graph = new DockPanel();
        Grid.SetColumn(graph, 1);
        body.Children.Add(graph);
        DockPanel.SetDock(view, Dock.Top);
        graph.Children.Add(view);
        graph.Children.Add(plot);
        view.SelectionChanged += (_, _) => Draw();
        plot.PreviewMouseLeftButtonDown += (_, args) =>
        {
            if (data is null || plot.ActualModel is null) return;
            var position = args.GetPosition(plot);
            var axes = plot.ActualModel.Axes.First(axis => axis.Position == AxisPosition.Bottom);
            double field = axes.InverseTransform(position.X);
            SetCursor((int)Math.Round((field - data["field"][0]) / (data["field"][1] - data["field"][0])));
            plot.Focus();
            args.Handled = true;
        };
        plot.PreviewKeyDown += (_, args) =>
        {
            if (args.Key is Key.Left or Key.Right)
            {
                SetCursor(cursorIndex + (args.Key == Key.Right ? 1 : -1));
                args.Handled = true;
            }
        };
        BuildInputs();
        SetInputs(new(), null);
        current = Storage.Document(Simulation.Run(new()), new());
        SetDocument(current);
    }

    private static MenuItem MenuAction(string title, RoutedEventHandler action)
    {
        var item = new MenuItem { Header = title };
        item.Click += action;
        return item;
    }

    private static Button IconButton(string glyph, string title, bool withText = true)
    {
        var content = new StackPanel { Orientation = Orientation.Horizontal };
        content.Children.Add(new TextBlock { Text = glyph, FontFamily = new("Segoe MDL2 Assets"), VerticalAlignment = System.Windows.VerticalAlignment.Center });
        if (withText) content.Children.Add(new TextBlock { Text = title, Margin = new(8, 0, 0, 0) });
        return new Button { Content = content, ToolTip = title, VerticalAlignment = System.Windows.VerticalAlignment.Center, MinHeight = 32 };
    }

    private void BuildInputs()
    {
        var tabs = new TabControl();
        for (int index = 0; index < 2; index++)
        {
            var panel = new StackPanel { Margin = new(4) };
            panel.Children.Add(Row("g-value", gInputs[index]));
            var grid = groupGrids[index];
            grid.AutoGenerateColumns = false;
            grid.CanUserAddRows = true;
            grid.CanUserDeleteRows = true;
            grid.MinHeight = 125;
            grid.MaxHeight = 280;
            grid.HeadersVisibility = DataGridHeadersVisibility.Column;
            foreach (var (header, property) in new[] { ("Nuclei", "Count"), ("Coupling / G", "Coupling"), ("Spin", "Spin") })
                grid.Columns.Add(new DataGridTextColumn { Header = header, Width = new DataGridLength(1, DataGridLengthUnitType.Star),
                    Binding = new Binding(property) { UpdateSourceTrigger = UpdateSourceTrigger.LostFocus, ValidatesOnExceptions = true } });
            panel.Children.Add(grid);
            tabs.Items.Add(new TabItem { Header = $"Radical {index + 1}", Content = panel });
        }
        inputs.Children.Add(new Expander { Header = "Radicals", IsExpanded = true, Content = tabs });
        var labels = new Dictionary<string, string>
        {
            ["Identical"] = "Identical radicals", ["SecondOrder"] = "Second-order splitting", ["CentreResonance"] = "Centre on resonance",
            ["Centre"] = "Centre / G", ["Sweep"] = "Sweep / G", ["Frequency"] = "Frequency / GHz", ["Points"] = "Points",
            ["Pair"] = "Pair (S, T, F)", ["Linear"] = "Linear coefficient", ["Exponent"] = "Q exponent (optional)",
            ["ObservedRadical"] = "Observed radical", ["NetEffect"] = "ST-minus net effect", ["Normalize"] = "Normalize components",
            ["RpmWeight"] = "ST-zero RPM", ["StmWeight"] = "ST-minus RPM", ["TmWeight"] = "Triplet mechanism", ["TmRatio"] = "TM polarization ratio"
        };
        AddFields("Field and physics", typeof(Config), configInputs, labels.Where(entry => !entry.Key.EndsWith("Weight") && entry.Key != "TmRatio"));
        AddFields("Mechanism mixture", typeof(Config), configInputs, labels.Where(entry => entry.Key.EndsWith("Weight") || entry.Key == "TmRatio"));
        var shapeLabels = new Dictionary<string, string>
        {
            ["T1"] = "T1 / us", ["T2"] = "T2 / us", ["Microwave"] = "Microwave / rad MHz", ["Delay"] = "Delay / us",
            ["Interval"] = "Interval / us", ["Samples"] = "Samples", ["Points"] = "Lineshape points",
            ["ExchangeTime"] = "Exchange time / us", ["Polarization"] = "Exchange polarization"
        };
        inputs.Children.Add(Row("Calculation", mode));
        AddFields("Lineshape", typeof(Shape), shapeInputs, shapeLabels);
    }

    private static Grid Row(string label, FrameworkElement control)
    {
        var row = new Grid { Margin = new(0, 3, 0, 3) };
        row.ColumnDefinitions.Add(new() { Width = new(176) });
        row.ColumnDefinitions.Add(new() { Width = new(1, GridUnitType.Star) });
        row.Children.Add(new TextBlock { Text = label, VerticalAlignment = System.Windows.VerticalAlignment.Center, TextWrapping = TextWrapping.Wrap });
        Grid.SetColumn(control, 1);
        row.Children.Add(control);
        return row;
    }

    private void AddFields(string title, Type type, Dictionary<string, FrameworkElement> map, IEnumerable<KeyValuePair<string, string>> fields)
    {
        var panel = new StackPanel();
        foreach (var (name, label) in fields)
        {
            var property = type.GetProperty(name)!;
            FrameworkElement input = property.PropertyType == typeof(bool) ? new CheckBox()
                : name == "Pair" ? new ComboBox { ItemsSource = new[] { "S", "T", "F" } }
                : name == "ObservedRadical" ? new ComboBox { ItemsSource = new[] { 1, 2 } } : new TextBox();
            map[name] = input;
            panel.Children.Add(Row(label, input));
        }
        inputs.Children.Add(new Expander { Header = title, Content = panel, IsExpanded = title == "Mechanism mixture" });
    }

    private static void WriteFields(object source, Dictionary<string, FrameworkElement> controls)
    {
        foreach (var (name, control) in controls)
        {
            object? value = source.GetType().GetProperty(name)!.GetValue(source);
            if (control is CheckBox check) check.IsChecked = (bool)value!;
            else if (control is ComboBox combo) combo.SelectedItem = value;
            else ((TextBox)control).Text = Convert.ToString(value, CultureInfo.InvariantCulture) ?? "";
        }
    }

    private static void ReadFields(object destination, Dictionary<string, FrameworkElement> controls)
    {
        foreach (var (name, control) in controls)
        {
            PropertyInfo property = destination.GetType().GetProperty(name)!;
            object? value;
            if (control is CheckBox check) value = check.IsChecked == true;
            else if (control is ComboBox combo) value = combo.SelectedItem;
            else
            {
                string text = ((TextBox)control).Text.Trim();
                Type type = Nullable.GetUnderlyingType(property.PropertyType) ?? property.PropertyType;
                if (text == "" && Nullable.GetUnderlyingType(property.PropertyType) is not null) value = null;
                else value = Convert.ChangeType(text, type, CultureInfo.InvariantCulture);
            }
            property.SetValue(destination, value);
        }
    }

    private void SetInputs(Config config, Shape? shape)
    {
        WriteFields(config, configInputs);
        WriteFields(shape ?? new Shape { ExchangeTime = 5.5 }, shapeInputs);
        Radical[] radicals = [config.Radical1, config.Radical2];
        for (int index = 0; index < 2; index++)
        {
            gInputs[index].Text = radicals[index].G.ToString(CultureInfo.InvariantCulture);
            groupGrids[index].ItemsSource = new ObservableCollection<Group>(radicals[index].Groups.Select(group =>
                new Group { Count = group.Count, Coupling = group.Coupling, Spin = group.Spin }));
        }
        mode.SelectedIndex = shape is null ? 0 : shape.ExchangeTime is null ? 1 : 2;
        status.Text = "Defaults / parameters loaded";
    }

    private (Config Config, Shape? Shape) ReadInputs()
    {
        var config = new Config();
        ReadFields(config, configInputs);
        var radicals = new Radical[2];
        for (int index = 0; index < 2; index++)
        {
            var grid = groupGrids[index];
            if (!grid.CommitEdit(DataGridEditingUnit.Cell, true) || !grid.CommitEdit(DataGridEditingUnit.Row, true))
                throw new ArgumentException("Correct the nuclear group values before calculating");
            radicals[index] = new Radical { G = double.Parse(gInputs[index].Text, CultureInfo.InvariantCulture),
                Groups = ((ObservableCollection<Group>)grid.ItemsSource).Select(group =>
                    new Group { Count = group.Count, Coupling = group.Coupling, Spin = group.Spin }).ToList() };
        }
        config.Radical1 = radicals[0];
        config.Radical2 = radicals[1];
        config.Validate();
        Shape? shape = null;
        if (mode.SelectedIndex > 0)
        {
            shape = new();
            ReadFields(shape, shapeInputs);
            if (mode.SelectedIndex == 1) shape.ExchangeTime = null;
            else if (shape.ExchangeTime is null) throw new ArgumentException("An exchange time is required");
            shape.Validate();
        }
        return (config, shape);
    }

    private async Task Calculate(bool reportErrors = true)
    {
        if (busy) return;
        busy = true;
        calculate.IsEnabled = menu.IsEnabled = inputs.IsEnabled = false;
        status.Text = "Calculating...";
        try
        {
            var (config, shape) = ReadInputs();
            var result = await Task.Run(() =>
            {
                var result = Simulation.Run(config);
                return shape is null ? result : Simulation.Broaden(result, shape);
            });
            SetDocument(Storage.Document(result, config, shape));
            status.Text = $"{config.Points} points";
        }
        catch (Exception error)
        {
            status.Text = error.GetBaseException().Message;
            if (!reportErrors) throw;
        }
        finally
        {
            calculate.IsEnabled = menu.IsEnabled = inputs.IsEnabled = true;
            busy = false;
        }
    }

    private void SetDocument(SpectrumDocument document)
    {
        current = document;
        data = document.ToData();
        view.ItemsSource = data.Keys.Where(key => key != "field" && key != "total_degeneracy").ToList();
        view.SelectedItem = data.ContainsKey("shaped") ? "shaped" : "mixed";
        Draw();
    }

    private void Draw()
    {
        if (data is null || view.SelectedItem is not string key) return;
        var model = new PlotModel { Title = key.Replace('_', ' '), TitleFontSize = 17, DefaultFont = "Segoe UI",
            PlotAreaBorderColor = OxyColors.LightGray, Background = OxyColors.White };
        model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "Magnetic field / Gauss", MajorGridlineStyle = LineStyle.Dot });
        model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Signal", MajorGridlineStyle = LineStyle.Dot });
        var series = new LineSeries { Color = OxyColor.FromRgb(8, 127, 114), StrokeThickness = 1.3 };
        for (int index = 0; index < data["field"].Length; index++) series.Points.Add(new(data["field"][index], data[key][index]));
        model.Series.Add(series);
        cursor = new LineAnnotation { Type = LineAnnotationType.Vertical, Color = OxyColor.FromRgb(196, 63, 82), StrokeThickness = 1 };
        model.Annotations.Add(cursor);
        plot.Model = model;
        SetCursor(cursorIndex);
    }

    private void SetCursor(int index)
    {
        if (data is null || cursor is null || view.SelectedItem is not string key) return;
        cursorIndex = Math.Clamp(index, 0, data["field"].Length - 1);
        cursor.X = data["field"][cursorIndex];
        readout.Text = $"{cursorIndex + 1} / {data["field"].Length}     {cursor.X:F6} G     {data[key][cursorIndex]:G8}";
        plot.InvalidatePlot(false);
    }

    private void Open(bool legacy)
    {
        var dialog = new OpenFileDialog { Filter = legacy ? "Legacy spectra|*.*" : "JSON spectra|*.json" };
        if (dialog.ShowDialog(this) != true) return;
        try
        {
            var document = legacy ? Storage.ImportLegacy(dialog.FileName) : Storage.Load(dialog.FileName);
            if (document.Config is not null) SetInputs(document.Config, document.Shape);
            SetDocument(document);
            status.Text = Path.GetFileName(dialog.FileName);
        }
        catch (Exception error) { MessageBox.Show(this, error.Message, "Cannot open spectrum", MessageBoxButton.OK, MessageBoxImage.Error); }
    }

    private void Save()
    {
        if (current is null) return;
        var dialog = new SaveFileDialog { Filter = "JSON spectra|*.json", FileName = "spectrum.json", DefaultExt = ".json" };
        if (dialog.ShowDialog(this) != true) return;
        try { Storage.Save(dialog.FileName, current); status.Text = "JSON saved"; }
        catch (Exception error) { MessageBox.Show(this, error.Message, "Cannot save spectrum", MessageBoxButton.OK, MessageBoxImage.Error); }
    }

    private void Clear()
    {
        data = null;
        current = null;
        plot.Model = null;
        view.ItemsSource = null;
        readout.Text = "";
        status.Text = "Cleared";
    }

    public async Task SmokeTest()
    {
        await Calculate(false);
        if (data is null || data["field"].Length != 1024) throw new Exception("Initial calculation failed");
        SetCursor(100);
        plot.RaiseEvent(new KeyEventArgs(Keyboard.PrimaryDevice, PresentationSource.FromVisual(plot), 0, Key.Right)
            { RoutedEvent = Keyboard.PreviewKeyDownEvent });
        if (cursorIndex != 101) throw new Exception("Keyboard cursor did not move exactly one sample");
        mode.SelectedIndex = 1;
        ((TextBox)shapeInputs["Points"]).Text = "32";
        await Calculate(false);
        if (!data!.ContainsKey("shaped")) throw new Exception("Lineshape calculation failed");
        Directory.CreateDirectory("artifacts");
        foreach (var (width, height) in new[] { (1200, 820), (800, 620) })
        {
            Width = width;
            Height = height;
            await Dispatcher.InvokeAsync(UpdateLayout, System.Windows.Threading.DispatcherPriority.ApplicationIdle);
            plot.InvalidatePlot(true);
            await Dispatcher.InvokeAsync(UpdateLayout, System.Windows.Threading.DispatcherPriority.ApplicationIdle);
            var bitmap = new RenderTargetBitmap((int)ActualWidth, (int)ActualHeight, 96, 96, PixelFormats.Pbgra32);
            bitmap.Render(this);
            var pixels = new byte[bitmap.PixelWidth * bitmap.PixelHeight * 4];
            bitmap.CopyPixels(pixels, bitmap.PixelWidth * 4, 0);
            int colored = 0;
            for (int offset = 0; offset < pixels.Length; offset += 4)
                if (pixels[offset + 1] > pixels[offset + 2] + 25 && pixels[offset + 1] > 65) colored++;
            if (colored < 100) throw new Exception("Spectrum render is blank");
            var encoder = new PngBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(bitmap));
            using var stream = File.Create($"artifacts/wpf-{width}.png");
            encoder.Save(stream);
        }
        Storage.Save("artifacts/csharp-spectrum.json", current!);
        File.WriteAllText("artifacts/wpf-smoke.txt", "PASS: default and shaped calculations, cursor key, 1200/800 pixel checks, JSON export");
    }
}