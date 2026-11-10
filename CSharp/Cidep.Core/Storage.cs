using System.Text.Json;

namespace Cidep.Core;

public sealed class SpectrumDocument
{
    [System.Text.Json.Serialization.JsonRequired]
    public int Version { get; set; } = 1;
    [System.Text.Json.Serialization.JsonRequired]
    public Dictionary<string, string> Units { get; set; } = new()
        { ["field"] = "Gauss", ["frequency"] = "GHz", ["time"] = "microseconds" };
    public string Source { get; set; } = "simulation";
    public Config? Config { get; set; }
    public Shape? Shape { get; set; }
    public Dictionary<string, double> Metadata { get; set; } = [];
    [System.Text.Json.Serialization.JsonRequired]
    public Dictionary<string, JsonElement> Data { get; set; } = [];

    public Dictionary<string, double[]> ToData() => Data.ToDictionary(entry => entry.Key,
        entry => entry.Key == "total_degeneracy" ? new[] { entry.Value.GetDouble() }
                 : entry.Value.Deserialize<double[]>(Storage.Options)!);
}

public static class Storage
{
    public static readonly JsonSerializerOptions Options = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
        WriteIndented = true
    };
    private static readonly HashSet<string> ArrayNames = ["field", "equilibrium", "rpm", "stm", "tm", "mixed",
        "shaped", "shaped_polarized", "shaped_equilibrium"];

    public static SpectrumDocument Document(Dictionary<string, double[]> data, Config? config = null, Shape? shape = null,
        string source = "simulation", Dictionary<string, double>? metadata = null)
    {
        var document = new SpectrumDocument { Config = config, Shape = shape, Source = source, Metadata = metadata ?? [] };
        foreach (var (key, values) in data)
            document.Data[key] = key == "total_degeneracy" ? JsonSerializer.SerializeToElement(values.Single())
                                : JsonSerializer.SerializeToElement(values);
        Validate(document);
        return document;
    }

    public static void Validate(SpectrumDocument document)
    {
        if (document.Version != 1 || document.Units is null || document.Units.Count != 3
            || document.Units.GetValueOrDefault("field") != "Gauss"
            || document.Units.GetValueOrDefault("frequency") != "GHz"
            || document.Units.GetValueOrDefault("time") != "microseconds")
            throw new ArgumentException("Unsupported JSON version or units");
        if (document.Data is null || !document.Data.ContainsKey("field") || !document.Data.ContainsKey("mixed")
            || !document.Data.ContainsKey("equilibrium"))
            throw new ArgumentException("Missing spectrum arrays");
        if (document.Metadata is null || document.Metadata.Values.Any(value => !double.IsFinite(value)))
            throw new ArgumentException("Metadata must contain finite numeric values");
        if (document.Data.Keys.Any(key => !ArrayNames.Contains(key) && key != "total_degeneracy"))
            throw new ArgumentException("Unknown spectrum data key");
        var data = document.ToData();
        var fields = data["field"];
        if (fields.Length is < 2 or > 16384 || fields.Any(value => !double.IsFinite(value)))
            throw new ArgumentException("Invalid field array");
        double increment = fields[1] - fields[0];
        for (int index = 1; index < fields.Length; index++)
            if (fields[index] <= fields[index - 1]
                || Math.Abs(fields[index] - fields[index - 1] - increment) > 1e-9 + 1e-8 * Math.Abs(increment))
                throw new ArgumentException("Field axis must be increasing and uniformly spaced");
        foreach (var (key, values) in data)
            if (values.Any(value => !double.IsFinite(value))
                || (key == "total_degeneracy" ? values.Length != 1 || values[0] <= 0 : values.Length != fields.Length))
                throw new ArgumentException($"Invalid spectrum array: {key}");
        document.Config?.Validate();
        document.Shape?.Validate();
        if (document.Config is not null && document.Config.Points != fields.Length)
            throw new ArgumentException("Configuration point count does not match data");
    }

    public static void Save(string path, SpectrumDocument document)
    {
        Validate(document);
        File.WriteAllText(path, JsonSerializer.Serialize(document, Options));
    }

    public static SpectrumDocument Load(string path)
    {
        var document = JsonSerializer.Deserialize<SpectrumDocument>(File.ReadAllText(path), Options)
                       ?? throw new ArgumentException("Empty spectrum document");
        Validate(document);
        return document;
    }

    public static SpectrumDocument ImportLegacy(string path)
    {
        using var stream = File.OpenRead(path);
        using var reader = new BinaryReader(stream);
        if (stream.Length < 14) throw new ArgumentException("Truncated legacy header");
        int count = reader.ReadInt16();
        if (count is < 2 or > 2048) throw new ArgumentException("Invalid legacy sample count (expected int16)");
        double sweep = reader.ReadSingle(), centre = reader.ReadSingle(), frequency = reader.ReadSingle();
        if (new[] { sweep, centre, frequency }.Any(value => !double.IsFinite(value) || value <= 0))
            throw new ArgumentException("Invalid legacy field or frequency metadata");
        var data = new Dictionary<string, double[]>();
        var metadata = new Dictionary<string, double> { ["sweep"] = sweep, ["centre"] = centre, ["frequency"] = frequency };
        string source;
        double[] ReadValues() => Enumerable.Range(0, count).Select(_ => (double)reader.ReadSingle()).ToArray();
        if (stream.Length == 14 + 8 * count)
        {
            data["mixed"] = ReadValues();
            data["equilibrium"] = ReadValues();
            source = "legacy-stick-int16-real32";
        }
        else if (stream.Length == 28 + 4 * count)
        {
            double microwave = reader.ReadSingle(), delay = reader.ReadSingle(), interval = reader.ReadSingle();
            int samples = reader.ReadInt16();
            if (new[] { microwave, delay, interval }.Any(value => !double.IsFinite(value)) || samples < 1)
                throw new ArgumentException("Invalid legacy lineshape metadata");
            metadata["microwave"] = microwave;
            metadata["delay"] = delay;
            metadata["interval"] = interval;
            metadata["samples"] = samples;
            data["mixed"] = ReadValues();
            data["shaped"] = (double[])data["mixed"].Clone();
            data["equilibrium"] = new double[count];
            source = "legacy-shaped-int16-real32";
        }
        else throw new ArgumentException("Unsupported legacy layout or file length; record-marked files are not supported");
        data["field"] = Enumerable.Range(0, count).Select(index => centre - sweep / 2 + sweep * index / (count - 1)).ToArray();
        return Document(data, source: source, metadata: metadata);
    }
}