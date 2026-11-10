using System.Text.Json;
using Cidep.Core;
using Xunit;

namespace Cidep.Tests;

public class PersistenceTests
{
    [Fact]
    public void PythonFixturesMatchIndependentCSharpCalculations()
    {
        var fixtures = JsonDocument.Parse(File.ReadAllText(Path.Combine(AppContext.BaseDirectory, "fixtures.json")));
        foreach (var fixture in fixtures.RootElement.EnumerateArray())
        {
            var document = fixture.GetProperty("document").Deserialize<SpectrumDocument>(Storage.Options)!;
            Storage.Validate(document);
            var expected = document.ToData();
            var actual = Simulation.Run(document.Config!);
            if (document.Shape is not null) actual = Simulation.Broaden(actual, document.Shape);
            foreach (var (key, values) in expected)
                for (int index = 0; index < values.Length; index++)
                    Assert.True(Math.Abs(values[index] - actual[key][index]) <= 1e-7 + 1e-8 * Math.Abs(values[index]),
                        $"{fixture.GetProperty("name")} {key}[{index}]: Python={values[index]:G17}, C#={actual[key][index]:G17}");
            string path = Path.GetTempFileName();
            try
            {
                Storage.Save(path, Storage.Document(actual, document.Config, document.Shape));
                Assert.Equal(actual["mixed"], Storage.Load(path).ToData()["mixed"]);
            }
            finally { File.Delete(path); }
        }
    }

    [Fact]
    public void LegacyImportAndMalformedData()
    {
        string path = Path.GetTempFileName();
        try
        {
            using (var writer = new BinaryWriter(File.Create(path)))
            {
                writer.Write((short)4);
                writer.Write(10f);
                writer.Write(3483f);
                writer.Write(9.777f);
                foreach (float value in new float[] { 1, -2, 3, 0, 1, 2, 3, 4 }) writer.Write(value);
            }
            var document = Storage.ImportLegacy(path);
            Assert.Equal(new double[] { 1, -2, 3, 0 }, document.ToData()["mixed"]);
            Assert.Equal(9.777, document.Metadata["frequency"], 5);
            using (var writer = new BinaryWriter(File.Create(path)))
            {
                writer.Write((short)4);
                foreach (float value in new float[] { 10, 3483, 9.777f, 0.1f, 1, 0.01f }) writer.Write(value);
                writer.Write((short)150);
                foreach (float value in new float[] { 1, 2, 3, 4 }) writer.Write(value);
            }
            Assert.Equal(150, Storage.ImportLegacy(path).Metadata["samples"]);
            document.Data["field"] = JsonSerializer.SerializeToElement(new[] { 1, 1, 2, 3 });
            Assert.Throws<ArgumentException>(() => Storage.Validate(document));
            File.WriteAllBytes(path, [1, 2, 3]);
            Assert.Throws<ArgumentException>(() => Storage.ImportLegacy(path));
            File.WriteAllText(path, "{\"units\":{},\"data\":{}}");
            Assert.Throws<JsonException>(() => Storage.Load(path));
        }
        finally { File.Delete(path); }
    }
}