using Cidep.Core;
using Xunit;

namespace Cidep.Tests;

public class NumericalTests
{
    [Fact]
    public void ExactTriplet()
    {
        var states = Simulation.GroupStates(new Group { Count = 2, Coupling = 1 });
        Assert.Equal(new double[] { -1, 0, 1 }, states.Select(state => state.Projection));
        Assert.Equal(new double[] { 1, 2, 1 }, states.Select(state => state.Degeneracy));
    }

    [Fact]
    public void AngularMomentumConservesStates()
    {
        foreach (double spin in new double[] { 0, 0.5, 1, 1.5 })
            foreach (int count in new[] { 1, 2, 3, 6 })
                Assert.Equal(Math.Pow(2 * spin + 1, count), Simulation.GroupStates(new Group
                    { Count = count, Spin = spin }, true).Sum(state => state.Degeneracy));
    }

    [Fact]
    public void IndependentSpinHalfRpm()
    {
        var config = new Config { Radical1 = new() { G = 2, Groups = [new() { Coupling = 2 }] },
            SecondOrder = false, Normalize = false, Pair = "S", Points = 101, Sweep = 10 };
        var result = Simulation.Run(config);
        var nonzero = result["rpm"].Where(value => value != 0).ToArray();
        Assert.Equal(2, nonzero.Length);
        Assert.Equal(1.138871796 * Math.Sqrt(2) / 4, nonzero[0], 12);
        Assert.Equal(-nonzero[0], nonzero[1], 12);
        Assert.Equal(2, result["equilibrium"].Sum());
    }

    [Fact]
    public void DefaultAndConvolutionAreFinite()
    {
        var result = Simulation.Run(new());
        Assert.Equal(128, result["total_degeneracy"][0]);
        Assert.Equal(100, result["rpm"].Max() - result["rpm"].Min(), 10);
        foreach (double? exchange in new double?[] { null, 5.5 })
        {
            var shaped = Simulation.Broaden(result, new() { Points = 16, Samples = 3, ExchangeTime = exchange });
            Assert.All(shaped["shaped"], value => Assert.True(double.IsFinite(value)));
        }
    }

    [Fact]
    public void TimeZeroAndSymmetry()
    {
        var kernels = Simulation.Kernels(new() { Delay = 0, Samples = 1, Points = 16 }, 0.1);
        for (int index = 0; index < 16; index++)
        {
            Assert.Equal(0, kernels.Polarized[index], 10);
            Assert.Equal(0, kernels.Equilibrium[index], 10);
            Assert.Equal(kernels.Polarized[index], kernels.Polarized[15 - index]);
        }
    }

    [Fact]
    public void InvalidInputsAreRejected()
    {
        Assert.Throws<ArgumentException>(() => Simulation.Run(new() { Points = 1 }));
        Assert.Throws<ArgumentException>(() => Simulation.Run(new() { Sweep = 0 }));
        Assert.Throws<ArgumentException>(() => Simulation.Run(new() { Exponent = -1 }));
    }
}