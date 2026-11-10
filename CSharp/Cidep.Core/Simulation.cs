using System.Numerics;
using MathNet.Numerics;

namespace Cidep.Core;

public sealed class Group
{
    public int Count { get; set; } = 1;
    public double Coupling { get; set; } = 21.94;
    public double Spin { get; set; } = 0.5;
}

public sealed class Radical
{
    public double G { get; set; } = 2.0025;
    public List<Group> Groups { get; set; } = [new(), new() { Count = 6, Coupling = 24.69 }];
}

public sealed class Config
{
    public Radical Radical1 { get; set; } = new();
    public Radical Radical2 { get; set; } = new() { G = 2.0007, Groups = [new() { Coupling = 1 }] };
    public bool Identical { get; set; } = true;
    public bool SecondOrder { get; set; } = true;
    public bool CentreResonance { get; set; } = true;
    public double Centre { get; set; } = 3483;
    public double Sweep { get; set; } = 180;
    public double Frequency { get; set; } = 9.777;
    public int Points { get; set; } = 1024;
    public string Pair { get; set; } = "T";
    public double Linear { get; set; }
    public double? Exponent { get; set; }
    public int ObservedRadical { get; set; } = 1;
    public bool NetEffect { get; set; } = true;
    public bool Normalize { get; set; } = true;
    public double RpmWeight { get; set; } = 1;
    public double StmWeight { get; set; }
    public double TmWeight { get; set; }
    public double TmRatio { get; set; } = -50;

    public void Validate()
    {
        if (Points is < 2 or > 16384 || Pair is not ("S" or "T" or "F") || ObservedRadical is < 1 or > 2)
            throw new ArgumentException("Invalid points, pair type, or observed radical");
        if (new[] { Centre, Sweep, Frequency }.Any(value => !double.IsFinite(value) || value <= 0))
            throw new ArgumentException("Field and frequency parameters must be finite and positive");
        if (new[] { Linear, RpmWeight, StmWeight, TmWeight, TmRatio }.Any(value => !double.IsFinite(value))
            || Exponent is double exponent && (!double.IsFinite(exponent) || exponent <= 0))
            throw new ArgumentException("Invalid polarization parameters");
        foreach (var radical in new[] { Radical1, Radical2 })
        {
            if (radical is null || !double.IsFinite(radical.G) || radical.G <= 0 || radical.Groups is null || radical.Groups.Count > 20)
                throw new ArgumentException("Invalid radical");
            foreach (var group in radical.Groups)
                if (group is null || group.Count is < 1 or > 15 || !double.IsFinite(group.Coupling)
                    || !double.IsFinite(group.Spin) || group.Spin is < 0 or > 5
                    || Math.Abs(2 * group.Spin - Math.Round(2 * group.Spin)) > 1e-10)
                    throw new ArgumentException("Invalid nuclear group");
        }
    }
}

public sealed class Shape
{
    public double T1 { get; set; } = 5;
    public double T2 { get; set; } = 0.6;
    public double Microwave { get; set; } = 0.1;
    public double Delay { get; set; } = 1;
    public double Interval { get; set; } = 0.01;
    public int Samples { get; set; } = 150;
    public int Points { get; set; } = 512;
    public double? ExchangeTime { get; set; }
    public double Polarization { get; set; } = -10;

    public void Validate()
    {
        if (new[] { T1, T2, Microwave, Interval }.Any(value => !double.IsFinite(value) || value <= 0)
            || !double.IsFinite(Delay) || Delay < 0 || !double.IsFinite(Polarization)
            || Samples is < 1 or > 10000 || Points is < 4 or > 2048 || Points % 2 != 0
            || ExchangeTime is double exchange && (!double.IsFinite(exchange) || exchange <= 0))
            throw new ArgumentException("Invalid relaxation, sampling, or exchange parameters");
    }
}

public readonly record struct NuclearState(double Projection, double Angular, double Degeneracy);
public readonly record struct Line(double Energy, double Degeneracy, double Multiplet);

public static class Simulation
{
    public static List<NuclearState> GroupStates(Group group, bool secondOrder = false)
    {
        double[] coefficients = [1];
        int width = (int)Math.Round(2 * group.Spin) + 1;
        for (int nucleus = 0; nucleus < group.Count; nucleus++)
        {
            var next = new double[coefficients.Length + width - 1];
            for (int index = 0; index < coefficients.Length; index++)
                for (int offset = 0; offset < width; offset++)
                    next[index + offset] += coefficients[index];
            coefficients = next;
        }
        double maximum = group.Count * group.Spin;
        var states = new List<NuclearState>();
        if (!secondOrder)
            return coefficients.Select((weight, index) => new NuclearState(index - maximum, 0, weight)).ToList();
        double previous = 0;
        for (int index = 0; index <= Math.Floor(maximum); index++)
        {
            double angular = maximum - index;
            double multiplicity = coefficients[index] - previous;
            previous = coefficients[index];
            if (multiplicity == 0) continue;
            for (double projection = -angular; projection <= angular; projection++)
                states.Add(new(projection, angular, multiplicity));
        }
        return states;
    }

    public static List<Line> RadicalStates(Radical radical, double centre, bool secondOrder = false, bool stm = false)
    {
        var groups = radical.Groups.Select(group => GroupStates(group, secondOrder || stm)).ToArray();
        long count = 1;
        foreach (var states in groups)
        {
            count *= states.Count;
            if (count > 4096) throw new ArgumentException("More than 4096 nuclear states; reduce groups");
        }
        var lines = new List<Line> { new(0, 1, 0) };
        for (int index = 0; index < groups.Length; index++)
        {
            var group = radical.Groups[index];
            var next = new List<Line>();
            foreach (var line in lines)
                foreach (var state in groups[index])
                {
                    double energy = group.Coupling * state.Projection;
                    if (secondOrder || stm)
                        energy += group.Coupling * group.Coupling * (state.Angular * (state.Angular + 1)
                                  - state.Projection * state.Projection) / (2 * centre);
                    double multiplet = stm ? (state.Angular - state.Projection + 1) * (state.Angular + state.Projection)
                                             * group.Coupling * group.Coupling * state.Degeneracy : 0;
                    next.Add(new(line.Energy + energy, line.Degeneracy * state.Degeneracy, line.Multiplet + multiplet));
                }
            lines = next;
        }
        return lines;
    }

    public static double[] Sample(double[] fields, IEnumerable<double> positions, IEnumerable<double> amplitudes)
    {
        double increment = fields[1] - fields[0];
        var result = new double[fields.Length];
        foreach (var (position, amplitude) in positions.Zip(amplitudes))
        {
            int lower = (int)Math.Floor((position - fields[0]) / increment);
            for (int index = lower; index <= lower + 1; index++)
                if (index >= 0 && index < fields.Length && Math.Abs(fields[index] - position) <= increment * (0.5 + 1e-10))
                    result[index] += amplitude;
        }
        return result;
    }

    public static Dictionary<string, double[]> Run(Config config)
    {
        config.Validate();
        Radical[] radicals = [config.Radical1, config.Identical ? config.Radical1 : config.Radical2];
        double centre = config.CentreResonance ? 714.4844244 * config.Frequency / radicals[0].G : config.Centre;
        var fields = Enumerable.Range(0, config.Points).Select(index => centre - config.Sweep / 2
                          + config.Sweep * index / (config.Points - 1)).ToArray();
        var normal = radicals.Select(radical => RadicalStates(radical, config.Centre,
                        config.SecondOrder && radical.Groups.Any(group => Math.Abs(group.Coupling) > 10))).ToArray();
        var stm = radicals.Select(radical => RadicalStates(radical, config.Centre, stm: true)).ToArray();
        int own = config.ObservedRadical - 1;
        int other = 1 - own;
        var states = normal[own];
        var partner = normal[other];
        double difference = config.Centre * (radicals[own].G - radicals[other].G) / 2.00232;
        double scale = config.Pair == "S" ? 1.138871796 : config.Pair == "T" ? -1.138871796 : -0.281368326;
        double total = states.Sum(line => line.Degeneracy) * partner.Sum(line => line.Degeneracy);
        var polarization = new double[states.Count];
        Parallel.For(0, states.Count, index =>
        {
            double value = 0;
            foreach (var counter in partner)
            {
                double delta = states[index].Energy + difference - counter.Energy;
                double absolute = Math.Abs(delta);
                double dependence = config.Exponent is double exponent ? Math.Pow(absolute, exponent)
                                    : Math.Sqrt(absolute) - config.Linear * absolute;
                value += Math.Sign(delta) * dependence * counter.Degeneracy;
            }
            polarization[index] = value * scale * states[index].Degeneracy / total;
        });
        double resonance = 714.4844244 * config.Frequency / radicals[own].G;
        var positions = states.Select(line => resonance - line.Energy).ToArray();
        var equilibrium = Sample(fields, positions, states.Select(line => line.Degeneracy));
        var rpm = Sample(fields, positions, polarization);
        double net = config.NetEffect ? stm[other].Sum(line => line.Multiplet) / stm[other].Sum(line => line.Degeneracy) : 0;
        double degeneracy = stm[own].Sum(line => line.Degeneracy);
        var stminus = Sample(fields, stm[own].Select(line => resonance - line.Energy),
            stm[own].Select(line => -(line.Multiplet + net * line.Degeneracy) / degeneracy));
        var tm = equilibrium.Select(value => config.TmRatio * value).ToArray();
        if (config.Normalize)
            foreach (var component in new[] { rpm, stminus, tm })
            {
                double extent = component.Max() - component.Min();
                if (extent > 0)
                    for (int index = 0; index < component.Length; index++) component[index] *= 100 / extent;
            }
        var mixed = Enumerable.Range(0, fields.Length).Select(index => config.RpmWeight * rpm[index]
                      + config.StmWeight * stminus[index] + config.TmWeight * tm[index]).ToArray();
        return new() { ["field"] = fields, ["equilibrium"] = equilibrium, ["rpm"] = rpm,
            ["stm"] = stminus, ["tm"] = tm, ["mixed"] = mixed, ["total_degeneracy"] = [degeneracy] };
    }

    public static (double[] Polarized, double[] Equilibrium) Kernels(Shape shape, double increment, double? fraction = null)
    {
        shape.Validate();
        double rate1 = 1 / shape.T1, rate2 = 1 / shape.T2, microwave = shape.Microwave;
        var pol = new double[shape.Points];
        var eq = new double[shape.Points];
        Parallel.For(0, shape.Points / 2, index =>
        {
            double delta2 = Math.Pow(index * increment * 17.60844567, 2);
            double[] coefficients;
            double alpha = rate2, gamma2 = 0;
            if (fraction is null)
                coefficients = [rate1 * rate2 * rate2 + microwave * microwave * rate2 + delta2 * rate1,
                    rate2 * rate2 + 2 * rate1 * rate2 + microwave * microwave + delta2, rate1 + 2 * rate2, 1];
            else
            {
                double exchange = shape.ExchangeTime > 100 * shape.T1 ? 0 : 1 / shape.ExchangeTime!.Value;
                double transfer = (1 - fraction.Value) * exchange;
                alpha = rate2 + transfer;
                double beta = rate1 + transfer, gamma = rate1 + exchange - transfer;
                double sumrate = beta + gamma, product = beta * gamma - transfer * (exchange - transfer);
                gamma2 = rate1 + exchange;
                coefficients = [alpha * alpha * product + microwave * microwave * alpha * gamma + product * delta2,
                    2 * alpha * product + sumrate * alpha * alpha + microwave * microwave * (alpha + gamma) + sumrate * delta2,
                    product + alpha * alpha + 2 * alpha * sumrate + microwave * microwave + delta2, sumrate + 2 * alpha, 1];
            }
            var roots = FindRoots.Polynomial(coefficients);
            Complex polarized = 0, equilibrium = 0, steadySum = 0;
            foreach (var root in roots)
            {
                Complex derivative = 0;
                for (int degree = coefficients.Length - 1; degree > 0; degree--)
                    derivative = derivative * root + degree * coefficients[degree];
                if (derivative.Magnitude < 1e-12) throw new ArgumentException("Coincident relaxation roots; adjust parameters");
                Complex amplitude = microwave * (root + alpha) / derivative;
                if (fraction is not null) amplitude *= root + gamma2;
                Complex steady = rate1 * amplitude / root;
                Complex factors = Complex.Exp(shape.Delay * root);
                if (shape.Samples > 1)
                    factors *= Expm1(shape.Samples * shape.Interval * root) / Expm1(shape.Interval * root);
                polarized += amplitude * factors;
                equilibrium += steady * factors;
                steadySum += steady;
            }
            double polarizedValue = polarized.Real;
            double equilibriumValue = (equilibrium - shape.Samples * steadySum).Real;
            pol[shape.Points / 2 - 1 - index] = pol[shape.Points / 2 + index] = polarizedValue;
            eq[shape.Points / 2 - 1 - index] = eq[shape.Points / 2 + index] = equilibriumValue;
        });
        return (pol, eq);
    }

    private static Complex Expm1(Complex value)
    {
        if (value.Magnitude > 1e-4) return Complex.Exp(value) - 1;
        Complex term = value, sum = value;
        for (int degree = 2; degree <= 12; degree++)
        {
            term *= value / degree;
            sum += term;
        }
        return sum;
    }

    private static double[] Convolve(double[] sticks, double[] kernel)
    {
        var result = new double[sticks.Length];
        Parallel.For(0, result.Length, index =>
        {
            double sum = 0;
            int fullIndex = index + kernel.Length / 2;
            for (int offset = Math.Max(0, fullIndex - kernel.Length + 1); offset <= Math.Min(sticks.Length - 1, fullIndex); offset++)
                sum += sticks[offset] * kernel[fullIndex - offset];
            result[index] = sum;
        });
        return result;
    }

    public static Dictionary<string, double[]> Broaden(Dictionary<string, double[]> data, Shape shape)
    {
        shape.Validate();
        double increment = data["field"][1] - data["field"][0];
        double[] pol, eq;
        if (shape.ExchangeTime is null)
        {
            var kernels = Kernels(shape, increment);
            pol = Convolve(data["mixed"], kernels.Polarized);
            eq = Convolve(data["equilibrium"], kernels.Equilibrium);
        }
        else
        {
            if (!data.TryGetValue("total_degeneracy", out var degeneracy))
                throw new ArgumentException("Slow exchange requires a simulated total degeneracy");
            pol = new double[data["field"].Length];
            eq = new double[pol.Length];
            foreach (double weight in data["equilibrium"].Distinct().Where(value => value > 0))
            {
                var kernels = Kernels(shape, increment, weight / degeneracy[0]);
                var sticks = data["equilibrium"].Select(value => value == weight ? weight : 0).ToArray();
                var contributionPol = Convolve(sticks, kernels.Polarized);
                var contributionEq = Convolve(sticks, kernels.Equilibrium);
                for (int index = 0; index < pol.Length; index++)
                {
                    pol[index] += shape.Polarization * contributionPol[index];
                    eq[index] += contributionEq[index];
                }
            }
        }
        var result = new Dictionary<string, double[]>(data)
        {
            ["shaped_polarized"] = pol, ["shaped_equilibrium"] = eq,
            ["shaped"] = pol.Zip(eq, (polarized, equilibrium) => polarized + equilibrium).ToArray()
        };
        return result;
    }
}