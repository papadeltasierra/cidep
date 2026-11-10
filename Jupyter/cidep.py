from dataclasses import asdict, dataclass, field
from concurrent.futures import ThreadPoolExecutor
import itertools
import math

import numpy as np
from scipy.signal import convolve


@dataclass
class Group:
    count: int = 1
    coupling: float = 21.94
    spin: float = 0.5


@dataclass
class Radical:
    g: float = 2.0025
    groups: list[Group] = field(default_factory=lambda: [Group(), Group(6, 24.69)])


@dataclass
class Config:
    radical1: Radical = field(default_factory=Radical)
    radical2: Radical = field(default_factory=lambda: Radical(2.0007, [Group(1, 1.0)]))
    identical: bool = True
    second_order: bool = True
    centre_resonance: bool = True
    centre: float = 3483.0
    sweep: float = 180.0
    frequency: float = 9.777
    points: int = 1024
    pair: str = "T"
    linear: float = 0.0
    exponent: float | None = None
    observed_radical: int = 1
    net_effect: bool = True
    normalize: bool = True
    rpm_weight: float = 1.0
    stm_weight: float = 0.0
    tm_weight: float = 0.0
    tm_ratio: float = -50.0

    @classmethod
    def from_dict(cls, data):
        data = dict(data)
        for key in ("radical1", "radical2"):
            if key in data:
                radical = dict(data[key])
                radical["groups"] = [Group(**group) for group in radical["groups"]]
                data[key] = Radical(**radical)
        return cls(**data)

    def validate(self):
        if not isinstance(self.points, int) or not 2 <= self.points <= 16384:
            raise ValueError("Points must be an integer from 2 to 16384")
        if self.pair not in ("S", "T", "F") or self.observed_radical not in (1, 2):
            raise ValueError("Invalid pair type or observed radical")
        for name in ("centre", "sweep", "frequency"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("linear", "rpm_weight", "stm_weight", "tm_weight", "tm_ratio"):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        if self.exponent is not None and (not math.isfinite(self.exponent) or self.exponent <= 0):
            raise ValueError("Q exponent must be finite and positive")
        for radical in (self.radical1, self.radical2):
            if not math.isfinite(radical.g) or radical.g <= 0 or len(radical.groups) > 20:
                raise ValueError("Invalid radical g-value or group count")
            for group in radical.groups:
                if (not isinstance(group.count, int) or not 1 <= group.count <= 15
                        or not math.isfinite(group.spin) or not 0 <= group.spin <= 5
                        or abs(2 * group.spin - round(2 * group.spin)) > 1e-10
                        or not math.isfinite(group.coupling)):
                    raise ValueError("Invalid nuclear group")


def group_states(group, second_order=False):
    coefficients = np.array([1.0])
    for _ in range(group.count):
        coefficients = np.convolve(coefficients, np.ones(round(2 * group.spin) + 1))
    maximum = group.count * group.spin
    if not second_order:
        return [(index - maximum, 0.0, value) for index, value in enumerate(coefficients)]
    states = []
    previous = 0.0
    for index in range(math.floor(maximum) + 1):
        angular = maximum - index
        multiplicity = coefficients[index] - previous
        previous = coefficients[index]
        if multiplicity:
            for projection in np.arange(-angular, angular + 0.5, 1):
                states.append((float(projection), angular, float(multiplicity)))
    return states


def radical_states(radical, centre, second_order=False, stm=False):
    groups = [group_states(group, second_order or stm) for group in radical.groups]
    count = math.prod(len(states) for states in groups)
    if count > 4096:
        raise ValueError("More than 4096 nuclear states; reduce the number of groups")
    rows = []
    for combination in itertools.product(*groups):
        energy, degeneracy, multiplet = 0.0, 1.0, 0.0
        for group, (projection, angular, weight) in zip(radical.groups, combination):
            energy += group.coupling * projection
            if second_order or stm:
                energy += group.coupling ** 2 * (angular * (angular + 1) - projection ** 2) / (2 * centre)
            degeneracy *= weight
            if stm:
                multiplet += (angular - projection + 1) * (angular + projection) * group.coupling ** 2 * weight
        rows.append((energy, degeneracy, multiplet))
    return np.array(rows)


def normalize(values):
    extent = float(np.ptp(values))
    return values * (100 / extent) if extent > 0 else values.copy()


def sample(fields, positions, amplitudes):
    increment = fields[1] - fields[0]
    result = np.zeros(len(fields))
    for position, amplitude in zip(positions, amplitudes):
        coordinate = (position - fields[0]) / increment
        lower = math.floor(coordinate)
        for index in (lower, lower + 1):
            if 0 <= index < len(fields) and abs(fields[index] - position) <= increment * (0.5 + 1e-10):
                result[index] += amplitude
    return result


def simulate(config):
    config.validate()
    radicals = [config.radical1, config.radical1 if config.identical else config.radical2]
    centre = 714.4844244 * config.frequency / radicals[0].g if config.centre_resonance else config.centre
    fields = np.linspace(centre - config.sweep / 2, centre + config.sweep / 2, config.points)
    normal = [radical_states(radical, config.centre, config.second_order and any(
        abs(group.coupling) > 10 for group in radical.groups)) for radical in radicals]
    stm = [radical_states(radical, config.centre, stm=True) for radical in radicals]
    own = config.observed_radical - 1
    other = 1 - own
    states, partner = normal[own], normal[other]
    difference = config.centre * (radicals[own].g - radicals[other].g) / 2.00232
    polarization = np.zeros(len(states))
    for index, energy in enumerate(states[:, 0]):
        delta = energy + difference - partner[:, 0]
        absolute = np.abs(delta)
        dependence = (absolute ** config.exponent if config.exponent is not None
                      else np.sqrt(absolute) - config.linear * absolute)
        polarization[index] = np.dot(np.sign(delta) * dependence, partner[:, 1])
    scale = {"S": 1.138871796, "T": -1.138871796, "F": -0.281368326}[config.pair]
    polarization *= scale * states[:, 1] / (states[:, 1].sum() * partner[:, 1].sum())
    resonance = 714.4844244 * config.frequency / radicals[own].g
    equilibrium = sample(fields, resonance - states[:, 0], states[:, 1])
    rpm = sample(fields, resonance - states[:, 0], polarization)
    own_stm, other_stm = stm[own], stm[other]
    net = other_stm[:, 2].sum() / other_stm[:, 1].sum() if config.net_effect else 0.0
    enhancement = -(own_stm[:, 2] + net * own_stm[:, 1]) / own_stm[:, 1].sum()
    stminus = sample(fields, resonance - own_stm[:, 0], enhancement)
    tm = config.tm_ratio * equilibrium
    components = {"rpm": rpm, "stm": stminus, "tm": tm}
    if config.normalize:
        components = {key: normalize(value) for key, value in components.items()}
    mixed = (config.rpm_weight * components["rpm"] + config.stm_weight * components["stm"]
             + config.tm_weight * components["tm"])
    return {"field": fields, "equilibrium": equilibrium, **components, "mixed": mixed,
            "total_degeneracy": float(own_stm[:, 1].sum())}


@dataclass
class Shape:
    t1: float = 5.0
    t2: float = 0.6
    microwave: float = 0.1
    delay: float = 1.0
    interval: float = 0.01
    samples: int = 150
    points: int = 512
    exchange_time: float | None = None
    polarization: float = -10.0

    def validate(self):
        for name in ("t1", "t2", "microwave", "interval"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not math.isfinite(self.delay) or self.delay < 0 or not math.isfinite(self.polarization):
            raise ValueError("Invalid delay or polarization")
        if not isinstance(self.samples, int) or not 1 <= self.samples <= 10000:
            raise ValueError("Samples must be an integer from 1 to 10000")
        if not isinstance(self.points, int) or not 4 <= self.points <= 2048 or self.points % 2:
            raise ValueError("Lineshape points must be even, from 4 to 2048")
        if self.exchange_time is not None and (not math.isfinite(self.exchange_time) or self.exchange_time <= 0):
            raise ValueError("Exchange time must be finite and positive")


def kernels(shape, increment, fraction=None):
    shape.validate()
    rate1, rate2, microwave = 1 / shape.t1, 1 / shape.t2, shape.microwave
    def at_offset(index):
        delta2 = (index * increment * 17.60844567) ** 2
        if fraction is None:
            coefficients = [1, rate1 + 2 * rate2,
                            rate2 ** 2 + 2 * rate1 * rate2 + microwave ** 2 + delta2,
                            rate1 * rate2 ** 2 + microwave ** 2 * rate2 + delta2 * rate1]
            numerator = np.poly1d([microwave, microwave * rate2])
        else:
            exchange = 0.0 if shape.exchange_time > 100 * shape.t1 else 1 / shape.exchange_time
            transfer = (1 - fraction) * exchange
            alpha, beta = rate2 + transfer, rate1 + transfer
            gamma = rate1 + exchange - transfer
            sumrate = beta + gamma
            product = beta * gamma - transfer * (exchange - transfer)
            coefficients = [1, sumrate + 2 * alpha,
                            product + alpha ** 2 + 2 * alpha * sumrate + microwave ** 2 + delta2,
                            2 * alpha * product + sumrate * alpha ** 2 + microwave ** 2 * (alpha + gamma) + sumrate * delta2,
                            alpha ** 2 * product + microwave ** 2 * alpha * gamma + product * delta2]
            numerator = microwave * np.poly1d([1, alpha]) * np.poly1d([1, rate1 + exchange])
        roots = np.roots(coefficients).astype(complex)
        derivative = np.polyval(np.polyder(coefficients), roots)
        if np.any(np.abs(derivative) < 1e-12):
            raise ValueError("Coincident relaxation roots; adjust lineshape parameters")
        amplitudes = numerator(roots) / derivative
        steady = rate1 * amplitudes / roots
        factors = np.exp(shape.delay * roots)
        if shape.samples > 1:
            factors *= np.expm1(shape.samples * shape.interval * roots) / np.expm1(shape.interval * roots)
        return float(np.real(np.sum(amplitudes * factors))), float(np.real(np.sum(steady * factors) - shape.samples * np.sum(steady)))
    with ThreadPoolExecutor() as executor:
        pairs = np.array(list(executor.map(at_offset, range(shape.points // 2))))
    return tuple(np.concatenate((pairs[::-1, index], pairs[:, index])) for index in (0, 1))


def broaden(result, shape):
    shape.validate()
    increment = result["field"][1] - result["field"][0]
    start = shape.points // 2
    count = len(result["field"])
    if shape.exchange_time is None:
        polarized, equilibrium = kernels(shape, increment)
        pol = convolve(result["mixed"], polarized, mode="full")[start:start + count]
        eq = convolve(result["equilibrium"], equilibrium, mode="full")[start:start + count]
    else:
        pol, eq = np.zeros(count), np.zeros(count)
        for weight in np.unique(result["equilibrium"]):
            if weight <= 0:
                continue
            polarized, equilibrium = kernels(shape, increment, weight / result["total_degeneracy"])
            sticks = np.where(result["equilibrium"] == weight, weight, 0)
            pol += convolve(sticks, shape.polarization * polarized, mode="full")[start:start + count]
            eq += convolve(sticks, equilibrium, mode="full")[start:start + count]
    return {**result, "shaped_polarized": pol, "shaped_equilibrium": eq, "shaped": pol + eq}