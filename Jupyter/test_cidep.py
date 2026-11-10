import unittest
import json
from pathlib import Path
import struct
import tempfile
from types import SimpleNamespace

import numpy as np

from cidep import Config, Group, Radical, Shape, broaden, group_states, kernels, radical_states, simulate
from storage import import_legacy, load, save


class NumericalTests(unittest.TestCase):
    def test_committed_cross_port_fixtures(self):
        path = Path(__file__).resolve().parents[1] / "verification" / "fixtures.json"
        for fixture in json.loads(path.read_text(encoding="utf-8")):
            payload = fixture["document"]
            result = simulate(Config.from_dict(payload["config"]))
            if payload["shape"]:
                result = broaden(result, Shape(**payload["shape"]))
            for key, expected in payload["data"].items():
                np.testing.assert_allclose(result[key], expected, rtol=1e-8, atol=1e-7, err_msg=fixture["name"])

    def test_exact_triplet(self):
        self.assertEqual(group_states(Group(2, 1, 0.5)), [(-1, 0, 1), (0, 0, 2), (1, 0, 1)])

    def test_total_angular_momentum(self):
        for spin in (0, 0.5, 1, 1.5):
            for count in (1, 2, 3, 6):
                states = group_states(Group(count, 1, spin), True)
                self.assertEqual(sum(row[2] for row in states), (2 * spin + 1) ** count)

    def test_independent_spin_half_rpm(self):
        config = Config(radical1=Radical(2, [Group(1, 2)]), second_order=False,
                        normalize=False, pair="S", points=101, sweep=10)
        result = simulate(config)
        nonzero = result["rpm"][result["rpm"] != 0]
        np.testing.assert_allclose(nonzero, [1.138871796 * np.sqrt(2) / 4, -1.138871796 * np.sqrt(2) / 4])
        self.assertEqual(result["equilibrium"].sum(), 2)

    def test_second_order_shift(self):
        states = radical_states(Radical(2, [Group(1, 2)]), 100, True)
        np.testing.assert_allclose(states[:, 0], [-0.99, 1.01])

    def test_default_and_mixing(self):
        config = Config(stm_weight=0.2, tm_weight=-0.3)
        result = simulate(config)
        self.assertEqual(result["total_degeneracy"], 128)
        self.assertAlmostEqual(np.ptp(result["rpm"]), 100)
        np.testing.assert_allclose(result["mixed"], result["rpm"] + 0.2 * result["stm"] - 0.3 * result["tm"])
        self.assertTrue(all(np.all(np.isfinite(value)) for value in result.values()))

    def test_zero_coupling(self):
        result = simulate(Config(radical1=Radical(2, [])))
        self.assertTrue(np.all(result["rpm"] == 0))

    def test_lineshape_symmetry_and_time_zero(self):
        polarized, equilibrium = kernels(Shape(delay=0, samples=1, points=16), 0.1)
        np.testing.assert_allclose(polarized, polarized[::-1])
        np.testing.assert_allclose(equilibrium, 0, atol=1e-12)
        np.testing.assert_allclose(polarized, 0, atol=1e-12)

    def test_convolution(self):
        result = simulate(Config(points=128))
        for exchange in (None, 5.5):
            shaped = broaden(result, Shape(points=16, samples=3, exchange_time=exchange))
            self.assertEqual(len(shaped["shaped"]), 128)
            self.assertTrue(np.all(np.isfinite(shaped["shaped"])))

    def test_invalid_inputs(self):
        for config in (Config(points=1), Config(sweep=0), Config(exponent=-1),
                       Config(radical1=Radical(2, [Group(1, 2, 0.3)]))):
            with self.assertRaises(ValueError):
                simulate(config)


class InterfaceTests(unittest.TestCase):
    def test_controls_and_precise_cursor(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from cidep_ui import Workbench
        workbench = Workbench()
        self.assertEqual(workbench.inputs()[0], Config())
        workbench.click(SimpleNamespace(inaxes=workbench.axes, xdata=workbench.result["field"][100]))
        self.assertEqual(workbench.cursor_index, 100)
        workbench.key(SimpleNamespace(key="right"))
        self.assertEqual(workbench.cursor_index, 101)
        workbench.key(SimpleNamespace(key="left"))
        self.assertEqual(workbench.cursor_index, 100)
        plt.close(workbench.figure)


class PersistenceTests(unittest.TestCase):
    def test_roundtrip(self):
        config = Config(points=64)
        shape = Shape(points=16)
        result = broaden(simulate(config), shape)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectrum.json"
            save(path, result, config, shape)
            payload = load(path)
            np.testing.assert_array_equal(payload["data"]["shaped"], result["shaped"])
            self.assertEqual(Config.from_dict(payload["config"]), config)
            invalid = json.loads(path.read_text())
            invalid["data"]["field"][1] = invalid["data"]["field"][0]
            path.write_text(json.dumps(invalid))
            with self.assertRaises(ValueError):
                load(path)

    def test_legacy_layouts_and_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "legacy.bin"
            header = struct.pack("<h3f", 4, 10, 3483, 9.777)
            path.write_bytes(header + struct.pack("<8f", 1, -2, 3, 0, 1, 2, 3, 4))
            payload = import_legacy(path)
            self.assertEqual(payload["data"]["mixed"], [1, -2, 3, 0])
            self.assertAlmostEqual(payload["metadata"]["frequency"], 9.777, places=5)
            path.write_bytes(header + struct.pack("<3fh4f", 0.1, 1, 0.01, 150, 1, 2, 3, 4))
            shaped = import_legacy(path)
            self.assertIn("shaped", shaped["data"])
            self.assertEqual(shaped["metadata"]["samples"], 150)
            path.write_bytes(header + b"bad")
            with self.assertRaises(ValueError):
                import_legacy(path)


if __name__ == "__main__":
    unittest.main()