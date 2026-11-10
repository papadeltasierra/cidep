import argparse
import json
from pathlib import Path

import numpy as np

from cidep import Config, Group, Radical, Shape, broaden, simulate
from storage import document, load


def write_fixtures():
    cases = [
        ("default", Config(points=64), None),
        ("mixed", Config(points=64, stm_weight=0.25, tm_weight=-0.4), None),
        ("distinct", Config(points=64, identical=False, normalize=False, pair="S"), None),
        ("counter_radical", Config(points=64, identical=False, observed_radical=2, exponent=0.7), None),
        ("spin_one", Config(points=64, radical1=Radical(2.003, [Group(2, 15, 1)]), net_effect=False), None),
        ("relaxation", Config(points=64), Shape(points=32, samples=4)),
        ("exchange", Config(points=64), Shape(points=32, samples=4, exchange_time=5.5)),
    ]
    fixtures = []
    for name, config, shape in cases:
        result = simulate(config)
        if shape:
            result = broaden(result, shape)
        fixtures.append({"name": name, "document": document(result, config, shape)})
    path = Path(__file__).resolve().parents[1] / "verification" / "fixtures.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(fixtures, indent=2, allow_nan=False), encoding="utf-8")
    print(f"Generated {len(fixtures)} cross-port fixtures: {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-csharp", type=Path)
    arguments = parser.parse_args()
    if arguments.check_csharp:
        payload = load(arguments.check_csharp)
        config = Config.from_dict(payload["config"])
        result = simulate(config)
        if payload["shape"]:
            result = broaden(result, Shape(**payload["shape"]))
        for key, expected in payload["data"].items():
            np.testing.assert_allclose(result[key], expected, rtol=1e-8, atol=1e-7, err_msg=key)
        print("PASS: Python read and reproduced the C# JSON export")
    else:
        write_fixtures()