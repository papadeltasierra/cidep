from dataclasses import asdict
import json
from pathlib import Path
import struct

import numpy as np

from cidep import Config, Shape


UNITS = {"field": "Gauss", "frequency": "GHz", "time": "microseconds"}
ARRAYS = {"field", "equilibrium", "rpm", "stm", "tm", "mixed", "shaped",
          "shaped_polarized", "shaped_equilibrium"}


def document(result, config=None, shape=None, source="simulation", metadata=None):
    payload = {"version": 1, "units": UNITS, "source": source,
               "metadata": metadata or {},
               "config": asdict(config) if config else None,
               "shape": asdict(shape) if shape else None,
               "data": {key: value.tolist() if isinstance(value, np.ndarray) else value
                        for key, value in result.items()}}
    validate(payload)
    return payload


def validate(payload):
    if payload.get("version") != 1 or payload.get("units") != UNITS:
        raise ValueError("Unsupported JSON version or units")
    if any(not isinstance(value, (int, float)) or not np.isfinite(value)
           for value in payload.get("metadata", {}).values()):
        raise ValueError("Metadata must contain finite numeric values")
    data = payload.get("data", {})
    if "field" not in data or "mixed" not in data or "equilibrium" not in data:
        raise ValueError("Missing spectrum arrays")
    fields = np.asarray(data["field"], dtype=float)
    if fields.ndim != 1 or not 2 <= len(fields) <= 16384:
        raise ValueError("Invalid field array")
    steps = np.diff(fields)
    if not np.all(np.isfinite(fields)) or np.any(steps <= 0) or not np.allclose(steps, steps[0], rtol=1e-8, atol=1e-9):
        raise ValueError("Field axis must be finite, increasing, and uniformly spaced")
    for key, value in data.items():
        if key in ARRAYS:
            array = np.asarray(value, dtype=float)
            if array.shape != fields.shape or not np.all(np.isfinite(array)):
                raise ValueError(f"Invalid spectrum array: {key}")
        elif key != "total_degeneracy" or not np.isfinite(value) or value <= 0:
            raise ValueError(f"Invalid spectrum metadata: {key}")
    if payload.get("config") is not None:
        config = Config.from_dict(payload["config"])
        config.validate()
        if config.points != len(fields):
            raise ValueError("Configuration point count does not match data")
    if payload.get("shape") is not None:
        Shape(**payload["shape"]).validate()


def save(path, result, config=None, shape=None, source="simulation", metadata=None):
    payload = document(result, config, shape, source, metadata)
    Path(path).write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")


def load(path):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    validate(payload)
    payload["data"] = {key: np.asarray(value, dtype=float) if key in ARRAYS else value
                       for key, value in payload["data"].items()}
    return payload


def import_legacy(path):
    raw = path if isinstance(path, bytes) else Path(path).read_bytes()
    if len(raw) < 14:
        raise ValueError("Truncated legacy header")
    count = struct.unpack_from("<h", raw)[0]
    if not 2 <= count <= 2048:
        raise ValueError("Invalid legacy sample count (expected a 16-bit integer)")
    sweep, centre, frequency = struct.unpack_from("<3f", raw, 2)
    if not all(np.isfinite(value) and value > 0 for value in (sweep, centre, frequency)):
        raise ValueError("Invalid legacy field or frequency metadata")
    metadata = {"sweep": sweep, "centre": centre, "frequency": frequency}
    if len(raw) == 14 + 8 * count:
        values = np.frombuffer(raw, dtype="<f4", count=2 * count, offset=14).astype(float)
        result = {"mixed": values[:count], "equilibrium": values[count:]}
        source = "legacy-stick-int16-real32"
    elif len(raw) == 28 + 4 * count:
        parameters = struct.unpack_from("<3fh", raw, 14)
        if not all(np.isfinite(value) for value in parameters[:3]) or parameters[3] < 1:
            raise ValueError("Invalid legacy lineshape metadata")
        metadata.update(zip(("microwave", "delay", "interval", "samples"), parameters))
        values = np.frombuffer(raw, dtype="<f4", count=count, offset=28).astype(float)
        result = {"mixed": values, "shaped": values.copy(), "equilibrium": np.zeros(count)}
        source = "legacy-shaped-int16-real32"
    else:
        raise ValueError("Unsupported legacy layout or file length; record-marked files are not supported")
    result["field"] = np.linspace(centre - sweep / 2, centre + sweep / 2, count)
    return document(result, source=source, metadata=metadata)