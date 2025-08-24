"""JSON serialization utilities for work state persistence."""

import json
import time
from dataclasses import asdict, is_dataclass
from typing import Any


def serialize_to_json(obj: Any) -> str:
    """Serialize object to JSON string with recursive handling of dataclasses and objects."""
    def recursive_serialize(o: Any) -> Any:
        if o is None:
            return {}
        if is_dataclass(o):
            return {k: recursive_serialize(v) for k, v in asdict(o).items()}
        if hasattr(o, "to_dict"):
            return recursive_serialize(o.to_dict())
        if isinstance(o, dict):
            return {k: recursive_serialize(v) for k, v in o.items()}
        if isinstance(o, (list, tuple, set)):
            return [recursive_serialize(v) for v in o]
        if isinstance(o, (str, int, float, bool)):
            return o
        # Fallback: convert to string

        try:
            return str(o)
        except Exception:
            return "<unserializable>"

    result = None
    try:
        result = json.dumps(recursive_serialize(obj), ensure_ascii=False)
    except Exception:
        # Em último recurso, transforma em str sem deixar gerar exceção
        result = json.dumps(str(obj), ensure_ascii=False)
    return result


def add_timestamp_if_missing(data: dict[str, Any]) -> dict[str, Any]:
    """Add timestamp to data if not present."""
    if "timestamp" not in data:
        data["timestamp"] = time.time()
    return data


def extract_object_metadata(obj: Any) -> dict[str, Any]:
    """Extract metadata from object for persistence."""
    meta = {}

    # Try to capture simple fields
    for attr in [
        "n",
        "L",
        "alphabet",
        "noise",
        "filename",
        "query",
        "db",
        "retmax",
        "type",
    ]:
        if hasattr(obj, attr):
            value = getattr(obj, attr)
            if value is not None:
                meta[attr] = value

    return meta
