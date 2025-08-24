"""Config pytest: garante que diretório raiz esteja no sys.path para importar 'algorithms'."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Força auto-discovery dos algoritmos
try:
    import algorithms  # noqa: F401
except Exception:  # pragma: no cover
    pass
