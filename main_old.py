"""
Compat layer for legacy tests expecting `main_old.app`.

Re-exports the Typer `app` from the current main.py.
"""

from main import app  # noqa: F401

if __name__ == "__main__":
    # Optional: allow running directly, same as `python main.py`
    import sys
    from typer import Typer

    if isinstance(app, Typer):  # safety check
        app()
    else:
        # Fallback: nothing to execute
        sys.exit(0)
