"""Entry point for MASW processing.

This module exposes a CLI and reusable pipeline helpers to prepare data,
compute MASW dispersion, and generate plots. See README for usage examples.
"""
from __future__ import annotations

from masw.cli import main


if __name__ == "__main__":
    main()
