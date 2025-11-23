"""MASW processing package.

Lazy imports keep the CLI help usable even when optional heavy
numerical dependencies are missing. The full functionality still
requires installing the dependencies listed in ``requirements.txt``.
"""
from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

__all__ = [
    "PathConfig",
    "PreprocessConfig",
    "MASWConfig",
    "AveragingConfig",
    "read_waveforms",
    "compute_masw",
    "compute_average_dispersion",
    "plot_wavefield",
    "plot_dispersion",
    "prepare_wavefield",
    "run_pipeline",
]

if TYPE_CHECKING:  # pragma: no cover - for static analysis only
    from .config import PathConfig, PreprocessConfig, MASWConfig, AveragingConfig
    from .io_utils import read_waveforms
    from .processing import compute_masw, compute_average_dispersion
    from .plotting import plot_wavefield, plot_dispersion
    from .pipeline import prepare_wavefield, run_pipeline


def __getattr__(name: str):
    if name in {"PathConfig", "PreprocessConfig", "MASWConfig", "AveragingConfig"}:
        module = import_module(".config", __name__)
        return getattr(module, name)
    if name in {"read_waveforms"}:
        return getattr(import_module(".io_utils", __name__), name)
    if name in {"compute_masw", "compute_average_dispersion"}:
        return getattr(import_module(".processing", __name__), name)
    if name in {"plot_wavefield", "plot_dispersion"}:
        return getattr(import_module(".plotting", __name__), name)
    if name in {"prepare_wavefield", "run_pipeline"}:
        return getattr(import_module(".pipeline", __name__), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
