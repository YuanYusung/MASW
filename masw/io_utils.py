"""I/O helpers for MASW processing."""
from __future__ import annotations

import glob
import os
import pickle
from typing import Iterable

import numpy as np
from obspy import read
from obspy.core.stream import Stream

from .config import PathConfig, PreprocessConfig
from .processing import taper_trace


def _load_cached_waveforms(cache_path: str) -> Stream | None:
    if not os.path.isfile(cache_path):
        return None
    with open(cache_path, "rb") as handle:
        return pickle.load(handle)


def _save_waveforms(cache_path: str, stream: Stream) -> None:
    with open(cache_path, "wb") as handle:
        pickle.dump(stream, handle, protocol=pickle.HIGHEST_PROTOCOL)


def _iter_sac_files(data_dir: str) -> Iterable[str]:
    return glob.glob(os.path.join(data_dir, "*.sac"))


def read_waveforms(
    paths: PathConfig,
    preprocess: PreprocessConfig,
    use_cache: bool = True,
) -> Stream:
    """Load and preprocess SAC files.

    Parameters
    ----------
    paths:
        File paths for data and caches.
    preprocess:
        Preprocessing parameters including bandpass and taper bounds.
    use_cache:
        When True, reuse the cached waveform pickle if present.
    """
    if use_cache:
        cached = _load_cached_waveforms(paths.waveform_cache)
        if cached is not None:
            return cached

    sac_files = sorted(_iter_sac_files(paths.data_dir))
    if not sac_files:
        raise FileNotFoundError(f"No SAC files found under {paths.data_dir}")

    stream: Stream | None = None
    for filename in sac_files:
        if stream is None:
            stream = read(filename)
        else:
            stream += read(filename)

    if stream is None:
        raise RuntimeError("Unexpected empty stream after reading SAC files")

    vmin, vmax = preprocess.taper_velocity_bounds
    for i, trace in enumerate(stream):
        trace.data = (trace.data + trace.data[::-1]) / 2
        dist = trace.stats.sac.dist
        tmin = dist / vmax
        tmax = dist / vmin
        stream[i] = taper_trace(trace, tmin, tmax)

    fmin, fmax = preprocess.bandpass
    stream.filter("bandpass", freqmin=fmin, freqmax=fmax, zerophase=True)

    if use_cache:
        _save_waveforms(paths.waveform_cache, stream)

    return stream
