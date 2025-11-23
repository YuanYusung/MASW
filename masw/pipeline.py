"""High-level pipeline utilities for MASW."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from .config import AveragingConfig, MASWConfig, PathConfig, PreprocessConfig
from .io_utils import read_waveforms
from .plotting import plot_dispersion, plot_wavefield
from .processing import compute_average_dispersion

logger = logging.getLogger(__name__)


def _load_dispersion_cache(cache_path: str) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    path = Path(cache_path)
    if not path.is_file():
        return None

    data = np.load(path)
    if not {"velocities", "frequencies", "energy"}.issubset(data.files):
        return None
    return data["velocities"], data["frequencies"], data["energy"]


def _save_dispersion_cache(cache_path: str, velocities: np.ndarray, frequencies: np.ndarray, energy: np.ndarray) -> None:
    Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache_path, velocities=velocities, frequencies=frequencies, energy=energy)


def prepare_wavefield(
    paths: PathConfig,
    preprocess: PreprocessConfig,
    use_cache: bool = True,
):
    stream = read_waveforms(paths, preprocess, use_cache=use_cache)

    waveforms, distances = [], []
    for trace in stream:
        distance = trace.stats.sac.dist
        normalized = trace.data / np.max(np.abs(trace.data))
        waveforms.append(normalized)
        distances.append(distance)

    reference = stream[0]
    t0 = reference.stats.sac.b
    npts = reference.stats.npts
    delta = reference.stats.delta
    times = np.linspace(t0, t0 + (npts - 1) * delta, npts)

    return np.array(waveforms), np.array(distances), times, reference.stats.sampling_rate


def run_pipeline(
    paths: PathConfig = PathConfig(),
    preprocess: PreprocessConfig = PreprocessConfig(),
    masw_config: MASWConfig = MASWConfig(),
    averaging: AveragingConfig = AveragingConfig(),
    use_cache: bool = True,
    wavefield_plot_path: str | None = "Figure1_wavefield.png",
    dispersion_plot_path: str | None = "Figure2_masw_avg_dispersion.png",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    wavefield, distances, times, sample_rate = prepare_wavefield(
        paths, preprocess, use_cache=use_cache
    )

    cached_dispersion = _load_dispersion_cache(paths.masw_cache) if use_cache else None
    if cached_dispersion:
        logger.info("Loaded dispersion from cache %s", paths.masw_cache)
        velocities, frequencies, energy = cached_dispersion
    else:
        velocities, frequencies, energy = compute_average_dispersion(
            distances,
            wavefield,
            times,
            sample_rate,
            averaging,
            masw_config,
        )
        if use_cache:
            _save_dispersion_cache(paths.masw_cache, velocities, frequencies, energy)

    if wavefield_plot_path:
        logger.info("Saving wavefield plot to %s", wavefield_plot_path)
        plot_wavefield(times, distances[np.argsort(distances)], wavefield[np.argsort(distances), :], wavefield_plot_path)

    if dispersion_plot_path:
        logger.info("Saving dispersion plot to %s", dispersion_plot_path)
        plot_dispersion(frequencies, velocities, energy, dispersion_plot_path)

    return velocities, frequencies, energy
