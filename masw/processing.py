"""Processing routines for MASW analysis."""
from __future__ import annotations

import logging
from typing import Tuple

import numpy as np
from obspy.core.trace import Trace

from .config import AveragingConfig, MASWConfig

logger = logging.getLogger(__name__)


def taper_trace(trace: Trace, tmin: float, tmax: float) -> Trace:
    """Apply a taper to a trace using a symmetric time window."""
    npts = trace.stats.sac.npts
    t0 = trace.stats.sac.b
    delta = trace.stats.delta
    taxis = np.linspace(t0, t0 + npts * delta - delta, npts)
    ind_taper = (taxis >= tmin) & (taxis <= tmax)

    tapered = trace.copy()
    tapered.data = trace.data[ind_taper]
    tapered.taper(max_percentage=0.05)

    data = np.zeros_like(trace.data)
    data[ind_taper] = tapered.data
    tapered.data = data
    return tapered


def compute_masw(
    distances_km: np.ndarray,
    times: np.ndarray,
    sample_rate: float,
    wavefield: np.ndarray,
    config: MASWConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute MASW dispersion energy."""
    npts = wavefield.shape[1]
    zeropad = np.zeros([wavefield.shape[0], config.zero_padding_factor * npts])
    zeropad[:, 0:npts] = wavefield

    freq = np.fft.fftfreq(zeropad.shape[1], d=1.0 / sample_rate)
    fmin, fmax = config.freq_range
    freq_mask = (freq > fmin) & (freq < fmax)
    selected_freqs = freq[freq_mask]

    ufft = np.fft.fft(zeropad, axis=-1)[:, freq_mask]
    omega = selected_freqs * 2 * np.pi

    rnorm = ufft / np.abs(ufft)
    asum = np.zeros([len(omega), len(config.velocity_grid)])
    for idx, w in enumerate(omega):
        if idx % 100 == 0:
            logger.info("Processing frequency %d / %d", idx, len(omega))
        accum = np.zeros([len(distances_km), len(config.velocity_grid)], dtype=complex)
        for station_idx, dist in enumerate(distances_km):
            accum[station_idx, :] = np.exp(1j * w * (dist / config.velocity_grid)) * rnorm[
                station_idx, idx
            ]
        asum[idx, :] = np.abs(np.sum(accum, axis=0))
        peak = np.max(asum[idx, :])
        if peak > 0:
            asum[idx, :] /= peak
        asum[idx, :] /= np.max(asum[idx, :])

    return config.velocity_grid, selected_freqs, asum


def compute_average_dispersion(
    distances_km: np.ndarray,
    wavefield: np.ndarray,
    times: np.ndarray,
    sample_rate: float,
    averaging: AveragingConfig,
    masw_config: MASWConfig,
):
    mask = distances_km > averaging.min_distance_km
    sorted_indices = np.argsort(distances_km)
    valid_indices = sorted_indices[mask[sorted_indices]]

    return compute_masw(
        distances_km[valid_indices],
        times,
        sample_rate,
        wavefield[valid_indices, :],
        MASWConfig(
            velocity_grid=averaging.velocity_grid,
            freq_range=averaging.freq_range,
            zero_padding_factor=masw_config.zero_padding_factor,
        ),
    )
