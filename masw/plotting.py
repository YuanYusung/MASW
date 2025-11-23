"""Plotting helpers for MASW outputs."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_wavefield(
    times: np.ndarray,
    distances_km: np.ndarray,
    wavefield: np.ndarray,
    output_path: str,
) -> None:
    X, Y = np.meshgrid(times, distances_km)
    fig, ax = plt.subplots(figsize=(8, 6))
    pcm = ax.pcolormesh(
        X,
        Y,
        wavefield,
        cmap="coolwarm",
        rasterized=True,
        shading="auto",
        vmin=-0.8,
        vmax=0.8,
    )
    fig.colorbar(pcm, ax=ax, label="Normalized amplitude")

    ymin, ymax = ax.get_ylim()
    t_line = np.array([ymin / 3.0, ymax / 3.0])
    ax.plot(t_line, [ymin, ymax], "k--", lw=1.2)
    ax.annotate("3 km/s", xy=(1.5, ymin * 0.8 + ymax * 0.2), rotation=63, color="k", fontsize=18)

    # ax.set_xlim(times.min(), times.max())
    ax.set_xlim(-2, 8)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("Correlation time (s)", fontsize=18)
    ax.set_ylabel("Interstation distance (km)", fontsize=18)
    ax.tick_params(labelsize=12)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_dispersion(
    frequencies: np.ndarray,
    velocities: np.ndarray,
    energy: np.ndarray,
    output_path: str,
) -> None:
    X, Y = np.meshgrid(frequencies, velocities)
    fig, ax = plt.subplots(figsize=(8, 6))
    pcm = ax.pcolormesh(X, Y, energy, cmap="coolwarm", shading="auto", rasterized=True)
    fig.colorbar(pcm, ax=ax, label="Normalized energy")

    ax.set_ylim(velocities.min(), velocities.max())
    ax.set_xlim(freqencies.min(), frequencies.max())
    ax.set_xlabel("Frequency (Hz)", fontsize=18)
    ax.set_ylabel("Phase Velocity (km/s)", fontsize=18)
    ax.tick_params(labelsize=12)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
