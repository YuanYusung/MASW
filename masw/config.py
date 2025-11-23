"""Configuration models for the MASW pipeline."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import numpy as np


@dataclass
class PathConfig:
    data_dir: str = "Data_CCFs"
    waveform_cache: str = "All_sacdata_ZZ.pickle"
    masw_cache: str = "MASW_raw_data.npz"


@dataclass
class PreprocessConfig:
    bandpass: Tuple[float, float] = (0.3, 3.0)
    taper_velocity_bounds: Tuple[float, float] = (0.2, 5.0)


@dataclass
class MASWConfig:
    velocity_grid: np.ndarray = field(
        default_factory=lambda: np.linspace(0.2, 2.5, 231)
    )
    freq_range: Tuple[float, float] = (0.2, 5.0)
    zero_padding_factor: int = 5


@dataclass
class AveragingConfig:
    velocity_grid: np.ndarray = field(
        default_factory=lambda: np.linspace(2.0, 4.5, 231)
    )
    freq_range: Tuple[float, float] = (0.3, 3.0)
    min_distance_km: float = 0.6
