"""Command line interface for MASW processing."""
from __future__ import annotations

import argparse
import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # Only imported for type checking to avoid eager heavy imports
    from .config import AveragingConfig, MASWConfig, PathConfig, PreprocessConfig
    from .pipeline import prepare_wavefield, run_pipeline

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="MASW processing pipeline")
    parser.add_argument("subcommand", choices=["prepare", "masw", "plot"], help="Action to run")
    parser.add_argument("--data-dir", default="Data_CCFs", help="Directory containing SAC files")
    parser.add_argument("--waveform-cache", default="All_sacdata_ZZ.pickle", help="Waveform cache path")
    parser.add_argument("--masw-cache", default="MASW_raw_data.npz", help="MASW cache path")
    parser.add_argument("--bandpass-low", type=float, default=0.3, help="Lower bandpass frequency")
    parser.add_argument("--bandpass-high", type=float, default=3.0, help="Upper bandpass frequency")
    parser.add_argument("--taper-vmin", type=float, default=0.2, help="Minimum velocity for taper window")
    parser.add_argument("--taper-vmax", type=float, default=5.0, help="Maximum velocity for taper window")
    parser.add_argument(
        "--min-distance", type=float, default=0.6, help="Minimum station spacing to include (km)"
    )
    parser.add_argument(
        "--no-cache", action="store_true", help="Disable reading/writing waveform cache"
    )
    parser.add_argument(
        "--zero-padding", type=int, default=5, help="Zero padding factor for MASW calculation"
    )
    parser.add_argument(
        "--wavefield-figure", default="Figure1_wavefield.png", help="Output path for wavefield plot"
    )
    parser.add_argument(
        "--dispersion-figure", default="Figure2_masw_avg_dispersion.png", help="Output path for dispersion plot"
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    from .config import AveragingConfig, MASWConfig, PathConfig, PreprocessConfig
    from .pipeline import prepare_wavefield, run_pipeline

    paths = PathConfig(
        data_dir=args.data_dir,
        waveform_cache=args.waveform_cache,
        masw_cache=args.masw_cache,
    )
    preprocess = PreprocessConfig(
        bandpass=(args.bandpass_low, args.bandpass_high),
        taper_velocity_bounds=(args.taper_vmin, args.taper_vmax),
    )
    masw_config = MASWConfig(zero_padding_factor=args.zero_padding)
    averaging = AveragingConfig(min_distance_km=args.min_distance)

    if args.subcommand == "prepare":
        prepare_wavefield(paths, preprocess, use_cache=not args.no_cache)
        return

    if args.subcommand == "masw":
        velocities, frequencies, energy = run_pipeline(
            paths,
            preprocess,
            masw_config,
            averaging,
            use_cache=not args.no_cache,
            wavefield_plot_path=None,
            dispersion_plot_path=None,
        )
        logging.info(
            "Computed MASW with velocity grid from %.2f to %.2f km/s",
            velocities.min(),
            velocities.max(),
        )
        return

    if args.subcommand == "plot":
        velocities, frequencies, energy = run_pipeline(
            paths,
            preprocess,
            masw_config,
            averaging,
            use_cache=not args.no_cache,
            wavefield_plot_path=args.wavefield_figure,
            dispersion_plot_path=args.dispersion_figure,
        )
        logging.info("Saved plots to %s and %s", args.wavefield_figure, args.dispersion_figure)
        return


if __name__ == "__main__":
    main()
