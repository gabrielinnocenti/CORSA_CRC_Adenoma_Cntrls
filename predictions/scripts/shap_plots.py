import argparse
import os
import sys
import glob

import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union


import joblib
import matplotlib.pyplot as plt


import numpy as np
import pandas as pd

from matplotlib.colors import LinearSegmentedColormap
from dataHandler import Config, MetadataHandler, OligosHandler, FeatureManager
from helpers import plot_shap_values, generate_feature_importance_table





# =====================
# SHAP helpers
# =====================
def _list_files(file_dir: Union[str, Path], file_pattern: str) -> List[Path]:
    file_dir = Path(file_dir)
    files = sorted(file_dir.glob(file_pattern))
    if not files:
        raise FileNotFoundError(
            f"No files matched pattern '{file_pattern}' in '{file_dir}'."
        )
    return files


def _load_shap_df(fn: Path, shap_key: str) -> pd.DataFrame:
    obj = joblib.load(fn)
    if shap_key not in obj:
        raise KeyError(f"Key '{shap_key}' not found in joblib file: {fn}")
    df = obj[shap_key]
    # if not isinstance(df, pd.DataFrame):
    #    raise TypeError(f"Expected '{shap_key}' to be a pandas DataFrame in {fn}, got {type(df)}")
    return df


def mean_shap_across_files(files, shap_key: str) -> pd.DataFrame:
    ref_index = None
    ref_cols = None
    acc = None
    for fn in files:
        df = _load_shap_df(fn, shap_key)  # returns a DataFrame

        if ref_index is None:
            ref_index = df.index
            ref_cols = df.columns
            acc = np.zeros(df.shape, dtype=np.float32)
        else:
            # alignment checks (fast + safe)
            if not df.index.equals(ref_index) or not df.columns.equals(ref_cols):
                raise ValueError(f"SHAP matrices not aligned. Offending file: {fn}")

        # accumulate (convert each df once)
        acc += df.to_numpy(dtype=np.float32, copy=False)

    acc /= len(files)

    return pd.DataFrame(acc, index=ref_index, columns=ref_cols)


def run_shap_summary_and_feature_table(
    config_file,
    file_dir,
    file_pattern,
    output_dir,
    output_name="shap_values",
    cmap="viridis",
    max_display=30,
    figure_size=(6, 5),
    shap_key=None,
    shap_fontsize=None,
    legend_labels=None,
    label_groups=None,
):
    files = _list_files(file_dir, file_pattern)
    if shap_key is None:
        shap_key = (
            "test_shap_values"
            if "validation" in file_pattern.lower()
            else "train_shap_values"
        )
    shap_values = mean_shap_across_files(files, shap_key)

    # Load data
    config = Config(config_file)
    metadata_handler = MetadataHandler(config)
    oligos_handler = OligosHandler(config)
    feature_manager = FeatureManager(
        config,
        metadata_handler,
        oligos_handler,
        subgroup="all",
        with_oligos=True,
        with_additional_features=False,
        filter_by_entropy=False,
        prevalence_threshold_min=0,
        prevalence_threshold_max=100,
    )
    X_train, y_train = feature_manager.get_features_target()

    # Align X/y to SHAP df (rows and columns)
    # missing_rows = shap_values.index.difference(X_train.index)
    # missing_cols = shap_values.columns.difference(X_train.columns)
    # if len(missing_rows) or len(missing_cols):
    #     raise KeyError(
    #         "SHAP df contains samples/features not found in X_train.\n"
    #         f"Missing rows in X_train: {list(missing_rows[:10])}{' ...' if len(missing_rows) > 10 else ''}\n"
    #         f"Missing cols in X_train: {list(missing_cols[:10])}{' ...' if len(missing_cols) > 10 else ''}"
    #     )

    X_train = X_train.loc[shap_values.index, shap_values.columns]
    y_train = y_train.loc[shap_values.index]

    group_tests = label_groups if label_groups is not None else config.group_tests

    plot_shap_values(
        shap_values.values,
        X_train,
        cmap=cmap,
        max_display=max_display,
        group_tests=group_tests,
        filename_label=output_name,
        fontsize=shap_fontsize,
        figure_size=figure_size,
        legend_labels=legend_labels,
        save_fig=True,
        figures_dir=output_dir,
    )

    oligos_metadata = oligos_handler.get_oligos_metadata_df()
    keep_cols = ["Description", "species", "genus", "family", "order", "pos", "len_seq"]
    # oligos_metadata.set_index(oligos_metadata.columns[0], inplace=True)
    missing = [c for c in keep_cols if c not in oligos_metadata.columns]
    if missing:
        raise KeyError(f"Oligos metadata missing expected columns: {missing}")
    oligos_metadata = oligos_metadata[keep_cols]

    generate_feature_importance_table(
        shap_values.values,
        X_train,
        y_train,
        oligos_metadata,
        group_tests=group_tests,
        filename_label=output_name,
        figures_dir=output_dir,
    )



def parse_cmap_arg(cmap: List[str]):
    """
    If a single str is provided, treat it as a matplotlib cmap name.
    If multiple color str are provided, treat them as colors to build a LinearSegmentedColormap.
    """
    if len(cmap) == 1:
        return cmap[0]
    return LinearSegmentedColormap.from_list("", cmap)


def parse_args_shap() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate SHAP summary plot and feature table"
    )

    parser.add_argument("--config_file", type=str, required=True)
    parser.add_argument("--file_dir", type=str, required=True)
    parser.add_argument("--file_pattern", type=str, required=True)
    parser.add_argument("--output_name", type=str, default=True)
    parser.add_argument("--output_dir", type=str, required=True)

    parser.add_argument("--max_display", type=int, default=30)
    parser.add_argument(
        "--figure_size",
        nargs=2,
        type=float,
        default=[6, 5],
        metavar=("WIDTH", "HEIGHT"),
    )

    # Either pass a matplotlib cmap name (e.g., viridis) or a list of hex colors
    parser.add_argument(
        "--cmap_colors",
        nargs="+",
        default=["#6699CC", "#CC6677"],
        help="Either a single matplotlib colormap name (e.g. 'viridis') or a list of colors (e.g. '#6699CC #CC6677').",
    )

    parser.add_argument(
        "--shap_key",
        type=str,
        default=None,
        help="Optional override for SHAP key: e.g. 'train_shap_values' or 'test_shap_values'.",
    )

    # Font sizes for plot_shap_values (kept compatible with your dict)
    parser.add_argument("--xlabel_fontsize", type=int, default=13)
    parser.add_argument("--xticks_fontsize", type=int, default=12)
    parser.add_argument("--yticks_fontsize", type=int, default=12)
    parser.add_argument("--colorbar_fontsize", type=int, default=12)
    parser.add_argument("--legend_fontsize", type=int, default=10)

    parser.add_argument("--legend_labels", nargs="+", default=["0", "1"])
    parser.add_argument("--label_groups", nargs="+", default=None)

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args_shap()

    font_sizes = {
        "xlabel": args.xlabel_fontsize,
        "xticks": args.xticks_fontsize,
        "yticks": args.yticks_fontsize,
        "legend": args.legend_fontsize,
        "colorbar": args.colorbar_fontsize,
    }
    cmap = parse_cmap_arg(args.cmap_colors)

    run_shap_summary_and_feature_table(
        config_file=args.config_file,
        file_dir=args.file_dir,
        file_pattern=args.file_pattern,
        output_dir=args.output_dir,
        output_name=args.output_name,
        cmap=cmap,
        max_display=args.max_display,
        figure_size=tuple(args.figure_size),
        shap_key=args.shap_key,
        shap_fontsize=font_sizes,
        label_groups=args.label_groups,
        legend_labels=args.legend_labels,
    )
