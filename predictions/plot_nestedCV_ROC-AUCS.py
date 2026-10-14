#!/usr/bin/env python3
"""
CLI tool to generate and preserve ROC AUC plots from nested_predictions joblib files.
Usage example:
  python plot_auc_cli.py \
    --joblib-dir /path/to/Diagnostics/Cirrhosis_HCC/HCC \
    --group1 Cirrhosis \
    --size1 50 \
    --group2 HCC \
    --size2 46 \
    --out-dir /path/to/output/plots \
    --out-base ROC_Cirrhosis_HCC_HCC
"""
import os
# Limit OpenMP/BLAS threads to avoid thread creation errors
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import joblib
import pickle
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.family'] = 'DejaVu Sans'

import matplotlib.pyplot as plt
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
from typing import Dict, List
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 1.5
#plt.rcParams['font.family'] = 'Arial'

def _load_and_interp(fn: Path, fpr_grid: np.ndarray):
    d   = joblib.load(fn)
    roc = d.get('roc_metrics_train', d.get('roc_metrics_test'))  # choose available key
    fpr = np.array(roc['fpr'])
    tpr = np.array(roc['tpr'])
    auc = roc['auc']
    interp = np.interp(fpr_grid, fpr, tpr)
    interp[0], interp[-1] = 0.0, 1.0
    return interp, auc

def plot_auc_for_dir(
    joblib_dir: str,
    group1: str,
    size1: int,
    group2: str,
    size2: int,
    colors: Dict[str,str],
    out_dir: str,
    out_base: str,
    fpr_grid: np.ndarray = np.linspace(0,1,200)
) -> Dict[str, plt.Axes]:
    """
    Scan `joblib_dir` for nested_predictions_*.joblib, group by prefix,
    and for each prefix produce & save an ROC plot PDF and serialized fig.
    Returns a dict mapping prefix -> Axes object.
    """
    joblib_dir = Path(joblib_dir)
    out_dir    = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if group1 == "HCC":
        label_loc =  "upper left" #"lower right" #"upper left"
    else:
        label_loc = "lower right"
    # collect and group files by prefix
    files = sorted(joblib_dir.glob("nested_xgboost_*.joblib"))  # or "random-forest_*.joblib"
    prefix_re = re.compile(r"^nested_xgboost_(.+?)_\d+$")
    buckets: Dict[str,List[Path]] = {}
    for fn in files:
        m = prefix_re.match(fn.stem)
        if not m:
            continue
        prefix = m.group(1)
        buckets.setdefault(prefix, []).append(fn)

    axes: Dict[str, plt.Axes] = {}
    # generate plots
    for prefix, fns in buckets.items():
        n_runs = len(fns)
        # load + interp in parallel
        with ThreadPoolExecutor() as exe:
            futures = [exe.submit(_load_and_interp, fn, fpr_grid) for fn in fns]
            results = [f.result() for f in as_completed(futures)]
        tprs, aucs = zip(*results)
        arr = np.vstack(tprs)

        # compute statistics
        mean_tpr    = arr.mean(axis=0)
        low_tpr     = np.percentile(arr, 2.5, axis=0)
        high_tpr    = np.percentile(arr, 97.5, axis=0)
        mean_auc    = np.mean(aucs)
        auc_ci_low, auc_ci_high = np.percentile(aucs, [2.5,97.5])

        # plotting
        fig, ax = plt.subplots(figsize=(6,6))
        for row in arr:
            ax.plot(fpr_grid, row, color=colors['indiv'], alpha=0.65, lw=1)
        ax.plot(fpr_grid, mean_tpr,
                color=colors['mean'], lw=2.5,
                label=f"Mean AUC={mean_auc:.2f}")
        ax.fill_between(fpr_grid, low_tpr, high_tpr,
                        color=colors['ci'], alpha=0.45,
                        label=f"95% CI [{auc_ci_low:.2f},{auc_ci_high:.2f}]")
        ax.plot([0,1],[0,1],'--', color=colors['rand'], lw=1, label='Random Classifier')

        # labels & title
        ax.set_xlim(0,1);
        ax.set_ylim(0,1)
        ax.set_xlabel("False Positive Rate", fontsize=16)
        ax.set_ylabel("True Positive Rate", fontsize=16)
        ax.tick_params(labelsize=15)

        title = (
#            f"Sex/Age Features\n"
            f"{group1} (n={size1}) vs {group2} (n={size2})"
        )
        ax.set_title(title, fontsize=16)
        ax.legend(fontsize=15, loc=label_loc)
        fig.tight_layout()

        # save PDF
        pdf_path = out_dir / f"{out_base}_{prefix}.pdf"
        fig.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
        # save SVG
        svg_path = out_dir / f"{out_base}_{prefix}.svg"
        fig.savefig(svg_path, format='svg', bbox_inches='tight')

        # serialize figure for later use
        fig_path = out_dir / f"{out_base}_{prefix}.fig.pkl"
        with open(fig_path, 'wb') as f:
            pickle.dump(fig, f)
        plt.close(fig)

        axes[prefix] = ax

    return axes

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate ROC AUC plots from nested_predictions joblib files.'
    )
    parser.add_argument('--joblib-dir',  required=True,
                        help='Directory with nested_predictions_*.joblib')
    parser.add_argument('--group1',      required=True, help='Label for group1')
    parser.add_argument('--size1',       type=int, required=True, help='Sample size of group1')
    parser.add_argument('--group2',      required=True, help='Label for group2')
    parser.add_argument('--size2',       type=int, required=True, help='Sample size of group2')
    parser.add_argument('--out-dir',     required=True, help='Directory to save PDFs and figures')
    parser.add_argument('--out-base',    required=True, help='Base filename for outputs')
    parser.add_argument('--color-indiv', default='moccasin', help='Color for individual curves')
    parser.add_argument('--color-mean',  default='#d95f02', help='Color for mean curve')
    parser.add_argument('--color-ci',    default='moccasin', help='Color for CI band')
    parser.add_argument('--color-rand',  default='gray',    help='Color for random diagonal')
    args = parser.parse_args()

    colors = {
        'indiv': args.color_indiv,
        'mean':  args.color_mean,
        'ci':    args.color_ci,
        'rand':  args.color_rand
    }

    axes = plot_auc_for_dir(
        joblib_dir=args.joblib_dir,
        group1    =args.group1,
        size1     =args.size1,
        group2    =args.group2,
        size2     =args.size2,
        colors    =colors,
        out_dir   =args.out_dir,
        out_base  =args.out_base
    )
    print("Generated plots for prefixes:", list(axes.keys()))
