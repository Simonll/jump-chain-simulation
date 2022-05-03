import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT_dir: str = Path(os.path.dirname(__file__)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

from scripts import colorblind8


def generate_figure_2(
    input_under: str, input_center: str, input_over: str, output: str
) -> bool:

    try:

        df_knn_under = pd.read_csv(
            input_under,
            sep="\t",
        )["lambda_CpG"]
    except Exception as e:
        print("something wrong when recovering %s, %s" % ("STRIP1", str(e)))
        return False

    try:
        df_knn_center = pd.read_csv(
            input_center,
            sep="\t",
        )["lambda_CpG"]
    except Exception as e:
        print("something wrong when recovering %s, %s" % ("GPAM", str(e)))
        return False

    try:
        df_knn_over = pd.read_csv(
            input_over,
            sep="\t",
        )["lambda_CpG"]
    except Exception as e:
        print("something wrong when recovering %s, %s" % ("WDR91", str(e)))
        return False

    df_knn_under_qt = df_knn_under[df_knn_under > np.quantile(a=df_knn_under, q=0.025)]
    df_knn_center_qt = df_knn_center[
        df_knn_center > np.quantile(a=df_knn_center, q=0.025)
    ]
    df_knn_over_qt = df_knn_under[df_knn_over > np.quantile(a=df_knn_over, q=0.025)]

    df_knn_under_weights = np.ones_like(df_knn_under) / df_knn_under.shape[0]
    df_knn_center_weights = np.ones_like(df_knn_center) / df_knn_center.shape[0]
    df_knn_over_weights = np.ones_like(df_knn_over) / df_knn_over.shape[0]

    df_knn_under_ratio = str(
        round((sum(df_knn_under_qt > 1)) / df_knn_under_qt.shape[0], 2)
    )
    df_knn_center_ratio = str(
        round((sum(df_knn_center_qt > 1)) / df_knn_center_qt.shape[0], 2)
    )
    df_knn_over_ratio = str(
        round((sum(df_knn_over_qt > 1)) / df_knn_over_qt.shape[0], 2)
    )

    bins = np.histogram(np.hstack([df_knn_under, df_knn_center, df_knn_over]), bins=15)[
        1
    ]

    try:
        fig, axes = plt.subplots(
            nrows=1,
            ncols=1,
            figsize=(6, 4),
            sharex=True,
            sharey=True,
            dpi=80,
            facecolor="w",
            edgecolor="black",
        )

        axes.hist(
            [df_knn_under],
            bins=bins,
            color=[colorblind8[5]],
            alpha=0.5,
            weights=[df_knn_under_weights],
            label=[r"$\lambda=8: $" + df_knn_under_ratio],
        )

        axes.hist(
            [df_knn_center],
            bins=bins,
            color=[colorblind8[6]],
            alpha=0.5,
            weights=[df_knn_center_weights],
            label=[r"$\lambda=8: $" + df_knn_center_ratio],
        )

        axes.hist(
            [df_knn_over],
            bins=bins,
            color=[colorblind8[7]],
            alpha=0.5,
            weights=[df_knn_over_weights],
            label=[r"$\lambda=8: $" + df_knn_over_ratio],
        )

        axes.axvline(8, color="black", linestyle="dashed", linewidth=2)

        print(
            "under mean %.2f, quantile_inf %.2f, quantile_sup %.2f"
            % (
                round(np.mean(df_knn_under), 2),
                round(np.quantile(a=df_knn_under, q=0.025), 2),
                round(np.quantile(a=df_knn_under, q=0.975), 2),
            )
        )
        print(
            "center mean %.2f, quantile_inf %.2f, quantile_sup %.2f"
            % (
                round(np.mean(df_knn_center), 2),
                round(np.quantile(a=df_knn_center, q=0.025), 2),
                round(np.quantile(a=df_knn_center, q=0.975), 2),
            )
        )
        print(
            "over mean %.2f, quantile_inf %.2f, quantile_sup %.2f"
            % (
                round(np.mean(df_knn_over), 2),
                round(np.quantile(a=df_knn_over, q=0.025), 2),
                round(np.quantile(a=df_knn_over, q=0.975), 2),
            )
        )

        axes.axvline(8, color="black", linewidth=2, linestyle="dashed")

        axes.set_xlabel(r"$\lambda$")
        axes.set_yticks([])
        axes.set_ylabel("density")
        fig.tight_layout()
        for pos in ["right", "top", "left"]:  #'bottom'
            fig.gca().spines[pos].set_visible(False)
        fig.savefig(
            output,
            dpi=300,
            transparent=False,
            bbox_inches="tight",
        )
    except Exception as e:
        print("something wrong when compiling figure %s, %s" % (output, str(e)))
        return False

    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input_under",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--input_center",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--input_over",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    generate_figure_2(
        input_under=args.input_under,
        input_center=args.input_center,
        input_over=args.input_over,
        output=args.output,
    )
