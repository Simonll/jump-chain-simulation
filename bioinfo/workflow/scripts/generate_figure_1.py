import argparse
import os
import sys
from pathlib import Path
from typing import Dict
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT_dir: str = Path(os.path.dirname(__file__)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

from scripts.generate_codeml_sup import compute_likelihood_ratio_test
from scripts.generate_codeml_sup import wrapper_recover_params_M0HKY_codeml
from scripts.generate_codeml_sup import wrapper_recover_params_M7HKY_codeml
from scripts.generate_codeml_sup import wrapper_recover_params_M8HKY_codeml

from . import CPG
from . import LIST_OF_GENES
from . import OMEGA_MIX
from . import colorblind8


def generate_figure_1(
    input_dir_M0HKY_codeml: str,
    metadata_M0HKY_codeml: Dict[str, Union[str, float]],
    input_dir_M7HKY_codeml: str,
    metadata_M7HKY_codeml: Dict[str, Union[str, float]],
    input_dir_M8HKY_codeml: str,
    metadata_M8HKY_codeml: Dict[str, Union[str, float]],
    output: str,
) -> bool:
    fig, axes = plt.subplots(
        nrows=2,
        ncols=3,
        figsize=(6, 4),
        sharey=False,
        sharex=False,
        dpi=80,
        facecolor="w",
        edgecolor="black",
    )
    axes = np.ravel(axes)

    if os.path.exists(input_dir_M0HKY_codeml + "/df_of_params_M0HKY.pickle"):
        df_of_params_M0HKY = pd.read_pickle(
            input_dir_M0HKY_codeml + "/df_of_params_M0HKY.pickle"
        )
    else:
        df_of_params_M0HKY = wrapper_recover_params_M0HKY_codeml(
            input_dir=input_dir_M0HKY_codeml, metadata=metadata_M0HKY_codeml
        )
        df_of_params_M0HKY.to_pickle(
            input_dir_M0HKY_codeml + "/df_of_params_M0HKY.pickle"
        )

    if df_of_params_M0HKY is None:
        print(
            "something wrong with recovering codeml params %s"
            % wrapper_recover_params_M0HKY_codeml.__name__
        )
    else:
        try:
            k = 0
            labels = ["A", "B", "C", "D", "E", "F"]
            for g in LIST_OF_GENES:
                if g in ["WDR91"]:

                    for w in [0.2, 0.5, 0.8]:
                        print("geneID: %s, omega %f" % (g, w))
                        CpG_1 = df_of_params_M0HKY.loc[
                            (df_of_params_M0HKY["fixed"] == w)
                            & (df_of_params_M0HKY["CpG"] == 1)
                            & (df_of_params_M0HKY["geneID"] == g)
                        ]["omega"].to_numpy()
                        CpG_2 = df_of_params_M0HKY.loc[
                            (df_of_params_M0HKY["fixed"] == w)
                            & (df_of_params_M0HKY["CpG"] == 2)
                            & (df_of_params_M0HKY["geneID"] == g)
                        ]["omega"].to_numpy()
                        CpG_4 = df_of_params_M0HKY.loc[
                            (df_of_params_M0HKY["fixed"] == w)
                            & (df_of_params_M0HKY["CpG"] == 4)
                            & (df_of_params_M0HKY["geneID"] == g)
                        ]["omega"].to_numpy()
                        CpG_8 = df_of_params_M0HKY.loc[
                            (df_of_params_M0HKY["fixed"] == w)
                            & (df_of_params_M0HKY["CpG"] == 8)
                            & (df_of_params_M0HKY["geneID"] == g)
                        ]["omega"].to_numpy()

                        CpG_1_weights = np.ones_like(CpG_1) / CpG_1.shape[0]
                        # CpG_2_weights = np.ones_like(CpG_2) /  CpG_2.shape[0]
                        CpG_4_weights = np.ones_like(CpG_4) / CpG_4.shape[0]
                        CpG_8_weights = np.ones_like(CpG_8) / CpG_8.shape[0]

                        CpG_1_ratio = str(round((sum(CpG_1 > w)) / CpG_1.shape[0], 2))
                        # CpG_2_ratio = str(round((sum(CpG_2 > w))/CpG_2.shape[0],2))
                        CpG_4_ratio = str(round((sum(CpG_4 > w)) / CpG_4.shape[0], 2))
                        CpG_8_ratio = str(round((sum(CpG_8 > w)) / CpG_8.shape[0], 2))

                        bins = np.histogram(
                            np.hstack([CpG_1, CpG_2, CpG_4, CpG_8]), bins=15
                        )[1]
                        axes[k].set_title(labels[k], loc="left")
                        axes[k].hist(
                            [CpG_1],
                            bins=bins,
                            color=[colorblind8[5]],
                            alpha=0.5,
                            stacked=False,
                            weights=[CpG_1_weights],
                            label=[r"$\lambda=1:$" + CpG_1_ratio],
                        )
                        print(r"$\lambda=1:$" + CpG_1_ratio)
                        axes[k].hist(
                            [CpG_4],
                            bins=bins,
                            color=[colorblind8[6]],
                            alpha=0.5,
                            stacked=False,
                            weights=[CpG_4_weights],
                            label=[r"$\lambda=4:$" + CpG_4_ratio],
                        )
                        print(r"$\lambda=4:$" + CpG_4_ratio)
                        axes[k].hist(
                            [CpG_8],
                            bins=bins,
                            color=[colorblind8[7]],
                            alpha=0.5,
                            stacked=False,
                            weights=[CpG_8_weights],
                            label=[r"$\lambda=8:$" + CpG_8_ratio],
                        )
                        print(r"$\lambda=8:$" + CpG_8_ratio)
                        axes[k].axvline(
                            w, color="black", linestyle="dashed", linewidth=2
                        )
                        axes[k].set_xlabel(r"$\omega$")
                        axes[k].set_yticks([])

                        # if k == 0:
                        #     ticks = [0.10, 0.2, 0.3]
                        #     axes[k].set_xticks(ticks=ticks)
                        #     axes[k].set_xticklabels(labels=ticks)
                        # if k == 1:
                        #     ticks = [0.4, 0.5, 0.6]
                        #     axes[k].set_xticks(ticks=ticks)
                        #     axes[k].set_xticklabels(labels=ticks)
                        # if k == 2:
                        #     ticks = [0.7, 0.8, 0.9]
                        #     axes[k].set_xticks(ticks=ticks)
                        #     axes[k].set_xticklabels(labels=ticks)
                        k += 1
                print("%s, first row compiled" % output)

        except Exception as e:
            print(
                "something wrong when compiling figure %s, first row, %s"
                % (output, str(e))
            )
            return False

    if os.path.exists(input_dir_M7HKY_codeml + "/df_of_params_M7HKY.pickle"):
        df_of_params_M7HKY = pd.read_pickle(
            input_dir_M7HKY_codeml + "/df_of_params_M7HKY.pickle"
        )
    else:
        df_of_params_M7HKY = wrapper_recover_params_M7HKY_codeml(
            input_dir=input_dir_M7HKY_codeml, metadata=metadata_M7HKY_codeml
        )
        df_of_params_M7HKY.to_pickle(
            input_dir_M7HKY_codeml + "/df_of_params_M7HKY.pickle"
        )

    if df_of_params_M7HKY is None:
        print(
            "something wrong with recovering codeml params %s"
            % wrapper_recover_params_M7HKY_codeml.__name__
        )
        return False

    if os.path.exists(input_dir_M8HKY_codeml + "/df_of_params_M8HKY.pickle"):
        df_of_params_M8HKY = pd.read_pickle(
            input_dir_M8HKY_codeml + "/df_of_params_M8HKY.pickle"
        )
    else:
        df_of_params_M8HKY = wrapper_recover_params_M8HKY_codeml(
            input_dir=input_dir_M8HKY_codeml, metadata=metadata_M8HKY_codeml
        )
        df_of_params_M8HKY.to_pickle(
            input_dir_M8HKY_codeml + "/df_of_params_M8HKY.pickle"
        )

    if df_of_params_M8HKY is None:
        print(
            "something wrong with recovering codeml params %s"
            % wrapper_recover_params_M8HKY_codeml.__name__
        )
        return False

    try:
        df_of_lrt: pd.DataFrame = compute_likelihood_ratio_test(
            df_h0=df_of_params_M7HKY,
            df_h1=df_of_params_M8HKY,
            list_of_genes=LIST_OF_GENES,
            list_of_omega_mix=OMEGA_MIX,
            list_of_CpG=CPG,
        )
    except Exception as e:
        print(
            "something wrong when  %s, %s"
            % (compute_likelihood_ratio_test.__name__, str(e))
        )
        return False

    try:
        for g in LIST_OF_GENES:
            if g in ["WDR91", "GPAM", "STRIP1"]:
                mix_1 = df_of_lrt.loc[
                    (df_of_lrt["mix"] == "mix1") & (df_of_lrt["geneID"] == g)
                ][["alpha 0.05"]]
                mix_2 = df_of_lrt.loc[
                    (df_of_lrt["mix"] == "mix2") & (df_of_lrt["geneID"] == g)
                ][["alpha 0.05"]]
                mix_3 = df_of_lrt.loc[
                    (df_of_lrt["mix"] == "mix3") & (df_of_lrt["geneID"] == g)
                ][["alpha 0.05"]]
                mix_4 = df_of_lrt.loc[
                    (df_of_lrt["mix"] == "mix4") & (df_of_lrt["geneID"] == g)
                ][["alpha 0.05"]]
                mix_5 = df_of_lrt.loc[
                    (df_of_lrt["mix"] == "mix5") & (df_of_lrt["geneID"] == g)
                ][["alpha 0.05"]]

                axes[k].set_title(labels[k], loc="left")
                axes[k].scatter(
                    x=[1, 2, 3, 4],
                    y=mix_1,
                    label=r"mixture 1: $\omega$ = {0.1, 0.2, 0.3}",
                    marker="o",
                    color="black",
                    alpha=0.3,
                )
                axes[k].scatter(
                    x=[1, 2, 3, 4],
                    y=mix_2,
                    label=r"mixture 2: $\omega$ = {0.4, 0.5, 0.6}",
                    marker=(5, 1),
                    color="black",
                    alpha=0.3,
                )
                axes[k].scatter(
                    x=[1, 2, 3, 4],
                    y=mix_3,
                    label=r"mixture 3: $\omega$ = {0.7, 0.8, 0.9}",
                    marker=(5, 2),
                    color="black",
                    alpha=0.3,
                )
                axes[k].scatter(
                    x=[1, 2, 3, 4],
                    y=mix_4,
                    label=r"mixture 4: $\omega$ = {0.2, 0.5, 0.7}",
                    marker="^",
                    color="black",
                    alpha=0.3,
                )
                axes[k].scatter(
                    x=[1, 2, 3, 4],
                    y=mix_5,
                    label=r"mixture 5: $\omega$ = {0.5, 0.7, 0.9}",
                    marker="s",
                    color="black",
                    alpha=0.3,
                )
                if k == 3:
                    axes[k].set_ylabel("% significant LRT")
                    axes[k].set_yticks([0, 25, 50, 75, 100])
                    axes[k].set_yticklabels(labels=[0, 25, 50, 75, 100])
                else:
                    axes[k].set_yticks([0, 25, 50, 75, 100])
                    axes[k].set_yticklabels(labels=[])
                axes[k].set_xlabel(r"$\lambda$")
                axes[k].set_xticks(ticks=[1, 2, 3, 4])
                axes[k].set_xticklabels(labels=[1, 2, 4, 8])
                k += 1
            fig.subplots_adjust(top=0.99, bottom=0.01, hspace=0.5, wspace=0.25)
            # fig.tight_layout(pad=1, w_pad=1, h_pad=1.0)
            fig.savefig(
                output,
                dpi=300,
                transparent=False,
                bbox_inches="tight",
            )
    except Exception as e:
        print(
            "something wrong when compiling figure %s, second row, %s"
            % (output, str(e))
        )
        return False
    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input_M0HKY_codeml_dir",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--metadata_M0HKY_codeml_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--input_M7HKY_codeml_dir",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--metadata_M7HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--input_M8HKY_codeml_dir",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--metadata_M8HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    generate_figure_1(
        input_dir_M0HKY_codeml=args.input_dir_M0HKY_codeml,
        metadata_M0HKY_codeml=args.metadata_M0HKY_codeml_codeml,
        input_dir_M7HKY_codeml=args.input_dir_M7HKY_codeml,
        metadata_M7HKY_codeml=args.metadata_M7HKY_codeml,
        input_dir_M8HKY_codeml=args.input_dir_M8HKY_codeml,
        metadata_M8HKY_codeml=args.metadata_M8HKY_codeml,
        output=args.output,
    )
