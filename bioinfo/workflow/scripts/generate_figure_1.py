import glob
import os
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

colorblind8: List[str] = [
    "#999999",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
]

LIST_OF_GENES: List[str] = [
    "WDR91",
    "GPAM",
    "STRIP1",
    "CSRP2BP",
    "HSPA4",
    "MEP1A",
    "MYB",
    "PDE6A",
    "SEMA3D",
    "TNFAIP3",
]

OMEGA_MIX: List[str] = ["mix1", "mix2", "mix3", "mix4", "mix5"]
CPG: List[float] = [1, 2, 4, 8]


def quantile(
    arr: List[float], quantile_inf: float = 0.025, quantile_sup: float = 0.975
) -> Tuple[float, float]:
    quantile_inf_v: float = np.quantile(a=arr, q=quantile_inf)
    quantile_sup_v: float = np.quantile(a=arr, q=quantile_sup)
    return (round(quantile_inf_v, 3), round(quantile_sup_v, 3))


def compute_likelihood_ratio_test(
    df_h0: pd.DataFrame,
    df_h1: pd.DataFrame,
    list_of_genes: List[str],
    list_of_omega_mix: List[str],
    list_of_CpG: List[float],
) -> pd.DataFrame:
    # 10.1093/oxfordjournals.molbev.a025957
    # 2delta_lnL = 2(log(h1) - log(h0))
    def lrt(h0: pd.DataFrame, h1: pd.DataFrame):
        inter = set(h0["fastaID"].to_list()).intersection(set(h1["fastaID"].to_list()))
        h0_ = h0.loc[h0["fastaID"].isin(inter)]
        h1_ = h1.loc[h1["fastaID"].isin(inter)]
        # M7 versus M8 : 36-34 = 2 degrees of freedom
        # h1_.sort_values(by="fastaID")["fastaID"].to_csv(
        #         self.get_dir(type="local", dir="tables")  + "h1.csv", sep="\t"
        #     )
        # h0_.sort_values(by="fastaID")["fastaID"].to_csv(
        #         self.get_dir(type="local", dir="tables")  + "h0.csv", sep="\t"
        #     )
        return 2 * (
            h1_.sort_values(by="fastaID")["lnL"].to_numpy()
            - h0_.sort_values(by="fastaID")["lnL"].to_numpy()
        )

    dict_of_lrt: Dict[int, Dict[str, Union[str, float, int, Tuple[float, float]]]] = {}
    k = 0

    for g in list_of_genes:
        for m in list_of_omega_mix:
            for CpG in list_of_CpG:
                arr_of_lrt: np.array = lrt(
                    h0=df_h0.loc[
                        (df_h0["mix"] == m)
                        & (df_h0["CpG"] == CpG)
                        & (df_h0["geneID"] == g)
                    ],
                    h1=df_h1.loc[
                        (df_h1["mix"] == m)
                        & (df_h1["CpG"] == CpG)
                        & (df_h1["geneID"] == g)
                    ],
                )
                dict_of_lrt[k] = {
                    "geneID": g,
                    "mix": m,
                    "CpG": str(CpG),
                    # "lrt": arr_of_lrt,
                    "mean": np.mean(arr_of_lrt),
                    "quantile": quantile(arr_of_lrt),
                    "alpha 0.10": round(
                        np.where(arr_of_lrt > 4.605)[0].shape[0]
                        / arr_of_lrt.shape[0]
                        * 100
                        if arr_of_lrt.shape[0] > 0
                        else 0,
                        1,
                    ),
                    "alpha 0.05": round(
                        np.where(arr_of_lrt > 5.991)[0].shape[0]
                        / arr_of_lrt.shape[0]
                        * 100
                        if arr_of_lrt.shape[0] > 0
                        else 0,
                        1,
                    ),
                    "alpha 0.01": round(
                        np.where(arr_of_lrt > 9.210)[0].shape[0]
                        / arr_of_lrt.shape[0]
                        * 100
                        if arr_of_lrt.shape[0] > 0
                        else 0,
                        1,
                    ),
                    "count": arr_of_lrt.shape[0],
                }
                k += 1
    return pd.DataFrame.from_dict(data=dict_of_lrt, orient="index")


def extract_M0_fixed_HKY(
    file: str, metadata: Dict[str, Union[str, float]]
) -> Dict[str, Union[str, float]]:
    dict_of_params: Dict[str, Union[str, float]] = {}

    for k, v in metadata.items():
        dict_of_params[k] = v

    try:
        with open(file, "r") as fh:
            lines = iter(fh.readlines())
            for l in lines:
                if l.startswith("lnL"):
                    dict_of_params["lnL"] = float(l.split()[-2])
                if l.startswith("tree"):
                    dict_of_params["tree"] = float(l.split()[-1])
                if l.startswith("kappa"):
                    dict_of_params["kappa"] = float(l.split()[-1])
                if l.startswith("omega"):
                    dict_of_params["omega"] = float(l.split()[-1])
    except Exception as e:
        print("something wrong with %s, %s" % (file, str(e)))
        raise RuntimeError
    return dict_of_params


def extract_M7M8HKY(
    file: str, metadata: Dict[str, Union[str, float]], M7_bool: bool = True
) -> Dict[str, Union[str, float]]:
    dict_of_params: Dict[str, Union[str, float]] = {}

    for k, v in metadata.items():
        dict_of_params[k] = v
    with open(file, "r") as fh:
        lines = iter(fh.readlines())
        for l in lines:
            if l.startswith("lnL"):
                dict_of_params["lnL"] = float(l.split()[-2])
            if l.startswith("tree"):
                dict_of_params["tree"] = float(l.split()[-1])
            if l.startswith("kappa"):
                dict_of_params["kappa"] = float(l.split()[-1])

            if l.startswith("Parameters in M"):
                l = next(lines)
                tup = l.split()
                dict_of_params["p0"] = 1.0
                dict_of_params["p"] = float(tup[2])
                dict_of_params["q"] = float(tup[5])

            if l.startswith("p:"):
                tup = l.split()
                dict_of_params["p_k1"] = float(tup[1])
                dict_of_params["p_k2"] = float(tup[2])
                dict_of_params["p_k3"] = float(tup[3])
                dict_of_params["p_k4"] = float(tup[4])
                dict_of_params["p_k5"] = float(tup[5])
                dict_of_params["p_k6"] = float(tup[6])
                dict_of_params["p_k7"] = float(tup[7])
                dict_of_params["p_k8"] = float(tup[8])
                dict_of_params["p_k9"] = float(tup[9])
                dict_of_params["p_k10"] = float(tup[10])
                if not M7_bool:
                    dict_of_params["p1"] = float(tup[11])

            if l.startswith("w:"):
                tup = l.split()
                dict_of_params["w_k1"] = float(tup[1])
                dict_of_params["w_k2"] = float(tup[2])
                dict_of_params["w_k3"] = float(tup[3])
                dict_of_params["w_k4"] = float(tup[4])
                dict_of_params["w_k5"] = float(tup[5])
                dict_of_params["w_k6"] = float(tup[6])
                dict_of_params["w_k7"] = float(tup[7])
                dict_of_params["w_k8"] = float(tup[8])
                dict_of_params["w_k9"] = float(tup[9])
                dict_of_params["w_k10"] = float(tup[10])
                if not M7_bool:
                    dict_of_params["w_p1"] = float(tup[11])

    return dict_of_params


def recover_params_M7M8HKY(
    files: List[str],
    metadata: Dict[str, Union[str, float]],
    M7_bool: bool = True,
) -> pd.DataFrame:

    dict_run: Dict[int, Dict[str, Union[str, float]]] = {}
    k = 0
    for f in files:
        dict_run[k] = extract_M7M8HKY(file=f, metadata=metadata, M7_bool=M7_bool)
        k += 1
    return pd.DataFrame.from_dict(data=dict_run, orient="index").dropna()


def recover_params_M0HKY(
    files: List[str], metadata: Dict[str, Union[str, float]]
) -> pd.DataFrame:

    dict_run: Dict[int, Dict[str, Union[str, float]]] = {}
    k = 0
    for f in files:
        dict_run[k] = extract_M0_fixed_HKY(file=f, metadata=metadata)
        k += 1
    return pd.DataFrame.from_dict(data=dict_run, orient="index").dropna()


def wrapper_recover_params_M0HKY(
    input: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    df_of_params: pd.DataFrame = pd.DataFrame()
    try:
        files: List[str] = glob.glob(input)
        print("number of files recovered from M0HKY: %d" % len(files))
        df_of_params = recover_params_M0HKY(files=files, metadata=metadata)
        return df_of_params
    except Exception as e:
        print("something wrong when recovering M0HKY %s" % (str(e)))
        return None


def wrapper_recover_params_M7HKY(
    input: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    df_of_params: pd.DataFrame = pd.DataFrame()
    try:
        files: List[str] = glob.glob(input)
        print(" recovered from M7HKY: %d" % len(files))
        df_of_params = recover_params_M7M8HKY(
            files=files, metadata=metadata, M7_bool=True
        )
        return df_of_params
    except Exception as e:
        print("something wrong when recovering M7HKY %s" % (str(e)))
        return None


def wrapper_recover_params_M8HKY(
    input: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    df_of_params: pd.DataFrame = pd.DataFrame()
    try:
        files: List[str] = glob.glob(input)
        print(" recovered from M8HKY: %d" % len(files))
        df_of_params = recover_params_M7M8HKY(
            files=files, metadata=metadata, M7_bool=False
        )
        return df_of_params
    except Exception as e:
        print("something wrong when recovering M8HKY %s" % (str(e)))
        return None


def generate_figure_1(
    input_M0HKY_dir: str,
    metadata_M0HKY: Dict[str, Union[str, float]],
    input_M7HKY_dir: str,
    metadata_M7HKY: Dict[str, Union[str, float]],
    input_M8HKY_dir: str,
    metadata_M8HKY: Dict[str, Union[str, float]],
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

    if os.path.exists(input_M0HKY_dir + "/df_of_params_M0HKY.pickle"):
        df_of_params_M0HKY = pd.read_pickle(
            input_M0HKY_dir + "/df_of_params_M0HKY.pickle"
        )
    else:
        df_of_params_M0HKY = wrapper_recover_params_M0HKY(
            input=input_M0HKY_dir, metadata=metadata_M0HKY
        )
        df_of_params_M0HKY.to_pickle(input_M0HKY_dir + "/df_of_params_M0HKY.pickle")

    if df_of_params_M0HKY is None:
        print(
            "something wrong with recovering codeml params %s"
            % wrapper_recover_params_M0HKY.__name__
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

    if os.path.exists(input_M7HKY_dir + "/df_of_params_M7HKY.pickle"):
        df_of_params_M7HKY = pd.read_pickle(
            input_M7HKY_dir + "/df_of_params_M7HKY.pickle"
        )
    else:
        df_of_params_M7HKY = wrapper_recover_params_M7HKY(
            input=input_M7HKY_dir, metadata=metadata_M7HKY
        )
        df_of_params_M7HKY.to_pickle(input_M7HKY_dir + "/df_of_params_M7HKY.pickle")

    if df_of_params_M7HKY is None:
        print(
            "something wrong with recovering codeml params %s"
            % wrapper_recover_params_M7HKY.__name__
        )
        return False

    if os.path.exists(input_M8HKY_dir + "/df_of_params_M8HKY.pickle"):
        df_of_params_M8HKY = pd.read_pickle(
            input_M8HKY_dir + "/df_of_params_M8HKY.pickle"
        )
    else:
        df_of_params_M8HKY = wrapper_recover_params_M8HKY(
            input=input_M8HKY_dir, metadata=metadata_M8HKY
        )
        df_of_params_M8HKY.to_pickle(input_M8HKY_dir + "/df_of_params_M8HKY.pickle")

    if df_of_params_M8HKY is None:
        print(
            "something wrong with recovering codeml params %s"
            % wrapper_recover_params_M8HKY.__name__
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
