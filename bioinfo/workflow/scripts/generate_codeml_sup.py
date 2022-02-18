import argparse
import glob
import os
import sys
from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import numpy as np
import pandas as pd

ROOT_dir: str = Path(os.path.dirname(__file__)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

from scripts import CPG
from scripts import LIST_OF_GENES
from scripts import OMEGA_MIX
from scripts.utils import quantile


def generate_supplements_M0HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    df_of_params: Optional[pd.DataFrame] = wrapper_recover_params_M0HKY_codeml(
        input_dir=input_dir, metadata=metadata
    )
    if df_of_params is not None:
        try:
            df_of_params.groupby(by=["geneID", "CpG", "fixed"])[
                ["tree", "kappa", "omega"]
            ].agg([np.mean, quantile]).round(3).to_csv(output, sep="\t")
            print("table %s supplement materials saved" % output)
            return True
        except Exception as e:
            print("something wrong with %s" % str(e))
            return False
    else:
        return False


def generate_supplements_M7HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    df_of_params: Optional[pd.DataFrame] = wrapper_recover_params_M7HKY_codeml(
        input_dir=input_dir, metadata=metadata
    )
    if df_of_params is not None:
        try:
            df_of_params.groupby(by=["geneID", "CpG", "mix"])[
                ["tree", "kappa", "p", "q"]
                + ["p_k" + str(i) for i in range(1, 11, 1)]
                + ["w_k" + str(i) for i in range(1, 11, 1)]
            ].agg([np.mean, quantile]).round(3).to_csv(output, sep="\t")
            print("table %s supplement materials saved" % output)
            return True
        except Exception as e:
            print("something wrong with %s" % str(e))
            return False
    else:
        return False


def generate_supplements_M8HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    df_of_params: Optional[pd.DataFrame] = wrapper_recover_params_M8HKY_codeml(
        input_dir=input_dir, metadata=metadata
    )
    if df_of_params is not None:
        try:
            df_of_params.groupby(by=["geneID", "CpG", "mix"])[
                ["tree", "kappa", "p", "q", "p0", "p1", "w_p1"]
                + ["p_k" + str(i) for i in range(1, 11, 1)]
                + ["w_k" + str(i) for i in range(1, 11, 1)]
            ].agg([np.mean, quantile]).round(3).to_csv(output, sep="\t")
            print("table %s supplement materials saved" % output)
            return True
        except Exception as e:
            print("something wrong with %s" % str(e))
            return False
    else:
        return False


def generate_supplements_M7M8HKY_codeml_LRT(
    input_dir_M7HKY_codeml: str,
    input_dir_M8HKY_codeml: str,
    metadata: Dict[str, Union[str, float]],
    output: str,
) -> bool:
    df_of_params_M7HKY: Optional[pd.DataFrame] = wrapper_recover_params_M7HKY_codeml(
        input_dir=input_dir_M7HKY_codeml, metadata=metadata
    )
    df_of_params_M8HKY: Optional[pd.DataFrame] = wrapper_recover_params_M8HKY_codeml(
        input_dir=input_dir_M8HKY_codeml, metadata=metadata
    )
    if df_of_params_M7HKY is not None and df_of_params_M8HKY is not None:
        try:
            df_of_lrt: pd.DataFrame = compute_likelihood_ratio_test(
                df_h0=df_of_params_M7HKY,
                df_h1=df_of_params_M8HKY,
                list_of_genes=LIST_OF_GENES,
                list_of_omega_mix=OMEGA_MIX,
                list_of_CpG=CPG,
            )
        except Exception as e:
            print("something wrong when recovering likelihood ratio test %s" % (str(e)))
            return False
        try:
            df_of_lrt[
                [
                    "geneID",
                    "CpG",
                    "mix",
                    "mean",
                    "quantile",
                    "alpha 0.10",
                    "alpha 0.05",
                    "alpha 0.01",
                    "count",
                ]
            ].round(3).to_csv(output, sep="\t")
            print("figure %s supplement materials saved" % output)
            return True
        except Exception as e:
            print("something wrong with %s" % str(e))
            return False
    return True


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
    dict_of_params["run"] = file.split("/")[-2]
    if isinstance(dict_of_params["run"], str):
        dict_of_params["CpG"] = float(dict_of_params["run"].split("-")[3])
        dict_of_params["geneID"] = dict_of_params["run"].split("-")[0]
        dict_of_params["fastaID"] = "-".join(dict_of_params["run"].split("-")[:7])
        dict_of_params["fixed"] = float(
            dict_of_params["run"].split("-")[4].split("_")[-1]
        )
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

    dict_of_params["run"] = file.split("/")[-2]
    if isinstance(dict_of_params["run"], str):
        dict_of_params["CpG"] = float(dict_of_params["run"].split("-")[3])
        dict_of_params["geneID"] = dict_of_params["run"].split("-")[0]
        dict_of_params["fastaID"] = "-".join(dict_of_params["run"].split("-")[:7])
        dict_of_params["mix"] = float(
            dict_of_params["run"].split("-")[4].split("_")[-1]
        )

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


def wrapper_recover_params_M0HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    try:
        files: List[str] = glob.glob(input_dir)
        print("number of files recovered from M0HKY: %d" % len(files))
        df_of_params: pd.DataFrame = recover_params_M0HKY(
            files=files, metadata=metadata
        )
        return df_of_params
    except Exception as e:
        print("something wrong when recovering M0HKY %s" % (str(e)))
        return None


def wrapper_recover_params_M7HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    try:
        files: List[str] = glob.glob(input_dir)
        print(" recovered from M7HKY: %d" % len(files))
        df_of_params: pd.DataFrame = recover_params_M7M8HKY(
            files=files, metadata=metadata, M7_bool=True
        )
        return df_of_params
    except Exception as e:
        print("something wrong when recovering M7HKY %s" % (str(e)))
        return None


def wrapper_recover_params_M8HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    try:
        files: List[str] = glob.glob(input_dir)
        print(" recovered from M8HKY: %d" % len(files))
        df_of_params: pd.DataFrame = recover_params_M7M8HKY(
            files=files, metadata=metadata, M7_bool=False
        )
        return df_of_params
    except Exception as e:
        print("something wrong when recovering M8HKY %s" % (str(e)))
        return None


def generate_supplements_codeml(
    input_dir_M0HKY_codeml: str,
    input_dir_M7HKY_codeml: str,
    input_dir_M8HKY_codeml: str,
    metadata_M0HKY_codeml: str,
    metadata_M7HKY_codeml: str,
    metadata_M8HKY_codeml: str,
    metadata_M7M8_LRT: str,
    output_M0HKY_codeml: str,
    output_M7HKY_codeml: str,
    output_M8HKY_codeml: str,
    output_M7M8_LRT: str,
):

    generate_supplements_M0HKY_codeml(
        input_dir=input_dir_M0HKY_codeml, metadata={}, output=output_M0HKY_codeml
    )
    generate_supplements_M7HKY_codeml(
        input_dir=input_dir_M7HKY_codeml, metadata={}, output=output_M7HKY_codeml
    )
    generate_supplements_M8HKY_codeml(
        input_dir=input_dir_M8HKY_codeml, metadata={}, output=output_M8HKY_codeml
    )
    generate_supplements_M7M8HKY_codeml_LRT(
        input_dir_M7HKY_codeml=input_dir_M7HKY_codeml,
        input_dir_M8HKY_codeml=input_dir_M8HKY_codeml,
        metadata={},
        output=output_M7M8_LRT,
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input_dir_M0HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--input_dir_M7HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--input_dir_M8HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--metadata_M0HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--metadata_M7HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--metadata_M8HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--metadata_M7M8_LRT",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output_M0HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output_M7HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output_M8HKY_codeml",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output_M7M8_LRT",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    generate_supplements_codeml(
        input_dir_M0HKY_codeml=args.input_dir_M0HKY_codeml,
        input_dir_M7HKY_codeml=args.input_dir_M7HKY_codeml,
        input_dir_M8HKY_codeml=args.input_dir_M8HKY_codeml,
        metadata_M0HKY_codeml=args.metadata_M0HKY_codeml,
        metadata_M7HKY_codeml=args.metadata_M7HKY_codeml,
        metadata_M8HKY_codeml=args.metadata_M8HKY_codeml,
        metadata_M7M8_LRT=args.metadata_M7M8_LRT,
        output_M0HKY_codeml=args.output_M0HKY_codeml,
        output_M7HKY_codeml=args.output_M7HKY_codeml,
        output_M8HKY_codeml=args.output_M8HKY_codeml,
        output_M7M8_LRT=args.output_M7M8_LRT,
    )
