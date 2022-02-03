import glob
import os
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import ete3 as phylo
import numpy as np
import pandas as pd
from bintools.phylobayes.mcmc_parser import posterior_M0_GTR
from generate_figure_1 import LIST_OF_GENES
from generate_figure_1 import compute_likelihood_ratio_test
from generate_figure_1 import wrapper_recover_params_M0HKY
from generate_figure_1 import wrapper_recover_params_M7HKY
from generate_figure_1 import wrapper_recover_params_M8HKY


def generate_tables(
    input_M0HKY_pb_mpi_dir: str,
    metadata_M0HKY_pb_mpi: Dict[str, Union[str, float]],
    input_M0HKY_codeml_dir: str,
    metadata_M0HKY_codeml: Dict[str, Union[str, float]],
    input_M7HKY_codeml_dir: str,
    metadata_M7HKY_codeml: Dict[str, Union[str, float]],
    input_M8HKY_codeml_dir: str,
    metadata_M8HKY_codeml: Dict[str, Union[str, float]],
    output_M0HKY_pb_mpi: str,
    output_M0HKY_codeml: str,
    output_M7HKY_codeml: str,
    output_M8HKY_codeml: str,
) -> bool:

    try:
        generate_supplements_M0HKY(
            input_dir=input_M0HKY_pb_mpi_dir,
            metadata=metadata_M0HKY_pb_mpi,
            output=output_M0HKY_pb_mpi,
        )
    except Exception as e:
        print("something wrong with %s, %s" % (input_M0HKY_pb_mpi_dir, str(e)))
        return False
    try:
        generate_supplements_M0HKY_codeml(
            input_dir=input_M0HKY_codeml_dir,
            metadata=metadata_M0HKY_codeml,
            output=output_M0HKY_codeml,
        )
    except Exception as e:
        print("something wrong with %s, %s" % (input_M0HKY_codeml_dir, str(e)))
        return False
    try:
        generate_supplements_M0HKY_codeml(
            input_dir=input_M7HKY_codeml_dir,
            metadata=metadata_M7HKY_codeml,
            output=output_M7HKY_codeml,
        )
    except Exception as e:
        print("something wrong with %s, %s" % (input_M7HKY_codeml_dir, str(e)))
        return False
    try:
        generate_supplements_M7HKY_codeml(
            input_dir=input_M8HKY_codeml_dir,
            metadata=metadata_M8HKY_codeml,
            output=output_M8HKY_codeml,
        )
    except Exception as e:
        print("something wrong with %s, %s" % (input_M8HKY_codeml_dir, str(e)))
        return False
    return True


def quantile(
    arr: List[float], quantile_inf: float = 0.025, quantile_sup: float = 0.975
) -> Tuple[float, float]:
    quantile_inf_v: float = np.quantile(a=arr, q=quantile_inf)
    quantile_sup_v: float = np.quantile(a=arr, q=quantile_sup)
    return (round(quantile_inf_v, 3), round(quantile_sup_v, 3))


def params_to_dict(
    params: posterior_M0_GTR, geneID: str, mixID: str
) -> Dict[int, Dict[str, Any]]:
    dict_of_stats: Dict[int, Dict[str, Any]] = {}
    for i in range(params.__sizeof__()):

        dict_of_omega: Dict[float, int] = {}
        for j in params.list_of_omega[i]:
            if j not in dict_of_omega:
                dict_of_omega[j] = 1
            else:
                dict_of_omega[j] += 1

        dict_of_stats[i] = {
            "geneID": geneID,
            "mix": mixID,
            "tree": tl_compute(list_of_trees=[params.list_of_trees[i]])[0],
            "phi_A": params.list_of_phi[i][0],
            "phi_C": params.list_of_phi[i][1],
            "phi_G": params.list_of_phi[i][2],
            "phi_T": params.list_of_phi[i][3],
            "rho_AC": params.list_of_rho[i][0],
            "rho_AG": params.list_of_rho[i][1],
            "rho_AT": params.list_of_rho[i][2],
            "rho_CG": params.list_of_rho[i][3],
            "rho_CT": params.list_of_rho[i][4],
            "rho_TG": params.list_of_rho[i][5],
        }

        for k, v in dict_of_omega.items():
            dict_of_stats[i].update(
                {"omega_" + str(k): v / len(params.list_of_omega[i])}
            )

    return dict_of_stats


def generate_supplements_M0HKY(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    if os.path.exists(input_dir + "/df_of_params_M0HKY.pickle"):
        df_of_params_M0HKY = pd.read_pickle(input_dir + "/df_of_params_M0HKY.pickle")
    else:
        df_of_params_M0HKY = wrapper_recover_params_M0HKY(
            input=input_dir, metadata=metadata
        )
        df_of_params_M0HKY.to_pickle(input_dir + "/df_of_params_M0HKY.pickle")

    try:
        df_of_params_M0HKY.fillna(0).groupby(by=["geneID", "mix"])[
            [
                "tree",
                "phi_A",
                "phi_C",
                "phi_G",
                "phi_T",
                "rho_AC",
                "rho_AG",
                "rho_AT",
                "rho_CG",
                "rho_CT",
                "rho_TG",
                "omega_0.2",
                "omega_0.5",
                "omega_0.8",
            ]
        ].agg([np.mean, quantile]).round(3).to_csv(output, sep="\t")
        return True
    except Exception as e:
        print("something wrong with %s, %s" % (output, str(e)))
        return False


def wrapper_recover_params_CABC(
    input: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    df_of_params: pd.DataFrame = pd.DataFrame()
    try:
        files: List[str] = glob.glob(input)
        print("number of files recovered from M0HKY: %d" % len(files))
        df_of_params = recover_params_CABC(files=files, metadata=metadata)
        return df_of_params
    except Exception as e:
        print("something wrong when recovering M0HKY %s" % (str(e)))
        return None


def recover_params_CABC(
    files: List[str], metadata: Dict[str, Union[str, float]]
) -> pd.DataFrame:

    list_of_df: List[pd.DataFrame] = []
    k = 0
    for f in files:
        list_of_df += [extract_CABC(file=f, metadata=metadata)]
        k += 1
    return pd.concat(list_of_df, ignore_index=True)


def extract_CABC(
    file: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    try:
        df_knn_cur = pd.read_csv(
            file,
            sep="\t",
        )
        for k, v in metadata.items():
            df_knn_cur[k] = [v] * df_knn_cur.shape[0]
        return df_knn_cur

    except Exception as e:
        print("something wrong when recovering %s, %s" % (file, str(e)))
        return None


def generate_supplements_CABC(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:

    df_concat: pd.DataFrame = wrapper_recover_params_CABC(
        input=input_dir, metadata=metadata
    )
    df_concat.rename(
        inplace=True,
        columns={
            "lambda_CpG": "lambda",
            "nucsA": "phi_A",
            "nucsC": "phi_C",
            "nucsG": "phi_G",
            "nucsT": "phi_T",
            "nucrrAC": "rho_AC",
            "nucrrAG": "rho_AG",
        },
    )
    try:
        df_concat[[c for c in df_concat.columns.to_list() if c != "chainID"]].groupby(
            by=["geneID", "CpG", "omega"]
        ).agg([np.mean, quantile]).round(3).to_csv(output, sep="\t")
    except Exception as e:
        print("something wrong %s, %s" % (output, str(e)))
        return False

    return True


def tl_compute(list_of_trees) -> List[float]:
    list_of_tl: List[float] = []
    for t in list_of_trees:
        sum_ = 0
        t_ = phylo.Tree(t)
        for b in t_.traverse("postorder"):
            sum_ += b.dist
        list_of_tl += [sum_]
    return list_of_tl


def generate_supplements_M0HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    df_of_params: Optional[pd.DataFrame] = wrapper_recover_params_M0HKY(
        input=input_dir, metadata=metadata
    )
    if df_of_params is not None:
        try:
            df_of_params.groupby(by=["geneID", "CpG", "fixed"])[
                ["tree", "kappa", "omega"]
            ].agg([np.mean, quantile]).round(3).to_csv(output, sep="\t")
            print("figure %s supplement materials saved" % output)
            return True
        except Exception as e:
            print("something wrong with %s" % str(e))
            return False
    else:
        return False


def generate_supplements_M7HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    df_of_params: Optional[pd.DataFrame] = wrapper_recover_params_M7HKY(
        input=input_dir, metadata=metadata
    )
    if df_of_params is not None:
        try:
            df_of_params.groupby(by=["geneID", "CpG", "mix"])[
                ["tree", "kappa", "p", "q"]
                + ["p_k" + str(i) for i in range(1, 11, 1)]
                + ["w_k" + str(i) for i in range(1, 11, 1)]
            ].agg([np.mean, quantile]).round(3).to_csv(output, sep="\t")
            print("figure %s supplement materials saved" % output)
            return True
        except Exception as e:
            print("something wrong with %s" % str(e))
            return False
    else:
        return False


def generate_supplements_M8HKY_codeml(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    df_of_params: Optional[pd.DataFrame] = wrapper_recover_params_M8HKY(
        input=input_dir, metadata=metadata
    )
    if df_of_params is not None:
        try:
            df_of_params.groupby(by=["geneID", "CpG", "mix"])[
                ["tree", "kappa", "p", "q", "p0", "p1", "w_p1"]
                + ["p_k" + str(i) for i in range(1, 11, 1)]
                + ["w_k" + str(i) for i in range(1, 11, 1)]
            ].agg([np.mean, quantile]).round(3).to_csv(output, sep="\t")
            print("figure %s supplement materials saved" % output)
            return False
        except Exception as e:
            print("something wrong with %s" % str(e))
            return False
    else:
        return False


def generate_supplements_M7M8HKY_codeml_LRT(
    input_dir: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    df_of_params_M7HKY: Optional[pd.DataFrame] = wrapper_recover_params_M7HKY(
        input=input_dir, metadata=metadata
    )
    df_of_params_M8HKY: Optional[pd.DataFrame] = wrapper_recover_params_M8HKY(
        input=input_dir, metadata=metadata
    )
    if df_of_params_M7HKY is not None and df_of_params_M8HKY is not None:
        try:
            df_of_lrt: pd.DataFrame = compute_likelihood_ratio_test(
                df_h0=df_of_params_M7HKY,
                df_h1=df_of_params_M8HKY,
                list_of_genes=LIST_OF_GENES,
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
