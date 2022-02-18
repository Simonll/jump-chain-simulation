import argparse
import glob
import inspect
import os
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Union

import numpy as np
import pandas as pd
from bintools.phylobayes.mcmc_parser import posterior_M0_GTR
from utils import quantile
from utils import tl_compute


def extract_values(
    file: str, metadata: Dict[str, Union[str, float]]
) -> Dict[int, Dict[str, Any]]:
    params: posterior_M0_GTR = posterior_M0_GTR.read(input_file=Path(file))
    dict_of_params: Dict[int, Dict[str, Any]] = {}

    for i in range(params.__sizeof__()):
        dict_of_omega: Dict[float, int] = {}
        for j in params.list_of_omega[i]:
            if j not in dict_of_omega:
                dict_of_omega[j] = 1
            else:
                dict_of_omega[j] += 1

        dict_of_params[i] = {
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
            dict_of_params[i].update(
                {"omega_" + str(k): v / len(params.list_of_omega[i])}
            )

        for k_, v_ in metadata.items():
            dict_of_params[i].update({k_: v_})

    return dict_of_params


def recover_params(
    files: List[str], metadata: Dict[str, Union[str, float]]
) -> pd.DataFrame:

    dict_run: Dict[int, Dict[int, Any]] = {}
    k = 0
    for f in files:
        dict_run[k] = extract_values(file=f, metadata=metadata)
        k += 1
    return pd.DataFrame.from_dict(data=dict_run, orient="index").dropna()


def wrapper_recover_params(
    input_dir: str, pattern: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    df_of_params: pd.DataFrame = pd.DataFrame()
    try:
        files: List[str] = glob.glob(input_dir + pattern)
        print("number of files recovered from M0HKY: %d" % len(files))
        df_of_params = recover_params(files=files, metadata=metadata)
        return df_of_params
    except Exception as e:
        print(
            "something wrong when recovering %s, using %s, %s"
            % (input_dir + pattern, inspect.stack()[0][3], str(e))
        )
        return None


def generate_supplements_M0HKY_pb_modified(
    input_dir: str, pattern: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:
    if os.path.exists(input_dir + "/df_of_params_M0HKY_pb_modified.pickle"):
        df_of_params_M0HKY = pd.read_pickle(
            input_dir + "/df_of_params_M0HKY_pb_modified.pickle"
        )
    else:
        df_of_params_M0HKY = wrapper_recover_params(
            input_dir=input_dir, pattern=pattern, metadata=metadata
        )
        df_of_params_M0HKY.to_pickle(
            input_dir + "/df_of_params_M0HKY_pb_modified.pickle"
        )

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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input_dir",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--pattern",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--metadata",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    generate_supplements_M0HKY_pb_modified(
        input_dir=args.input_dir,
        pattern=args.pattern,
        metadata=args.metadata,
        output=args.output,
    )
