import argparse
import glob
import os
import sys
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Union

import numpy as np
import pandas as pd
from bintools.phylobayes.mcmc_parser import posterior_M0_GTR

ROOT_dir: str = Path(os.path.dirname(__file__)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)


from scripts.utils import quantile
from scripts.utils import tl_compute


def extract_values(
    file: str, metadata: Dict[str, Union[str, float]]
) -> Dict[int, Dict[str, Any]]:
    params: posterior_M0_GTR = posterior_M0_GTR.parse_mcmc(
        mcmc_path=Path(file), burnin=600
    )
    dict_of_params: Dict[int, Dict[str, Any]] = {}
    for i in range(params.burnin, params.__sizeof__(), 1):
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
            "omega": params.list_of_omega[i],
        }
        for k, v in metadata.items():
            dict_of_params[i].update({k: v})
        geneID: str = file.split("/")[-1].split("-")[0]
        dict_of_params[i].update({"geneID": geneID})
    return dict_of_params


def recover_params(
    files: List[str], metadata: Dict[str, Union[str, float]]
) -> pd.DataFrame:
    list_of_df: List[pd.DataFrame] = []
    for f in files:
        list_of_df += [
            pd.DataFrame.from_dict(
                data=extract_values(file=f, metadata=metadata), orient="index"
            )
        ]
    return pd.concat(list_of_df, ignore_index=True, axis=0)


def wrapper_recover_params(
    input_dir: str, pattern: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    try:
        files: List[str] = glob.glob(input_dir + pattern)
        assert len(files) > 0, "Nothing to glob %s" % input_dir + pattern
        print("number of files recovered from M0GTR: %d" % len(files))
        df_of_params: pd.DataFrame = recover_params(files=files, metadata=metadata)
        return df_of_params
    except Exception as e:
        print("something wrong when recovering M0GTR %s" % (str(e)))
        return None


def generate_supplements_M0GTR_pb_mpi(
    input_dir: str, pattern: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:

    df_of_params_M0HKY = wrapper_recover_params(
        input_dir=input_dir, pattern=pattern, metadata=metadata
    )
    if df_of_params_M0HKY is None:
        print("something wrong with %s" % input_dir + pattern)
        return False
    try:
        df_of_params_M0HKY.fillna(0).groupby(by=["geneID"])[
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
                "omega",
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
        default="",
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

    generate_supplements_M0GTR_pb_mpi(
        input_dir=args.input_dir,
        pattern=args.pattern,
        metadata=args.metadata,
        output=args.output,
    )
