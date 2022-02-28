import argparse
import glob
import os
import sys
from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional
from typing import Union

import numpy as np
import pandas as pd

ROOT_dir: str = Path(os.path.dirname(__file__)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

from scripts.utils import quantile


def wrapper_recover_params(
    input_dir: str, pattern: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    try:
        files: List[str] = glob.glob(input_dir + "/" + pattern)
        print("number of files recovered from M0HKY: %d" % len(files))
        df_of_params: pd.DataFrame = recover_params(files=files, metadata=metadata)
        return df_of_params
    except Exception as e:
        print("something wrong when recovering cabc %s" % (str(e)))
        return None


def recover_params(
    files: List[str], metadata: Dict[str, Union[str, float]]
) -> pd.DataFrame:

    list_of_df: List[pd.DataFrame] = []
    k = 0
    for f in files:
        list_of_df += [extract_values(file=f, metadata=metadata)]
        k += 1
    return pd.concat(list_of_df, ignore_index=True, axis=0)


def extract_values(
    file: str, metadata: Dict[str, Union[str, float]]
) -> Optional[pd.DataFrame]:
    try:
        tup = file.split("/")[-1].split("-")
        geneID = tup[0]
        modelID = tup[1]
        CpG = tup[3]
        omega = tup[4].split("_")[1]
        df_knn_cur = pd.read_csv(
            file,
            sep="\t",
        )
        for k, v in metadata.items():
            df_knn_cur[k] = [v] * df_knn_cur.shape[0]
        df_knn_cur["geneID"] = [geneID] * df_knn_cur.shape[0]
        df_knn_cur["modelID"] = [modelID] * df_knn_cur.shape[0]
        df_knn_cur["CpG"] = [CpG] * df_knn_cur.shape[0]
        df_knn_cur["omega"] = [omega] * df_knn_cur.shape[0]
        return df_knn_cur

    except Exception as e:
        print("something wrong when recovering %s, %s" % (file, str(e)))
        return None


def generate_supplements_CABC(
    input_dir: str, pattern: str, metadata: Dict[str, Union[str, float]], output: str
) -> bool:

    df_concat: pd.DataFrame = wrapper_recover_params(
        input_dir=input_dir, pattern=pattern, metadata=metadata
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
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    generate_supplements_CABC(
        input_dir=args.input_dir,
        pattern=args.pattern,
        metadata=args.metadata,
        output=args.output,
    )
