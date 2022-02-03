import argparse
import os
import pickle
from typing import List
from typing import Optional

import pandas as pd
from sklearn.base import BaseEstimator


def recover_r_abc_results(
    knn_file: str, model_preprocessing_params_file: str, output: str
) -> bool:
    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )
    try:
        model_preprocessing_params: Optional[BaseEstimator] = None
        with open(model_preprocessing_params_file, "rb") as fh:
            model_preprocessing_params = pickle.load(fh)
    except Exception as e:
        print("something wrong with %s" % model_preprocessing_params_file)
        return False
    try:
        df_knn: pd.DataFrame = pd.read_feather(knn_file)
    except Exception as e:
        print("something wrong with %s" % knn_file)
        return False
    try:
        columns: List[str] = [c for c in df_knn.columns.to_list() if c != "chainID"]
        df_knn_: pd.DataFrame = pd.DataFrame(
            data=model_preprocessing_params.inverse_transform(df_knn[columns]),
        )

        df_knn_["chainID"] = df_knn["chainID"]
        df_knn_["phi_sum"] = df_knn_[["nucsA", "nucsC", "nucsG", "nucsT"]].sum(axis=1)
        if all(
            [
                True if i in columns else False
                for i in [
                    "nucrrAC",
                    "nucrrAG",
                    "nucrrAT",
                    "nucrrCG",
                    "nucrrCT",
                    "nucrrGT",
                ]
            ]
        ):
            df_knn_["rho_sum"] = df_knn_[
                ["nucrrAC", "nucrrAG", "nucrrAT", "nucrrCG", "nucrrCT", "nucrrGT"]
            ].sum(axis=1)
            for rho_i in [
                "nucrrAC",
                "nucrrAG",
                "nucrrAT",
                "nucrrCG",
                "nucrrCT",
                "nucrrGT",
            ]:
                df_knn_[rho_i] = df_knn_[rho_i].div(df_knn_["rho_sum"], axis=0)
        else:
            df_knn_["rho_sum"] = df_knn_[["nucrrAC", "nucrrAG"]].sum(axis=1)
            for rho_i in ["nucrrAC", "nucrrAG"]:
                df_knn_[rho_i] = df_knn_[rho_i].div(df_knn_["rho_sum"], axis=0)
        for phi_i in ["nucsA", "nucsC", "nucsG", "nucsT"]:
            df_knn_[phi_i] = df_knn_[phi_i].div(df_knn_["phi_sum"], axis=0)
        df_knn_[columns + ["chainID"]].to_csv(output, sep="\t", index=False)
        return True
    except Exception as e:
        print("something wrong with %s, %s" % (output, str(e)))
        return False


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--knn_file",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--model_preprocessing_params_file",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    recover_r_abc_results(
        knn_file=args.knn_file,
        model_preprocessing_params_file=args.model_preprocessing_params_file,
        output=args.output,
    )
