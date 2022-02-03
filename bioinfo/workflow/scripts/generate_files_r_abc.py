import argparse
import os
import pickle
from typing import List
from typing import Optional

import numpy as np
import pandas as pd
from scipy.spatial import distance
from sklearn.base import BaseEstimator
from sklearn.preprocessing import FunctionTransformer
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PowerTransformer
from sklearn.preprocessing import QuantileTransformer
from sklearn.preprocessing import StandardScaler


def transform(trans_name: str, df: pd.DataFrame) -> Optional[BaseEstimator]:
    pre_process: Optional[BaseEstimator] = None
    if trans_name == "qt":
        pre_process = QuantileTransformer()
        pre_process.fit(df)
    if trans_name == "boxcox":
        pre_process = PowerTransformer(method="box-cox", standardize=False)
        pre_process.fit(df)
    if trans_name == "minmax":
        pre_process = MinMaxScaler()
        pre_process.fit(df)
    elif trans_name == "standard":
        pre_process = StandardScaler()
        pre_process.fit(df)
    elif trans_name == "log":
        pre_process = FunctionTransformer(func=np.log, inverse_func=np.exp)
        pre_process.fit(df)
    elif trans_name == "log2":
        pre_process = FunctionTransformer(func=np.log2, inverse_func=np.exp2)
        pre_process.fit(df)
    elif trans_name == "log10":
        pre_process = FunctionTransformer(
            func=np.log10, inverse_func=lambda x: np.power(x, 10)
        )
        pre_process.fit(df)
    elif trans_name == "none":
        pre_process = FunctionTransformer(func=lambda x: x, inverse_func=lambda x: x)
        pre_process.fit(df)
    else:
        print("something wrong with transformation name: %s" % trans_name)

    return pre_process


def generate_files_r_abc(
    simu_space_file: str,
    knn: int,
    params: str,
    ss: str,
    preprocessing_ss: str,
    preprocessing_params: str,
    true_ss_file: str,
    output_dir: str,
    prefix: str,
) -> bool:

    list_of_params: List[str] = list(params.split())
    list_of_ss: List[str] = list(ss.split())

    os.makedirs(
        os.path.dirname(output_dir) + "/",
        exist_ok=True,
    )

    print("Here the list of params and summary statistics")
    print(list_of_params)
    print(list_of_ss)

    df_simu_space: pd.DataFrame = pd.read_csv(
        simu_space_file, sep="\t", index_col=False
    )

    print("simulation space dimensions", np.shape(df_simu_space))

    col: List[str] = list(df_simu_space.columns)
    df_simu_space.drop(
        [i for i in col if i not in list_of_params + list_of_ss + ["chainID"]],
        axis=1,
        inplace=True,
    )
    print("shape of reference table", np.shape(df_simu_space))
    df_simu_space.dropna(inplace=True)
    print("shape of reference table after droping na", np.shape(df_simu_space))

    # reading realdata
    df_true_ss: pd.DataFrame = pd.read_csv(true_ss_file, sep="\t", index_col=False)
    df_true_ss.drop(
        [i for i in list(df_true_ss.columns) if (i not in list_of_ss)],
        axis=1,
        inplace=True,
    )
    random_row: int = np.random.randint(0, df_true_ss.shape[0] - 1)
    df_true_ss = df_true_ss.iloc[random_row : random_row + 1, :]
    print("realdata table ", np.shape(df_true_ss))

    model_preprocessing_ss: Optional[BaseEstimator] = transform(
        trans_name=preprocessing_ss, df=df_simu_space[list_of_ss]
    )

    if model_preprocessing_ss is None:
        print("Something wrong with pre-processing of summary statistics")
        raise RuntimeError

    model_preprocessing_params: Optional[BaseEstimator] = transform(
        trans_name=preprocessing_params, df=df_simu_space[list_of_params]
    )

    if model_preprocessing_params is None:
        print("Something wrong with pre-processing of params")
        raise RuntimeError

    try:
        with open(output_dir + prefix + "-model_preprocessing_ss.pickle", "wb") as fh:
            pickle.dump(model_preprocessing_ss, fh, pickle.HIGHEST_PROTOCOL)
    except Exception as e:
        print(
            "something wrong while dumping %s on %s"
            % (prefix + "model_preprocessing_ss", output_dir)
        )

    try:
        with open(
            output_dir + prefix + "-model_preprocessing_params.pickle", "wb"
        ) as fh:
            pickle.dump(model_preprocessing_params, fh, pickle.HIGHEST_PROTOCOL)
    except Exception as e:
        print(
            "something wrong while dumping %s on %s"
            % (prefix + "model_preprocessing_params", output_dir)
        )

    assert df_true_ss.shape[1] == df_simu_space[list_of_ss].shape[1]

    df_simu_space["metric"] = distance.cdist(
        model_preprocessing_ss.transform(
            df_simu_space[list_of_ss],
        ),
        model_preprocessing_ss.transform(df_true_ss),
        lambda u, v: ((u - v) ** 2).sum(),
    )

    df_simu_space.sort_values(by="metric", inplace=True)
    df_simu_space.reset_index(inplace=True)
    df_simu_space = df_simu_space.iloc[:knn, :]
    assert df_simu_space.shape[0] == knn, "something wrong with %d %d" % (
        df_simu_space.shape[0],
        knn,
    )
    try:
        pd.DataFrame(
            data=df_simu_space["chainID"],
            columns=["chainID"],
        ).to_feather(output_dir + prefix + "-df_simu_space_knn_chainID.feather")
        pd.DataFrame(
            data=model_preprocessing_params.transform(df_simu_space[list_of_params]),
            columns=list_of_params,
        ).to_feather(output_dir + prefix + "-df_simu_space_knn_params.feather")
        pd.DataFrame(
            data=model_preprocessing_ss.transform(df_simu_space[list_of_ss]),
            columns=list_of_ss,
        ).to_feather(output_dir + prefix + "-df_simu_space_knn_ss.feather")
        pd.DataFrame(
            data=model_preprocessing_ss.transform(df_true_ss[list_of_ss]),
            columns=list_of_ss,
        ).reset_index()[list_of_ss].to_feather(
            output_dir + prefix + "-df_true_ss.feather"
        )
        # http://onnx.ai/sklearn-onnx/ for better portability and durability

    except Exception as e:
        print("Something wrong saving files for abc r package %s" % str(e))

        return False

    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--simu",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--knn",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--params",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--ss",
        type=str,
        required=True,
    )

    parser.add_argument("--output_dir", type=str, required=True)

    parser.add_argument(
        "--trans_fct_ss",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--trans_fct_params",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--true_ss_file",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--prefix",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    generate_files_r_abc(
        simu_space_file=args.simu,
        knn=args.knn,
        params=args.params,
        ss=args.ss,
        preprocessing_ss=args.trans_fct_ss,
        preprocessing_params=args.trans_fct_params,
        true_ss_file=args.true_ss_file,
        prefix=args.prefix,
        output_dir=args.output_dir,
    )
