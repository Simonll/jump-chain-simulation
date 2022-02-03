import argparse
import os
import re
from typing import Dict

from bintools.run.run import generate_abc_cmd


def generate_mapping(local, docker, input) -> str:
    s: str = docker + re.sub(local, "", input)
    return s


def wrapper_generate_cabc_cmd(
    df_true_ss,
    df_simu_space_knn_params,
    df_simu_space_knn_ss,
    df_simu_space_knn_chainID,
    output,
    local,
    docker,
) -> str:
    kwargs: Dict[str, str] = {
        "--reg_model": "loclinear",
        "--hcorr": "TRUE",
        "--kernel": "epanechnikov",
        "--transf": "none",
        "--output": generate_mapping(local=local, docker=docker, input="$1"),
        "--df_true_ss": generate_mapping(local=local, docker=docker, input=df_true_ss),
        "--df_simu_space_knn_params": generate_mapping(
            local=local, docker=docker, input=df_simu_space_knn_params
        ),
        "--df_simu_space_knn_ss": generate_mapping(
            local=local, docker=docker, input=df_simu_space_knn_ss
        ),
        "--df_simu_space_knn_chainID": generate_mapping(
            local=local, docker=docker, input=df_simu_space_knn_chainID
        ),
    }

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )

    cmd: str = generate_abc_cmd(method="abc", mapping=local + ":" + docker, **kwargs)

    with open(output, "w") as fh:
        fh.write(cmd)

    return cmd


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--df_true_ss",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--df_simu_space_knn_params",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--df_simu_space_knn_ss",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--df_simu_space_knn_chainID",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--local",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--docker",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    wrapper_generate_cabc_cmd(
        df_true_ss=args.df_true_ss,
        df_simu_space_knn_params=args.df_simu_space_knn_params,
        df_simu_space_knn_ss=args.df_simu_space_knn_ss,
        df_simu_space_knn_chainID=args.df_simu_space_knn_chainID,
        output=args.output,
        local=args.local,
        docker=args.docker,
    )
