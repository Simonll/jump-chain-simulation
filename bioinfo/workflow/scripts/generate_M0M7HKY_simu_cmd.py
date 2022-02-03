import argparse
import os
import re
from typing import Dict

from bintools.cabc.cabc import generate_cabc_conf
from bintools.run.run import generate_simu_cmd


def generate_mapping(local, docker, input) -> str:
    s: str = docker + re.sub(local, "", input)
    return s


def wrapper_generate_M0M7HKY_simu_cmd(
    phylip, mcmc, params, ss, map, nsimu, localparams, model, output, local, docker
) -> str:
    kwargs: Dict[str, str] = {
        "align": generate_mapping(local=local, docker=docker, input=phylip),
        "chainname": generate_mapping(local=local, docker=docker, input=mcmc)[:-6],
        "param": params,
        "ss": ss,
        "map": map,
        "sampling": "600 4 999",
        "nrun": " ".join(["1", nsimu]),
        "nthreads": "1",
        "output": generate_mapping(local=local, docker=docker, input=output[:-3]),
        "localparam": localparams,
    }
    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )

    with open(
        output[:-3] + ".conf",
        "w",
    ) as fh:
        for l in generate_cabc_conf(
            method=model, mapping=local + ":" + docker, **kwargs
        ):
            fh.write(l + "\n")

    kwargs = {"-m": model}

    logger: str = " ".join(["2>", output[:-3] + ".log"])

    cmd: str = generate_simu_cmd(
        method=model,
        mapping=local + ":" + docker,
        config_filename=generate_mapping(local=local, docker=docker, input=output[:-3])
        + ".conf",
        logger=logger,
        **kwargs
    )

    with open(output, "w") as fh:
        fh.write(cmd)

    return cmd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--phylip",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--mcmc",
        type=str,
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
    parser.add_argument(
        "--map",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--localparams",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--model",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--nsimu",
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
    wrapper_generate_M0M7HKY_simu_cmd(
        phylip=args.phylip,
        mcmc=args.mcmc,
        ss=args.ss,
        params=args.params,
        map=args.map,
        localparams=args.localparams,
        model=args.model,
        nsimu=args.nsimu,
        output=args.output,
        local=args.local,
        docker=args.docker,
    )
