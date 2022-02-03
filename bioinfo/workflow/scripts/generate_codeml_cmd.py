import argparse
import copy
import os
import re
from typing import Any
from typing import Dict
from typing import Optional

from bintools.codeml.codeml import generate_codeml_conf
from bintools.run.run import generate_codeml_cmd
from bintools.utils.utils import get_yaml_config


def generate_mapping(local, docker, input) -> str:
    s: str = docker + re.sub(local, "", input)
    return s


def wrapper_generate_codeml_cmd(
    fasta, tree, yaml, modelID, output, local, docker
) -> str:

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )

    step_config: Optional[Dict[Any, Any]] = get_yaml_config(yaml_file=yaml)

    if step_config is None:
        print("something wrong with %s" % yaml)
        raise RuntimeError

    kwargs: Dict[str, str] = copy.deepcopy(step_config["codeml"][modelID])
    kwargs.update(
        {
            "seqfile": generate_mapping(local=local, docker=docker, input=fasta),
            "outfile": generate_mapping(local=local, docker=docker, input=output[:-3]),
            "treefile": generate_mapping(local=local, docker=docker, input=tree),
        }
    )

    with open(
        output[:-3] + ".conf",
        "w",
    ) as fh:
        for l in generate_codeml_conf(**kwargs):
            fh.write(l)

    logger: str = " ".join(
        [
            "2>",
            output[:-3] + ".log",
        ]
    )

    cmd: str = generate_codeml_cmd(
        method=modelID[0:2],
        mapping=local + ":" + docker,
        config_filename=generate_mapping(local=local, docker=docker, input=output[:-3])
        + ".conf",
        logger=logger,
    )
    with open(output, "w") as fh:
        fh.write(cmd)

    return cmd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--tree",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--yaml",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--modelID",
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
    wrapper_generate_codeml_cmd(
        fasta=args.fasta,
        tree=args.tree,
        modelID=args.modelID,
        yaml=args.yaml,
        output=args.output,
        local=args.local,
        docker=args.docker,
    )
