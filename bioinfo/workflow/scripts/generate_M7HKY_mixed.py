import argparse
import os
from pathlib import Path
from typing import Dict

from bintools.align.align import read_phylip
from bintools.phylobayes.mcmc_parser import posterior_M0_GTR
from bintools.phylobayes.priors import prior_M7M8_mix

OMEGA_MIX = {
    "mix1": "0.1 0.2 0.3",
    "mix2": "0.4 0.5 0.6",
    "mix3": "0.7 0.8 0.9",
    "mix4": "0.2 0.5 0.7",
    "mix5": "0.5 0.7 0.9",
}


def generate_M7HKY_mixed(phylip, mcmc, burnin, omega_mix, output) -> bool:

    os.makedirs(
        os.path.dirname(output) + "/",
        exist_ok=True,
    )
    omega_mix = list(OMEGA_MIX[omega_mix].split())
    Nsite: int = -1
    with open(phylip) as stream:
        Nsite = int(read_phylip(fh=stream).get_n_site() / 3)
        pM0GTR: posterior_M0_GTR = posterior_M0_GTR.parse_mcmc(
            mcmc_path=Path(mcmc), burnin=burnin
        )
        params: Dict[str, str] = pM0GTR.sample_hky()
        params.update({"omega_site": prior_M7M8_mix(N=Nsite, mixture=omega_mix)})
        if os.path.exists(output):
            append_write = "a"  # append if already exists
        else:
            append_write = "w"  # make a new file if not

        with open(
            output,
            append_write,
        ) as fl:
            pM0GTR.write_values(dict_of_params=params, file_handler=fl)
    return True


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
        "--burnin",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--omega_mix",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    generate_M7HKY_mixed(
        phylip=args.phylip,
        mcmc=args.mcmc,
        burnin=args.burnin,
        omega_mix=args.omega_mix,
        output=args.output,
    )
