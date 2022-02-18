import os
import sys
from pathlib import Path

ROOT_dir: Path = Path(os.path.dirname(workflow.snakefile)).__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"
GENEID = config["geneID"]
REPID = ["A"]#["A","B"]
DRAW = [i for i in range(0, 100)][0:1]
OMEGA = [0.2, 0.5, 0.8]
OMEGA_MIX = {"mix1":"0.1 0.2 0.3","mix2":"0.4 0.5 0.6","mix3":"0.7 0.8 0.9","mix4":"0.2 0.5 0.7","mix5":"0.5 0.7 0.9"}
CPG = [1, 2, 4, 8]
rule all:
    input:
        # PB_MPI M0GTR
        expand(ROOT_dir+"/outputs/data/{geneID}.phylip", geneID=GENEID),
        expand(ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.sh",geneID=GENEID, repID=REPID),
        # expand(ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.chain",geneID=GENEID, repID=REPID),
        # # SIMU M0HKY
        # expand(ROOT_dir+"/outputs/modified_chain_m0hky/{geneID}-M0HKY-{omega}-{draw}-{repID}.chain", geneID=GENEID, omega=OMEGA, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.sh",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.conf",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.simu",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # # SIMU M7HKY
        # expand(ROOT_dir+"/outputs/modified_chain_m7hky/{geneID}-M7HKY-{omega_mix}-{draw}-{repID}.chain",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.sh",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.conf",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.simu",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0.fasta",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),

include: ROOT_dir + "/rules/simu.smk"
