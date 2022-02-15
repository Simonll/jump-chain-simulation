import os
import sys
from pathlib import Path

ROOT_dir: Path = Path(os.path.dirname(workflow.snakefile)).__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"
GENEID = ["CSRP2BP"]#config["geneID"]
REPID = ["A"]#["A","B"]
DRAW = [i for i in range(0, 100)][0:1]
OMEGA = [0.2, 0.5, 0.8]
OMEGA_MIX = {"mix1":"0.1 0.2 0.3","mix2":"0.4 0.5 0.6","mix3":"0.7 0.8 0.9","mix4":"0.2 0.5 0.7","mix5":"0.5 0.7 0.9"}
CPG = [1, 2, 4, 8]





rule all:
    input:
        CABC M0HKY
        expand(ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.sh",geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.conf",geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.simu",geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_chainID.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_params.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_ss.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_true_ss.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-model_preprocessing_params.pickle",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-cabc.sh",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-knn.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-post.tsv",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID)

include: ROOT_dir + "/rules/cabc.smk"
