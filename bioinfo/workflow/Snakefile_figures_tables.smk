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

include: ROOT_dir + "/rules/figures_tables.smk"

rule all:
    input:
        ROOT_dir+"/reports/figure_1.png",
        ROOT_dir+"/reports/figure_2.png",
        ROOT_dir+"/reports/tableSup1.csv",
        ROOT_dir+"/reports/tableSup2.csv",
        ROOT_dir+"/reports/tableSup3.csv",
        ROOT_dir+"/reports/tableSup4.csv",
        ROOT_dir+"/reports/tableSup5.csv",
        ROOT_dir+"/reports/tableSup6.csv",
        ROOT_dir+"/reports/tableSup7.csv",
        ROOT_dir+"/reports/tableSup8.csv",
