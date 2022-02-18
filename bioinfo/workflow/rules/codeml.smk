import os
import sys
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"

ruleorder: generate_codeml_M0HKY_cmd > generate_codeml_M7HKY_cmd > generate_codeml_M8HKY_cmd > run_codeml_M0HKY_cmd > run_codeml_M7HKY_cmd > run_codeml_M8HKY_cmd

rule generate_codeml_M0HKY_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_codeml_cmd.py",
        fasta=ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        yaml=ROOT_dir+"/configs/codeml.yaml"
    output:
        sh=ROOT_dir+"/outputs/codeml/m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml/m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
    params:
        modelID="M0HKY",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --fasta={input.fasta} \
        --tree={input.tree} \
        --yaml={input.yaml} \
        --modelID={params.modelID} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker} \
        """


rule generate_codeml_M7HKY_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_codeml_cmd.py",
        fasta=ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        yaml=ROOT_dir+"/configs/codeml.yaml"
    output:
        sh=ROOT_dir+"/outputs/codeml/m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml/m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
    params:
        modelID="M7HKY",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --fasta={input.fasta} \
        --tree={input.tree} \
        --yaml={input.yaml} \
        --modelID={params.modelID} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker} \
        """



rule generate_codeml_M8HKY_cmd:
    input:
        scirpt=ROOT_dir+"/scripts/generate_codeml_cmd.py",
        fasta=ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        yaml=ROOT_dir+"/configs/codeml.yaml"
    output:
        sh=ROOT_dir+"/outputs/codeml/m8hky/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml/m8hky/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
    params:
        modelID="M8HKY",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --fasta={input.fasta} \
        --tree={input.tree} \
        --yaml={input.yaml} \
        --modelID={params.modelID} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker} \
        """


rule run_codeml_M0HKY_cmd:
    input:
        ROOT_dir+"/outputs/codeml/m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:
        ROOT_dir+"/outputs/codeml/m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule run_codeml_M7HKY_cmd:
    input:
        ROOT_dir+"/outputs/codeml/m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:
        ROOT_dir+"/outputs/codeml/m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule run_codeml_M8HKY_cmd:
    input:
        ROOT_dir+"/outputs/codeml/m8hky/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:
        ROOT_dir+"/outputs/codeml/m8hky/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
    shell:
        """
        echo bash {input}
        bash {input}
        """
