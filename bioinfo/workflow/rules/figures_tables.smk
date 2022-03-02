import os
import sys
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"


rule generate_figure_1:
    input:
        script=ROOT_dir+"/scripts/generate_figure_1.py",
        input_M0HKY_codeml_dir=ROOT_dir+"/outputs/codeml/m0hky/",
        input_M7HKY_codeml_dir=ROOT_dir+"/outputs/codeml/m7hky/",
        input_M8HKY_codeml_dir=ROOT_dir+"/outputs/codeml/m8hky/",
    output:
        ROOT_dir+"/reports/figure_1.png"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_M0HKY_dir={input.input_M0HKY_codeml_dir} \
        --input_M7HKY_dir={input.input_M7HKY_codeml_dir} \
        --input_M8HKY_dir={input.input_M8HKY_codeml_dir} \
        --output={output} \
        """


rule generate_figure_2:
    """
    STRIP1 was used to show under-estimation of lambda parameter
    GPAM was used to show close centered estimation of lambda parameter
    WDR91 was used to show over-estimation of lambda parameter
    """

    input:
        script=ROOT_dir+"/scripts/generate_figure_2.py",
        input_under=ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-post.tsv,
        input_center=ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-post.tsv,
        input_over=ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-post.tsv,
    output:
        ROOT_dir+"/reports/figure_2.png"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_under={input.input_under} \
        --input_center={input.input_center} \
        --input_over={input.input_over} \
        --outout={output}
        """


rule generate_tableSup1:
    input:
        script=ROOT_dir+"/scripts/generate_M0HKY_pb_modified_sup.py",
        input_dir=ROOT_dir+"/outputs/modified_chain_m0hky/"
    params:
        pattern=""
        metadata="{}"
    output:
        ROOT_dir+"/reports/tableSup1.csv"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_dir={input.input_dir} \
        --pattern={params.pattern} \
        --metadata={params.metadata} \
        --outout={output}
        """

rule generate_tableSup2:
    input:
        script=ROOT_dir+"/scripts/generate_M0HKY_pb_modified_sup.py",
        input_dir=ROOT_dir+"/outputs/modified_chain_m7hky/"
    params:
        pattern=""
        metadata="{}"
    output:
        ROOT_dir+"/reports/tableSup2.csv"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_dir={input.input_dir} \
        --pattern={params.pattern} \
        --metadata={params.metadata} \
        --outout={output}
        """

rule generate_tableSup3_6:

    input:
        script=ROOT_dir+"/scripts/generate_codeml_sup.py",
        input_M0HKY_codeml_dir=ROOT_dir+"/outputs/",
        input_M7HKY_codeml_dir=ROOT_dir+"/outputs/",
        input_M8HKY_codeml_dir=ROOT_dir+"/outputs/",
    output:
        output_M0HKY_codeml=ROOT_dir+"/reports/tableSup1.csv",
        output_M7HKY_codeml=ROOT_dir+"/reports/tableSup2.csv",
        output_M8HKY_codeml=ROOT_dir+"/reports/tableSup3.csv",
        output_M7M8_LRT=ROOT_dir+"/reports/tableSup4.csv",
    params:
        metadata_M0HKY_codeml="",
        metadata_M7HKY_codeml="",
        metadata_M8HKY_codeml="",
        metadata_M7M8_LRT="",
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_M0HKY_codeml_dir={input.input_M0HKY_codeml_dir} \
        --metadata_M0HKY_codeml={params.metadata_M0HKY_codeml} \
        --input_M7HKY_codeml_dir={input.input_M7HKY_codeml_dir} \
        --metadata_M7HKY_codeml={params.metadata_M7HKY_codeml} \
        --input_M8HKY_codeml_dir={input.input_M8HKY_codeml_dir} \
        --metadata_M8HKY_codeml={params.metadata_M8HKY_codeml} \
        --metadata_M7M8_LRT={param.metadata_M7M8_LRT} \
        --output_M0HKY_codeml={output.output_M0HKY_codeml} \
        --output_M7HKY_codeml={output.output_M7HKY_codeml} \
        --output_M8HKY_codeml={output.output_M8HKY_codeml} \
        --output_M7M8_LRT={output.output_M7M8_LRT}
        """



rule generate_supplements_CABC:
    input:
        script=ROOT_dir+"/scripts/generate_cabc_sup.py",
        input_dir=ROOT_dir+"/outputs/cabc/step_3/"
    params:
        pattern="*-post.tsv",
        metadata="{}",
    output:
        ROOT_dir + "/reports/tableSup7.csv",
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_dir={input.input_dir} \
        --pattern={params.pattern} \
        --metadata={params.metadata} \
        --outout={output}
        """


rule generate_supplements_M0GTR_pb_mpi:
    input:
        script=ROOT_dir+"/scripts/generate_M0GTR_pb_mpi_sup.py",
        input_dir=ROOT_dir+"/outputs/pbmpi_m0gtr/"
    params:
        pattern="*-chain",
        metadata="{}",
    output:
        ROOT_dir + "/reports/tableSup8.csv",
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_dir={input.input_dir} \
        --pattern={params.pattern} \
        --metadata={params.metadata} \
        --outout={output}
        """
