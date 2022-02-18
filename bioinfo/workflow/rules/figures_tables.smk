import os
import sys
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"



rule generate_tableSup1_4:

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



rule generate_figure_1:
    input:
        script=ROOT_dir+"/scripts/generate_figure_1.py",
        input_M0HKY_dir=ROOT_dir+"/outputs/",
        metadata_M0HKY=ROOT_dir+"/outputs/",
        input_M7HKY_dir=ROOT_dir+"/outputs/",
        metadata_M7HKY=ROOT_dir+"/outputs/",
        input_M8HKY_dir=ROOT_dir+"/outputs/",
        metadata_M8HKY=ROOT_dir+"/outputs/",
        output=ROOT_dir+"/reports/"
    output:
        ROOT_dir+"/reports/figure_1.png"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_M0HKY_dir={input.input_M0HKY_dir} \
        --metadata_M0HKY={input.metadata_M0HKY} \
        --input_M7HKY_dir= \
        --input_M8HKY_dir= \
        --input_M8HKY_dir= \
        --metadata_M8HKY= \
        --output={output} \
        """


rule generate_figure_2:
    input:
        script=ROOT_dir+"/scripts/generate_figure_2.py",
        input_under=ROOT_dir+"/outputs/",
        input_center=ROOT_dir+"/outputs/",
        input_over=ROOT_dir+"/outputs/",
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
