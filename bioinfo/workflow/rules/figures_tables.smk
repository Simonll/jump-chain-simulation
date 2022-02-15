import os
import sys
from pathlib import Path

ROOT_dir: Path = Path(os.path.dirname(workflow.snakefile)).__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"



rule generate_tables:

    input:
        script=ROOT_dir+"/scripts/generate_tables.py",
        input_M0HKY_pb_mpi_dir=ROOT_dir+"/outputs/",
        input_M0HKY_codeml_dir=ROOT_dir+"/outputs/",
        input_M7HKY_codeml_dir=ROOT_dir+"/outputs/",
        input_M8HKY_codeml_dir=ROOT_dir+"/outputs/",
    output:
        output_M0HKY_pb_mpi=ROOT_dir+"/reports/",
        output_M0HKY_codeml=ROOT_dir+"/reports/",
        output_M7HKY_codeml=ROOT_dir+"/reports/",
        output_M8HKY_codeml=ROOT_dir+"/reports/",
    params:
        metadata_M0HKY_pb_mpi="",
        metadata_M0HKY_codeml="",
        metadata_M7HKY_codeml="",
        metadata_M8HKY_codeml="",
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input_M0HKY_pb_mpi_dir={input.input_M0HKY_pb_mpi_dir} \
        --metadata_M0HKY_codeml={params.metadata_M0HKY_codeml} \
        --input_M0HKY_codeml_dir={input.input_M0HKY_codeml_dir} \
        --metadata_M0HKY_codeml={params.metadata_M0HKY_codeml} \
        --input_M7HKY_codeml_dir={input.input_M7HKY_codeml_dir} \
        --metadata_M7HKY_codeml={params.metadata_M7HKY_codeml} \
        --input_M8HKY_codeml_dir={input.input_M8HKY_codeml_dir} \
        --metadata_M8HKY_codeml={params.metadata_M8HKY_codeml} \
        --output_M0HKY_pb_mpi={output.output_M0HKY_pb_mpi} \
        --output_M0HKY_codeml={output.output_M0HKY_codeml} \
        --output_M7HKY_codeml={output.output_M7HKY_codeml} \
        --output_M8HKY_codeml={output.output_M8HKY_codeml}
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
    output:ROOT_dir+"/reports/figure_1.png"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python4 {input.script} \
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
    output: ROOT_dir+"/reports/figure_2.png"
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
