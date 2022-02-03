import os
import sys
from pathlib import Path

ROOT_dir: Path = Path(os.path.dirname(workflow.snakefile)).__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"





rule fasta2phylip:
    input: 
        script=ROOT_dir+"/scripts/data_prep.py",
        fasta=ROOT_dir+"/data/{geneID}.fasta"
    output: ROOT_dir+"/outputs/data/{geneID}.phylip"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --fasta={input.fasta} \
        --phylip={output}
        """


rule generate_pb_mpi_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/data/{geneID}.phylip",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    output:ROOT_dir+"/outputs/pb_mpi_m0gtr/{geneID}-M0GTR-{repID}.sh"
    params:
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={input.phylip} \
        --tree={input.tree} \
        --output={output} \
        --local={params.local} \
        --docker={params.docker}
        """

rule run_pb_mpi_cmd:
    input: ROOT_dir+"/outputs/pb_mpi_m0gtr/{geneID}-M0GTR-{repID}.sh"
    output: ROOT_dir+"/outputs/pb_mpi_m0gtr/{geneID}-M0GTR-{repID}.chain"
    shell:
        """
        echo bash {input}
        bash {input}
        """


rule generate_M0HKY_fixed:
    input:
        script=ROOT_dir+"/scripts/generate_M0HKY_fixed.py",
        phylip=ROOT_dir+"/outputs/data_prep/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/pb_mpi_m0gtr/{geneID}-M0GTR-{repID}.chain",
    output:ROOT_dir+"/outputs/modified_chain_m0hky/{geneID}-M0HKY-{omega}-{draw}-{repID}.chain"
    params:
        burnin=600,
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --burnin={params.burnin} \
        --omega={wildcards.omega} \
        --output={output}
        """
rule generate_M0HKY_simu_with_fixed_values_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_M0M7HKY_simu_cmd.py",
        phylip=ROOT_dir+"/outputs/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/modified_chain_m0hky/{geneID}-M0HKY-{omega}-{draw}-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.sh",
        conf=ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.conf"
    params:
        ss="'pwAC pwAG pwAT pwCG pwCT pwGT dinuc31CG dinuc31TG dinuc31CA nuc3A nuc3C nuc3G nuc3T pwAA'",
        params="'chainID root lambda_CpG lambda_TBL lambda_omega nucsA nucsC nucsG nucsT nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT'",
        map="'Nsub Nsynsub dinucCGCA dinucCGTG dinucNSCGCA dinucNSCGTG dinucSCGCA dinucSCGTG gtnrAA gtnrAC gtnrAG gtnrAT gtnrCA gtnrCC gtnrCG gtnrCT gtnrGA gtnrGC gtnrGG gtnrGT gtnrTA gtnrTC gtnrTG gtnrTT'",
        localparams="'-iscodon -code Universal -freeroot Echinops Procavia -rootlength 100 -fixlambdaCpG {CpG} -tofasta'",
        model="M7",
        nsimu="1",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --ss={params.ss} \
        --map={params.map} \
        --params={params.params} \
        --localparams={params.localparams} \
        --model={params.model} \
        --nsimu={params.nsimu} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker}
        """

rule run_M0HKY_simu_with_fixed_values_cmd:
    input: ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.sh"
    output:
        simu=ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.simu",
        fasta=ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta"
    shell:
        """
        echo bash {input}
        bash {input}
        """


rule generate_M7HKY_mixed:
    input:
        script=ROOT_dir+"/scripts/generate_M7HKY_mixed.py",
        phylip=ROOT_dir+"/outputs/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/pb_mpi_m0gtr/{geneID}-M0GTR-{repID}.chain",
    output:ROOT_dir+"/outputs/modified_chain_m7hky/{geneID}-M7HKY-{omega_mix}-{draw}-{repID}.chain"
    params:
        burnin=600,
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3  {input.script} \
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --burnin={params.burnin} \
        --omega_mix={wildcards.omega_mix} \
        --output={output}
        """
rule generate_M7HKY_simu_with_mixed_values_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_M0M7HKY_simu_cmd.py",
        phylip=ROOT_dir+"/outputs/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/modified_chain_m7hky/{geneID}-M7HKY-{omega_mix}-{draw}-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.sh",
        conf=ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.conf"
    params:
        ss="'pwAC pwAG pwAT pwCG pwCT pwGT dinuc31CG dinuc31TG dinuc31CA nuc3A nuc3C nuc3G nuc3T pwAA'",
        params="'chainID root lambda_CpG lambda_TBL lambda_omega nucsA nucsC nucsG nucsT nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT'",
        map="'Nsub Nsynsub dinucCGCA dinucCGTG dinucNSCGCA dinucNSCGTG dinucSCGCA dinucSCGTG gtnrAA gtnrAC gtnrAG gtnrAT gtnrCA gtnrCC gtnrCG gtnrCT gtnrGA gtnrGC gtnrGG gtnrGT gtnrTA gtnrTC gtnrTG gtnrTT'",
        localparams="'-iscodon -code Universal -freeroot Echinops Procavia -rootlength 100 -fixlambdaCpG {CpG} -tofasta'",
        model="M7",
        nsimu="1",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --ss={params.ss} \
        --map={params.map} \
        --params={params.params} \
        --localparams={params.localparams} \
        --model={params.model} \
        --nsimu={params.nsimu} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker}
        """

rule run_M7HKY_simu_with_mixed_values_cmd:
    input: ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.sh"
    output:
        simu=ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.simu",
        fasta=ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0.fasta"
    shell:
        """
        echo bash {input}
        bash {input}
        """


rule generate_codeml_M0HKY_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_codeml_cmd.py",
        fasta=ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        yaml=ROOT_dir+"/configs/codeml.yaml"
    output:
        sh=ROOT_dir+"/outputs/codeml_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
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
        sh=ROOT_dir+"/outputs/codeml_m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml_m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
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
        sh=ROOT_dir+"/outputs/codeml_m8hky/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml_m8hky/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
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
    input:ROOT_dir+"/outputs/codeml_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:ROOT_dir+"/outputs/codeml_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule run_codeml_M7HKY_cmd:
    input:ROOT_dir+"/outputs/codeml_m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:ROOT_dir+"/outputs/codeml_m7hky/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule run_codeml_M8HKY_cmd:
    input:ROOT_dir+"/outputs/codeml_m8hky/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:ROOT_dir+"/outputs/codeml_m8hky/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
    shell:
        """
        echo bash {input}
        bash {input}
        """



rule generate_M0HKY_simu_cmd:
    input:
        phylip=ROOT_dir+"/outputs/data_prep/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/step_1/{geneID}-M0GTR-{repID}.chain",
    output:
        sh=ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.sh",
        conf=ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.conf"
    params:
        ss="'pwAC pwAG pwAT pwCG pwCT pwGT dinuc31CG dinuc31TG dinuc31CA nuc3A nuc3C nuc3G nuc3T pwAA'",
        params="'chainID root lambda_CpG lambda_TBL lambda_omega nucsA nucsC nucsG nucsT nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT'",
        map="'Nsub Nsynsub dinucCGCA dinucCGTG dinucNSCGCA dinucNSCGTG dinucSCGCA dinucSCGTG gtnrAA gtnrAC gtnrAG gtnrAT gtnrCA gtnrCC gtnrCG gtnrCT gtnrGA gtnrGC gtnrGG gtnrGT gtnrTA gtnrTC gtnrTG gtnrTT'",
        localparams="'-iscodon -code Universal -freeroot Echinops Procavia -rootlength 100 -freelambdaCpG -freelambdaTBL -freehky -freelambdaomega'",
        model="CodonMutSelFinite",
        nsimu="100000",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/generate_M0M7HKY_simu_cmd.py \
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --ss={params.ss} \
        --map={params.map} \
        --params={params.params} \
        --localparams={params.localparams} \
        --model={params.model} \
        --nsimu={params.nsimu} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker}
        """

rule run_M0HKY_simu_cmd:
    input: ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.sh"
    output: ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.simu"
    shell:
        """
        echo bash {input}
        bash {input}
        """
