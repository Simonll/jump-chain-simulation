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
        #PB_MPI M0GTR
        expand(ROOT_dir+"/outputs/data_prep/{geneID}.phylip", geneID=GENEID),
        expand(ROOT_dir+"/outputs/step_1/{geneID}-M0GTR-{repID}.sh",geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/step_1/{geneID}-M0GTR-{repID}.chain",geneID=GENEID, repID=REPID),
        # SIMU M0HKY
        expand(ROOT_dir+"/outputs/step_2/{geneID}-M0HKY-{omega}-{draw}-{repID}.chain", geneID=GENEID, omega=OMEGA, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.sh",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.conf",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.simu",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # CODEML M0HKY
        expand(ROOT_dir+"/outputs/codeml_M0HKY/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/codeml_M0HKY/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/codeml_M0HKY/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # SIMU M7HKY
        expand(ROOT_dir+"/outputs/step_4/{geneID}-M7HKY-{omega_mix}-{draw}-{repID}.chain",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.sh",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.conf",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.simu",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0.fasta",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        # CODEML M7HKY
        expand(ROOT_dir+"/outputs/codeml_M7HKY/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0/codeml.sh",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/codeml_M7HKY/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0/codeml.sh",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/codeml_M7HKY/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0/codeml",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        # CODEML M8HKY
        expand(ROOT_dir+"/outputs/codeml_M8HKY/{geneID}-M8HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0/codeml.sh",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/codeml_M8HKY/{geneID}-M8HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0/codeml.conf",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/codeml_M8HKY/{geneID}-M8HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0/codeml",geneID=GENEID, omega_mix=[i for i,j in OMEGA_MIX.items()], CpG=CPG, draw=DRAW, repID=REPID),
        # CABC M0HKY
        # expand(ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.sh",geneID=GENEID, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.conf",geneID=GENEID, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.simu",geneID=GENEID, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_chainID.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_params.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_ss.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_true_ss.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-model_preprocessing_params.pickle",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-cabc.sh",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-knn.feather",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-post.tsv",geneID=GENEID, omega=OMEGA, CpG=CPG, draw=DRAW, repID=REPID)

rule fasta2phylip:
    input: ROOT_dir+"/data/{geneID}.fasta"
    output: ROOT_dir+"/outputs/data_prep/{geneID}.phylip"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/data_prep.py \
        --fasta={input} \
        --phylip={output}
        """


rule generate_pb_mpi_cmd:
    input:
        phylip=ROOT_dir+"/outputs/data_prep/{geneID}.phylip",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    output:ROOT_dir+"/outputs/step_1/{geneID}-M0GTR-{repID}.sh"
    params:
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/generate_pb_mpi_cmd.py \
        --phylip={input.phylip} \
        --tree={input.tree} \
        --output={output} \
        --local={params.local} \
        --docker={params.docker}
        """

rule run_pb_mpi_cmd:
    input: ROOT_dir+"/outputs/step_1/{geneID}-M0GTR-{repID}.sh"
    output: ROOT_dir+"/outputs/step_1/{geneID}-M0GTR-{repID}.chain"
    shell:
        """
        echo bash {input}
        bash {input}
        """


rule generate_M0HKY_fixed:
    input:
        phylip=ROOT_dir+"/outputs/data_prep/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/step_1/{geneID}-M0GTR-{repID}.chain",
    output:ROOT_dir+"/outputs/step_2/{geneID}-M0HKY-{omega}-{draw}-{repID}.chain"
    params:
        burnin=600,
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/generate_M0HKY_fixed.py \
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --burnin={params.burnin} \
        --omega={wildcards.omega} \
        --output={output}
        """
rule generate_M0HKY_simu_with_fixed_values_cmd:
    input:
        phylip=ROOT_dir+"/outputs/data_prep/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/step_2/{geneID}-M0HKY-{omega}-{draw}-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.sh",
        conf=ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.conf"
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

rule run_M0HKY_simu_with_fixed_values_cmd:
    input: ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.sh"
    output:
        simu=ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.simu",
        fasta=ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta"
    shell:
        """
        echo bash {input}
        bash {input}
        """


rule generate_M7HKY_mixed:
    input:
        phylip=ROOT_dir+"/outputs/data_prep/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/step_1/{geneID}-M0GTR-{repID}.chain",
    output:ROOT_dir+"/outputs/step_4/{geneID}-M7HKY-{omega_mix}-{draw}-{repID}.chain"
    params:
        burnin=600,
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/generate_M7HKY_mixed.py \
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --burnin={params.burnin} \
        --omega_mix={wildcards.omega_mix} \
        --output={output}
        """
rule generate_M7HKY_simu_with_mixed_values_cmd:
    input:
        phylip=ROOT_dir+"/outputs/data_prep/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/step_4/{geneID}-M7HKY-{omega_mix}-{draw}-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.sh",
        conf=ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.conf"
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

rule run_M7HKY_simu_with_mixed_values_cmd:
    input: ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.sh"
    output:
        simu=ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.simu",
        fasta=ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0.fasta"
    shell:
        """
        echo bash {input}
        bash {input}
        """


rule generate_codeml_M0HKY_cmd:
    input:
        fasta=ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        yaml=ROOT_dir+"/configs/codeml.yaml"
    output:
        sh=ROOT_dir+"/outputs/codeml_M0HKY/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml_M0HKY/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
    params:
        modelID="M0HKY",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/generate_codeml_cmd.py \
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
        fasta=ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        yaml=ROOT_dir+"/configs/codeml.yaml"
    output:
        sh=ROOT_dir+"/outputs/codeml_M7HKY/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml_M7HKY/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
    params:
        modelID="M7HKY",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/generate_codeml_cmd.py \
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
        fasta=ROOT_dir+"/outputs/step_5/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0.fasta",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        yaml=ROOT_dir+"/configs/codeml.yaml"
    output:
        sh=ROOT_dir+"/outputs/codeml_M8HKY/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh",
        conf=ROOT_dir+"/outputs/codeml_M8HKY/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.conf"
    params:
        modelID="M8HKY",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/generate_codeml_cmd.py \
        --fasta={input.fasta} \
        --tree={input.tree} \
        --yaml={input.yaml} \
        --modelID={params.modelID} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker} \
        """


rule run_codeml_M0HKY_cmd:
    input:ROOT_dir+"/outputs/codeml_M0HKY/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:ROOT_dir+"/outputs/codeml_M0HKY/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule run_codeml_M7HKY_cmd:
    input:ROOT_dir+"/outputs/codeml_M7HKY/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:ROOT_dir+"/outputs/codeml_M7HKY/{geneID}-M7HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule run_codeml_M8HKY_cmd:
    input:ROOT_dir+"/outputs/codeml_M8HKY/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml.sh"
    output:ROOT_dir+"/outputs/codeml_M8HKY/{geneID}-M8HKY-{omega}-{CpG}-{draw}-{repID}-1_0/codeml"
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



# rule generate_files_r_abc:
#     input:
#         simu=ROOT_dir+"/outputs/cabc_step_1/{geneID}-M0HKY-{repID}.simu",
#         true_ss=ROOT_dir+"/outputs/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.simu"
#     output:
#         chainID=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_chainID.feather",
#         preprocessed_params=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_params.feather",
#         preprocessed_ss=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_ss.feather",
#         true_ss=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{draw}-{CpG}-{repID}-df_true_ss.feather",
#         model_preprocessing_params_file=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-model_preprocessing_params.pickle",
#     params:
#         preprocessing_ss="log2",
#         preprocessing_params="log",
#         knn=1000,
#         params="'root lambda_CpG lambda_TBL lambda_omega nucsA nucsC nucsG nucsT nucrrAC nucrrAG'",
#         ss="'pwAC pwAG pwAT pwCG pwCT pwGT pwAA nuc3A nuc3C nuc3G nuc3T dinuc31CA dinuc31CG dinuc31TG'",
#         output_dir=ROOT_dir+"/outputs/cabc_step_2/",
#         prefix="{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}",
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 scripts/generate_files_r_abc.py \
#         --simu={input.simu} \
#         --true_ss_file={input.true_ss} \
#         --knn={params.knn} \
#         --params={params.params} \
#         --ss={params.ss} \
#         --output_dir={params.output_dir} \
#         --trans_fct_params={params.preprocessing_params} \
#         --trans_fct_ss={params.preprocessing_ss}  \
#         --prefix={params.prefix}
#         """


# rule generate_cabc_cmd:
#     input:
#         chainID=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_chainID.feather",
#         preprocessed_params=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_params.feather",
#         preprocessed_ss=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_ss.feather",
#         true_ss=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_true_ss.feather",
#     output: ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-cabc.sh"
#     params:
#         local=ROOT_dir,
#         docker="/data"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 scripts/generate_cabc_cmd.py \
#         --df_true_ss={input.true_ss} \
#         --df_simu_space_knn_params={input.preprocessed_params} \
#         --df_simu_space_knn_ss={input.preprocessed_ss} \
#         --df_simu_space_knn_chainID={input.chainID} \
#         --output={output} \
#         --local={params.local} \
#         --docker={params.docker}
#         """
# rule run_abc_cmd:
#     input: ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-cabc.sh"
#     output: ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-knn.feather"
#     params: "/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-knn.feather"
#     shell:
#         """
#         echo bash {input} {params}
#         bash {input} {params}
#         """

# rule generate_final_r_abc:
#     input:
#         knn_file=ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-knn.feather",
#         model_preprocessing_params_file=ROOT_dir+"/outputs/cabc_step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-model_preprocessing_params.pickle"
#     output: ROOT_dir+"/outputs/cabc_step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-post.tsv"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 scripts/generate_final_r_abc.py \
#         --knn_file={input.knn_file} \
#         --model_preprocessing_params_file={input.model_preprocessing_params_file} \
#         --output={output} \
#         """

# rule generate_figure_1:
# rule generate_figure_2:
# rule generate_tables:
