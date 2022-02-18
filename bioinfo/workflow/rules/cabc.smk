import os
import sys
from pathlib import Path

# ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).__str__()
# workdir: ROOT_dir
# if ROOT_dir not in sys.path:
#     sys.path.append(ROOT_dir)

# configfile: ROOT_dir + "/configs/configs.yaml"

ruleorder: generate_M0HKY_simu_cmd > run_M0HKY_simu_cmd > generate_files_r_abc > generate_cabc_cmd > run_abc_cmd > generate_final_r_abc

rule generate_M0HKY_simu_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_M0M7HKY_simu_cmd.py",
        phylip=ROOT_dir+"/outputs/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.chain",
    output:
        sh=ROOT_dir+"/outputs/cabc/step_1/{geneID}-M0HKY-{repID}.sh",
        conf=ROOT_dir+"/outputs/cabc/step_1/{geneID}-M0HKY-{repID}.conf"
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

rule run_M0HKY_simu_cmd:
    input:
        ROOT_dir+"/outputs/cabc/step_1/{geneID}-M0HKY-{repID}.sh"
    output:
        ROOT_dir+"/outputs/cabc/step_1/{geneID}-M0HKY-{repID}.simu"
    shell:
        """
        echo bash {input}
        bash {input}
        """



rule generate_files_r_abc:
    input:
        script=ROOT_dir+"/scripts/generate_files_r_abc.py",
        simu=ROOT_dir+"/outputs/cabc/step_1/{geneID}-M0HKY-{repID}.simu",
        true_ss=ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.simu"
    output:
        chainID=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_chainID.feather",
        preprocessed_params=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_params.feather",
        preprocessed_ss=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_ss.feather",
        true_ss=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{draw}-{CpG}-{repID}-df_true_ss.feather",
        model_preprocessing_params_file=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-model_preprocessing_params.pickle",
    params:
        preprocessing_ss="log2",
        preprocessing_params="log",
        knn=1000,
        params="'root lambda_CpG lambda_TBL lambda_omega nucsA nucsC nucsG nucsT nucrrAC nucrrAG'",
        ss="'pwAC pwAG pwAT pwCG pwCT pwGT pwAA nuc3A nuc3C nuc3G nuc3T dinuc31CA dinuc31CG dinuc31TG'",
        output_dir=ROOT_dir+"/outputs/cabc/step_2/",
        prefix="{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}",
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --simu={input.simu} \
        --true_ss_file={input.true_ss} \
        --knn={params.knn} \
        --params={params.params} \
        --ss={params.ss} \
        --output_dir={params.output_dir} \
        --trans_fct_params={params.preprocessing_params} \
        --trans_fct_ss={params.preprocessing_ss}  \
        --prefix={params.prefix}
        """


rule generate_cabc_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_cabc_cmd.py",
        chainID=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_chainID.feather",
        preprocessed_params=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_params.feather",
        preprocessed_ss=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_simu_space_knn_ss.feather",
        true_ss=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-df_true_ss.feather",
    output:
        ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-cabc.sh"
    params:
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --df_true_ss={input.true_ss} \
        --df_simu_space_knn_params={input.preprocessed_params} \
        --df_simu_space_knn_ss={input.preprocessed_ss} \
        --df_simu_space_knn_chainID={input.chainID} \
        --output={output} \
        --local={params.local} \
        --docker={params.docker}
        """

rule run_abc_cmd:
    input:
        ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-cabc.sh"
    output:
        ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-knn.feather"
    params:
        "/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-knn.feather"
    shell:
        """
        echo bash {input} {params}
        bash {input} {params}
        """

rule generate_final_r_abc:
    input:
        script=ROOT_dir+"/scripts/generate_final_r_abc.py",
        knn_file=ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-knn.feather",
        model_preprocessing_params_file=ROOT_dir+"/outputs/cabc/step_2/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-model_preprocessing_params.pickle"
    output:
        ROOT_dir+"/outputs/cabc/step_3/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}-post.tsv"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --knn_file={input.knn_file} \
        --model_preprocessing_params_file={input.model_preprocessing_params_file} \
        --output={output} \
        """
