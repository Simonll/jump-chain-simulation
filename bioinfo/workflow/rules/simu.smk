import os
import sys
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).parent.__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"

ruleorder: fasta2phylip > generate_pb_mpi_cmd  > generate_M0HKY_fixed > generate_M0HKY_simu_with_fixed_values_cmd > run_M0HKY_simu_with_fixed_values_cmd > generate_M7HKY_mixed > generate_M7HKY_simu_with_mixed_values_cmd > run_M7HKY_simu_with_mixed_values_cmd

rule fasta2phylip:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        fasta=ROOT_dir+"/data/{geneID}.fasta"
    output:
        ROOT_dir+"/outputs/data/{geneID}.phylip"
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
    output:
        docker=ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.sh",
        cluster=ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}-cluster.sh"
    params:
        local=ROOT_dir,
        docker="/data",
        phylip="/outputs/data/{geneID}.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        cluster="/home/sll/jump-chain-simulation/bioinfo/workflow",
        chainname="{geneID}-M0GTR-{repID}",
        np="4",
        output="{geneID}-M0GTR-{repID}-cluster.sh",
        model="'-mutsel -catfix uniform -freeomega'"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={input.phylip} \
        --tree={input.tree} \
        --model={params.model} \
        --output={output.docker} \
        --local={params.local} \
        --docker={params.docker}
        echo #####
        echo "#qsub -cwd -pe mpi {params.np} {params.output}" >> {output.cluster}
        echo "LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH" >> {output.cluster}
        echo "CXX=`which g++`" >> {output.cluster}
        echo "CC=`which gcc`" >> {output.cluster}
        echo "export LD_LIBRARY_PATH" >> {output.cluster}
        echo "export CC" >> {output.cluster}
        echo "export CXX" >> {output.cluster}
        echo mpirun -np {params.np} /home/sll/pbmpi/data/pb_mpi -d {params.cluster}{params.phylip} -T {params.cluster}{params.tree} {params.model} -s -x 1 2000 {params.chainname} >> {output.cluster}
        """


# rule run_pb_mpi_cmd:
#     input: ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.sh"
#     output: ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.chain"
#     shell:
#         """
#         echo bash {input}
#         bash {input}
#         """


rule generate_M0HKY_fixed:
    input:
        script=ROOT_dir+"/scripts/generate_M0HKY_fixed.py",
        phylip=ROOT_dir+"/outputs/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.chain",
    output:
        ROOT_dir+"/outputs/modified_chain_m0hky/{geneID}-M0HKY-{omega}-{draw}-{repID}.chain"
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
    input:
        ROOT_dir+"/outputs/simu_m0hky/{geneID}-M0HKY-{omega}-{CpG}-{draw}-{repID}.sh"
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
        mcmc=ROOT_dir+"/outputs/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.chain",
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
    input:
        ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.sh"
    output:
        simu=ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}.simu",
        fasta=ROOT_dir+"/outputs/simu_m7hky/{geneID}-M7HKY-{omega_mix}-{CpG}-{draw}-{repID}-1_0.fasta"
    shell:
        """
        echo bash {input}
        bash {input}
        """
