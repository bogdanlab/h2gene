import numpy as np
import pandas as pd
import os
from os.path import join
import itertools

# constants
# contain the plink and also the pre-computed mean + std
CHR_I = 1
PLINK_PATH = f"/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/genotype/raw/chr{CHR_I}"
LD_PATH = "/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/ld/"
PARTITION_PATH = f"/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/partition/fourier_ls-chr{CHR_I}.bed"
LDSCORE_PATH = f"/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/ldscore/chr{CHR_I}.l2.ldscore.gz"
GENE_LIST_PATH = "/u/project/pasaniuc/kangchen/h2gene/data/gene_list.bed"

N_SIM=30

SRC_DIR = "~/project-pasaniuc/h2gene/src"

with open(PARTITION_PATH, 'r') as f:
    N_PARTITION = len(f.readlines()) - 1
    
# {n_causal_gene}_{background_causal_prob}_{n_causal_snps_body}_{n_causal_snps_tss}
sim_prefix_list = ['_'.join([str(p) for p in param]) for param in itertools.product([20, 50, 100], [0.001, 0.01], [5], [3])]
sim_prefix_list.extend(['_'.join([str(p) for p in param]) for param in itertools.product([20, 50, 100], [0.001, 0.01], [10], [6])])
# sim_prefix_list = ["20_0.001_5_3"]
rule all:
    input:
        "out/data/gene_list.bed",
        expand("out/sim/gwas/{sim_prefix}/beta_hat.npy", sim_prefix=sim_prefix_list),
        expand("out/sim/gwas/{sim_prefix}/sim_0_par_0.tsv.gz", sim_prefix=sim_prefix_list),
        # expand("out/estimate/{sim_prefix}/sim_{sim_i}_par_{par_i}.rds", sim_prefix=sim_prefix_list, sim_i=np.arange(N_SIM), par_i=np.arange(N_PARTITION)),
        # expand("out/estimate/{sim_prefix}/summary.tsv", sim_prefix=sim_prefix_list)
        expand("out/estimate_hess/{sim_prefix}/summary.tsv", sim_prefix=sim_prefix_list)
            
# prepare genotype / protein coding gene list
rule prepare_data:
    resources:
        mem_gb=16,
        time_min=60
    input:
    output:
        gene_list="out/data/gene_list.bed",
        partition="out/data/partition.bed"
    params:
        log=lambda wildcards: "log/prepare_data.log",
        min_gene_num_snps=25
    run:
        # process partition file
        partition = pd.read_csv(PARTITION_PATH, delim_whitespace=True)
        partition.columns = ['CHR', 'START', 'STOP']
        partition['CHR'] = partition['CHR'].apply(lambda x : int(x[3:]))
        partition.to_csv(output.partition, sep='\t', index=False)
        
        # process gene_list
        bim = pd.read_csv(PLINK_PATH + '.bim', header=None, delim_whitespace=True, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
        gene_list = pd.read_csv(GENE_LIST_PATH, header=None, 
                                names=['CHR', 'START', 'STOP', 'STRAND', 'NAME'],
                                delim_whitespace=True)
        gene_list = gene_list[gene_list.CHR == str(CHR_I)].reset_index(drop=True)
        
        # for each gene, filter out genes overlapping with boundary / no SNPs
        filtered_index = []
        for par_i, par in partition.iterrows():
            par_gene_list = gene_list.loc[(par.START <= gene_list.START) & (gene_list.STOP < par.STOP), :]
            for gene_i, gene in par_gene_list.iterrows():
                if sum((gene.START <= bim.BP) & (bim.BP < gene.STOP)) > params.min_gene_num_snps:
                    filtered_index.append(gene_i)
        
        gene_list.iloc[filtered_index, :].to_csv(output.gene_list, sep='\t', index=False)
        
        
rule simulate:
    resources:
        mem_gb=16,
        time_min=600
    input:
        bed=PLINK_PATH + ".bed",
        gene_list="out/data/gene_list.bed",
        mean_std=PLINK_PATH + "_mean_std.txt"
    output:
        beta_hat="out/sim/gwas/{sim_prefix}/beta_hat.npy",
    params:
        log=lambda wildcards: "log/simulate.log",
        input_bfile=lambda wildcards, input: input.bed[:input.bed.rfind('.')],
        seed=1234,
        total_h2=0.05,
        body_h2=0.03,
        tss_h2=0.01,
    run:
        num_causal_genes, background_causal_prob, num_causal_snps_body, num_causal_snps_tss = wildcards.sim_prefix.split('_')
        out_dir = os.path.dirname(output.beta_hat)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        from simulate import simulate, simulate2
        from pysnptools.snpreader import Bed
        
        geno = Bed(params.input_bfile, count_A1=False)
        bim = pd.read_csv(params.input_bfile + '.bim', header=None, 
                delim_whitespace=True, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
        mean_std = pd.read_csv(input.mean_std, delim_whitespace=True)
        gene_list = pd.read_csv(input.gene_list, delim_whitespace=True)
        
        np.random.seed(params.seed)
        config = {
            'total_h2': params.total_h2,
            'body_h2': params.body_h2,
            'tss_h2': params.tss_h2,
            'num_causal_genes': int(num_causal_genes),
            'background_causal_prob': float(background_causal_prob),
            'num_causal_snps_body': int(num_causal_snps_body),
            'num_causal_snps_tss': int(num_causal_snps_tss)
        }
        print(config)
        data = simulate2(geno=geno, bim=bim, mean_std=mean_std, gene_list=gene_list, config=config)
        
        for d in ["beta", "phe", "beta_hat"]:
            np.save(os.path.join(out_dir, '{}.npy'.format(d)), data[d])
        data["causal_genes"].to_csv(os.path.join(out_dir, 'causal_genes.csv'), index=False)
            

"""
Convert npy to GWAS format
"""
rule format_gwas:
    resources:
        mem_gb=12,
        time_min=20
    input:
        bed=PLINK_PATH + ".bed",
        mean_std=PLINK_PATH + "_mean_std.txt",
        beta_hat="out/sim/gwas/{sim_prefix}/beta_hat.npy",
        ldscore=LDSCORE_PATH,
        partition="out/data/partition.bed"
    output:
        "out/sim/gwas/{sim_prefix}/sim_0_par_0.tsv.gz",
    params:
        log="log/misc.log",
        input_bfile=lambda wildcards, input: input.bed[:input.bed.rfind('.')],
    run:
        from pysnptools.snpreader import Bed
        
        partition = pd.read_csv(input.partition, delim_whitespace=True)
        beta_hat = np.load(input.beta_hat)
        n_sim = beta_hat.shape[1]
        
        out_dir = os.path.dirname(input.beta_hat)
        
        geno = Bed(params.input_bfile, count_A1=False)
        bim = pd.read_csv(params.input_bfile + '.bim', header=None, 
                delim_whitespace=True, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
        ldscore = pd.read_csv(input.ldscore, delim_whitespace=True)
        mean_std = pd.read_csv(input.mean_std, delim_whitespace=True)
        
        n_indiv = geno.shape[0]
        bim["N"] = n_indiv
        bim["MAF"] = 1. - mean_std['mean'].values / 2
        bim["LDSCORE"] = ldscore['L2']
        for sim_i in range(n_sim):
            assoc = bim.copy()
            assoc["Z"] = np.sqrt(n_indiv) * beta_hat[:, sim_i]
            
            for par_i, par in partition.iterrows():
                par_snps = np.where((par.CHR == assoc.CHR.values) & \
                                    (par.START <= assoc.BP.values) & \
                                    (assoc.BP.values < par.STOP))[0]
                filename=join(out_dir, f"sim_{sim_i}_par_{par_i}.tsv.gz")
                assoc.iloc[par_snps, ].to_csv(filename, sep='\t', index=False, float_format='%.6f')

# estimate for each LD region
rule estimate:
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 3 + 4,
        time_min=lambda wildcards, attempt: attempt * 20 if attempt <= 2 else attempt * 60
    input:
        gene_list="out/data/gene_list.bed",
        sumstats="out/sim/gwas/{sim_prefix}/sim_{sim_i}_par_{par_i}.tsv.gz"
    output:
        estimate="out/estimate/{sim_prefix}/sim_{sim_i}_par_{par_i}.rds",
    params:
        Rscript="~/project-pasaniuc/software/anaconda3/envs/r_env/bin/Rscript",
        log="log/misc.log",
        out_dir=lambda wildcards, output: output.estimate[:output.estimate.rfind('/')]
    shell:
        """

        mkdir -p {params.out_dir}

        . /u/local/Modules/default/init/modules.sh
        module load gcc/6.3.0
        
        {params.Rscript} {SRC_DIR}/h2gene_cli.R \
            --ld_prefix {LD_PATH}/{CHR_I}/par_{wildcards.par_i} \
            --gene_list {input.gene_list} \
            --sumstats {input.sumstats} \
            --out {output.estimate}
        """

rule summary:
    resources:
        mem_gb=16,
        time_min=lambda wildcards, attempt: attempt * 60
    input:
        gene_list="out/data/gene_list.bed",
        partition="out/data/partition.bed"
    output:
        summary="out/estimate/{sim_prefix}/summary.tsv"
    params:
        log="log/summary.log",
        Rscript="~/project-pasaniuc/software/anaconda3/envs/r_env/bin/Rscript",
        estimate_dir="out/estimate",
        sim_dir="out/sim/gwas",
        out_prefix=lambda wildcards, output: output.summary[:output.summary.rfind('.')]
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/6.3.0
        {params.Rscript} sim_summary.R \
            --sim_dir {params.sim_dir}/{wildcards.sim_prefix} \
            --ld_dir {LD_PATH}/{CHR_I} \
            --bim {PLINK_PATH}.bim \
            --estimate_dir {params.estimate_dir}/{wildcards.sim_prefix} \
            --partition {input.partition} \
            --gene_list {input.gene_list} \
            --out_prefix {params.out_prefix}
        """

# # benchmark other ways of estimating gene-level heritability
rule estimate_hess:
    resources:
        mem_gb=12,
        time_min=lambda wildcards, attempt: attempt * 30
    input:
        bed=PLINK_PATH + ".bed",
        gene_list="out/data/gene_list.bed",
        partition="out/data/partition.bed",
        ld=expand(join(LD_PATH, str(CHR_I), "par_{par_i}.npy"), par_i=np.arange(N_PARTITION))
    output:
        estimate="out/estimate_hess/{sim_prefix}/summary.tsv"
    params:
        python="~/project-pasaniuc/software/anaconda3/bin/python",
        bfile = lambda wildcards, input: input.bed[:input.bed.rfind('.')],
        log="log/sim_estimate.log",
        ld_dir=lambda wildcards, input: input.ld[0][:input.ld[0].rfind('/')],
        out_dir=lambda wildcards, output: output.estimate[:output.estimate.rfind('/')]
    shell:
        """
        mkdir -p {params.out_dir}
        
        . /u/local/Modules/default/init/modules.sh
        module load gcc/6.3.0

        {params.python} estimate_compare.py hess_cli \
            --ld_dir {params.ld_dir} \
            --gene_list {input.gene_list} \
            --partition {input.partition} \
            --sumstats_dir out/sim/gwas/{wildcards.sim_prefix} \
            --n_sim {N_SIM} \
            --out {output.estimate}
        """

