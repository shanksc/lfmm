import sys
import os 
from glob import glob


#CHANGE TO INTERGRATE SLURM DIRECTLY 
'''
Create input for LFMM2 and run 
'''

'''
Config example 

pop_labels:
    pop1: eur_0_athletes
    pop2: eur_0_non_athletes
    ... add more 

scripts_dir

tsv: ../working/
#make sure pop labels match tsv locations ie working/pop_label/temp/

#make same as env
name: eur_0

working: eur_0
'''

'''
OVERVIEW

-Split variants into chunks (10k variants or so), or we can just partition each chr into N chunks 
-Run latent factor model on each chromosome, save object and pass to each respective chunk 
'''

print(config)

#since we're io limited , split into groups of variants, can we save fitted latent factors? 
N_VARIANTS=10000
SUFFIX_LENGTH=3 

#wildcard_constraints:
    #chrom = 

rule all:
    input:
        config['bcf'],
        config['working'],
        f"{config['working']}/{config['name']}/env.txt",
    
        #expand("{prefix}/{name}/{chrom}.pruned.lfmm", prefix=config['working'], name=config['name'], chrom=range=(1,23))
        #f"{config['working']}/{config['name']}/sample_order.txt"
        
        #run this alone prior to required inputs below
        #expand("{prefix}/{name}/{chrom}.ready.txt", prefix=config['working'], name=config['name'], chrom=range(1,23))
        
        #"{prefix}/{name}/chunks/{chrom}", prefix=config['working'], name=config['name'], chrom=range(1,23)),
        #"{prefix}/{name}/{chrom}.pruned.lfmm", prefix=config['working'], name=config['name'], chrom=range(1,23)),
        #"{prefix}/{name}/{chrom}.rds", prefix=config['working'], name=config['name'], chrom=range(1,23)),
        
        #expand("{prefix}/{name}/{chrom}.ready.txt", prefix=config['working'], name=config['name'], chrom=range(1,23))
        #expand("{prefix}/{name}/{chrom}.rds", prefix=config['working'], name=config['name'], chrom=range=(1,23))
        
	expand("{prefix}/{name}/{chrom}.done.txt", prefix=config['working'], name=config['name'], chrom=range(1,23))    
        #expand("{prefix}/{name}/{chrom}.rds", prefix=config['working'], name=config['name'], chrom=range(1,23))
rule get_sample_order:
    input:
        config['bcf']
    output:
	    f"{config['working']}/{config['name']}/sample_order.txt"
    shell:
	    """
	    bcftools query -l {config[bcf]} > {config[working]}/{config[name]}/sample_order.txt 
	    """

#create info file for each chunk of snps 
checkpoint create_chunk_info:
    input:
        config['bcf']
    output:
        #directory("{prefix}/{name}/chunks/{chrom}"),
        #directory("{prefix}/{name}/lfmm/{chrom}")
        directory(f"{config['working']}/{config['name']}/chunks/" + "{chrom, [0-9]+}"),
        directory(f"{config['working']}/{config['name']}/lfmm/" + "{chrom, [0-9]+}"),
        directory(f"{config['working']}/{config['name']}/pvals/" + "{chrom, [0-9]+}")
        #directory("{prefix}/{name}/")
        #directory("{prefix}/{name}/chunks/{chrom}")
    shell:
        """
        mkdir {config[working]}/{config[name]}/chunks/{wildcards.chrom}
        mkdir {config[working]}/{config[name]}/lfmm/{wildcards.chrom}
        mkdir {config[working]}/{config[name]}/pvals/{wildcards.chrom}
        bcftools query -r chr{wildcards.chrom} -f'%CHROM\t%POS\n' {config[bcf]} | split -l {N_VARIANTS} --additional-suffix .txt -d -a {SUFFIX_LENGTH} - {config[working]}/{config[name]}/chunks/{wildcards.chrom}/
        """

#then convert each chunk to lfmm
'''
rule convert_chunks_to_lfmm:
'''
rule convert_pruned_bcf_to_lfmm:
    input:
        config['pruned_bcf']
        #f"{config['working']}/{config['name']}/sample_order.txt.missing"
    output:
        "{prefix}/{name}/{chrom}.pruned.lfmm"
    shell:
        """
        bash {config[scripts_dir]}/convert_bcf_to_lfmm.sh {config[pruned_bcf]} {wildcards.chrom} {wildcards.prefix}/{wildcards.name}/{wildcards.chrom}.pruned.lfmm {config[scripts_dir]}
        """

#use pruned lfmm files to fit latent factors for each chrom 
#pass aditional k since we'll need to run for 3-4 different k's most likely 
rule fit_lfmm_per_chrom:
    input:
        "{prefix}/{name}/{chrom}.pruned.lfmm"
    output:
        "{prefix}/{name}/{chrom}.rds"
    shell:
        "Rscript fit_lfmm.R {wildcards.prefix}/{wildcards.name}/{wildcards.chrom}.pruned.lfmm {config[working]}/{config[name]}/env.txt {output}"


rule convert_bcf_to_lfmm:
    input:
        #"{prefix}/{name}/chunks/{chrom}/{chunk}.txt"
        #f"{config}['working']/{config['name']}/" + "{chrom}.ready.txt",
        f"{config['working']}/{config['name']}/chunks/" + "{chrom}/{chunk}.txt"
    output:
        #lfmm works and not lfmms? 
        f"{config['working']}/{config['name']}/lfmm/" + "{chrom}/{chunk}.lfmm"
    shell:
        """
        bash {config[scripts_dir]}/convert_bcf_to_lfmm_chunk.sh {config[bcf]} {input} {output} {config[scripts_dir]}
        """

#calc pvals and combine with chunk.txt into a tsv file containng coordinates and pvals 
rule calc_pvals:
    input:
        #f"{config['working']}/{config['name']}/" + "{chrom}.rds",
        #f"{config}['working']/{config['name']}/" + "{chrom}.ready.txt" #should this be here? how does this work with the checkpoint
        f"{config['working']}/{config['name']}/lfmm/" + "{chrom}/{chunk}.lfmm"
    output:
        f"{config['working']}/{config['name']}/pvals/" + "{chrom}/{chunk, [0-9]+}.tsv"
    shell:
        """
        Rscript run_lfmm.R {config[working]}/{config[name]}/{wildcards.chrom}.rds {input} {config[working]}/{config[name]}/env.txt {config[working]}/{config[name]}/lfmm/{wildcards.chrom}/{wildcards.chunk}.pvals
        paste {config[working]}/{config[name]}/chunks/{wildcards.chrom}/{wildcards.chunk}.txt {config[working]}/{config[name]}/lfmm/{wildcards.chrom}/{wildcards.chunk}.pvals > {output}
        """

def aggregate_input(wildcards):
    #print(wildcards)
    checkpoint_output = checkpoints.create_chunk_info.get(**wildcards).output[0]
    #print(checkpoint_output)
    chunk_files = [fname.split('/')[-1].split('.')[0] for fname in glob(f"{checkpoint_output}/*.txt")]
    #chunk_files = [fname.split('')]
    #print(checkpoint_output, chunk_files[0])
    #sys.exit()
    #print(list(expand("{prefix}/{name}/lfmm/{chrom}/{chunk}.lfmm", prefix=wildcards.prefix, name=wildcards.name, chrom=wildcards.chrom, chunk=chunk_files))[:3])
    #return expand("{prefix}/{name}/lfmm/{chrom}/{chunk}.lfmm", prefix=wildcards.prefix, name=wildcards.name, chrom=wildcards.chrom, chunk=chunk_files)
    #print(expand(f"{config['working']}/{config['name']}/lfmm/" + "{chrom}/{chunk}.lfmm", chrom=wildcards.chrom, chunk=chunk_files)[:3])
    #return expand(f"{config['working']}/{config['name']}/lfmm/" + "{chrom}/{chunk}.lfmm", chrom=wildcards.chrom, chunk=chunk_files)
    print(expand(f"{config['working']}/{config['name']}/pvals/" + "{chrom}/{chunk}.tsv", chrom=wildcards.chrom, chunk=chunk_files)[:3])
    return expand(f"{config['working']}/{config['name']}/pvals/" + "{chrom}/{chunk}.tsv", chrom=wildcards.chrom, chunk=chunk_files)


rule get_chunks:
    input:
        "{prefix}/{name}/chunks/{chrom}",
        "{prefix}/{name}/{chrom}.pruned.lfmm",
        "{prefix}/{name}/{chrom}.rds",
    output:
        "{prefix}/{name}/{chrom}.ready.txt"
    shell:
        "touch {output}"

rule get_pvals:
    input:
        "{prefix}/{name}/{chrom}.ready.txt",
        aggregate_input
    output:
        "{prefix}/{name}/{chrom}.done.txt"
    shell:
        "touch {output}"

#rule calc_all_pvals:
#    input:
        

#make this more effecient by globbing all the chunks and then applying a rule to each chunk
''' 
rule convert_chunks_to_lfmm_by_chrom:
    input:
        chunks = "{prefix}/{name}/chunks/{chrom}",
    output:
        "{prefix}/{name}/{chrom}.ready.txt"
    shell:
        """
        for filename in {input.chunks}/*.txt; do
            output_path="$filename"
            bash {config[scripts_dir]}/convert_bcf_to_lfmm_chunk.sh {config[bcf]} $filename $output_path.lfmm {config[scripts_dir]}
        done
        touch {output}
        """
'''




#rule calc_p_vals:
#    input:

"""
rule get_lfmm_format_per_chrom:
    input:
        aggregate_input
    output:
        "{prefix}/{name}/{chrom}.ready.txt"
"""       
    

'''
rule run_lfmm_chunks_per_chrom:
    input:
        chunks = "{prefix}/{name}/chunks/{chrom}",
        "{prefix}/{name}/chunks/{chrom}.ready.txt",
        rds = "{prefix}/{name}/{chrom}.rds"
    output:
        "{prefix}/{name}/{chrom}.lfmm.done.txt"
    shell:
        """
        for filename in {input.chunks}/*.txt.lfmm; do
            output_path="$filename"
            Rscript run_lfmm.R {input.rds} $filename {config[working]}/{config[name]}/env.txt $output_path.lfmm.pvals 
        done
        touch {output}
        """
'''

#for each chrom, load the lfmm object and compute p-values 
'''
rule run_lfmm_per_chunk:
    input:

    output:
    
    shell:
'''

##################################################
# Just use smartpca for tracy-widom test for now #
##################################################
#run sMNF to determine K latent factors
#rule determine_k:

#rule run_tracy_widom:
    #input:

rule create_lfmm_env:
    input:
        config['bcf'],
        config['metadata']
    output:
        f"{config['working']}/{config['name']}/env.txt",
        #f"{config['working']}/{config['name']}/sample_order.txt.missing"
    shell:
        """
        bcftools query -l {config[bcf]} > {config[working]}/{config[name]}/sample_order.txt
        python3 {config[scripts_dir]}/create_lfmm_env.py {config[metadata]} {config[working]}/{config[name]}/sample_order.txt > {config[working]}/{config[name]}/env.txt
        """

#make a rule so that we run this several times with different K latent factors
'''
rule run_lfmm_per_chrom:
    input:
        f"{config['working']}/{config['name']}/env.txt",
        "{prefix}/{name}/{chrom}.lfmm",
        "{prefix}/{name}/{chrom}.pruned.lfmm"
    output:
        "{prefix}/{name}/{chrom}_{part}_p_values.txt"
    shell:
        """
        Rscript run_lfmm.R {wildcards.prefix}/{wildcards.name}/{wildcards.chrom}.lfmm {wildcards.prefix}/{wildcards.name}/{wildcards.chrom}.pruned.lfmm {config[working]}/{config[name]}/env.txt {wildcards.prefix}/{wildcards.name}/{wildcards.chrom}_p_values.txt
        """
'''

'''
rule create_tsv:
    input:
        expand("{prefix}/{name}/{chrom}_{part}_p_values.txt", prefix=config['working'], name=config['name'], chrom=range(1,23))
    output:
        f"{config['working']}/{config['name']}/lfmm.tsv"
    shell:
        """
        touch {output}
        for chrom in $(seq 1 22); do
		paste {config[working]}/{config[name]}/{{$chrom}}_snp_info.tsv {config[working]}/{config[name]}/{{$chrom}}_p_values.txt >> {output}
        done
        """
'''
