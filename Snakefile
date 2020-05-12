import os
import shutil
import gzip
import yaml
from re import match
from copy import deepcopy
import subprocess
import pandas as pd
from snakemake.utils import validate

def dict_merge(a, b):
    """
    Deep merge 2 dicts together
    """
    if not isinstance(b, dict):
        return b
    result = deepcopy(a)
    for k, v in b.items():
        if k in result and isinstance(result[k], dict):
            result[k] = dict_merge(result[k], v)
        else:
            result[k] = deepcopy(v)
    return result
	
def get_sample_perBatch(wildcards,prefix,suffix):
    return prefix+samples.loc[samples['batch']==wildcards.batch, "sample"].unique()+suffix

def getRawFile(wildcards,prefix):
    return prefix+"/"+samples.loc[(samples['batch']==wildcards.batch) & (samples['sample']==wildcards.sample), "file"].unique()

# default executable for snakemake
shell.executable("bash")

# include configuration file
configfile: srcdir("config/config.default.yaml")

SRC_dir = srcdir("workflow/scripts/")
ENV_dir = srcdir("workflow/envs/")

# custom configuration file
CUSTOM_CONFIG_PATH = os.environ.get("CONFIGFILE","../nothing/here")

# merge 2 configurations files together
if os.path.exists(CUSTOM_CONFIG_PATH):
    print("read configuration from "+CUSTOM_CONFIG_PATH)
    with open(CUSTOM_CONFIG_PATH, 'r') as rhandle:
        data = yaml.load(rhandle)
        config = dict_merge(config, data)
#validate(config,schema="schemas/config.schema.yaml")

EMAIL = config['email']
if EMAIL != "":
    if not re.fullmatch(r"^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$", EMAIL):
        EMAIL = ""
        print("Your email address is not valid, you will not receive notifications.")

if os.path.isabs(os.path.expandvars(config['sample_table'])):
    SAM_PATH = os.path.expandvars(config['sample_table'])
else:
    SAM_PATH = os.getcwd() + "/" + os.path.expandvars(config['sample_table'])
try:
    samples = pd.read_table(SAM_PATH)
except:
    print("Sample table was not found. Please enter the absolute path and file name in the config file.")
    raise

if 'batch' in samples.columns:
    samples['batch'] = samples['batch'].astype(str)
    for r in samples['batch']:
        if not re.match(r"[0-9a-zA-Z]",r):
            raise Exception('please start batch names with a letter or number')
else:
    samples['batch'] = ["batch1"] * samples.shape[0]
    print("adding column with batch info")
samples = samples.set_index(["sample","batch"],drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels]) 
samples['sample'] = samples['sample'].astype(str)
for lib in samples['sample']:
    if not re.match(r"[0-9a-zA-Z]",lib):
        raise Exception('please start sample names with a letter or number')

if 'file' not in samples.columns:
    raise Exception("You haven't provided file names - column should be named file.")


if os.path.isabs(os.path.expandvars(config['raw_directory'])):
    RAW = os.path.expandvars(config['raw_directory'])
else:
    RAW = os.getcwd() + "/" + os.path.expandvars(config['raw_directory'])

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))

OUTPUTDIR = os.path.expandvars(config['outputdir'])

# temporary directory will be stored inside the OUTPUTDIR directory
# unless an absolute path is set
TMPDIR = os.path.expandvars(config['tmp_dir'])
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.abspath(os.path.join(OUTPUTDIR, TMPDIR))
elif not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

workdir:
    OUTPUTDIR

f = open('full.config.yaml', 'w+')
yaml.dump(config, f, allow_unicode=True,default_flow_style=False)


if EMAIL == "":
    onsuccess:
        shell("mkdir -p job.errs.outs; mv slurm* job.errs.outs || touch job.errs.outs; mv snakejob.* job.errs.outs || touch job.errs.outs; mv *log job.errs.outs || touch job.errs.outs; mv *logfile job.errs.outs || touch job.errs.outs")
else:
    onsuccess:
        shell('mkdir -p job.errs.outs; mv slurm* job.errs.outs || touch job.errs.outs; mv snakejob.* job.errs.outs || touch job.errs.outs; mv *log job.errs.outs || touch job.errs.outs; mv *logfile job.errs.outs || touch job.errs.outs; echo "$(date) {config[sessionName]}" | mail -s "dadasnake finished" {EMAIL} ')
    onerror:
        shell('echo "$(date) {config[sessionName]}" | mail -s "dadasnake exited with error" {EMAIL} ')
    onstart:
        shell('echo "$(date) {config[sessionName]}" | mail -s "dadasnake started" {EMAIL} ')


localrules: ALL, SamplesPrint, fastqc_all

# master command
rule ALL:
    input:
        "status/fastqc.done",
        "sample_table.tsv",
        "miRNA.workspace.RDS"
    output:
        touch('workflow.done')
    threads: 1


rule SamplesPrint:
    input:
        SAM_PATH
    output:
        "sample_table.tsv"
    threads: 1
    params:
        runtime="00:10:00",
        mem="8G"
    run:
        samples.to_csv(path_or_buf=output[0],sep="\t",index=False,index_label=False)

rule fastqc_all:
    input:
        expand("fastqc/{samples.batch}/{samples.sample}/{step}_fastqc.html",samples=samples.itertuples(),step=["raw","trim","trim_trim2","trim_trim2_filt"])
    output:
        touch("status/fastqc.done")
    threads: 1


# fastqc 1: raw
rule fastqcReads_raw:
    input:
        files = lambda wildcards: getRawFile(wildcards,RAW)
    output:
        "fastqc/{batch}/{sample}/raw_fastqc.html"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00",
        outputdir="fastqc"
    conda: ENV_dir +"fastqc_env.yml"
    message: "fastqc reads in {input}."
    log: "logs/fastqcReads_raw.{batch}.{sample}.log"
    shell:
        """
        mkdir -p {params.outputdir}/{wildcards.batch}/{wildcards.sample}
        cp {input} {params.outputdir}/{wildcards.batch}/{wildcards.sample}/raw.fq
        fastqc -o {params.outputdir}/{wildcards.batch}/{wildcards.sample} --extract -f fastq {params.outputdir}/{wildcards.batch}/{wildcards.sample}/raw.fq -t {threads} -d {TMPDIR} > {log} 2>&1
        rm {params.outputdir}/{wildcards.batch}/{wildcards.sample}/raw.fq
        """
        
# fastqc 2: clipped
rule fastqcReads_clipped:
    input:
        "cln.adapt/{batch}/{sample}.trim.fq"
    output:
        "fastqc/{batch}/{sample}/trim_fastqc.html"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00",
        outputdir="fastqc"
    conda: ENV_dir +"fastqc_env.yml"
    message: "fastqc reads in {input}."
    log: "logs/fastqcReads_clipped.{batch}.{sample}.log"
    shell:
        """
        mkdir -p {params.outputdir}/{wildcards.batch}/{wildcards.sample}
        cp {input} {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim.fq
        fastqc -o {params.outputdir}/{wildcards.batch}/{wildcards.sample} --extract -f fastq {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim.fq -t {threads} -d {TMPDIR} > {log} 2>&1
        rm {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim.fq
    """
        
# fastqc 3: trimmed
rule fastqcReads_trimmed:
    input:
        "cln.trim/{batch}/{sample}.trim_trim2.fq"
    output:
        "fastqc/{batch}/{sample}/trim_trim2_fastqc.html"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00",
        outputdir="fastqc"
    conda: ENV_dir +"fastqc_env.yml"
    message: "fastqc reads in {input}."
    log: "logs/fastqcReads_trimmed.{batch}.{sample}.log"
    shell:
        """
        mkdir -p {params.outputdir}/{wildcards.batch}/{wildcards.sample}
        cp {input} {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim_trim2.fq
        fastqc -o {params.outputdir}/{wildcards.batch}/{wildcards.sample} --extract -f fastq {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim_trim2.fq -t {threads} -d {TMPDIR} > {log} 2>&1
        rm {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim_trim2.fq
        """
        
# fastqc 4: final
rule fastqcReads_final:
    input:
        "cln.filt/{batch}/{sample}.trim_trim2_filt.fq"
    output:
        "fastqc/{batch}/{sample}/trim_trim2_filt_fastqc.html"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00",
        outputdir="fastqc"
    conda: ENV_dir +"fastqc_env.yml"
    message: "fastqc reads in {input}."
    log: "logs/fastqcReads_final.{batch}.{sample}.log"
    shell:
        """
        mkdir -p {params.outputdir}/{wildcards.batch}/{wildcards.sample}
        cp {input} {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim_trim2_filt.fq
        fastqc -o {params.outputdir}/{wildcards.batch}/{wildcards.sample} --extract -f fastq {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim_trim2_filt.fq -t {threads} -d {TMPDIR} > {log} 2>&1
        rm {params.outputdir}/{wildcards.batch}/{wildcards.sample}/trim_trim2_filt.fq
        """

# cutadapt
rule cutReads:
    input:
        files = lambda wildcards: getRawFile(wildcards,RAW)
    output:
        "cln.adapt/{batch}/{sample}.trim.fq"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00",
        five_prime=config['adapters']['five_prime'],
        three_prime=config['adapters']['three_prime'],
        length=config['read_length']
    conda: ENV_dir +"cutadapt_env.yml"
    message: "cutting adapters from reads in {input}."
    log: "logs/cutReads.{batch}.{sample}.log"
    shell:
        """
        cutadapt -a {params.five_prime}...{params.three_prime} -O 10 -m {params.length} -o {output} {input} &> {log}
        """

# fastq_quality_trimmer
rule trimReads:
    input:
        "cln.adapt/{batch}/{sample}.trim.fq"
    output:
        "cln.trim/{batch}/{sample}.trim_trim2.fq"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00",
        phred=config['phred_format'],
        length=config['read_length']
    conda: ENV_dir +"fastx_env.yml"
    message: "quality trimming reads in {input}."
    log: "logs/trimReads.{batch}.{sample}.log"
    shell:
        """
        fastq_quality_trimmer -v -t 25 {params.phred} -l {params.length} -i {input} -o {output} &> {log}
        """

# fastq_quality_filter
rule filterReads:
    input:
        "cln.trim/{batch}/{sample}.trim_trim2.fq"
    output:
        "cln.filt/{batch}/{sample}.trim_trim2_filt.fq"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00",
        phred=config['phred_format']
    conda: ENV_dir +"fastx_env.yml"
    message: "quality filtering reads in {input}."
    log: "logs/filterReads.{batch}.{sample}.log"
    shell:
        """
        fastq_quality_filter -v -q 25 -p 100 {params.phred} -i {input} -o {output} &> {log}
        """

# fastx_collapser
rule collapseReads:
    input:
        "cln.filt/{batch}/{sample}.trim_trim2_filt.fq"
    output:
        "count_fas/{batch}/{sample}.trim_trim2_filt_counts.fa"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00"
    conda: ENV_dir +"fastx_env.yml"
    message: "collapsing unique reads in {input}."
    log: "logs/collapseReads.{batch}.{sample}.log"
    shell:
        """
        fastx_collapser -i {input} -o {output} &> {log}
        """

# convertCounts: counts2tabLengthFilter.pl
rule tableCounts:
    input:
        "count_fas/{batch}/{sample}.trim_trim2_filt_counts.fa"
    output:
        "count_tabs/{batch}/{sample}.countTab.txt"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00",
        min=config['tag_min_length'],
        max=config['tag_max_length']
    conda: ENV_dir +"perl_env.yml"
    message: "converting collapsed reads to table for {input}."
    log: "logs/tableCounts.{batch}.{sample}.log"
    shell:
        """
        perl {SRC_dir}/counts2tabLengthFilter.pl {input} {params.min} {params.max} >> {output}
        """

# prepare mature DB
rule prepDB:
    input:
        config['miR_DB']
    output:
        "mature.txt"
    threads: 1
    params:
        mem="8G",
        runtime="2:00:00"
    message: "Preparing mature miR DB."
    log: "logs/prepDB.log"
    shell:
        """
        cp {input} {output} &> {log}
        """

# isomiRSEA
rule isomirSEA:
    input:
        "count_tabs/{batch}/{sample}.countTab.txt",
        "mature.txt"
    output:
        "isomiR_SEA/{batch}/{sample}/out_result_mature_unique.txt",
        "isomiR_SEA/{batch}/{sample}/out_result_mature_ambigue.txt",
        "isomiR_SEA/{batch}/{sample}/summary.txt"
    threads: 1
    params:
        mem="8G",
        runtime="24:00:00",
        isomirSEA_binary=config['path_to_isomirSEA'],
        length=config['miR_min_length'],
        outdir="isomiR_SEA"
    message: "Running isomiR-SEA on {wildcards.sample}."
    log: "logs/isomirSEA.{batch}.{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}/{wildcards.batch}/{wildcards.sample}
        cd {params.outdir}/{wildcards.batch}/{wildcards.sample}
        cp ../../../{input[1]} .
        cp ../../../{input[0]} .
        FILE={input[0]}
        IFILE=${{FILE##*/}}
        IN=${{IFILE%.txt}}
        {params.isomirSEA_binary} -s hsa -l {params.length} -ss 6 -h 11 -m mature -b 4 -i . -p {params.outdir} -t $IN >> ../../../{output[2]} 2> ../../../{log}
        rm mature.txt
        rm mature_db_group_sorted.txt
        rm mature_db_group.txt
        rm $IFILE 
        """

# R to combine per-sample results
rule combineResults:
    input:
        uni=expand("isomiR_SEA/{samples.batch}/{samples.sample}/out_result_mature_unique.txt", samples=samples.itertuples()),
        ambi=expand("isomiR_SEA/{samples.batch}/{samples.sample}/out_result_mature_ambigue.txt", samples=samples.itertuples())
    output:
        "miRNA.workspace.RDS"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00"
    conda: ENV_dir +"R_env.yml"
    message: "combining sample results in one R workspace."
    log: "logs/combineResults.log"
    script:
        SRC_dir+"combine_results.R"



