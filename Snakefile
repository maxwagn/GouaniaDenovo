import pandas as pd
import numpy as np
import os
jn = os.path.join

#snakemake -p --profile qsub_hs_test --keep-going --immediate-submit --use-conda
# ${progWriteSoapConfig} -insLength ${libList} -r1 ${forwardReads} -r2 ${reverseReads} -max ${maxReadLength} -ru 3 -rank -o ${soapConf}\n" >> ${fileout}


## Config
configfile: "metaAndconfig/config.yaml"

CHROMOSOMES = [str(i) for i in range(1,25) if i not in [23] ] # +  ['MT']

sample_mt = pd.read_csv(config["sample_mt"],
                                dtype=str, sep='\t').set_index("sequence_id", drop=False)
sample_mt.drop_duplicates(inplace=True)
n_samples = len(sample_mt)

readnames = ["1P", "2P"]
id_list = sample_mt["sequence_id"].tolist()

#SOAPdenovo2 has two commands, SOAPdenovo-63mer and SOAPdenovo-127mer. 
#The first one is suitable for assembly with k-mer values less than 63 bp, requires less memory and runs faster. The latter one works for k-mer values less than 127 bp.
K_mers = [65, 75, 85, 95, 105, 115]
soapout_ext = [".newContigIndex",".links",".scaf_gap",".scaf",".gapSeq",".scafSeq",".contigPosInscaff",".bubbleInScaff",".scafStatistics"]


rule all:
    input:
        #expand("reports/raw/{id}_{read}_fastqc.html", id = id_list, read = readnames), 
        #expand("reports/trim/{id}_{read}_trim_fastqc.html", id = id_list, read = readnames)
        #"reports/raw/multiqc_report.html",
        #"reports/trim/multiqc_report.html"
        #expand("metaAndconfig/soapconfig/{id}_soapconfig", id = id_list)
        ###### SOAP DENOVO ######
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.newContigIndex", id = id_list[4],  kmer = K_mers[0]),
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.links", id = id_list[4],  kmer = K_mers[0]),
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.scaf_gap", id = id_list[4],  kmer = K_mers[0]),
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.scaf", id = id_list[4],  kmer = K_mers[0]),
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.gapSeq", id = id_list[4],  kmer = K_mers[0]),
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.scafSeq", id = id_list[4],  kmer = K_mers[0]),
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.contigPosInscaff", id = id_list[4],  kmer = K_mers[0]),
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.bubbleInScaff", id = id_list[4],  kmer = K_mers[0]),
        #expand("assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.scafStatistics", id = id_list[4],  kmer = K_mers[0]),
        ###### ABySS ######
        #expand("assemblies/{id}_ABySS/{kmer}KSize/{id}_abyss_K{kmer}-contigs.fa", id = id_list[4],  kmer = K_mers[0])
        ###### SPAdes #####
        expand("assemblies/{id}_SPAdes/{kmer}KSize/{id}_spades_K{kmer}/scaffolds.fasta", id = id_list[4],  kmer = K_mers[0])


rule fastqc_raw:
    input:
        fastq = "../../rawdata/{id}_{read}.fastq"
    output:
        "reports/raw/{id}_{read}_fastqc.html",
        "reports/raw/{id}_{read}_fastqc.zip"
    shell: 
        "fastqc -o ./reports/raw/ {input.fastq}"

rule multiqc_raw:
    input:
        expand("reports/raw/{id}_{read}_fastqc.html", id = id_list, read = readnames)
    output:
        "reports/raw/multiqc_report.html"    
    shell:
        "multiqc ./reports/raw/ -o ./reports/raw/"

rule trimmomatic:
    input:
        f = "../../rawdata/{id}_1P.fastq",
        r = "../../rawdata/{id}_2P.fastq"
    output:
        fout = "trimmed/{id}_1P_trim.fastq",
        funp = "trimmed/{id}_1P_unpaired.fastq",
        rout = "trimmed/{id}_2P_trim.fastq",
        runp = "trimmed/{id}_2P_unpaired.fastq"
    shell:
        "trimmomatic PE -threads 8 {input.f} {input.r} {output.fout} {output.funp} {output.rout} {output.runp} ILLUMINACLIP:adapterseq/Adapters_PE.fa:2:30:10: LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:80"

rule fastqc_trim:
    input:
        fastq = "trimmed/{id}_{read}_trim.fastq"
    output:
        "reports/trim/{id}_{read}_trim_fastqc.html",
        "reports/trim/{id}_{read}_trim_fastqc.zip"
    shell:
        "fastqc -o reports/trim/ {input.fastq}"

rule multiqc_trim:
    input:
        expand("reports/trim/{id}_{read}_trim_fastqc.html", id = id_list, read = readnames)
    output:
        "reports/trim/multiqc_report.html"
    shell:
        "multiqc reports/trim/ -o reports/trim/"


#https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/de_novo_assembly_tools/soapdenovo2/

rule soap_config:
    input:
        fRead = "trimmed/{id}_1P_trim.fastq",
        rRead = "trimmed/{id}_2P_trim.fastqi"
    output:
        soapconfig = "metaAndconfig/soapconfig/{id}_soapconfig"
    shell:
        "python3 scripts/writesoapconfig.py {output.soapconfig} {input.fRead} {input.rRead}"

rule soapdenovo:
    input:
        config = "metaAndconfig/soapconfig/{id}_soapconfig"
    output:
        # get info about all outpufiles here: https://www.animalgenome.org/bioinfo/resources/manuals/SOAP.html
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.newContigIndex",
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.links",
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.scaf_gap", 
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.scaf", 
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.gapSeq",
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.scafSeq",
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.contigPosInscaff",
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.bubbleInScaff", 
        "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}.scafStatistics",
        outdir = "assemblies/{id}_SOAPDENOVO/{kmer}KSize/{id}_soap_K{kmer}"
    params:
        ksize = "{kmer}"
    conda:
        "envs/SoapDenovo.yml"
    resources:
        mem_mb = 80,
        walltime = 72
    shell:
        """
        if (({params.ksize} <= 63)); 
        then
        SOAPdenovo-63mer all -s {input.config} -K {params.ksize} -R -o {output.outdir}
        elif ((63 < {params.ksize} <= 127));
        then
        SOAPdenovo-127mer all -s {input.config} -K {params.ksize} -R -o {output.outdir}
        else
        echo "K-mer size has to be set between 1 and 127"
        fi
        """
rule ABySS:
    input:
        fRead = "trimmed/{id}_1P_trim.fastq",
        rRead = "trimmed/{id}_2P_trim.fastq"
    output:
        "assemblies/{id}_ABySS/{kmer}KSize/{id}_abyss_K{kmer}-contigs.fa",
        outdir = "assemblies/{id}_ABySS/{kmer}KSize/{id}_abyss_K{kmer}"
    params:
        ksize = "{kmer}"
    conda:
        "envs/ABySS.yml"
    shell:
        "abyss-pe name={output.outdir} k={params.ksize} in='{input.fRead} {input.rRead}'"

rule SPADES:
    input:
        fRead = "trimmed/{id}_1P_trim.fastq",
        rRead = "trimmed/{id}_2P_trim.fastq"
    output:
        "assemblies/{id}_SPAdes/{kmer}KSize/{id}_spades_K{kmer}/scaffolds.fasta",
        outdir = "assemblies/{id}_SPAdes/{kmer}KSize/{id}_spades_K{kmer}/"
    params:
        ksize = "{kmer}"
    conda:
        "envs/SPADES.yml"
    shell:
        "spades.py -k {params.ksize}-1 {input.fRead} -2 {input.rRead} --careful --covcutiff -o {output.outdir}"

