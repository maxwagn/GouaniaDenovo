import pandas as pd
import numpy as np
import os
jn = os.path.join

#snakemake -p --profile qsub_hs_test --keep-going --immediate-submit --use-conda
#snakemake -np --dag | dot -Tsvg > dag.svg # draw dag graph


## Config
configfile: "metaAndconfig/config.yaml"

CHROMOSOMES = [str(i) for i in range(1,25) if i not in [23] ] # +  ['MT']

sample_mt = pd.read_csv(config["sample_mt"],
                                dtype=str, sep='\t').set_index("sequence_id", drop=False)
sample_mt.drop_duplicates(inplace=True)
n_samples = len(sample_mt)
id_list = sample_mt["sequence_id"].tolist()
#id_list = ["TESTFILE"]
readnames = ["1P", "2P"]
kmer_ext = ["lower", "upper", "optimal"]

rule all:
    input:
        #expand("reports/trim/{id}_{read}_trim_fastqc.html", id = id_list, read = readnames)
        #"reports/raw/multiqc_report.html",
        #"reports/trim/multiqc_report.html"
        #expand("metaAndconfig/soapconfig/{id}_soapconfig", id = id_list)
        ##### KmerGenie #####
        #expand("metaAndconfig/KmerGenie/{id}_KmerGenie.config", id = id_list)
        #expand("reports/KmerGenie/{id}_report.html", id = id_list)
        #expand("metaAndconfig/KmerGenie/{id}_lowerK.config", id = id_list),
        #expand("metaAndconfig/KmerGenie/{id}_upperK.config", id = id_list),
        #expand("metaAndconfig/KmerGenie/{id}_optimalK.config", id = id_list)
        ###### SOAP DENOVO ######
        #expand("assemblies/{id}_SOAPDENOVO/{id}_soap_{kmer}K.scafStatistics", id = id_list,  kmer = kmer_ext)
        ###### ABySS ######
        #expand("assemblies/{id}_ABySS/{id}_abyss_{kmer}K-stats.csv", id = id_list[1:],  kmer = kmer_ext)
        #expand("assemblies/{id}_ABySS/{id}_abyss_{kmer}K-scaffolds.fa", id = id_list[1:],  kmer = kmer_ext)
        ###### QUAST #####
        #expand("metaAndconfig/bestN50/{id}_DISCOVAR_highestN50.txt", id = id_list)

rule fastqc_raw:
    input:
        fastq = "../../rawdata/{id}_{read}.fastq"
    output:
        "reports/raw/{id}_{read}_fastqc.html",
        "reports/raw/{id}_{read}_fastqc.zip"
    conda:
        "envs/QCenv.yml"
    shell: 
        "fastqc -o ./reports/raw/ {input.fastq}"

rule multiqc_raw:
    input:
        expand("reports/raw/{id}_{read}_fastqc.html", id = id_list, read = readnames)
    output:
        "reports/raw/multiqc_report.html"    
    conda:
        "envs/SoapDenovo.yml"
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
    threads: 8
    shell:
        "trimmomatic PE -threads {threads} {input.f} {input.r} {output.fout} {output.funp} {output.rout} {output.runp} ILLUMINACLIP:adapterseq/Adapters_PE.fa:2:30:10: LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:80"

rule fastqc_trim:
    input:
        fastq = "trimmed/{id}_{read}_trim.fastq"
    output:
        "reports/trim/{id}_{read}_trim_fastqc.html",
        "reports/trim/{id}_{read}_trim_fastqc.zip"
    conda:
        "envs/QCenv.yml"
    shell:
        "fastqc -o reports/trim/ {input.fastq}"

rule multiqc_trim:
    input:
        expand("reports/trim/{id}_{read}_trim_fastqc.html", id = id_list, read = readnames)
    output:
        "reports/trim/multiqc_report.html"
    conda:
        "envs/SoapDenovo.yml"
    shell:
        "multiqc reports/trim/ -o reports/trim/"


#https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/de_novo_assembly_tools/soapdenovo2/

rule soap_config:
    input:
        fRead = "trimmed/{id}_1P_trim.fastq",
        rRead = "trimmed/{id}_2P_trim.fastq"
    output:
        soapconfig = "metaAndconfig/soapconfig/{id}_soapconfig"
    shell:
        "python3 scripts/writesoapconfig.py {output.soapconfig} {input.fRead} {input.rRead}"

rule DiscoVarPrep:
    input:
        fRead = "trimmed/{id}_1P_trim.fastq",
        rRead = "trimmed/{id}_2P_trim.fastq",
    output:
        bam = "trimmed/bam/{id}_trim.bam"
    threads: 24
    resources:
        mem_gb = 150,
        walltime = 12
    conda:
        "envs/bamprep.yml"
    shell:
        "picard FastqToSam F1={input.fRead} F2={input.rRead} O={output.bam} SM={wildcards.id}" # RG=rg0013"

rule KmerGenie_prep:
    input:
        f = "trimmed/{id}_1P_trim.fastq", 
        r = "trimmed/{id}_2P_trim.fastq"
    output:
        "metaAndconfig/KmerGenie/{id}_KmerGenie.config"
    shell:
        "ls -1 {input.f} {input.r} > {output}"

rule KmerGenie:
    #https://onestopdataanalysis.com/estimate-genome-size-best-k-mer-size-for-assembly/
    input:
        readlist = "metaAndconfig/KmerGenie/{id}_KmerGenie.config" 
    output:
        "reports/KmerGenie/{id}_report.html",
        logfile = "reports/KmerGenie/{id}_kmergenie.log"
    params:
        outdir = "reports/KmerGenie/{id}",
        #logfile = "reports/KmerGenie/{id}_kmergenie.log"
    threads: 14
    resources:
        mem_gb = 400,
        walltime = 48
    conda:
        "envs/KmerGenie.yml"
    shell:
        "kmergenie {input.readlist} -l 21 -k 121 -s 6 -o {params.outdir} -t {threads} --diploid | tee > {output.logfile}"

rule KmerInputLowerUpperBest:
    input: 
        "reports/KmerGenie/{id}_kmergenie.log"
    output: 
        lower = "metaAndconfig/KmerGenie/{id}_lowerK.config",
        upper = "metaAndconfig/KmerGenie/{id}_upperK.config",
        optimal = "metaAndconfig/KmerGenie/{id}_optimalK.config"
    shell:
        "python3 scripts/kmerextractor.py {input} {output.lower} {output.upper} {output.optimal}"
    
rule soapdenovo:
    input:
        ksize = "metaAndconfig/KmerGenie/{id}_{kmer}K.config",
        config = "metaAndconfig/soapconfig/{id}_soapconfig"
    output:
        # get info about all outpufiles here: https://www.animalgenome.org/bioinfo/resources/manuals/SOAP.html
        "assemblies/{id}_SOAPDENOVO/{id}_soap_{kmer}K.scafSeq",
        "assemblies/{id}_SOAPDENOVO/{id}_soap_{kmer}K.scafStatistics"
    threads: 26 
    params:
        outdir = "assemblies/{id}_SOAPDENOVO/{id}_soap_{kmer}K",
        #config = "metaAndconfig/soapconfig/{id}_soapconfig"
    conda:
        "envs/SoapDenovo.yml"
    resources:
        mem_gb = 200,
        walltime = 48 
    shell:
        """
        file={input.ksize}
        kSIZE=$(cat "$file")
        echo $kSIZE
        if (($kSIZE <= 63)); 
        then
        SOAPdenovo-63mer all -p {threads} -s {input.config} -K $kSIZE -R -o {params.outdir}
        elif ((63 < $kSIZE <= 127));
        then
        SOAPdenovo-127mer all -p {threads} -s {input.config} -K $kSIZE -R -o {params.outdir}
        else
        echo "K-mer size has to be set between 1 and 127"
        fi
        """

rule ABySS:
    input:
        ksize = "metaAndconfig/KmerGenie/{id}_{kmer}K.config",
        fRead = "trimmed/{id}_1P_trim.fastq",
        rRead = "trimmed/{id}_2P_trim.fastq"
    output:
        "assemblies/{id}_ABySS/{id}_abyss_{kmer}K-scaffolds.fa",
        #"assemblies/{id}_ABySS/{kmer}KSize/{id}_abyss_{kmer}K-stats.csv" 
    params:
        pwd = os.getcwd(),
        outdir = "assemblies/{id}_ABySS/",
        name = "{id}_abyss_{kmer}K",
        logfile = "assemblies/{id}_ABySS/{id}_abyss_{kmer}K.log"
    threads: 26 
    resources:
        mem_gb = 200,
        walltime = 48 
    conda:
        "envs/ABySS.yml"
    shell:
        """
        file={input.ksize}
        kSIZE=$(cat "$file")
        echo $kSIZE
        export TMPDIR={params.pwd}/{params.name}-tmp
        abyss-pe np={threads} -C {params.pwd}/{params.outdir} name={params.name} k=$kSIZE in='{params.pwd}/{input.fRead} {params.pwd}/{input.rRead}' | tee {params.logfile}
        """

rule DiscovarDenovo:
    input:
        bam = "trimmed/bam/{id}_trim.bam"
    output:
        "assemblies/{id}_DISCOVAR/a.final/a.fasta"
    params:
        outdir = "assemblies/{id}_DISCOVAR/",
        logfile = "assemblies/{id}_DISCOVAR/{id}_DISCOVAR.log"
    threads: 24
    resources:
        mem_gb = 210,
        walltime = 48
    conda:
        "envs/discovar-denovo.yml"
    shell:
        "DiscovarDeNovo READS={input.bam} OUT_DIR={params.outdir} NUM_THREADS={threads} MAX_MEM_GB={resources.mem_gb} TEE={params.logfile}"

rule renameDIscovar:
    input:
        "assemblies/{id}_DISCOVAR/a.final/a.fasta"
    output:
        "assemblies/{id}_DISCOVAR/a.final/{id}_discovar.fasta"
    shell:
        "cp {input} {output}"

rule quast_contigs_summary:
    input:
        soapdenovo = expand("assemblies/{id}_SOAPDENOVO/{id}_soap_{kmer}K.scafSeq", id = id_list,  kmer = kmer_ext),
        abyss = expand("assemblies/{id}_ABySS/{id}_abyss_{kmer}K-scaffolds.fa", id = id_list,  kmer = kmer_ext),
        discovar = expand("assemblies/{id}_DISCOVAR/a.final/{id}_discovar.fasta", id = id_list)
    output:
        "reports/QUASTcontigs_finalTest/transposed_report.tsv"
    params:
        outdir = "reports/QUASTcontigs_finalTest/"
    threads: 12
    conda:
        "envs/QCenv.yml"
    shell:
        "quast.py -o {params.outdir} {input.soapdenovo} {input.abyss} {input.discovar} --threads {threads}" #-R {input.ref} -G {input.gff}"

rule extract_best_N50_contigs_for_BUSCO:
    input:
        quast_report = "reports/QUASTcontigs_finalTest/transposed_report.tsv"
    output:
        SOAP = "metaAndconfig/bestN50/{id}_SOAP_highestN50.txt",
        ABYSS = "metaAndconfig/bestN50/{id}_ABYSS_highestN50.txt",
        DISCOVAR = "metaAndconfig/bestN50/{id}_DISCOVAR_highestN50.txt"
    params:
        #pwd = os.getcwd(),
        metadata = config["sample_mt"],
        outdir = "metaAndconfig/bestN50/"
    shell:
        "python scripts/gethighestN50.py wilcards.id {params.metadata} {input.quast_reports} {params.outdir}"

#rule BUSCO_contigs:
    #input:




