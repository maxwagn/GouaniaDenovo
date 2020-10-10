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
#id_list = sample_mt["sequence_id"].tolist()
id_list = ["TESTFILE"]
readnames = ["1P", "2P"]
kmer_ext = ["lower", "upper", "optimal"]
assembly = ["SOAP", "ABYSS", "DISCOVAR"]

rule all:
    input:
        "reports/raw/multiqc_report.html",
        "reports/trim/multiqc_report.html",
        "reports/QUAST_contigs/transposed_report.tsv",
        expand("reports/BUSCO_contigs/{{id}}_{{assembler}}_busco/{{id}}_{{assembler}}_busco/run_{}/short_summary.txt".format(config["BUSCO"]["lineage"]), id = id_list, assembler = assembly),
        "reports/final/BUSCOandQUAST_summary_final.tsv"

rule fastqc_raw:
    input:
        fastq = "{}/{{id}}_{{read}}.fastq".format(config["raw_dir"])
    output:
        "reports/raw/{id}_{read}_fastqc.html",
        "reports/raw/{id}_{read}_fastqc.zip"
    conda:
        "envs/QCenv.yml"
    threads: config["resources"]["fastqc"]["threads"]
    resources:
        mem_gb = config["resources"]["fastqc"]["mem_gb"],
        walltime = config["resources"]["fastqc"]["walltime"]
    shell: 
        "fastqc -o ./reports/raw/ {input.fastq} -t {threads}"

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
        f = "{}/{{id}}_1P.fastq".format(config["raw_dir"]),
        r = "{}/{{id}}_2P.fastq".format(config["raw_dir"])
    output:
        fout = "trimmed/{id}_1P_trim.fastq",
        funp = "trimmed/{id}_1P_unpaired.fastq",
        rout = "trimmed/{id}_2P_trim.fastq",
        runp = "trimmed/{id}_2P_unpaired.fastq"
    threads: config["resources"]["trimmomatic"]["threads"]
    resources:
        mem_gb = config["resources"]["trimmomatic"]["mem_gb"],
        walltime = config["resources"]["trimmomatic"]["walltime"]
    conda:
        "envs/SoapDenovo.yml"
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
    threads: config["resources"]["fastqc"]["threads"]
    resources:
        mem_gb = config["resources"]["fastqc"]["mem_gb"],
        walltime = config["resources"]["fastqc"]["walltime"]
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
        rRead = "trimmed/{id}_2P_trim.fastq"
    output:
        bam = "trimmed/bam/{id}_trim.bam"
    threads: config["resources"]["DiscoVar"]["threads"]
    resources:
        mem_gb = config["resources"]["DiscoVar"]["mem_gb"],
        walltime = config["resources"]["DiscoVar"]["walltime"]
    conda:
        "envs/bamprep.yml"
    shell:
        "picard FastqToSam F1={input.fRead} F2={input.rRead} O={output.bam} SM={wildcards.id}"

rule KmerGenie_prep:
    input:
        f = "trimmed/{id}_1P_trim.fastq", 
        r = "trimmed/{id}_2P_trim.fastq"
    output:
        "metaAndconfig/KmerGenie/{id}_KmerGenie.config"
    shell:
        "ls -1 {input.f} {input.r} > {output}"

rule KmerGenie:
    input:
        readlist = "metaAndconfig/KmerGenie/{id}_KmerGenie.config" 
    output:
        "reports/KmerGenie/{id}_report.html",
        logfile = "reports/KmerGenie/{id}_kmergenie.log"
    params:
        outdir = "reports/KmerGenie/{id}",
    threads: config["resources"]["KmerGenie"]["threads"]
    resources:
        mem_gb = config["resources"]["KmerGenie"]["mem_gb"],
        walltime = config["resources"]["KmerGenie"]["walltime"]
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
        "assemblies/{id}_SOAPDENOVO/{id}_soap_{kmer}K.scafSeq",
        "assemblies/{id}_SOAPDENOVO/{id}_soap_{kmer}K.scafStatistics"
    params:
        outdir = "assemblies/{id}_SOAPDENOVO/{id}_soap_{kmer}K",
    conda:
        "envs/SoapDenovo.yml"
    threads: config["resources"]["SoapDenovo"]["threads"]
    resources:
        mem_gb = config["resources"]["SoapDenovo"]["mem_gb"],
        walltime = config["resources"]["SoapDenovo"]["walltime"] 
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
    threads: config["resources"]["ABySS"]["threads"]
    resources:
        mem_gb = config["resources"]["ABySS"]["mem_gb"],
        walltime = config["resources"]["ABySS"]["walltime"]
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
    threads: config["resources"]["DiscoVar"]["threads"]
    resources:
        mem_gb = config["resources"]["DiscoVar"]["mem_gb"],
        walltime = config["resources"]["DiscoVar"]["walltime"]
    conda:
        "envs/discovar-denovo.yml"
    shell:
        "DiscovarDeNovo READS={input.bam} OUT_DIR={params.outdir} NUM_THREADS={threads} MAX_MEM_GB={resources.mem_gb} TEE={params.logfile}"

rule renameDiscovar:
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
        "reports/QUAST_contigs/transposed_report.tsv"
    params:
        outdir = "reports/QUAST_contigs/"
    threads: config["resources"]["QUAST"]["threads"]
    conda:
        "envs/QCenv.yml"
    shell:
        "quast.py -o {params.outdir} {input.soapdenovo} {input.abyss} {input.discovar} --threads {threads}" #-R {input.ref} -G {input.gff}"

rule extract_best_N50_contigs_for_BUSCO:
    input:
        "reports/QUAST_contigs/transposed_report.tsv"
    output:
        "metaAndconfig/bestN50/{id}_{assembler}_highestN50.txt"
    params:
        metadata = config["sample_mt"],
        outdir = "metaAndconfig/bestN50/"
    shell:
        "python scripts/gethighestN50.py {wildcards.id} {params.metadata} {input} {params.outdir}"

rule BUSCO_contigs:
    input:
        "metaAndconfig/bestN50/{id}_{assembler}_highestN50.txt"
    output:
        "reports/BUSCO_contigs/{{id}}_{{assembler}}_busco/{{id}}_{{assembler}}_busco/run_{}/short_summary.txt".format(config["BUSCO"]["lineage"])
    params:
        outdir = "reports/BUSCO_contigs/",
        mode = config["BUSCO"]["mode"],
        lineage = config["BUSCO"]["lineage"],
        ids = "{id}_{assembler}_busco"
    threads: config["resources"]["BUSCO"]["threads"]
    resources:
        mem_gb = config["resources"]["BUSCO"]["mem_gb"],
        walltime = config["resources"]["BUSCO"]["walltime"]
    conda:
        "envs/BUSCO.yml"
    shell:
        """
        contigfile=$(cat {input})
        echo Making BUSCO search on $contigfile
        busco -i $contigfile -l {params.lineage} -o {params.ids} --out_path {params.outdir} -m {params.mode} -c {threads} -f
        """

rule RagTag:
    input:
        query = "metaAndconfig/bestN50/{id}_{assembler}_highestN50.txt",
        ref = config["reference"]["fasta"]
    output:
        "RagTag/{id}_{assembler}_ragtag/{id}_{assembler}_ragtag.scaffolds.fasta"
    params:
        id = {id},
        outdir = "RagTag/{id}_{assembler}_ragtag/"
    threads: config["resources"]["RagTag"]["threads"]
    resources:
        mem_gb = config["resources"]["RagTag"]["mem_gb"],
        walltime = config["resources"]["RagTag"]["walltime"]
    conda:
        "envs/RagTag.yml"
    shell:
        """
        contigfile=$(cat {input.query})
        echo Conducting RagTag correction on $contigfile
        ragtag.py correct {input.ref} $contigfile -o {params.outdir} -t {threads}
        filename=$(basename $(cat {input.query}) | awk -F "." '{{print $1".corrected.fasta"}}')
        ragtag.py scaffold {input.ref} {params.outdir}$filename -o {params.outdir} -t {threads}
        mv {params.outdir}ragtag.scaffolds.fasta {output}
        """

rule quast_RagTag:
    input:
        expand("RagTag/{id}_{assembler}_ragtag/{id}_{assembler}_ragtag.scaffolds.fasta", id = id_list, assembler = assembly)
    output:
        "reports/QUAST_RAGTAG/transposed_report.tsv"
    params:
        outdir = "reports/QUAST_RAGTAG/"
    threads: config["resources"]["QUAST"]["threads"] 
    conda:
        "envs/QCenv.yml"
    shell:
        "quast.py -o {params.outdir} {input} --threads {threads}"


rule BUSCO_RagTag_scaffolds:
    input:
        "RagTag/{id}_{assembler}_ragtag/{id}_{assembler}_ragtag.scaffolds.fasta"
    output:
        "reports/BUSCO_RagTag_scaffolds/{{id}}_{{assembler}}_RagTag_busco/run_{}/short_summary.txt".format(config["BUSCO"]["lineage"])
    params:
        outdir = "reports/BUSCO_RagTag_scaffolds/",
        mode = config["BUSCO"]["mode"],
        lineage = config["BUSCO"]["lineage"],
        ids = "{id}_{assembler}_RagTag_busco"
    threads: config["resources"]["BUSCO"]["threads"]
    resources:
        mem_gb = config["resources"]["BUSCO"]["mem_gb"],
        walltime = config["resources"]["BUSCO"]["walltime"]
    conda:
        "envs/BUSCO.yml"
    shell:
        "busco -i {input} -l {params.lineage} -o {params.ids} --out_path {params.outdir} -m {params.mode} -c {threads} -f"

rule finalReportPrep:
    input:
        expand("reports/BUSCO_RagTag_scaffolds/{{id}}_{{assembler}}_RagTag_busco/run_{}/short_summary.txt".format(config["BUSCO"]["lineage"]), id = id_list, assembler = assembly)
    output:
        "reports/BUSCO_RagTag_scaffolds/BUSCO_report_AllinOne.txt"
    shell:
        "cat {input} >> {output}"

rule finalReport:
    input:
        busco = "reports/BUSCO_RagTag_scaffolds/BUSCO_report_AllinOne.txt",
        quast = "reports/QUAST_RAGTAG/transposed_report.tsv"
    output:
        "reports/final/QUAST_summary_PLOT.pdf",
        "reports/final/BUSCO_summary_PLOT.pdf",
        "reports/final/BUSCOandQUAST_summary_final.tsv"
    params:
        metadata = config["sample_mt"]
    shell:
        "python scripts/mergeBUSCO_QUAST.py {input.busco} {input.quast} {params.metadata}"

