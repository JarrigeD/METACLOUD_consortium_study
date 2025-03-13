#####################################################################################################################
# Workflow purpose :                                                                                                #
#     Treats metatranscriptomic datasets from their raw state through quality check, cleaning and filtering,        #
#           then mapping to reference genomes of the reads with STAR and feature counting.                          #
# Author : Domitille Jarrige                                                                                        #
# Date   : 2025-03-13                                                                                               #
# ------------------------------------------------------------------------------------------------------------------#
#                     CHANGE PARAMETERS AND DIRECTORIES IN THE CONFIG FILE BEFORE RUNNING                           #
# ------------------------------------------------------------------------------------------------------------------#
#                                                                                                                   #
# to run workflow:   snakemake --cores={YOURCHOICE} --snakefile {SNAKEFILE} --configfile {CONFIG_FILE}              #
#                    --use-conda --conda-frontend conda                                                             #
#                          Do not forget to build your project architecture before running!                         #
#####################################################################################################################
import os

(RUNS,SAMPLES,FILES,) = glob_wildcards(config["data_directory"] + "Run_{run}/{sample}_R{file}.fastq.gz")


RUNS = list(set(RUNS))
RUNS.sort()
SAMPLES = list(set(SAMPLES))
SAMPLES.sort()
FILES = list(set(FILES))
FILES.sort()
file1=FILES[0]
file2=FILES[1]

print(SAMPLES)
print(RUNS)
print(FILES)

rule all:
    input:
        expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_fastqc.html", run=RUNS, sample=SAMPLES, file = FILES
        ),
        expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_fastqc.zip", run=RUNS, sample=SAMPLES, file = FILES
        ),
        expand(
            config["cleaning_directory"] + "Run_{run}/{sample}_R{file}_cleaned.fastq.gz", run=RUNS, sample=SAMPLES, file=FILES
        ),
        expand(
            config["filtering_directory"] + "Run_{run}/{sample}_R{file}_filtered.fastq.gz", run=RUNS, sample=SAMPLES, file=FILES
        ),
        expand(
            config["rrna_directory"] + "Run_{run}/{sample}_R{file}_rRNA.fastq.gz", run=RUNS, sample=SAMPLES, file=FILES
        ),
        expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_cleaned-fastqc.html", run=RUNS, sample=SAMPLES, file=FILES
        ),
        expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_cleaned-fastqc.zip", run=RUNS, sample=SAMPLES, file=FILES
        ),
        expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_filtered-fastqc.html", run=RUNS, sample=SAMPLES, file=FILES
        ),
        expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_filtered-fastqc.zip", run=RUNS, sample=SAMPLES, file=FILES
        ),
        config["quality_control_directory"] + "multiqc_report.html",
        config["genome_index_directory"] + "community_PDD/Genome",
        expand(
            config["mapping_directory"] + "community_PDD/{sample}_Aligned.out.bam", sample=SAMPLES
        ),
        expand(
            config["mapping_directory"] + "community_PDD/{sample}_rRNA-Aligned.out.bam", sample=SAMPLES
        ),
        expand(
            config["mapping_directory"] + "community_PDD/{sample}_Log.final.out", sample=SAMPLES
        ),
        expand(
            config["mapping_directory"] + "community_PDD/{sample}_rRNA-Log.final.out", sample=SAMPLES
        ),
        expand(
            config["mapping_directory"] + "community_PDD/{sample}_Aligned.out.bam.bai", sample=SAMPLES
        ),
        config["counts_directory"] + "community_PDD_counts.csv",
        config["counts_directory"] + "community_PDD_rRNA-counts.csv",
        config["mapping_directory"] + "summary.txt",




rule fastqc:
    input:
        config["data_directory"] + "Run_{run}/{sample}_R{file}.fastq.gz",
    output:
        html=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_fastqc.html",
        zip=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_fastqc.zip",
        flag=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_qc_done.tmp",
    priority: 12
    threads: config["parameters"]["threads"]
    params:
        outdir=config["quality_control_directory"] + "Run_{run}/",
    conda:
        "envs/envBash.yml"
    shell:
        """
        mkdir -p {params.outdir}
        rm -f {output.flag}
        srun -p fast -c {threads} fastqc {input} --threads={threads} -o {params.outdir} -q
        touch {output.flag}
        """


rule fastp:
    input:
        in1=config["data_directory"] + "Run_{run}/{sample}_R" + file1 + ".fastq.gz",
        in2=config["data_directory"] + "Run_{run}/{sample}_R" + file2 + ".fastq.gz",
    output:
        out1=config["cleaning_directory"] + "Run_{run}/{sample,[^r]}_R" + file1 + "_cleaned.fastq.gz",
        out2=config["cleaning_directory"] + "Run_{run}/{sample,[^r]}_R" + file2 + "_cleaned.fastq.gz",
        flag1=temporary(touch(config["cleaning_directory"] + "Run_{run}/{sample}_R" + file1 + "_clean.tmp")),
        flag2=temporary(touch(config["cleaning_directory"] + "Run_{run}/{sample}_R" + file2 + "_clean.tmp")),
    priority: 11
    threads: config["parameters"]["threads"]
    params: 
        cleaning_dir=config["cleaning_directory"] + "Run_{run}/",
        qc_dir=config["quality_control_directory"],
        out_html= lambda wildcards, output : os.path.basename(output[0][:-20]) + "_fastp.html",
    conda:
        "envs/envBash.yml"
    shell:
        """
        rm -f {output.flag1} {output.flag2}
        srun -p fast -c {threads} fastp --thread {threads} --in1 {input.in1} --in2 {input.in2} --out1 {output.out1} --out2 {output.out2} --html {params.qc_dir} --json /dev/null --qualified_quality_phred 25 --unqualified_percent_limit 15 --cut_mean_quality 25 --cut_window_size 8 --cut_front --cut_tail --length_required 130 --correction
        touch {output.flag1} {output.flag2}
        """



rule fastqc_post_cleaning:
    input:
        infile=config["cleaning_directory"] + "Run_{run}/{sample}_R{file}_cleaned.fastq.gz",
        flag=config["cleaning_directory"] + "Run_{run}/{sample}_R{file}_clean.tmp",
    output:
        html=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_cleaned-fastqc.html",
        zip=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_cleaned-fastqc.zip",
        flag=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_cleaned-qc_done.tmp",
    threads: config["parameters"]["threads"]
    priority: 10
    params:
        outdir=config["quality_control_directory"] + "Run_{run}/",
        tmp_html=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_cleaned_fastqc.html",
        tmp_zip=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_cleaned_fastqc.zip",
    conda:
        "envs/envBash.yml"
    shell:
        """
        mkdir -p {params.outdir}
        rm -f {output.flag}
        srun -p fast -c {threads} fastqc {input.infile} --threads={threads} -o {params.outdir} -q
        mv {params.tmp_html} {output.html}
        mv {params.tmp_zip} {output.zip}
        touch {output.flag}
        """
       

rule remove_rRNA:
    input:
        ref=config["rRNA_database"],
        reads_1=config["cleaning_directory"] + "Run_{run}/{sample}_R" + file1 + "_cleaned.fastq.gz",
        reads_2=config["cleaning_directory"] + "Run_{run}/{sample}_R" + file2 + "_cleaned.fastq.gz",
    output:
        rrna1=config["rrna_directory"] + "Run_{run}/{sample}_R" + file1 + "_rRNA.fastq.gz",
        rrna2=config["rrna_directory"] + "Run_{run}/{sample}_R" + file2 + "_rRNA.fastq.gz",
        out1=config["filtering_directory"] + "Run_{run}/{sample}_R" + file1 + "_filtered.fastq.gz",
        out2=config["filtering_directory"] + "Run_{run}/{sample}_R" + file2 + "_filtered.fastq.gz",
    threads: config["parameters"]["threads"]
    params: 
        work_dir=config["filtering_directory"] + "Run_{run}/{sample}/",
        out_dir=config["filtering_directory"] + "Run_{run}/",
        rrna1=config["filtering_directory"] + "Run_{run}/{sample}/rRNA",
        rrna2=config["filtering_directory"] + "Run_{run}/{sample}/rRNA",
        out1=config["filtering_directory"] + "Run_{run}/{sample}/filtered",
        out2=config["filtering_directory"] + "Run_{run}/{sample}/filtered",
    priority: 9
    conda:
        "envs/envBash.yml"
    shell:
        """
        rm -drf {params.work_dir}
        mkdir -p {params.work_dir}
        mkdir -p {params.out_dir} ##
        srun -p fast -c {threads} sortmerna --ref {input.ref} --reads {input.reads_1} --reads {input.reads_2} --fastx --aligned {params.rrna1} --aligned {params.rrna2} --other {params.out1} --other {params.out2} --out2 --threads {threads} --workdir {params.work_dir} -v False --paired_in True
        mv {params.work_dir}filtered_fwd.fq.gz {output.out1}
        mv {params.work_dir}filtered_rev.fq.gz {output.out2}
        mv {params.work_dir}rRNA_fwd.fq.gz {output.rrna1}
        mv {params.work_dir}rRNA_rev.fq.gz {output.rrna2}
        rm -drf {params.work_dir}
        """
        
        
rule fastqc_post_rRNA_filtering:
    input:
        infile=config["filtering_directory"] + "Run_{run}/{sample}_R{file}_filtered.fastq.gz",
    output:
        html=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_filtered-fastqc.html",
        zip=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_filtered-fastqc.zip",
        flag=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_filtered-qc_done.tmp",
    threads: config["parameters"]["threads"]
    priority: 8
    params:
        outdir=config["quality_control_directory"] + "Run_{run}/",
        tmp_html=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_filtered_fastqc.html",
        tmp_zip=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_filtered_fastqc.zip",
    conda:
        "envs/envBash.yml"
    shell:
        """
        mkdir -p {params.outdir}
        rm -f {output.flag}
        srun -p fast -c {threads} fastqc {input.infile} --threads={threads} -o {params.outdir} -q
        mv {params.tmp_html} {output.html}
        mv {params.tmp_zip} {output.zip}
        touch {output.flag}
        """
        
        
        
rule fastqc_rRNA:
    input:
        infile=config["rrna_directory"] + "Run_{run}/{sample}_R{file}_rRNA.fastq.gz",
    output:
        html=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_rRNA-fastqc.html",
        zip=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_rRNA-fastqc.zip",
        flag=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_rRNA-qc_done.tmp",
    threads: config["parameters"]["threads"]
    priority: 8
    params:
        outdir=config["quality_control_directory"] + "Run_{run}/",
        tmp_html=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_rRNA_fastqc.html",
        tmp_zip=config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_rRNA_fastqc.zip",
    conda:
        "envs/envBash.yml"
    shell:
        """
        mkdir -p {params.outdir}
        rm -f {output.flag}
        srun -p fast -c {threads} fastqc {input.infile} --threads={threads} -o {params.outdir} -q
        mv {params.tmp_html} {output.html}
        mv {params.tmp_zip} {output.zip}
        touch {output.flag}
        """
        
        
rule multiqc:
    input:
        indir=config["quality_control_directory"],
        flag=expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_qc_done.tmp", run=RUNS, sample=SAMPLES, file=FILES
            ),
        flag2=expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_cleaned-qc_done.tmp", run=RUNS, sample=SAMPLES, file=FILES
             ),
        flag3=expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_filtered-qc_done.tmp", run=RUNS, sample=SAMPLES, file=FILES
            ),
        flag4=expand(
            config["quality_control_directory"] + "Run_{run}/{sample}_R{file}_rRNA-qc_done.tmp", run=RUNS, sample=SAMPLES, file=FILES
            ),
    output:
        config["quality_control_directory"] + "multiqc_report.html",
    priority: 7
    params:
        outdir=config["quality_control_directory"][:-1],
    conda:
        "envs/envBash.yml"
    shell:
        """
        multiqc -d {input.indir} -o {params.outdir} -f --interactive
        """        
        
   
        
rule star_indexing:
    input:
        fasta=config["genome_fasta_directory"] + "community.fasta.gz",
        gff=config["genome_gff_directory"] + "community.gff",
    output:
        config["genome_index_directory"] + "community_PDD/Genome",
    threads: config["parameters"]["threads"]
    params:
        out_dir=config["genome_index_directory"] + "community_PDD/",
    priority: 5
    conda:
        "envs/envBash.yml"
    shell:
        """
        mkdir -p "{params.out_dir}"
        srun -p fast -c {threads} zcat {input.fasta} > "{input.fasta}.tmp"
        srun -p fast -c {threads} STAR --runMode genomeGenerate --genomeDir {params.out_dir} --runThreadN {threads} --genomeFastaFiles "{input.fasta}.tmp" --sjdbGTFfile {input.gff} --sjdbOverhang 100 --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 11 
        rm -f "{input.fasta}.tmp"
        """
         

rule star_mapping:
     input:
         in1_1=config["filtering_directory"] + "Run_1/{sample}_R" + file1 + "_filtered.fastq.gz",
         in1_2=config["filtering_directory"] + "Run_1/{sample}_R" + file2 + "_filtered.fastq.gz",
         in2_1=config["filtering_directory"] + "Run_2/{sample}_R" + file1 + "_filtered.fastq.gz",
         in2_2=config["filtering_directory"] + "Run_2/{sample}_R" + file2 + "_filtered.fastq.gz",
         #in1_1=config["filtering_directory"] + "Run_1/run1_{wildcards.short_sample}_R" + file1 + "_filtered.fastq.gz",
         #in1_2=config["filtering_directory"] + "Run_1/run1_{wildcards.short_sample}_R" + file2 + "_filtered.fastq.gz",
         #in2_1=config["filtering_directory"] + "Run_2/run2_{wildcards.short_sample}_R" + file1 + "_filtered.fastq.gz",
         #in2_2=config["filtering_directory"] + "Run_2/run2_{wildcards.short_sample}_R" + file2 + "_filtered.fastq.gz",
         index=config["genome_index_directory"] + "community_PDD/Genome",
     output:
         out=config["mapping_directory"] + "community_PDD/{sample}_Aligned.out.bam",
         log=config["mapping_directory"] + "community_PDD/{sample}_Log.final.out", ###### added 2024-06-18
         #out=config["mapping_directory"] + "community_PDD/{wildcards.sample}_Aligned.out.bam",
     threads: config["parameters"]["threads"] // 2
     params:
         index_dir=config["genome_index_directory"] + "community_PDD",
         gen_dir=config["mapping_directory"] + "community_PDD",
         out_dir=config["mapping_directory"] + "community_PDD/",
         out_prefix=config["mapping_directory"] + "community_PDD/{sample}_",
         #out_prefix=config["mapping_directory"] + "community_PDD/{wildcards.sample}_",
         old_summary=config["mapping_directory"] + "summary.txt",
     priority: 4
     conda:
         "envs/envBash.yml"
     shell:
         """
         touch {output.log}
         mkdir -p {params.out_dir}
         rm -drf {params.out_prefix}tmp
         mkdir -p {params.gen_dir}
         rm -f {params.old_summary}
         srun -p fast -c {threads} STARlong --runThreadN {threads} --readFilesIn {input.in1_1},{input.in2_1} {input.in1_2},{input.in2_2} --readFilesCommand zcat --genomeDir {params.index_dir} --outSAMattrRGline ID:run1 , ID:run2 --outTmpDir "{params.out_prefix}tmp" --outFileNamePrefix {params.out_prefix} --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0.02 --outFilterMatchNminOverLread 0.02 --outFilterMatchNmin 0 --outFilterMismatchNmax 10 --alignEndsProtrude 10 ConcordantPair --peOverlapNbasesMin 1 --outFilterMultimapNmax 10 --outSAMprimaryFlag AllBestScore --outReadsUnmapped Fastx
         mv {params.out_prefix}Aligned.sortedByCoord.out.bam {output.out}
         rm -drf "{params.out_prefix}tmp"
         ls {output.log}
         """


rule star_mapping_rRNA:
     input:
         in1_1=config["rrna_directory"] + "Run_1/{sample}_R" + file1 + "_rRNA.fastq.gz",
         in1_2=config["rrna_directory"] + "Run_1/{sample}_R" + file2 + "_rRNA.fastq.gz",
         in2_1=config["rrna_directory"] + "Run_2/{sample}_R" + file1 + "_rRNA.fastq.gz",
         in2_2=config["rrna_directory"] + "Run_2/{sample}_R" + file2 + "_rRNA.fastq.gz",
         #in1_1=config["filtering_directory"] + "Run_1/run1_{wildcards.#}_R" + file1 + "_rRNA.fastq.gz",
         #in1_2=config["filtering_directory"] + "Run_1/run1_{wildcards.short_sample}_R" + file2 + "_rRNA.fastq.gz",
         #in2_1=config["filtering_directory"] + "Run_2/run2_{wildcards.short_sample}_R" + file1 + "_rRNA.fastq.gz",
         #in2_2=config["filtering_directory"] + "Run_2/run2_{wildcards.short_sample}_R" + file2 + "_rRNA.fastq.gz",
         index=config["genome_index_directory"] + "community_PDD/Genome",
     output:
         out=config["mapping_directory"] + "community_PDD/{sample}_rRNA-Aligned.out.bam",
         log=config["mapping_directory"] + "community_PDD/{sample}_rRNA-Log.final.out",
         #out=config["mapping_directory"] + "community_PDD/{wildcards.sample}_rRNA_Aligned.out.bam",
     threads: config["parameters"]["threads"] // 2
     params:
         index_dir=config["genome_index_directory"] + "community_PDD",
         gen_dir=config["mapping_directory"] + "community_PDD",
         out_dir=config["mapping_directory"] + "community_PDD/",
         out_prefix=config["mapping_directory"] + "community_PDD/{sample}_rRNA_",
         #out_prefix=config["mapping_directory"] + "community_PDD/{sample}_rRNA_",
         #add_param=bacterium_eval("{genome}"),
         tmp_out=config["mapping_directory"] + "community_PDD/{sample}_rRNA_Aligned.out.bam",
         tmp_log=config["mapping_directory"] + "community_PDD/{sample}_rRNA_Log.final.out",
     priority: 4
     conda:
         "envs/envBash.yml"
     shell:
         """
         mkdir -p {params.out_dir}
         rm -drf "{params.out_prefix}tmp"
         mkdir -p {params.gen_dir}
         srun -p fast -c {threads} STARlong --runThreadN {threads} --readFilesIn {input.in1_1},{input.in2_1} {input.in1_2},{input.in2_2} --readFilesCommand zcat --genomeDir {params.index_dir} --outSAMattrRGline ID:run1 , ID:run2 --outTmpDir "{params.out_prefix}tmp" --outFileNamePrefix {params.out_prefix} --outSAMtype BAM Unsorted --outFilterScoreMinOverLread 0.02 --outFilterMatchNminOverLread 0.02 --outFilterMatchNmin 0 --outFilterMismatchNmax 10 --alignEndsProtrude 10 ConcordantPair --peOverlapNbasesMin 1 --outFilterMultimapNmax 10 --outSAMprimaryFlag AllBestScore
         mv {params.tmp_out} {output.out}
         mv {params.tmp_log} {output.log}
         rm -drf "{params.out_prefix}tmp"
         """
         
         
rule feature_counts:
    input:
        gff=config["genome_gff_directory"] + "community.gff",
        bam=[config["mapping_directory"] + "community_PDD/" + samp for samp in expand("{sample}_Aligned.out.bam", sample=SAMPLES)],
    output:
        config["counts_directory"] + "community_PDD_counts.csv",
    threads: config["parameters"]["threads"]
    priority: 3
    params:
        dir=config["mapping_directory"] + "community_PDD/",
    conda:
        "envs/envBash.yml"
    shell:
        """
        srun -p fast -c {threads} featureCounts -a {input.gff} -o {output} {input.bam} -g "ID" -T {threads} -M -O --fraction -p -t "exon,CDS,misc_RNA,ncRNA,rRNA,tRNA,tmRNA,mRNA" --byReadGroup --extraAttributes "gene,product,product_name,Name,locus_tag,transcriptId"
        """
        

rule feature_counts_rRNA:
    input:
        gff=config["genome_gff_directory"] + "community.gff",
        bam=[config["mapping_directory"] + "community_PDD/" + samp for samp in expand("{sample}_rRNA-Aligned.out.bam", sample=SAMPLES)],
    output:
        config["counts_directory"] + "community_PDD_rRNA-counts.csv",
    threads: config["parameters"]["threads"]
    priority: 3
    params:
        dir=config["mapping_directory"] + "community_PDD/",
    conda:
        "envs/envBash.yml"
    shell:
        """
        srun -p fast -c {threads} featureCounts -a {input.gff} -o {output} {input.bam} -g "ID" -T {threads} -M -O --fraction -p -t "exon,CDS,misc_RNA,ncRNA,rRNA,tRNA,tmRNA,mRNA" --byReadGroup --extraAttributes "gene,product,product_name,Name,locus_tag,transcriptId"
        """


rule mapping_done:
   input:
       expand(config["mapping_directory"] + "community_PDD/{sample}_Aligned.out.bam", sample=SAMPLES),
   output:
       config["mapping_directory"] + "community_PDD/mapping_done",
   priority: 3
   shell:
       """
       echo "All samples have been mapped."
       touch {output}
       """
        

rule mapping_summary:
    input: 
        dir=config["mapping_directory"] + "community_PDD/",
        file=expand(
            config["mapping_directory"] + "community_PDD/{sample}_Log.final.out", sample=SAMPLES
            )[0],
        flag=config["mapping_directory"] + "community_PDD/mapping_done",
    output: config["mapping_directory"] + "summary.txt"
    threads: config["parameters"]["threads"]
    priority: 2
    params:
        dir="\/".join(config["mapping_directory"].split("/")), 
        # this operation add backslashes to escape the slashes in the sed expression
    shell:
        """
        echo "sample" > {output}
        cut -f 1 {input.file} | grep "|" | sed s'/|//'g | awk '{{$1=$1;print}}' >> {output}
        for samp in {input.dir}/*Log.final.out; do
            echo $samp | sed 's/{params.dir}community_PDD\///g' | sed 's/_Log.final.out//g'> tmp
            grep "|" $samp | cut -f 2 >> tmp
            paste -d "\t" {output} tmp > tmp2
            mv -f tmp2 {output}
            rm -f tmp
            rm -f {input.flag}
        done
        """