############### git ###############
# Make sure you're using SSH URL (you've already set this)
git remote -v  # Should show git@github.com:yemingx/xieyeming1.github.io.git
# Stage your changes
git add --all .
# Commit changes
git commit -m "update"
# Push to remote
git push origin master
git push origin main

# local sync repo by command line
git clone git@github.com:genometube/MICC-seq.git
# local update repo
git fetch --all
git pull origin main

# local sync repo by github desktop app
click "file", then click "clone repository"
username/repository: genometube/MICC-seq
local path: C:\Users\xieyeming1\Documents\GitHub
click "Clone"
# local sync repo by github desktop app
click "Fetch origin"
click "pull origin"
############### rmd ##############
---
author: "yemingxie@gmail.com"
date: "`r format(Sys.time(), '%a %b %d %H:%M:%S %Y')`"
output: 
  github_document:
    html_preview: false
    fig_width: 7
    fig_height: 5
    dev: "png"
---

```{r}
knitr::opts_chunk$set(echo = TRUE)


```

plot_name=RPII_enhancer_micc_violin
outdir=/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/enhancer_Micc_noMicc/
source /mnt/software/anaconda3/bin/activate R4_4
Rscript -e "rmarkdown::render('${plot_name}.Rmd', output_dir = '${outdir}')"
sed -i "s|${outdir}||g" ${outdir}${plot_name}.md

############## hpc ##############
### bio format
sample="SCA_hek_wt_b1"
# ${sample}_R1.fq.gz    ${sample}_R2.fq.gz    ${sample}_fragment.bed.gz

chunk="0"
${sample}_${chunk}
# file name: ${sample}_${chunk}_sorted.bam ${sample}_${chunk}_merged.paf
# readID="C${chunk}R1000"

# meta_table
# file name: ${sample}.meta_table
# column name: sample        seq      cell_line treatment batch experiment_date
#              SCA_hek_wt_b1 sca-seq  hek293t   wt        b1    20231017


### tool
conda="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/bin/conda"
source ~/conda_init.sh
picard='/home/yeming/miniconda3/envs/bioinfo/bin/picard'
cutadapt="/research/zhangchen/software/anaconda3/bin/cutadapt"
cutadapt="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/envs/cutadapt/bin/cutadapt"
bt2='/share/app/bowtie2-2.2.5/bowtie2'
bt2_build='/share/app/bowtie2-2.2.5/bowtie2-build'
bt2='/hwfssz4/BC_PUB/Software/03.Soft_ALL/bowtie2-2.3.4.3/bowtie2'
bt2_build='/hwfssz4/BC_PUB/Software/03.Soft_ALL/bowtie2-2.3.4.3/bowtie2-build'
python='/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/bin/python'
java="/ifswh7/BC_PS/xieyeming/software/java/jre1.8.0_45/bin/java"
Rscript='/mnt/Software/anaconda3/envs/Bioinfo/bin/Rscript'
Rscript="/mnt/software/anaconda3/envs/Bioinfo/bin/Rscript"
Rscript="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/envs/R/bin/Rscript"
Rscript="/research/xieyeming1/software/Miniconda/envs/fyt/bin/Rscript"
Rscript="/mnt/software/anaconda3/envs/R4_4/bin/Rscript"
seqtk='/ifswh7/BC_PS/xieyeming/software/seqtk-master/seqtk'
fastp='/ifswh7/BC_PS/xieyeming/software/fastp/fastp'
fastp='/hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/miniconda3/bin/fastp'
fastqc='/share/app/fastqc/0.11.9/fastqc'
RSeQC='/hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/bioinfo/RSeQC-2.6.4'
featureCounts="/hwfswh2/BC_PUB/Software/07.User-defined/02.Research/chenweitian/subread-2.0.1-source/bin/featureCounts"
featureCounts="/hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/bioinfo/subread-2.0.1-Linux-x86_64/bin/featureCounts"
samtools="/share/app/samtools-1.2/bin/samtools"
samtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/samtools"
bedtools="/share/app/bedtools/2.29.2/bin/bedtools"
bedtools="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/bin/bedtools"
bedtools="/mnt/Software/anaconda3/envs/Bioinfo/bin/bedtools"
bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools"
bedGraphToBigWig='/research/zhangchen/software/bedGraphToBigWig'
bedGraphToBigWig='/ifswh7/BC_PS/zhangchen/software/UCSC_application/bedGraphToBigWig'
bamCoverage="/research/zhangchen/software/anaconda3/bin/bamCoverage"
/share/app/bedtools/2.29.2/bin/multiIntersectBed
bcftools="/zfssz5/BC_PUB/Software/03.Soft_ALL/bcftools-1.5/bcftools"
trinity='/share/app/trinityrnaseq-2.0.6/Trinity'
trinity_rsem="/share/app/trinityrnaseq-2.0.6/util/align_and_estimate_abundance.pl"
porechop="/hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/miniconda3/bin/porechop"
bwa='/share/app/bwa/0.7.17/bwa'
bwa="/research/zhangchen/software/anaconda3/envs/nanopolish/bin/bwa"
STAR='/share/app/STAR-2.40/STAR'
hisat='/share/app/hisat/2.2.1/hisat2'
rsem_cal_exp="/share/app/rsem/1.3.3/bin/rsem-calculate-expression"
conda='/hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/miniconda3/bin/conda'
conda='/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/bin/conda'
in_house_py="/hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/bioinfo/in_house_py"
nanopolish='/hwfssz4/BC_PUB/Software/07.User-defined/02.Research/zhangchen3/anaconda3/envs/nanopolish/bin/nanopolish'
macs='/zfssz4/BC_RD_P1/PROJECT/P19Z12200N0059_liwenqing/P19Z12200N0059_liwenqing/Software/anaconda3/bin/macs2'
macs='/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/envs/macs2/bin/macs2'
macs='/mnt/software/anaconda3/envs/fyt/bin/macs2'
cellranger='/hwfssz4/BC_PUB/pipeline/RNA/10x/cellranger-3.1.0/cellranger'
stringtie='/ifswh1/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2018a/software'
fast5_subset='/ifswh1/BC_RD/RD_COM/USER/chenweitian/software/python/envs/py38/bin/fast5_subset'
fast5_subset="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/bin/fast5_subset"
fast5_subset="/research/zhangchen/software/anaconda3/envs/megalodon/bin/fast5_subset"
fastq_dump='/research/xieyeming1/software/sra_toolkit/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump.3.0.0'
minimap2="/ifswh7/BC_PS/xieyeming/software/minimap2/minimap2"
minimap2="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/envs/bioinfo_py3.3/bin/minimap2"
minimap2="/research/zhangchen/software/anaconda3/bin/minimap2"
k8="/ifswh7/BC_PS/xieyeming/software/minimap2/misc/k8-0.2.4/k8-Linux"
paftool="/ifswh7/BC_PS/xieyeming/software/minimap2/misc/paftools.js"
gatk="/ifswh7/BC_PS/xieyeming/software/gatk/gatk-4.2.6.0/gatk"
greenhill="/research/xieyeming1/software/GreenHill/src"
fastq_dump="/mnt/Software/anaconda3/envs/Bioinfo/bin/fastq-dump"
cooler="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/bin/cooler"
pairix="/ifswh7/BC_PS/xieyeming/software/pairix-master/bin/pairix"
bgzip="/ifswh7/BC_PS/xieyeming/software/pairix-master/bin/bgzip"
pairtools="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/envs/scHiCExplorer/bin/pairtools"
Rscript="/mnt/software/anaconda3/envs/cwf-hicrange/bin/Rscript"
juicer_tools="/research/liwenqing/Software/juicer/juicer_tools.1.9.9_jcuda.0.8.jar"
juicer_tools="/ifswh7/BC_PS/xieyeming/software/juicer/juicer_tools.1.9.9_jcuda.0.8.jar"

hic_pro="/research/xieyeming1/software/hic_pro/HiC-Pro-master/hic_pro"
cooltools="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/cooltools"
cooler="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/cooler"

source /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/zhangchen3/anaconda3/bin/activate
source /ifswh7/BC_PBG_R/fangyitong/Software/anaconda3/bin/activate
source /home/chenwt/software/conda/bin/activate megalodon

# output
############### check monitor error ###############
/ifswh1/BC_PUB/biosoft/pipeline/Package/pymonitor-1.1/monitor.cpu stat

############### nanopolish ##############
#new aliyun
export HDF5_PLUGIN_PATH=/mnt/chenweitian/nanopolish_cpg_gpc/plugin
nanopolish="/mnt/chenweitian/nanopolish_cpg_gpc/nanopolish/nanopolish"

############### conda ###############
conda info --envs
conda env list
conda create --name myenv
conda env remove --name myenv
conda create -n myenv -c conda-forge python=3.6 scipy=0.15.0 astroid babel
conda create -n R -c conda-forge r-base
# check files tar.gz for specific package version
conda create -n linear_algebra python=3.7
pip install vtk==8.1.2
pip install mayavi ipyevents
conda install matplotlib numpy pandas scipy sympy basemap ipython ipympl scikit-learn
conda install freetype
conda install -c conda-forge r-devtools
conda create --name hapcut2 -c bioconda htslib hapcut2
conda create -n hint -c bioconda hint
conda install -c conda-forge --strict-channel-priority r-arrow
conda install r-xml
gChain ComplexHeatmap, VariantAnnotation
install.packages()
"XML","restfulr"
BiocManager::install("biomaRt")
"biomaRt","rtracklayer","BSgenome","GenomicFeatures","VariantAnnotation"
# bioinfo package    conda install -c CHANNEL_NAME PACKAGE_NAME
conda install -c bioconda pysam pybedtools samtools=1.8 sra-tools macs2 bedtools

conda install -c conda-forge glob2 imagemagick scikit-learn jupyterlab

pip install -U seaborn
conda config --add channels conda-forge

# conda config --append channels CHANNEL_NAME

# Install the conda environment from the environment.yml file:
conda env create -f scVI-3D_conda_environment.yml

# For quick python module installation
pip install -r python-requirements.txt

#config module file:
cat /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/module_file/bedtools/2.26.0
#%Module1.0
prepend-path PATH /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/bioinfo/longranger-2.2.2/bedtools/v2.26.0/

# use module
module use --append /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/module_file/
module load bedtools/2.26.0

# conda windows
set PATH=%PATH%;C:\Anaconda3;C:\Anaconda3\Scripts\
conda install -c conda-forge jupyterlab scikit-learn numpy pandas seaborn matplotlib fastcluster

conda install -c anaconda ipykernel
python -m ipykernel install --user --name=<YOUR_CONDA_ENV>
1、多个conda环境下,只需要装一个jupyter notebook,环境的切换是通过切换 kernel实现的；
2、不同环境下,都需要安装 ipykernel,conda install ipykernel
3、不同环境下,需要生成内核 kernel ,你切换到你的环境,然后 python -m ipykernel install --user --name 环境名称 --display-name "在notebook中显示的环境名"
注意：1个环境,对应2个kernel.json, 一个位于package里面,一个在share路径下面。上面 #3 生成内核实际就是生成kernel.json文件。

megalodon                /research/zhangchen/software/anaconda3/envs/megalodon
nanopolish               /research/zhangchen/software/anaconda3/envs/nanopolish
pore-c-snakemake         /research/zhangchen/software/anaconda3/envs/pore-c-snakemake

# /research/xieyeming1/proj_2022/lm_hia5_20220628/dimelo_data_v2_t2t/mega_test_newAliyun.sh
# fix env
conda install -c anaconda libstdcxx-ng

# glibc install https://serverkurma.com/linux/how-to-update-glibc-newer-version-on-centos-6-x/
############### qsub ###############
qsub -cwd -l vf=1G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q -e log -o log [bash.sh]
-cwd run curent working dir; -l resource=value; -q queue list; -e standard error; -o standard output
qstat -j 5595966
qdel 5595966

############### check if PyTorch with CUDA support is installed successfully ###############
import torch
print(torch.__version__)
print(torch.cuda.is_available())
print(torch.version.cuda)

############### fastp fastqc ##############
/mnt/software/anaconda3/envs/Bioinfo/bin/fastp -i Diploid_11_TACGCTGC_GAGCCTTA_R1.fastq.gz test.fq
/mnt/software/anaconda3/envs/Bioinfo/bin/fastqc --noextract --nogroup -o fastqc raw_data/V350212326_L01_read_1.fq.gz

############### samtools ###############
# convert bam to sam
${samtools} view in.bam > out.sam

## check flag
# convert sam to bam
${samtools} view in.sam -o out.bam -O BAM
${samtools} view -bT reference.fa test.sam > test.bam
# add header
${samtools} view -h chunk_53.bam.test > chunk_53.sam.w_header.test

# check flag
${samtools} flagstat in.bam > out.txt
samtools view -F 256 input.bam # only primary alignments
samtools view -f 256 input.bam  # only secondary alignments

# properly paired 
${samtools} view -f 0x0002 
# forward strand
samtools view -F 0x10

# get mapped reads, keep header
samtools view -h -F 4
# get unmapped reads, keep header
samtools view -h -f 4
## remove dup
#sort by read_id(query)
${samtools} sort -n -o sample_Sorted_names.bam -O BAM sample.sam
# fill in mate coordinates and insert size fields
${samtools} fixmate -m sample_Sorted_names.bam sample_Fixmate.bam
# explanation: This will sort based on chromosome number and coordinates
${samtools} sort -o sample_Sorted.bam sample_Fixmate.bam
# explanation: This will remove all the duplicates and also print some basic stats about the result file
${samtools} markdup -r -s sample_Sorted.bam sample_Final_File.bam

# get genome locus sequence
${samtools} faidx genome.fa chr2:98665000-98665100

# get genome locus bam
samtools view input.bam "Chr10:18000-45500" > output.bam

# get depth
samtools depth input.bam > output.depth
samtools depth input.bam | awk '$1 == "chr10" && $2 > 72295300 && $2 < 72297375' 

# view bam file in IGV
samtools sort file_name.bam -O BAM -o sorted_file_name.bam
samtools index sorted_file_name.bam
bamCoverage -b hek293_rnaseq_S1_rsem.transcript.bam -o hek293_rnaseq_S1_rsem.transcript.bw
${bedGraphToBigWig} hek_nome3h_45M_chr7_GpC_methy.10.bedgraph ${ref_fai} hek_nome3h_45M_chr7_GpC_methy.10.bw

# bed to fasta
samtools="/share/app/samtools/1.11/bin/samtools"
ref='/ifswh7/BC_PBG_R/xieyeming/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa'
site=($(cat _50cell_noF_ED2IP1_DE.bed|sed 's/\t/\:/'|sed 's/\t/\-/'))
for i in ${site[@]}; do
  ${samtools} faidx ${ref} ${i} >> _50cell_noF_ED2IP1_DE.bed.fa
done

############### cutadapt ###############
cutadapt="/research/zhangchen/software/anaconda3/bin/cutadapt"
cutadapt="/ifswh7/BC_PS/xieyeming/software/conda/anaconda3/envs/cutadapt/bin/cutadapt"
paste -d $'' R1.fq R2.fq|sed '3~4 s/++/+/g'  > combine.fq
${cutadapt} --rc -j 1 -b "GTTCGTCTTCTGCCGTATGCTCGAGAAGGCTAC;o=33" -o out_combine.fq combine.fq

${cutadapt} --rc -j 1 -B "GTTCGTCTTCTGCCGTATGCTCGAGAAGGCTAC;o=33" -o out_R1.fq -p out_R2.fq R1.fq R2.fq


############### bedtools ###############
${bedtools} intersect -a chr7.fa_200 -b hek_nome3h_45M_chr7_merged.bed -c > hek_nome3h_45M_chr7_merged.200.bedgraph
${bedtools} bamtobed -bed12 -i sorted.bam > sorted.bed
sort -T tmp -k1,1 -k2,2n
${bedtools} intersect -a chr7.fa_${bin_size} -b ${sample}_${type}_site.bed -c -nobuf -sorted> ${sample}_${type}_site.${bin_size}.bedgraph

${bedtools} genomecov -i hek_nome3h_45M_hg19_mt_merged.bed -g hg19mt_genomecov.xls > hek_nome3h_45M_hg19_mt_merged.bed.genomecov

cat /research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa.fai |cut -f1,2 > hg19mt_genomecov.xls
# mean cov
/mnt/Software/anaconda3/envs/Bioinfo/bin/bedtools coverage -mean -a hg19mt_genomecov.xls.bed -b NPC_1.bed > NPC_1.bed.genomecov
# point cov  
/share/app/bedtools/2.29.2/bin/bedtools genomecov -i ${sample}_sorted.bed -g chr7.genome > chr7.cov 

awk '{print $3"\t"$2"\t"$2+1}' # convert site to bed

############### bedtools ###############
/mnt/software/anaconda3/envs/fyt/bin/macs2 bdgpeakcall -i ${sample}.bedgraph -o ${sample}.narrowpeak

############### gatk ###############
${samtools} sort -T ${sample} -o ${sample}.sorted.bam ${sample}.clean.bam
${samtools} index ${sample}.sorted.bam
${java} -jar ${picard} MarkDuplicates INPUT=${sample}.sorted.bam METRICS_FILE=${sample}.sorted.bam.metrics OUTPUT=${sample}.dedup.bam REMOVE_DUPLICATES=true
${java} -jar ${picard} CollectInsertSizeMetrics INPUT=${sample}.dedup.bam OUTPUT=${sample}.insert_size_metrics.txt HISTOGRAM_FILE=${sample}_insert_size_hist.pdf
${gatk} MarkDuplicates --INPUT ${sample}_R1.sorted.bam --METRICS_FILE ${sample}_R1.sorted.bam.metrics --OUTPUT ${sample}_R1.dedup.bam  --REMOVE_DUPLICATES true
############### seqtk ###############
# random 1000 reads
${seqtk} sample -s100 in.fq 1000 > out.fq

sed 's#$#\/1#g' read.id > read_1.id
${seqtk} subseq R1.fq read_1.id > subseq_R1.fq

${seqtk} trimfq -L [retain at most INT bp from the 5-end,int] in.fq > out.fq
${seqtk} trimfq -b [trim INT bp from left,int] -e [trim INT bp from right,int] in.fq > out.fq
seqkit stats *.fa
############### fastp ###############
${fastp} -i R1.fq -I R2.fq -o R1_trim.fq -O R2_trim.fq -h path/to/html_output -j path/to/json_output
${fastp} -M 30 -5 -3 -W 20
${fastp} -u 5 -q 30 

############### blastn ###############
# blast+
/usr/bin/makeblastdb -in ref.fa -dbtype nucl -out blastn_db/ref.fa
/usr/bin/blastn -query query.fa -db blastn_db/ref.fa -dust no -outfmt 6 -evalue 1e-5 -num_alignments 1 -num_descriptions 1 -num_threads 2 -out query.fa.blast_out
# output
# query_seqid subject_seqid pct_ident align_length mismatch gap_open q_start q_end s_start s_end evalue bitscore

############### bowtie2 ###############
${bt2} -x ${ref_dir}/${ref} -1 R1.fq -2 R2.fq -S out.sam
${bt2} -x ${ref_dir}/${ref} -U R1.fq -S out.sam
${bt2} -x genomeindex -U out.fq|${samtools} view -bS - > out.bam
${bt2_build} -q ${ref_dir}/${ref} ${ref_dir}/${ref}
# --no-unal --no-hd --no-sq --score-min C,0,0 -N 0 --end-to-end --norc -L 20

############### featureCounts ###############
${featureCounts} -a ${gtf} -g transcript_id -o count.txt in.sam -F GTF
${featureCounts} -a ${gtf} -g transcript_id -R CORE -o in.sam.CORE in.sam -F GTF

############### bwa ###############
${bwa} mem -t 20 ${ref} R1.fq R2.fq > out.sam
${bwa} index ref.fa

############### minimap ###############
${minimap2} -t ${thread} -a -x map-ont ${ref} ${fastq_pass}/${chip_id}_${k}.fastq > chunk/${i}/${i}.sam

############### STAR ###############
# make ref
${STAR} --runThreadN 5 --runMode genomeGenerate --genomeDir STARIndex \
--genomeFastaFiles genome.fa --sjdbOverhang 100 --sjdbGTFfile genes.gtf
# align
${STAR} --runThreadN 5 --genomeDir STARIndex \
--sjdbGTFfile genes.gtf --sjdbOverhang 100 \
--readFilesIn R1.fastq.gz R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix path/to/output/prefix.

############### hisat ###############
${hisat} --sensitive -p ${thread} -x ${hisat_index} -U ${read_2_dir}/${read_2} -S ${sam} 

############### rsem ###############
/share/app/rsem/1.3.3/bin/rsem-prepare-reference --bowtie2 --bowtie2-path /share/app/bowtie2-2.2.5/ hg19_mRNA.fa hg19_mRNA.fa
 ${rsem_cal_exp} --bowtie2 --bowtie2-path /share/app/bowtie2-2.2.5/ -p ${thread} --paired-end ${sample}_R1.clean.fq ${sample}_R2.clean.fq ${ref} ${sample}_rsem

############### cgat ###############
cgat gtf2gtf

############### trinity ###############
additional_trinity_parameter="--bypass_java_version_check --min_kmer_cov 3"
${trinity} ${additional_trinity_parameter} --seqType fq --max_memory 1G \
       --left R1.fq \
       --right R2.fq \
       --CPU 1 --output ${out_trinity_cov_dir}/trinity_${j}

${trinity} ${additional_trinity_parameter} --seqType fq --max_memory 1G \
       --single R1.fq \
       --CPU 1 --output ${out_trinity_cov_dir}/trinity_${j}
# needs permission to fastq_dir, --max_memory required, output folder name must start with "trinity"
${trinity_rsem} \
  --transcripts ${out_trinity_cov_dir}/trinity_${j}/Trinity.fasta \
  --left R1.fq \
  --right R2.fq \
  --seqType fq --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference \
  --output_dir ${out_trinity_cov_dir}/trinity_${j}/rsem

############### nanopolish ###############
git clone --recursive https://github.com/jts/nanopolish.git
sudo make
#sudo apt-get install zlib1g-dev
mv ~/Downloads/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/ /usr/local/
https://github.com/nanoporetech/vbz_compression/releases
export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin

############### sam2tsv ################
sam2tsv='/research/xieyeming1/software/sam2tsv/jvarkit/dist/sam2tsv.jar'
 ${picard} CreateSequenceDictionary -R ${ref}
 ${samtools} faidx ${ref}
java -jar ${sam2tsv} -R ${ref} ${_02_mut_quant}/${lane}.bam | cut -f -4,6-| awk '{if ($9 ~/M|=|X/ ) print $0}' > ${_02_mut_quant}/${lane}_both_strand_all.tsv

${java} -jar ${sam2tsv} -R ${ref} ctcfBMporeC.sorted.chunk_0.bam --regions hek293_atac_chr7.bed |awk '$10 == "M"'|cut -f1,2,4- |head
 ${samtools} view -M -L in.bed -O BAM in.bam | java -jar dist/sam2tsv.jar (...)

############### fastq-dump ###############
fastq-dump -I --split-files SRR390728

############### pycharm ###############
# To select all occurrences in the file, press Ctrl+Alt+Shift+J
# ctrl + /, comment a block of code
############### docker ###############
# check image
docker images

# remove image
docker rmi [image_id]

# download image
sudo docker pull bdgenomics/rhapsody:1.5.1

# view files or scripts in docker images
docker run -it 3a59c950d406 ls -lh
docker run -it 3a59c950d406 find . -type f -name *.py

sudo service docker restart    # Fix internet error
docker image prune    # Remove unused images
docker system prune
docker system df

############### aws ###############
aws configure # set Access Key Id and Secret Access Key
aws s3api list-buckets --query "Buckets[].Name" # list bucket name
aws s3 ls s3://yanlab.data.share --recursive --human-readable --summarize # list files in bucket
aws s3 cp --recursive s3://yanlab.data.share/ . # download files in bucket

############### aspera ##############
wget https://download.asperasoft.com/download/sw/connect/3.9.9/ibm-aspera-connect-3.9.9.177872-linux-g2.12-64.tar.gz
ascp -i [私钥] -T -K 1 -l [最大传输速度] [下载地址及SPA数据编号] [下载输出位置]

############### wget ##############
/mnt/software/anaconda3/envs/chm/bin/axel -n 12

############### .parquet ###############
conda install -c conda-forge parquet-tools # /research/zhangchen/software/anaconda3/envs/pore-c-snakemake/bin/parquet-tools
parquet-tools csv NlaIII_run02_draft1_unphased.contacts.parquet|head
parquet-tools csv -n 3 NlaIII_run02_draft1_unphased.contacts.parquet

############### BGI_pipe ###############
# RNA-seq
/ifswh1/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a_test/example

############### in_house_bash ###############
# qsub -cwd -l vf=1G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q 
# compress all .fq and remove .fq under given DIR
sh /ifswh7/BC_PBG_R/xieyeming/scripts/bash/compress.sh INDIR suffix

# whole genome .pairs to .mcool
/research/xieyeming1/proj_2022/zhichao_smHiC_20230313/chia_pet_hic/bedpe_to_cooler_CTCF/cooler_newAliyun.sh

# hires reformat
/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hires/clean_RNA_Vala/clean_RNA_Vala_single.sh

# cellranger RNA
/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hires/10x_RNA_format/test_10x_pipe/reformat/reformat.sh
/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hires/10x_RNA_format/hires_10x_pipe/cellranger_count.sh

# dipC hickit


### down sample ###
/ifswh7/BC_PBG_R/xieyeming/proj_2023/pubData_ngsHiC/juicer/juicer_single_hek293_ngsHiC_3M.sh

############### in_house_tools ###############
/ifswh7/BC_PBG_R/xieyeming/scripts/

/hwfssz4/BC_PUB/Software/07.User-defined/02.Research/xieyeming1/bioinfo/in_house_py/dist.py
/ifswh7/BC_PBG_R/xieyeming/proj_2023/SCA_seq_demo_2023/03.adjacent_contact_2D/hek_nome3h_45M_hg19mt/multi_pairwise.sh
/research/xieyeming1/proj_2022/zhichao_smHiC_20230116/dedup_combine_analysis/cartesian_join_V1.py
/research/xieyeming1/proj_2022/fy_MATdamNGS_20230313/fy_MATdamNGS_20230313/03.frag_bed_count_/merge_overlap_range_by_read.py
/research/xieyeming1/proj_2022/fy_MATdam_20230711/subset.py
/research/xieyeming1/proj_2023/fy_laminDamNGS_20230829/02.frag_bed_count/subset_V1.py
/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/dip_C/raw_data/download.sh
/research/xieyeming1/db/genome/hg19/Annotation/feature_bed/motif_hg19/feature_bed_hg19.py

# R2 single barcode demultiplex
/ifswh7/BC_PBG_R/xieyeming/proj_2022/yaning_293tNGSatac_20220802/demultiplex/cut_tag/test_run.sh
# R2 i5 + i7 barcode demultiplex
/research/xieyeming1/proj_2023/zhichao_porecHiC_20230724/demultiplex/L01/demultiplex_i5_i7.sh
# R1 10x + R2 i5 + i7 barcode demultiplex
/research/xieyeming1/proj_2022/zhichao_smHiC_20230313/zhichao_smHiC_20230313/01.clean_fq/L01/clean_fq.sh

# bedpe to cooler and .hic chr7
/ifswh7/BC_PBG_R/xieyeming/proj_2023/guomei_circHiC_20231113/chiapet_v3/bedpe_to_cooler_hic_guomei_L02/cooler_hic.sh

############### in_house_R ###############
# create personal library and install package
#hic pipeline
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install("GenomicInteractions")


############### useful resources ###############
https://www.encodeproject.org/files/
https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/

############### archive ###############
/hwfssz4/BC_COM_P1/RD_P1/P19Z12200N0089_xym/P19Z12200N0089_xym
/zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_xieyeming1
/zfssz4/BC_RD_P1/PROJECT/P18Z12200N0396_xieyeming1
db='/zfssz4/BC_RD_P1/PROJECT/fangyitong/shared_database'
attribute="/zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_xieyeming1/archive/genome/hg19/Annotation/attribute_table"


############### bgi pipe ###############
/ifswh4/BC_PUB_T1/Pipeline/RNA/DrTom_bigRNA/Pipe_for_RNA/example/
/ifswh4/BC_PUB_T1/Pipeline/EPI/ATAC-seq/2017a/example/run.sh
/ifswh4/BC_PUB_T2/BioSysDB/BGI/v2201/Index/9606/NCBI/GCF_000001405.25_GRCh37.p13

perl /ifswh4/BC_PUB_T1/Pipeline/EPI/ATAC-seq/2017a/ATAC-seq_2017a.pl -i config


