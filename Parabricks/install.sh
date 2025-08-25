# https://docs.nvidia.com/clara/parabricks/latest/gettingstarted/gettingthesoftware.html
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1

docker run \
      --gpus all \
      --rm \
      --volume $(pwd):/workdir \
      --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1 \
    pbrun fq2bam \
      --ref /workdir/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta \
      --in-fq /workdir/parabricks_sample/Data/sample_1.fq.gz /workdir/parabricks_sample/Data/sample_2.fq.gz \
      --out-bam /outputdir/fq2bam_output.bam

# wget -O parabricks_sample.tar.gz \
# "https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz"

# docker run \
#       --gpus all \
#       --rm \
#       --volume $(pwd):/workdir \
#       --volume $(pwd):/outputdir \
#     nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1 \
#     pbrun fq2bam \
#       --ref /workdir/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta \
#       --in-fq /workdir/parabricks_sample/Data/sample_1.fq.gz /workdir/parabricks_sample/Data/sample_2.fq.gz \
#       --out-bam /outputdir/fq2bam_output.bam --low-memory


yum-config-manager --add-repo https://nvidia.github.io/nvidia-docker/centos7/x86_64/nvidia-docker.repo
distribution=$(. /etc/os-release;echo $ID$VERSION_ID) \
   && curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.repo | sudo tee /etc/yum.repos.d/nvidia-docker.repo
sudo yum install -y nvidia-container-toolkit
sudo systemctl restart docker


# /research/xieyeming1/proj_2025/wj_chip4micc_250527/align/align_single.sh
ref="/research/xieyeming1/db/genome/hg19/Sequence/bt2_index_hg19mt/genome.fa"
sample=65
R1="/research/xieyeming1/proj_2025/wj_chip4micc_250527/raw_data/250523_M022_V350344234_L01_CWHPE25050299-${sample}/V350344234_L01_${sample}_1.fq.gz"
R2="/research/xieyeming1/proj_2025/wj_chip4micc_250527/raw_data/250523_M022_V350344234_L01_CWHPE25050299-${sample}/V350344234_L01_${sample}_2.fq.gz"
mkdir -p tmp

docker run \
      --gpus all \
      --rm \
      --volume $(pwd):/workdir \
      --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1 \
    pbrun fq2bam \
      --ref ${ref} \
      --in-fq ${R1} ${R2} \
      --out-bam /outputdir/fq2bam_output.bam --low-memory --tmp-dir tmp

