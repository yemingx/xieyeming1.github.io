# https://docs.nvidia.com/clara/parabricks/latest/gettingstarted/gettingthesoftware.html
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1

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


yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
distribution=$(. /etc/os-release;echo $ID$VERSION_ID) \
   && curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.repo | sudo tee /etc/yum.repos.d/nvidia-docker.repo
sudo yum install -y nvidia-container-toolkit
sudo systemctl restart docker

sudo rpm -ivh libnvidia-container-tools-1.13.5-1.x86_64.rpm
sudo rpm -ivh libnvidia-container1-1.13.5-1.x86_64.rpm
sudo rpm -ivh nvidia-container-toolkit-1.13.5-1.x86_64.rpm
sudo rpm -ivh nvidia-container-toolkit-base-1.13.5-1.x86_64.rpm


# change yum repo
# https://blog.csdn.net/2302_77114273/article/details/142769901?spm=wolai.workspace.0.0.92e02944vpD6Lu
mv /etc/yum.repos.d/CentOS-Base.repo /etc/yum.repos.d/CentOS-Base.repo.bak
wget -O /etc/yum.repos.d/CentOS-Base.repo http://mirrors.aliyun.com/repo/Centos-7.repo
yum clean all

###add below to  /etc/yum.repos.d/CentOS-Base.repo###
[centos-sclo-rh]
name=CentOS-7 - SCLo rh
baseurl=http://vault.centos.org/centos/7/sclo/$basearch/rh/
gpgcheck=1
enabled=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo
 
[centos-sclo-sclo]
name=CentOS-7 - SCLo sclo
baseurl=http://vault.centos.org/centos/7/sclo/$basearch/sclo/
gpgcheck=1
enabled=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo
###add above to  /etc/yum.repos.d/CentOS-Base.repo###

yum makecache


# install docker
# https://developer.aliyun.com/article/1551022
sudo yum update -y
sudo yum install -y yum-utils device-mapper-persistent-data lvm2
sudo yum-config-manager --add-repo http://mirrors.aliyun.com/docker-ce/linux/centos/docker-ce.repo
sudo yum makecache fast
sudo yum install -y docker-ce docker-ce-cli containerd.io
sudo systemctl start docker
sudo systemctl enable docker

sudo mkdir -p /etc/docker
sudo tee /etc/docker/daemon.json <<-'EOF'
{
  "registry-mirrors": ["https://w0pc1i5g.mirror.aliyuncs.com"]
}
EOF

sudo systemctl daemon-reload
sudo systemctl restart docker

docker -v
sudo docker run --rm hello-world

# install nvidia-container-runtime 
# https://blog.51cto.com/kusorz/13135792
tar -zxvf nvidia-container-runtime.tar.gz
cd nvidia-container-runtime
rpm -Uvh --force --nodeps *.rpm
systemctl restart docker
whereis nvidia-container-runtime




docker pull nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1

# /research/xieyeming1/proj_2025/wj_chip4micc_250527/align/align_single.sh
# ref="/research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa"
ref="/research/xieyeming1/software/parabricks/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta"
sample=65
R1="/research/xieyeming1/proj_2025/wj_chip4micc_250527/raw_data/250523_M022_V350344234_L01_CWHPE25050299-${sample}/V350344234_L01_${sample}_1.fq.gz"
R2="/research/xieyeming1/proj_2025/wj_chip4micc_250527/raw_data/250523_M022_V350344234_L01_CWHPE25050299-${sample}/V350344234_L01_${sample}_2.fq.gz"
outdir=/research/xieyeming1/self_study/git/xieyeming1.github.io/Parabricks/out
mkdir -p tmp
docker run \
      --gpus all \
      --rm \
      --volume $(pwd) \
      --volume $(pwd) \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
    pbrun fq2bam \
      --ref /research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa \
      --in-fq ${R1} ${R2} \
      --out-bam outputdir/fq2bam_output.bam --low-memory --tmp-dir tmp



docker run --gpus all nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 pbrun fq2bam -h