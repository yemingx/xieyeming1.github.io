# 在阿里云CentOS 7.9服务器上搭建A10 GPU环境

## 修改yum源
参考链接: https://blog.csdn.net/2302_77114273/article/details/142769901

mv /etc/yum.repos.d/CentOS-Base.repo /etc/yum.repos.d/CentOS-Base.repo.bak
wget -O /etc/yum.repos.d/CentOS-Base.repo http://mirrors.aliyun.com/repo/Centos-7.repo
yum clean all

###add this section to  /etc/yum.repos.d/CentOS-Base.repo###

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

###add this section to  /etc/yum.repos.d/CentOS-Base.repo###

yum makecache


## install docker
参考链接:  https://developer.aliyun.com/article/1551022

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

## install nvidia-container-runtime 
参考链接:  https://blog.51cto.com/kusorz/13135792

tar -zxvf nvidia-container-runtime.tar.gz
cd nvidia-container-runtime
rpm -Uvh --force --nodeps *.rpm
systemctl restart docker
whereis nvidia-container-runtime

# run Parabricks demo
参考链接:  https://docs.nvidia.com/clara/parabricks/4.0.1/tutorials/fq2bam_tutorial.html

wget -O parabricks_sample.tar.gz \
"https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz"

tar xvf parabricks_sample.tar.gz

## test with 10 cpu threads
ref=parabricks_sample/Ref/Homo_sapiens_assembly38.fasta 
fq1=parabricks_sample/Data/sample_1.fq.gz 
fq2=parabricks_sample/Data/sample_2.fq.gz 
date
/usr/bin/bwa mem -t 10 ${ref} ${fq1} ${fq2} > aligned2.sam
samtools view --threads 10 -Sb aligned2.sam > aligned2.bam
date

## 10 cpu threads out log, used 15min wall clock time
Tue Aug 26 11:41:08 CST 2025
[main] Version: 0.7.17-r1188
[main] CMD: /usr/bin/bwa mem -t 10 parabricks_sample/Ref/Homo_sapiens_assembly38.fasta parabricks_sample/Data/sample_1.fq.gz parabricks_sample/Data/sample_2.fq.gz
[main] Real time: 791.180 sec; CPU: 7934.364 sec
Tue Aug 26 11:56:05 CST 2025

## test with one A10 GPU
systemctl is-active docker
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1
mkdir -p tmp
docker run \
      --gpus all \
      --rm \
      --volume $(pwd):/workdir \
      --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
    pbrun fq2bam \
      --ref /workdir/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta \
      --in-fq /workdir/parabricks_sample/Data/sample_1.fq.gz /workdir/parabricks_sample/Data/sample_2.fq.gz \
      --out-bam /outputdir/fq2bam_output.bam --tmp-dir tmp


## one A10 GPU out log, used 3min wall clock time
[Parabricks Options Mesg]: Checking argument compatibility
[Parabricks Options Mesg]: Automatically generating ID prefix
[Parabricks Options Mesg]: Read group created for /workdir/parabricks_sample/Data/sample_1.fq.gz and
/workdir/parabricks_sample/Data/sample_2.fq.gz
[Parabricks Options Mesg]: @RG\tID:HK3TJBCX2.1\tLB:lib1\tPL:bar\tSM:sample\tPU:HK3TJBCX2.1
[PB Info 2025-Aug-26 03:32:20] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:32:20] ||                 Parabricks accelerated Genomics Pipeline                 ||
[PB Info 2025-Aug-26 03:32:20] ||                              Version 4.0.0-1                             ||
[PB Info 2025-Aug-26 03:32:20] ||                       GPU-BWA mem, Sorting Phase-I                       ||
[PB Info 2025-Aug-26 03:32:20] ------------------------------------------------------------------------------
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[PB Info 2025-Aug-26 03:32:23] GPU-BWA mem
[PB Info 2025-Aug-26 03:32:23] ProgressMeter    Reads           Base Pairs Aligned
[PB Info 2025-Aug-26 03:32:42] 5043564          580000000
[PB Info 2025-Aug-26 03:32:59] 10087128 1160000000
[PB Info 2025-Aug-26 03:33:15] 15130692 1740000000
[PB Info 2025-Aug-26 03:33:32] 20174256 2320000000
[PB Info 2025-Aug-26 03:33:48] 25217820 2900000000
[PB Info 2025-Aug-26 03:34:04] 30261384 3480000000
[PB Info 2025-Aug-26 03:34:21] 35304948 4060000000
[PB Info 2025-Aug-26 03:34:37] 40348512 4640000000
[PB Info 2025-Aug-26 03:34:54] 45392076 5220000000
[PB Info 2025-Aug-26 03:35:10] 50435640 5800000000
[PB Info 2025-Aug-26 03:35:24] 
GPU-BWA Mem time: 180.801377 seconds
[PB Info 2025-Aug-26 03:35:24] GPU-BWA Mem is finished.


[main] CMD: /usr/local/parabricks/binaries//bin/bwa mem -Z ./pbOpts.txt /workdir/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta /workdir/parabricks_sample/Data/sample_1.fq.gz /workdir/parabricks_sample/Data/sample_2.fq.gz @RG\tID:HK3TJBCX2.1\tLB:lib1\tPL:bar\tSM:sample\tPU:HK3TJBCX2.1
[main] Real time: 184.391 sec; CPU: 3267.032 sec
[PB Info 2025-Aug-26 03:35:24] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:24] ||        Program:                      GPU-BWA mem, Sorting Phase-I        ||
[PB Info 2025-Aug-26 03:35:24] ||        Version:                                           4.0.0-1        ||
[PB Info 2025-Aug-26 03:35:24] ||        Start Time:                       Tue Aug 26 03:32:20 2025        ||
[PB Info 2025-Aug-26 03:35:24] ||        End Time:                         Tue Aug 26 03:35:24 2025        ||
[PB Info 2025-Aug-26 03:35:24] ||        Total Time:                            3 minutes 4 seconds        ||
[PB Info 2025-Aug-26 03:35:24] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:24] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:24] ||                 Parabricks accelerated Genomics Pipeline                 ||
[PB Info 2025-Aug-26 03:35:24] ||                              Version 4.0.0-1                             ||
[PB Info 2025-Aug-26 03:35:24] ||                             Sorting Phase-II                             ||
[PB Info 2025-Aug-26 03:35:24] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:24] progressMeter - Percentage
[PB Info 2025-Aug-26 03:35:24] 0.0       0.00 GB
[PB Info 2025-Aug-26 03:35:34] Sorting and Marking: 10.000 seconds
[PB Info 2025-Aug-26 03:35:34] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:34] ||        Program:                                  Sorting Phase-II        ||
[PB Info 2025-Aug-26 03:35:34] ||        Version:                                           4.0.0-1        ||
[PB Info 2025-Aug-26 03:35:34] ||        Start Time:                       Tue Aug 26 03:35:24 2025        ||
[PB Info 2025-Aug-26 03:35:34] ||        End Time:                         Tue Aug 26 03:35:34 2025        ||
[PB Info 2025-Aug-26 03:35:34] ||        Total Time:                                     10 seconds        ||
[PB Info 2025-Aug-26 03:35:34] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:34] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:34] ||                 Parabricks accelerated Genomics Pipeline                 ||
[PB Info 2025-Aug-26 03:35:34] ||                              Version 4.0.0-1                             ||
[PB Info 2025-Aug-26 03:35:34] ||                         Marking Duplicates, BQSR                         ||
[PB Info 2025-Aug-26 03:35:34] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:34] progressMeter -  Percentage
[PB Info 2025-Aug-26 03:35:44] 45.5      8.65 GB
[PB Info 2025-Aug-26 03:35:54] 100.0     0.00 GB
[PB Info 2025-Aug-26 03:35:59] BQSR and writing final BAM:  25.067 seconds
[PB Info 2025-Aug-26 03:35:59] ------------------------------------------------------------------------------
[PB Info 2025-Aug-26 03:35:59] ||        Program:                          Marking Duplicates, BQSR        ||
[PB Info 2025-Aug-26 03:35:59] ||        Version:                                           4.0.0-1        ||
[PB Info 2025-Aug-26 03:35:59] ||        Start Time:                       Tue Aug 26 03:35:34 2025        ||
[PB Info 2025-Aug-26 03:35:59] ||        End Time:                         Tue Aug 26 03:35:59 2025        ||
[PB Info 2025-Aug-26 03:35:59] ||        Total Time:                                     25 seconds        ||
[PB Info 2025-Aug-26 03:35:59] ------------------------------------------------------------------------------
Please visit https://docs.nvidia.com/clara/#parabricks for detailed documentation


