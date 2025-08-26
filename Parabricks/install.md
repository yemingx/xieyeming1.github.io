# change yum repo
refernece link: https://blog.csdn.net/2302_77114273/article/details/142769901?spm=wolai.workspace.0.0.92e02944vpD6Lu

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


# install docker
refernece link: https://developer.aliyun.com/article/1551022

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
refernece link: https://blog.51cto.com/kusorz/13135792

tar -zxvf nvidia-container-runtime.tar.gz
cd nvidia-container-runtime
rpm -Uvh --force --nodeps *.rpm
systemctl restart docker
whereis nvidia-container-runtime
