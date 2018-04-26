#Configure, softwares and databases

# Configure



# Softwares

## BioConda

#https://bioconda.github.io/
cd ~/Downloads
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
#确认许可协议，默认安装目录为conda，不添加环境变量
bash Miniconda3-latest-Linux-x86_64.sh # yes /conda no 
# 手动设置3个常用命令到环境，添加目录会替换为Python3环境
ln /conda/bin/conda /usr/local/bin/
ln /conda/bin/activate /usr/local/bin/
ln /conda/bin/deactivate /usr/local/bin/
#添加bioconda，越靠后优先越高
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ 
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ 
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ 
conda config --set show_channel_urls yes
conda config --add channels r # Optional
#显示已有的通道
conda config --get channels

#安装bwa
conda install bwa
#建立环境并安装系列工具
conda create -n aligners bwa bowtie hisat2 star
source activate aligners
source deactivate
#py27
conda create -n py27 python=2.7
source activate py27
conda install deeptools RseQC 
source deactivate

### QIIME1.9.1
#http://qiime.org/install/install.html
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose # -c bioconda
source activate qiime1
print_qiime_config.py -t
source deactivate
conda remove --name qiime1 --all






## R

## Python2.7

## Python3.5

## Other




# Databases
