#1. 创建conda环境
##1.1 创建新的环境
conda create -n gcf python=3
##1.2 激活新创建的环境
envac gcf
##1.3 把pip的默认源改为国内(例如阿里云)
pip3 config set global.index-url https://mirrors.aliyun.com/pypi/simple/

#2. 安装配置bigslice
##2.1 安装bigslice
pip3 install bigslice
##2.2 安装hmmer
conda install hmmer
##2.3 检测是否已下载HMM database
if [ ! -d /dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/envs/gcf/bin/bigslice-models ]
then
download_bigslice_hmmdb
fi
##2.4 测试bigslice是否配置完成
bigslice --version

#3. 为新环境安装R和R包
##3.1 安装R
conda install -c skyblued r-base
##3.2 安装stringi
conda install -c conda-forge r-stringi
##3.3 在R中安装DBI, circlize, iNEXT, RSQLite
install.packages("DBI")
install.packages("RSQLite")
install.packages("circlize")
install.packages("iNEXT")