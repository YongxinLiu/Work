# ath -- Arabidopsis thaliana




# rice -- Oryza Sativa


# wheat -- Triticum aestivum L. 

# 华大基因测序数据
## 2017/10/16 景美苜蓿真菌
nohup wget --user 20171013F17FTSNCKF3414 --password LIBkboR ftp://cdts-wh.genomics.cn/F17FTSNCKF3414_LIBkboR/Raw/CWHPEPI00001583/CWHPEPI00001583_1.fq.gz &bg
nohup wget --user 20171013F17FTSNCKF3414 --password LIBkboR ftp://cdts-wh.genomics.cn/F17FTSNCKF3414_LIBkboR/Raw/CWHPEPI00001583/CWHPEPI00001583_2.fq.gz &bg # 无效，不支持wget


## 171121华大两个lane数据
cd ~/seq
mkdir 171121
cd 171121
#nohup wget -dcr --user  20171013F17FTSNCKF3414 --password LIBpfjR ftp://cdts-wh.genomics.cn/F17FTSNCKF3414_LIBpfjR/Clean/ &bg # 530 Login incorrect
for file in `find ./ -name *.gz`; do
	fastqc $file &
done
# 查找姜婷Index1，正向，反向分别尝试；参考序列~/ref/culture/IlluminaIndex48.txt
ls|grep 'CGTGAT' # result无结果
ls|grep 'ATCACG' # 找到结果 
mv *L2* ~/seq/171122.wheatNP.Aq/
cd ~/seq/171122.wheatNP.Aq
for file in `find | grep 'L2'| grep 'fq.gz'`; do
	mv $file ./
done
rename 's/FCHY55LBCXY_L2_CWHPEPI00001611_Index-//' *.gz
fastqc *.fq.gz -t 96 &
multiqc . # 汇总评估报告



## 171212华大两个lane数据
# 筛选第二个lane的数据rice minicore到171213
cd ~/seq/171212
mkdir -p ~/seq/171213riceMinicore/
for file in `find | grep 'L2'| grep 'fq.gz'`; do
	mv $file ~/seq/171213riceMinicore/
done
cd ~/seq/171213riceMinicore/
rename 's/FCH3L7FBCX2_L2_CWHPEI17110010-//' *.gz
# 将index取反向互补
format_seq_revcom.pl -i lane4.index
format_seq_revcom.pl -i lane5.index


## 171212 华大lane7
cd ~/seq/171220.lane7.GuoXX.JiangT
# 简化文件名
rename 's/171214_I188_FCH3MHMBCX2_L1_CWHPEPI00001645/lane7/g' *.gz
# 制作反向互补的index与测序结果一致
format_seq_revcom.pl -i index.txt
# 拆分文库: grep处匹配需要提示less查看文件头的具体格式
for index in `cat index.txt.out`; do
	echo ${index}
	zcat lane7_1.fq.gz|grep -A 3 "#${index}AT_TCTTTCCC/"|grep -v -P '^--$' |gzip > ${index}_1.fq.gz &bg
	zcat lane7_2.fq.gz|grep -A 3 "#${index}AT_TCTTTCCC/"|grep -v -P '^--$' |gzip > ${index}_2.fq.gz &bg
done
fastqc *.gz -t 96 # 质控所有数据
multiqc . # 汇总评估报告
cut -f 1,16 multiqc_data/multiqc_fastqc.txt|grep '_1'|sed 's/_1//;s/\.0//' > multiqc_data/multiqc_fastqc.count
awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $0,a[$1]}' multiqc_data/multiqc_fastqc.count index.txt.ou


## 171221 华大lane6 重测
cd ~/seq/171221.lane6.reseq/raw/lane6
rename 's/171217_I191_FCH3N7TBCX2_L2_CWHPEPI00001634/lane/g' *.gz
md5sum *.gz > md5.txt # 获得文件的md5植
md5sum -c md5.txt # 检查原始md5，确定文件传输是否正确
fastqc *.gz -t 2 &
cat >index.txt # 粘贴index列表
format_seq_revcom.pl -i index.txt # 翻转序列
for index in `cat index.txt.out`; do
	echo ${index}
	zcat lane_1.fq.gz|grep -A 3 "#${index}AT_TCTTTCCC/"|grep -v -P '^--$' |gzip > ${index}_1.fq.gz &bg
	zcat lane_2.fq.gz|grep -A 3 "#${index}AT_TCTTTCCC/"|grep -v -P '^--$' |gzip > ${index}_2.fq.gz &bg
done


## 180108 华大lane8
cd ~/seq/180108.lane8.rice/lane8
md5sum *.gz
format_seq_revcom.pl -i index.txt
# 拆分文库: grep处匹配需要提示less查看文件头的具体格式
for index in `cat index.txt.out`; do
	echo ${index}
	zcat lane_1.fq.gz|grep -A 3 "#${index}AT_TCTTTCCC/"|grep -v -P '^--$' |gzip > ${index}_1.fq.gz &bg
	zcat lane_2.fq.gz|grep -A 3 "#${index}AT_TCTTTCCC/"|grep -v -P '^--$' |gzip > ${index}_2.fq.gz &bg
done
fastqc *.gz -t 96 # 质控所有数据
multiqc . # 汇总评估报告
cut -f 1,8 multiqc_data/multiqc_fastqc.txt|grep '_1'|sed 's/_1//;s/\.0//' > multiqc_data/multiqc_fastqc.count
awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $0,a[$1]}' multiqc_data/multiqc_fastqc.count index.txt.out


## 180210.lane9.ath3T
cd ~/seq/180210.lane9.ath3T/Clean/CWHPEPI00001683/
rename 's/FCHCMCCBCX2_L1_CWHPEPI00001683/lane/g' *.gz
md5sum *.gz > md5.txt # 获得文件的md5植
fastqc *.gz -t 2 &
cat >index.txt # 粘贴index列表
format_seq_revcom.pl -i index.txt # 翻转序列
for index in `cat index.txt.out`; do
	echo ${index}
	zcat lane_1.fq.gz|grep -A 3 "#${index}ATC_TCTTTCCC/"|grep -v -P '^--$' |gzip > ${index}_1.fq.gz &bg
	zcat lane_2.fq.gz|grep -A 3 "#${index}ATC_TCTTTCCC/"|grep -v -P '^--$' |gzip > ${index}_2.fq.gz &bg
done

# 生成3个库的测试文件
zless lane_1.fq.gz|head -n 10000|gzip > lane1_1.fq.gz
zless lane_2.fq.gz|head -n 10000|gzip > lane1_2.fq.gz


## 180301 华大lane10
cd ~/seq/180301.lane10/Clean/CWHPEPI00001684
rename 's/FCH7L5VBCX2_L1_CWHPEPI00001684/lane/g' *.gz
md5sum *.gz > md5.txt & # 获得文件的md5植
fastqc *.gz -t 2 &
cp /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt index.txt # 筛选为本实验使用的
nohup time parallel -j 32 "zcat lane_1.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > {1}_1.fq" ::: `cat index.txt` &bg
nohup time parallel -j 32 "zcat lane_2.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > {1}_2.fq" ::: `cat index.txt` &bg

## 2018/5/28 华大lane11
cd ~/seq
mv F18FTSNCKF1459_LIBulyR/ 180528.lane11
cd 180528.lane11/Clean/CWHPEPI00001823/
rename 's/FCH7GG7BCX2_L2_CWHPEPI00001823/lane_64/' *.gz
# 格式转换
time fastp -i lane_64_1.fq.gz -I lane_64_2.fq.gz -o lane_1.fq.gz -O lane_2.fq.gz -6 -A -G -Q -L -w 9
cut -f 1,3 /mnt/bai/xiaoning/seq_raw_data/180528.lane11/Clean/CWHPEPI00001823/index.txt | sed '1 i LibraryID\tIndexRC' > library.txt
# 并行拆库
mkdir -p seq
parallel --xapply -j 32 "zcat lane_1.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > seq/{1}_1.fq" ::: `tail -n+2 library.txt | cut -f 2`
parallel --xapply -j 32 "zcat lane_2.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > seq/{1}_2.fq" ::: `tail -n+2 library.txt | cut -f 2`

# 180601.ath.WGS

## 诺合数据下载工具下载和安装
	# https://helpcdn.aliyun.com/document_detail/50452.html
	cd ~/bin
	wget http://docs-aliyun.cn-hangzhou.oss.aliyun-inc.com/assets/attach/50452/cn_zh/1524643963683/ossutil64?spm=a2c4g.11186623.2.6.TfLBVN -O ossutil
	chmod +x ossutil

## 数据下载指定目录

	# 进入screen环境防断网
	screen -d -r seq
	# 设置下载目录，每次不同数据请修改此处
	wd=~/seq/180601.ath.WGS/
	mkdir -p $wd
	# 创建配置文件，不同用户名请修改用户名和密码
	ossutil config -e oss.aliyuncs.com -i LTAIhcpcybDqnbxI -k c3oJcIIDzc9z3H7f8CENdC4taSqUoG
	# 下载文件
	ossutil cp oss://novo-data-nj/customer-sMIXvUoG/ $wd -r -f --jobs 3 --parallel 2



# 180719华大lane12

	cd /mnt/bai/yongxin/seq/180719.nrt1.1a.lane12/
	tar xvzf upload.tar.gz
	# 使用远程桌面查看分析报告数据量、质量评估结果
	cd /mnt/bai/yongxin/seq/180719.nrt1.1a.lane12/Clean/AC/
	# 简单文件名
	rename 's/FCH7F2JBCX2_L1_CWHPE18070021-//g' *.gz
	# 检查数据质量格式：33可以，64还需转换为33
	determine_phred-score.pl AAAATG_1.fq.gz

## 180816 lane13 sd1
    
    # 设置批次编号
    wd=L180816
    cd ~/seq/$wd
    # :附录：查看质量报告和评估文件完整性
    # 统一lane文件为当前目录的Lane_1/2.fq.gz
    mv Clean/D20170718/FCHFLH2BCX2_L2_CWHPE18070298_* ./
    rename 's/FCHFLH2BCX2_L2_CWHPE18070298/lane/' *.fq.gz
	# 附录：检查数据质量，如果64跳到附录2
    # 附录：准备Index列表拆分，需要IndexRC列，可由ID或Index列检索到
    tail -n+2 ~/rice/zn.sd1/b2/doc/library.txt |awk '{print "Index"$5"\t"$3"\t"$2}' > index.txt
    # 附录：并行拆库，注意必须有wd变量，且indexRC位于第三列
    # 附录：质控与压缩

## 181009.lane14

    # 移动全部文件至根目录
    cd ~/seq/181009.lane14
    for f in `find . -maxdepth 3 | grep 'fq.gz'`; do
        mv $f ./
    done
    # 改名
    rename 's/FCHFLMTBCX2_L1_Index-/L181009_/' *.gz


## 181023WangGD

    # 简化目录和文件名
    mv Clean/CCPM/* ./
    rename 's/FCHL5W2BCX2_L2_CWHPE18090036/lane/' *.fq.gz
    # 并行拆库
    prefix=L181023
    

## lane16 181024WangChao

    # 简化目录和文件名
    mv Clean/DYMSLR18/FCHFLMNBCX2_L2_CWHPE18100007_* ./
    rename 's/FCHFLMNBCX2_L2_CWHPE18100007/lane/' *.fq.gz
    # 并行拆库
    prefix=L181024


## 181025.lane15

    # 简化目录和文件名
    mv Clean/SDRS20180912/FCHL72LBCX2_L2_CWHPEPI00001955_* ./
    rename 's/FCHL72LBCX2_L2_CWHPEPI00001955/lane/' *.fq.gz
    # 按附录代码操作

    ## 郭晓璇的Index13没有结果，沟通后尝试Index11的RC GGCTAC在lane中检索发现较多，修改index列表最后一列，重新分析
    rm L181025_AGTCAA* # 删除错误Index相关文件
    idx=GGCTAC
    prefix=L181025


## 181122.lane17
    
    wd=L181122
    cd ~/seq
    mv F18FTSNCKF3152_LIBagiR $wd
    cd $wd
    # 简化目录和文件名
    mv Clean/G5WH_WWAT/FCHHH57BCX2_L2_CWHPEPI00001961_* ./
    rename 's/FCHHH57BCX2_L2_CWHPEPI00001961/lane/' *.fq.gz
    # 按附录代码操作



## 扩增子lane处理通用代码

### 1. 华大lane返回通用处理代码

    # 查看质量报告和评估文件完整性
    tar xvzf upload.tar.gz # 解压报告
    firefox upload/index.html # firefox、下载或远程桌面查看
    md5sum -c md5.txt # 检查下载数据准确性，一致输出OK

	# 检查数据质量，如果64跳到附录2
	determine_phred-score.pl lane_1.fq.gz

    # 准备Index列表拆分，需要IndexRC列，可由ID或Index列检索到
	# 情况1：如保存ID列为index.ID，获取index和IndexRC，并补充到实验设计
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' \
        /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt index.ID > index.txt
    # 情况2：如保存Index序列，保存为index，获取indexID和IndexRC，并补充到实验设计
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$2]=$0} NR>FNR {print a[$1]}' \
        /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt index > index.txt
    cat index.txt
    # 保存至库文件和登记 SeqLibraryList.xlsx

    # 并行拆库，按目录名见附件编号A18xxxx，j可按文库数量调整，推荐默认24
    parallel -j 24 "zcat lane_1.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > ${wd}_{1}_1.fq" ::: `cut -f 3 index.txt | grep  -v '^$'` 
    parallel -j 24 "zcat lane_2.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > ${wd}_{1}_2.fq" ::: `cut -f 3 index.txt | grep  -v '^$'`

    # 质控和报告汇总
    fastqc -t 48 L*.fq.gz
    multiqc ./ # 详见multiqc_report.html
    # 提取各样品数据量，不同批数据量的列会有变化，可能要更改列的数值5、6等
    l=`head -n1 multiqc_data/multiqc_fastqc.txt|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep 'Total Sequences'|cut -f 1`
    cut -f 1,${l} multiqc_data/multiqc_fastqc.txt | sed 's/_.\t/\t/' | uniq > datasize.txt
    cat datasize.txt

    # 压缩空间及上传
    pigz *.fq
    # 分双端统计md5值
    md5sum L*_1.fq.gz > /tmp/md5sum1.txt
    md5sum L*_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt


### 2. lane文件格式64转33标准化

	# 如果为64，改原始数据为33
	rename 's/lane/lane_33/' lane_*
	# 关闭质量控制，主要目的是格式转换64至33，而不是过滤序列，否则扩增子必须成对不满足则usearch无法合并
	time fastp -i lane_64_1.fq.gz -I lane_64_2.fq.gz \
		-o lane_1.fq.gz -O lane_2.fq.gz -6 -A -G -Q -L -w 9
	# 1lane 80GB, 2 threads, 102min



# meta 

    meta@meta:/mnt/m2/data/meta/

## 2017/12/11 ath 2.5T
    
    cd /mnt/m2/data/meta/ath/2.5T
    # 只有去宿主后的clean_data 137G，而且重复了两份clean_data(tar包含有fq, fastqc, log和out有总体宿主或剩余数量)和remove_host_fq，没有去宿主前和统计信息
    # 样本已经为单样一个文件，直接改名
    mkdir -p seq
    awk 'BEGIN{OFS=FS="\t"}{system("mv remove_host_fq/"$3"/"$3".rmhost.clean.1.fq.gz seq/"$1"_1.fq.gz ");system("mv remove_host_fq/"$3"/"$3".rmhost.clean.2.fq.gz seq/"$1"_2.fq.gz ");}' <(tail -n+2 metadata.txt)
    # 统计样本md5值和数据量见附录代码


## 2017/12/27 rice nrt1.1b
    
    cd /mnt/m2/data/meta/rice/nrt1.1b
    # 只有去宿主后的clean_data 137G，而且重复了两份clean_data(tar包含有fq, fastqc, log和out有总体宿主或剩余数量)和remove_host_fq，没有去宿主前和统计信息
    # 样本已经为单样一个文件，直接改名
    mkdir -p seq
	# 按实验设计按文件夹批量合并再改名，需要输入文件每个样本一个目录
	p=30
    for i in `tail -n+2 metadata.txt|cut -f3`; do
        zcat `find 1.clean_reads/${i}/ -name *.gz | grep '_1.fq'` | pigz -p ${p} > seq/${i}_1.fq.gz &
        zcat `find 1.clean_reads/${i}/ -name *.gz | grep '_2.fq'` | pigz -p ${p} > seq/${i}_2.fq.gz &
    done
    awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$3"_1.fq.gz seq/"$1"_1.fq.gz ");system("mv seq/"$3"_2.fq.gz seq/"$1"_2.fq.gz ");}' <(tail -n+2 metadata.txt)
    # 统计样本md5值和数据量见附录代码


## 2018/11/19 ath 3T
    
    cd /mnt/m2/data/meta/ath/3T
    # 整理代码见 /mnt/m1/yongxin/ath/3T/manual.md
    for i in `tail -n+2 design.txt|cut -f3`; do
        zcat `find 01.filter/${i}/ -name *.gz | grep '_1.fq'` | pigz -p ${p} > seq/${i}_1.fq.gz &
        zcat `find 01.filter/${i}/ -name *.gz | grep '_2.fq'` | pigz -p ${p} > seq/${i}_2.fq.gz &
    done
    awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$3"_1.fq.gz seq/"$1"_1.fq.gz ");system("mv seq/"$3"_2.fq.gz seq/"$1"_2.fq.gz ");}' <(tail -n+2 design.txt)


## 2018/11/21 rice miniCore

    cd /mnt/m2/data/meta/rice/miniCore
    # 了解文件格式
    less -S 00.rawdata1/data.info # 样本和180个原始数据对应表
    find 01.cleandata/*/* -name *.gz | wc -l # 整理好的180个对应文件，己去除宿主
    # 整理原始序列按样本合并、上传、改名
    find 00.rawdata*/*/* -name *.gz|wc # 725文件，单数吗？
    find 00.rawdata*/*/* -name *.gz|cut -f 3 -d '/'|sort|uniq|wc -l # 720，有5个重复，应该是重复多次
    # 非冗余整理到同一文件夹中
    mkdir -p temp
    for i in `find 00.rawdata*/*/* -name *.gz`; do
        ln -sf `pwd`/$i temp/
    done
    # 制作实验设计批量两文件合并并重命名
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$2"_1.fq.gz temp/"$3"_1.fq.gz | pigz -p 30 > seq/"$1"_1.fq.gz")}' <(tail -n+2 metadata.txt|cut -f 1,9,10)
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$2"_2.fq.gz temp/"$3"_2.fq.gz | pigz -p 30 > seq/"$1"_2.fq.gz")}' <(tail -n+2 metadata.txt|cut -f 1,9,10)


## 2018/11/28 rice wetdry
    
    # 制作实验设计改名并移动至seq目录中
    cd /mnt/m2/data/meta/rice/wetdry
    mkdir -p seq
    awk 'BEGIN{OFS=FS="\t"}{system("mv 2.cleandata/"$2"/"$2"_1.clean.fq.gz seq/"$1"_1.fq.gz");system("mv 2.cleandata/"$2"/"$2"_2.clean.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 design.txt)


## 附录: 宏基因组通用操作代码
    
    # 质控和报告汇总
    cd seq
    fastqc -t 48 *.fq.gz
    multiqc ./ # 详见multiqc_report.html
    # 提取各样品数据量，不同批数据量的列会有变化，可能要更改列的数值5、6等
    l=`head -n1 multiqc_data/multiqc_fastqc.txt|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep 'Total Sequences'|cut -f 1`
    cut -f 1,${l} multiqc_data/multiqc_fastqc.txt | sed 's/_.\t/\t/' | uniq > datasize.txt
    cat datasize.txt

    # 分双端统计md5值
    md5sum L*_1.fq.gz > /tmp/md5sum1.txt
    md5sum L*_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt



# rna

    meta@meta:/mnt/m2/data/rna/

## 2017/11/27 王鑫拟南芥无菌苗 /mnt/m2/data/rna/171127sterility/ 诺禾致源测序


## 2018/11/20 王鑫拟南芥时间序列 /mnt/m2/data/rna/181120TimeCourse 美格基因测序
    
    # 数据位于meta服务器meta用户的/mnt/m2/data/rna/181120TimeCourse目录[meta@meta:/mnt/m2/data/rna/181120TimeCourse]$
    sed -i 's/\r//' md5sum.txt # 文件尾有windows换行符，替换后才能运行md5sum
    # 整理文件名与样本名对照表
    rename 's/_\w+_R/_/' *.gz
    rename 's/_001.fastq/.fq/' *.gz # 两次删除和替换与库名一致
    # FTP上传GSA/181120RNASeq


## 附录. RNA-Seq通用操作代码

    # 测试文件完整性
    md5sum -c md5sum.txt
    # 根据实验设计批量改名
    mkdir -p seq
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$9"_1.fq.gz temp/"$10"_1.fq.gz | pigz -p 16 > seq/"$1"_1.fq.gz")}' <(tail -n+2 metadata.txt) &
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$9"_2.fq.gz temp/"$10"_2.fq.gz | pigz -p 16 > seq/"$1"_2.fq.gz")}' <(tail -n+2 metadata.txt) &

    # 质控和报告汇总
    fastqc -t 48 *.fq.gz
    multiqc ./ # 详见multiqc_report.html
    # 提取各样品数据量，不同批数据量的列会有变化，可能要更改列的数值5、6等
    l=`head -n1 multiqc_data/multiqc_fastqc.txt|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep 'Total Sequences'|cut -f 1`
    cut -f 1,${l} multiqc_data/multiqc_fastqc.txt | sed 's/_.\t/\t/' | uniq > datasize.txt
    cat datasize.txt

    # 分双端统计md5值
    md5sum L*_1.fq.gz > /tmp/md5sum1.txt
    md5sum L*_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt
