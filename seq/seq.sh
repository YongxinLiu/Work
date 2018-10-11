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
	# Qaulity access 
	fastqc *.gz -t 99
	# Merge all fastqc report, result in multiqc_report.html
	multiqc .



# 181009.lane14

    # 准备index列表
    cp /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt ~/seq/
    sed -i '1 i IndexID\tIndex\tIndexRC' ~/seq/IlluminaIndex48.txt

    # 移动全部文件至根目录
    cd ~/seq/181009.lane14
    for f in `find . -maxdepth 3 | grep 'fq.gz'`; do
        mv $f ./
    done
    # 改名
    rename 's/FCHFLMTBCX2_L1_Index-/L181009_/' *.gz
    # 获取ID中编号
    ls *.fq.gz|cut -f 2 -d '_'|uniq|sed '1 i IndexRC'>index.txt
    # 添加index的编号
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$3]=$0} NR>FNR {print a[$1]}' ~/seq/IlluminaIndex48.txt index.txt|sed 's/Index//'> index.tsv
    # md5sum
    # 分双端统计md5值
    md5sum *_1.fq.gz > md5sum.txt
    md5sum *_2.fq.gz >> md5sum.txt
    cat md5sum.txt
    # md5值校验
    md5sum -c md5sum.txt > md5sum.check
    cat md5sum.check