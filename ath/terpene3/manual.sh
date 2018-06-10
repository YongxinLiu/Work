#!/bin/bash

# 16S扩增子操作手册 16S Amplicon pipeline manual, v1.4 2017/12/23

# 按如下命令执行，并修改makefile中的参数
# 进入工作目录
cd /mnt/bai/yongxin/ath/jt.HuangAC/batch3all


# 0. 配置实验设计参数

# 复制最新版流程配置文件+操作手册至工作目录
cp ~/rice/timecourse/ma* ./

# 修改工作目录wd
pwd

# 初始化工作目录
make init

# 人工编写实验设计，按design.xlsx中要求每个表另存为txt格式上传于doc/文件夹中
# 编辑并上传实验设计至doc/，必选(mappingfiles/summary/material/library/group_compare), 可选(group_venn, group_tern)

# 如文库较多，可将所有库保存于design_raw.txt中，使用split_design.pl命令生成按库分割的mappingfile
dos2unix doc/*
split_design.pl -i doc/design_raw.txt
# 基于mappingfile生成实验设计，需要指定一个存在文件提取header，默认为L1
cat <(head -n1 doc/L1.txt|sed 's/#//g') doc/L* |grep -v '#' > doc/design.txt
# 检查实验设计首行是否名称唯一
wc -l doc/design.txt && cut -f 1 doc/design.txt | sort | uniq | wc -l

# 上传原始数据*.fq.gz至clean_data/
for i in `tail -n+2 doc/library.txt|cut -f 2`;do
	ln ~/seq/180108.lane8.rice/lane8/${i}_1.fq.gz clean_data/${i}_1.fq.gz
	ln ~/seq/180108.lane8.rice/lane8/${i}_2.fq.gz clean_data/${i}_2.fq.gz
done

# 简化文件名(可选) 
rename 's/Med-//g;s/_HLY73BCXY_L1//g;s/clean\.//g' *.gz 
# 原始数据文件名重命名为L1/2/3_1/2.fq.gz，并格式化文本信息
make rename
# 修改文库列表list
# ll clean_data/*.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|sort -u|tr "\n" " "
tail -n+2 doc/library.txt |cut -f 1|tr "\n" " " # 另一种方法，但需要有library.txt
# 显示实验设计,修改实验的主要和次要分组类型g1/g2，如groupID, genotype, compartment, soiltype，没有g2可使用默认batch
head -n 3 doc/design.txt
# 获取主、次、第三分组信息：具体根据实验设计进一步筛选
tail -n+2 doc/design.txt|cut -f 5|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
# 最简单的实验只需一组信息，第二组填batch和1；复杂实验需要次、三分组信息
# 仔细全文阅读配置文件，修改所有与本项目相关参数



# 1. 数据质控
make qc # 原始数据质控
make merge # 双端数据合并再质控
make extract_barcodes # 提取样本的barcodes
make split # 拆分样品
make cutadapt # 切除双端接头
make stat # 统计文库各步操作stat each library

make multiqc # 汇总数据合并前后的评估报告

grep -v -c '#' doc/L* | sed 's/doc\/L//;s/.txt:/\t/' | sort -k1,1n # 统计各库中样品数据
# 展示库数据量，数据列可变，有可能为13，18，19，20
cut -f 1,20 clean_data/multiqc_data/multiqc_fastqc.txt|grep '_1'|sed 's/L//;s/_1//;s/\.0//'|sort -k1,1n 

# 2. 标准流程
make merge_library # 合并所有文库，并简单统计绘图






# 从这里开始：手动合并两批样品
cp ../batch2/doc/design.txt doc/designb2.txt
cp ../batch3/doc/design.txt doc/designb3.txt
# 检查标题行，不对应，需要修改design2
head -n1 doc/designb2.txt doc/designb3.txt
# 修改为v1版本
dos2unix doc/* # 去除windows字符
head doc/*v1.txt # 检查表头是否一致
cat doc/designb2_v1.txt <(tail -n+2 doc/designb3_v1.txt) > doc/design.txt # 生成新design


cat ../batch2/temp/seqs_usearch.fa ../batch3/temp/seqs_usearch.fa > temp/seqs_usearch.fa
make derep # 数据去冗余
make unoise # 鉴定生物特征序列
make rm_chimeras # (可选)去嵌合体，执行完Unoise后，不执行此步可直接跳过
make otu_table # 生成OTU表
make assign_tax # 基于gg13.8注释物种
make tree # 构建进化树，用于多样性分析
make alpha # 计算alpha多样性
make alpha_usearch_rare # (可选)，采用usearch结果绘制
make alpha_qiime_rare # (可选)，采用qiime进行稀释曲线分析，大量抽样比较费时；可继续运行其它步骤
make beta # 计算beta多样性
make graphlan # 绘制高丰度物种树
make ggtree # 绘制高丰度进化树
make diversity # 多样性绘制箱线和散点

cat result/qc.sum # 显示数据量统计


# 3. 个性分析
make filter # 过滤高丰度菌分析，默认0.05%。如lefse只有1%，ruben文章常用0.5-0.1%
make rediv # 新表重新多样性分析
make draw_div # 多样性绘图
make draw_tax # 物种门-属级别比较
make draw_otu # OTUs水平差异比较
make draw_ter # OTUs三元图比较

make alpha_usearch_rare # usearch rarefraction
make alpha_qiime_rare # alpha rarefraction

# 功能差异分析
make picrust_calc
make picrust_draw # 结果差异不明显，没有之前2半萜例子明显？

make rmd # 生成可读性报告 write report 



# 4. 可选高级分析
make culture_graphlan # 绘制高丰度菌可培养比例物种树
# 选择97%的结果位于result/otu_cultured.txt，原始比对结果在result/rep_seqs.blastn，OTU1/2无差异很大没有培养，统计结果filter_""_k1/culture.txt

cat filter_""_k1/culture.txt # 查看统计信息

# 4.1 添加土和三萜中差异的菌
batch=""
#for batch in Root; do
sample=/mnt/bai/yongxin/ath/jt.HuangAC/batch3all/filter_${batch}_k1
# 筛选样品中对应的库中OTU和相似度
awk '$3*$13>=9700' ${sample}/rep_seqs.blastn|cut -f 1-3 |sed 's/\tOTU_/\t/g' > ${sample}/otu_cultured.txt
# 添加自然样品丰度，存在一个菌对应多个OTU，按丰度从小到大排序，读取为保留最大的
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {printf ($0"\t%.4f\n",(a[$1]*100))}' ${sample}/otu_table.txt.mean ${sample}/otu_cultured.txt | sort -k4,4n | sed '1 s/^/OTU\tStock\tSimilar\t\%\n/' > ${sample}/otu_cultured.info

awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {printf a[$1]"\n"}' ${sample}/../result/rep_seqs_tax.txt ${sample}/otu_cultured.info | cut -f 2-8 |sed "1 s/^/Kindom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/" > ${sample}/temp.a

# 写每个样品对应的结果，然后删除

## 添加b3Colvsb3BS
#awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$4"\t"$5"\t"$6} NR>FNR {print a[$1]}' ${sample}/../result_k1-c/otu_b3Colvsb3BS ${sample}/otu_cultured.info | sed '1 s/$/Log2FC\tPvalue\tb3Colvsb3BS/' > ${sample}/temp.b3Colvsb3BS
#
## 引用分组变量
#i=b3Colvsb3BS
#awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$4"\t"$5"\t"$6} NR>FNR {print a[$1]}' ${sample}/../result_k1-c/otu_${i} ${sample}/otu_cultured.info | sed "1 s/$/Log2FC\tPvalue\t${i}/" > ${sample}/temp.${i}

# 添加每个基因型比较的差异
# 循环group_compare.txt
for i in `sed 's/\t/vs/g' doc/group_compare.txt`; do
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$4"\t"$5"\t"$6} NR>FNR {print a[$1]}' ${sample}/../result_k1-c/otu_${i} ${sample}/otu_cultured.info | sed "1 s/$/Log2FC\tPvalue\t${i}/" > ${sample}/temp.${i}
done
paste ${sample}/otu_cultured.info ${sample}/temp.* > ${sample}/otu_cultured.xls

# 添加最新版本4基因差异注释
batch=b3_3
sample=/mnt/bai/yongxin/ath/jt.HuangAC/batch3all/filter_${batch}_k1
mkdir -p ${sample}
cp filter__k1/otu_cultured.info filter_${batch}_k1/
# 循环group_compare.txt
for i in `sed 's/\t/vs/g' doc/b3_3/group_compare.txt`; do
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$4"\t"$5"\t"$6} NR>FNR {print a[$1]}' ${sample}/../result_k1-c/otu_${i} ${sample}/otu_cultured.info | sed "1 s/$/Log2FC\tPvalue\t${i}/" > ${sample}/temp.${i}
done
paste ${sample}/otu_cultured.info ${sample}/temp.* > ${sample}/otu_cultured_4gene.xls

# 更新菌保资源为Root最新200kb库
# 添加最新版本4基因差异注释
batch=Root
sample=/mnt/bai/yongxin/ath/jt.HuangAC/batch3all/filter_${batch}_k1
mkdir -p ${sample}
rm ${sample}/temp.*
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {printf ($0"\t%.4f\n",(a[$1]*100))}' ${sample}/otu_table.txt.mean ${sample}/otu_cultured.txt | sort -k4,4n | sed '1 s/^/OTU\tStock\tSimilar\t\%\n/' > ${sample}/otu_cultured.info
# 循环group_compare.txt
for i in `sed 's/\t/vs/g' doc/b3_5/group_compare.txt`; do
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} NR>FNR {print a[$1]}' ${sample}/../result_k1-c/otu_${i} ${sample}/otu_cultured.info | sed "1 s/$/A_mean\tB_mean\tLog2FC\tPvalue\tFDR\t${i}/" > ${sample}/temp.${i}
done
paste ${sample}/otu_cultured.info ${sample}/temp.* > ${sample}/otu_cultured_4gene.xls


# 4.2 检查物种注释准确性
# OTU_56/62/63/188序列只差1个碱基，但物种连门都不同。doc/实验菌注释A144 A215.xlsx
grep 'OTU_56\t' result/rep_seqs_tax.txt


# 4.3 手动比较组 compare_edgeR.r 新增wilcon秩和检验，结果在*_wilcon.txt
mkdir -p compare
Rscript compare_edgeR.r -i result_k1-c/otu_table.txt -d doc/design.txt -n groupID -c b3ThasKO2-b3Col

# 与edgeR负二项分布比较 pvalue共有的比较
# 添加p值至差异OTUs
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print $0,a[$1]}' \
	<(cut -f 1-4 compare/edgeR_b3ThasKO2-b3Col_wilcon.txt) \
	<(cut -f 1-8 compare/edgeR_b3ThasKO2-b3Col_sig.txt) \
	> compare/edgeR_b3ThasKO2-b3Col_lrt_wilcon.txt

# 维恩图比较ThasKO2 vs Col下两种方法的OTUs一致性。
# edgeR: http://210.75.224.110/report/16s/AC_b3_all_b3 thasko2 p&fdr<0.05, lfc>1.3, 151,73;
# wilcox: http://210.75.224.110/report/16s/AC_b3_all_b3_2/ thasko2 p<0.05, fdr<0.2, lfc>1.3, 95,150
# Venny中比较，图片保存于comapre中。
# 上调中共有62个，占wilcon的65%；其中OTU_1只有wilcon中特有；
# 下调中共有40个，占edgeR中54%，下调重合更少，表明：不稳定
# 使用edgeR输出pvalue和FDR时，与原始只有Pvalue的比较
# 之前的edgeR中使用的fdr和fold change分别为0.05和1.3倍，现改为0.05，与wilcon一致比较
# http://210.75.224.110/report/16s/AC_b3_all_b3_2_lrt p<0.05, fdr<0.2, lfc>1.3, 333,251
# 上调中共有91个，占wilcon的95.7%；其中OTU_1只有wilcon中特有；
# 下调中共有114个，占ewilcon的76%，下调重合更少，表明：不稳定



# 4.4 检查OTU_1是否在其它菌库存在

# 查找某菌在拟南芥中是否分离
mkdir -p culture
grep -A 5 ">OTU_1$" result/rep_seqs.fa > culture/OTU_1.fa # 筛选OTU1，多行fasta，V5-V7通常为5行

# 德国菌库：最近也只有94%
blastn -query culture/OTU_1.fa -db /mnt/bai/yongxin/other/baiyang/16S_nature/Root_isolates_WGS.fasta -out culture/OTU_1.blastn -num_alignments 3 -num_descriptions 3 -evalue 1e-10 -num_threads 9 
cat culture/OTU_1.blastn

# 查找某菌在所有菌库中是否分离:也只有94%
for i in ath medicago rice tomato wheat; do
	blastn -query culture/OTU_1.fa -db /mnt/bai/yongxin/culture/${i}/result/culture_select.fa -outfmt 6 \
	-out culture/OTU_1.${i}.blastn -num_alignments 3 -evalue 1e-10 -num_threads 9 
	echo culture/OTU_1.${i}.blastn
	cat culture/OTU_1.${i}.blastn
done


# 4.5 快速查询一组OTU是否可培养和相对丰度
cat >temp/temp # 保存OTU至文件
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' filter__k1/otu_cultured.info temp/temp



# 4.6 多组共有上下调OTU热图展示看规律 heatmap_multip_group.r

# 4.6.1 三组差异OTU共有绘图c("b3Col","b3ThasKO2","b3ACT2KO")
cut -f 1,7 <(cat result_k1-c/otu_b3ThasKO2vsb3Col result_k1-c/otu_b3ACT2KOvsb3Col)|grep -P 'enriched|depleted'|cut -f 1|sort|uniq|sed '1 s/^/otu\n/'> scripts/list_3.txt
# 默认聚类OTU和分组结果均可以，但对数2/10结果观察无规律，关闭行聚类更好些

# 4.6.2 五组差异OTU共有绘图 c("b3Col","b3ThasKO1","b3ThasKO2","b3ACT2KO","b3ACT2CR")
cut -f 1,7 <(cat result_k1-c/otu_b3ThasKO2vsb3Col result_k1-c/otu_b3ThasKO1vsb3Col result_k1-c/otu_b3ACT2KOvsb3Col result_k1-c/otu_b3ACT2CRvsb3Col)|grep -P 'enriched|depleted'|cut -f 1|sort|uniq|sed '1 s/^/otu\n/'> scripts/list_5.txt
wc -l scripts/list_5.txt # 609个OTUs
cut -f 1,7 <(cat result_k1-c/otu_b3ThasKO2vsb3Col result_k1-c/otu_b3ThasKO1vsb3Col result_k1-c/otu_b3ACT2KOvsb3Col result_k1-c/otu_b3ACT2CRvsb3Col)|grep -P 'enriched|depleted'|sort|uniq|sed '1 s/^/otu\tlevel\n/'> scripts/list_5.txt
wc -l scripts/list_5.txt # 612个OTUs，有3个不同类型
cut -f 1  scripts/list_5.txt |sort|uniq -d # 查看不稳定的OTU，居然有5/123/1269
awk '{a[$1]=$2} END {for (i in a) print i"\t"a[i]}' scripts/list_5.txt | sed '1 s/^/otu\tlevel\n/'> scripts/list_5_unique.txt # 随机去掉不稳定的类型

# 4.6.3 五组差异OTU共有绘图 c("b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR")
cut -f 1,7 <(cat result_k1-c/otu_b3ThasKO2vsb3Col result_k1-c/otu_b3ThahKOvsb3Col result_k1-c/otu_b3ThadKOvsb3Col result_k1-c/otu_b3ACT2CRvsb3Col)|grep -P 'enriched|depleted'|sort|uniq|sed '1 s/^/otu\tlevel\n/'> scripts/list_4.txt
wc -l scripts/list_4.txt # 633个OTUs，有3个不同类型
cut -f 1  scripts/list_4.txt |sort|uniq -d > temp/confilt.txt # 查看不稳定的OTU，居然有 1283 274 450 86
cut -f 1  scripts/list_4.txt |sort|uniq -u > temp/unique.txt # 查看不稳定的OTU，居然有 1283 274 450 86
# awk '{a[$1]=$2} END {for (i in a) print i"\t"a[i]}' scripts/list_4.txt | sed '1 s/^/otu\tlevel\n/'> scripts/list_4_unique.txt # 随机去掉不稳定的类型
awk 'NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' scripts/list_4.txt temp/unique.txt>scripts/list4_unique.txt

# 4.7 以为Thas had act2 thah 共6个基因型, pcao, tax, man, heatmap, venny，版本b3_4


# 4.7.1 去除ThasKO1(弱，安诚一直用KO2)和ACT2KO(强，但ACT2为下游终产物，强弱下结果弱更合理)；更改ACT2CR为KO
# 添加KO2差异表达结果中有各组均值
sed -i '1 s/^/otu\t/' result_k1-c/database.txt # 补齐列名
head -n1 result_k1-c/database.txt|tr '\t' '\n'|awk '{print NR"\t"$0}' # 查看列号
cut -f 1,2,9-13 result_k1-c/database.txt > result_k1-c/database.txt2 # 筛选ID和平均及各组平均
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result_k1-c/database.txt2 result_k1-c/otu_b3ThasKO2vsb3Col_enriched.xls > culture/otu_b3ThasKO2vsb3Col_enriched.txt
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result_k1-c/database.txt2 result_k1-c/otu_b3ThasKO2vsb3Col_depleted.xls >  culture/otu_b3ThasKO2vsb3Col_depleted.txt


# 4.7.2 修正volcano图下边空的问题，修改代码默认Y轴最小值为5

head -n1 result_k1-c/man_otu_b3ACT2KOvsb3Col.txt|tr '\t' '\n'|awk '{print NR"\t"$0}'



# 2018/5/11 所有OTU各种差异列表，添加可培养注释

# 所有差异OTUs合并
for i in `sed 's/\t/vs/g' doc/b3_5/group_compare.txt`; do
	cut -f 1-7 result_k1-c/otu_${i} | sed "1 s/A_mean\tB_mean/${i}/" | sed 's/vs/\t/' | less -S > temp/temp.${i}
done
paste temp/temp.* > manual/DA_otu_all.xls




# 文章图表整理 fig/

## Fig3. 菌群差异分析

### a CPCoA 整体
beta_cpcoa.sh -i `pwd`/result/otu_table.txt -m '"bray","jaccard"' \
	-d `pwd`/doc/design.txt -A groupID -B '"b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO"' -E TRUE \
	-o `pwd`/fig/1/ -w 5 -h 3 -s 7
# 基此基础上修改script/beta_cpcoa_batch.R，添修改design.txt中添加batch2/3两种分组方案，其中batch3的结果更好。

### b. stackplot of taxonomy 门纲
修改为6组，重复计算taxonom并保存文件"b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO","b3BS"，修改sub中的b3_5，为b3_2为4组
rm draw_tax ; make draw_tax
查找result_k1-c/tax_stack_phylumpro_abc.pdf

### c manhattan plot 目
原始数据位于 http://210.75.224.110/report/16s/AC_b3_all_b3_3_lrt 新版放大在http://bailab.genetics.ac.cn/report/16s/AC_b3_all_KO_b3_5，颜色与b不一致

http://bailab.genetics.ac.cn/report/16s/AC_b3_all_KO_b3_2 # 5组，选择 phylum或order，建议order

### d. 科水平维恩图

需要提供数据重画
	cd ~/ath/jt.HuangAC/batch3all/fig/3
	batch_venn.pl -i venn.txt -d family-venn.txt

### e. OTU水平共有

http://bailab.genetics.ac.cn/report/16s/AC_b3_all_b3_3_lrt_soil # 维恩部分

	# 共有部分再添加饼形图，注释是否要根中特异的
	# 1. 四组共同上调
	cd ~/ath/jt.HuangAC/batch3all
	tail -n 28 AC_b3_all_b3_3_lrt_soil/result_k1-c/otu.txt.vennb3ThasKO2vsb3Col_enrichedb3ThahKOvsb3Col_enrichedb3ThadKOvsb3Col_enrichedb3ACT2KOvsb3Col_enriched.xls.xls > fig/3/venn_pie/enriched.txt # 共有上调OTU列表，无表头

	tail -n 11 AC_b3_all_b3_3_lrt_soil/result_k1-c/otu.txt.vennb3ThasKO2vsb3Col_depletedb3ThahKOvsb3Col_depletedb3ThadKOvsb3Col_depletedb3ACT2KOvsb3Col_depleted.xls.xls> fig/3/venn_pie/depleted.txt # 共有上调OTU列表，无表头

	cat AC_b3_all_b3_3_lrt_soil/result_k1-c/otu_b3Colvsb3BS_enriched.xls <(tail -n+2 AC_b3_all_b3_3_lrt_soil/result_k1-c/otu_b3Colvsb3BS_depleted.xls) > fig/3/venn_pie/colvsroot.txt # 标注类型的
	


## Fig4. 抑菌实验

	cd ~/ath/jt.HuangAC/batch3all/fig/4
	format_seq2fasta.pl -i "*.seq" -o candidate.fa
	makeblastdb -in ~/ath/jt.HuangAC/batch3all/result_k1-c/rep_seqs.fa -dbtype nucl
	blastn -query candidate.fa -db ~/ath/jt.HuangAC/batch3all/result_k1-c/rep_seqs.fa -out candidate.blastn -outfmt 6 -num_alignments 30 -evalue 1 -num_threads 9
	awk '$3>97' candidate.blastn > candidate.txt
	cut -f 2 candidate.txt|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
	# 根据候选前3绘制箱线图 result_k1-c/otu_table.txt百分比， result/otu_table.txt 1万抽样
	cd ~/ath/jt.HuangAC/batch3all
	alpha_boxplot.sh -i `pwd`/result/otu_table_norm.txt -m '"OTU_311","OTU_103","OTU_231","OTU_501","OTU_69","OTU_584","OTU_422","OTU_324","OTU_172","OTU_70","OTU_239","OTU_365","OTU_21","OTU_179","OTU_16","OTU_345","OTU_177","OTU_125","OTU_561","OTU_112","OTU_25","OTU_479","OTU_323","OTU_232","OTU_66","OTU_556","OTU_166","OTU_205","OTU_163","OTU_80","OTU_46","OTU_200","OTU_578","OTU_229","OTU_804","OTU_264","OTU_116","OTU_1073","OTU_262","OTU_24","OTU_521","OTU_215","OTU_775","OTU_44","OTU_1532"' \
	-d `pwd`/doc/design.txt -A groupID -B '"b3BS","b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO"' \
	-o `pwd`/fig/4/box  -w 5 -h 3 -s 7 -t TRUE -n FALSE
	# 2018/6/4 绘制OTU_136
	alpha_boxplot.sh -i `pwd`/result/otu_table_norm.txt -m '"OTU_136"' \
	-d `pwd`/doc/design.txt -A groupID -B '"b3BS","b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO"' \
	-o `pwd`/fig/4/box  -w 5 -h 3 -s 7 -t TRUE -n FALSE
	# 绘制OTU_136所有基因型
	alpha_boxplot.sh -i `pwd`/result/otu_table_norm.txt -m '"OTU_136"' \
	-d `pwd`/doc/design.txt -A groupID -B '"b3ACT1KD","b3ACT2CR","b3ACT2KO","b3ACT2ThahDK","b3ACT3GK","b3ACT3KO","b3ALDHKO","b3Col","b3ThadKO","b3ThahKO","b3ThasKO1","b3ThasKO2","b3TL1KO","b3BS"' \
	-o `pwd`/fig/4/box_all  -w 10 -h 6 -s 7 -t TRUE -n FALSE
	



### FigS1.a/b Alpha boxplot shannon and observed_otus

alpha_boxplot.sh -i `pwd`/result/alpha.txt -m '"shannon","chao1","observed_otus"' \
	-d `pwd`/doc/design.txt -A groupID -B '"b3BS","b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO"' \
	-o `pwd`/fig/S1/ -w 5 -h 3 -s 7

### FigS2.a/b/c/d Heatmap of DA OTUs 发现差异数量少时会出现热图黑框，整体发黑；用AI查看差异数量少的结果有灰线框；调小图尺寸可去除；

初步结果位于 http://bailab.genetics.ac.cn/report/16s/AC_b3_all_b3_3_lrt_cultureNew ，保存于 D:\work\ath\jt.HuangAC\batch3all\fig\S2\

修改图片为4x2.5，更新结果位于 http://bailab.genetics.ac.cn/report/16s/AC_b3_all_KO_b3_5 ，保存替换D:\work\ath\jt.HuangAC\batch3all\fig\S2\

# 尝试手绘热图，发现聚类没有原来的好，有被错误聚类的样品
cat result_k1-c/otu_b3ACT2KOvsb3Col_enriched.txt <(tail -n+2 result_k1-c/otu_b3ACT2KOvsb3Col_depleted.txt) > result_k1-c/otu_b3ACT2KOvsb3Col_sig.txt
plot_heatmap.sh -i result_k1-c/otu_b3ACT2KOvsb3Col_sig.txt -o fig/S2/ -w 5 -h 3 -l 12


