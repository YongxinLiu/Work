<!-- TOC -->

- [1. 标准流程 Standard pipeline](#1-标准流程-standard-pipeline)
    - [1.1. 按实验设计拆分lane为文库](#11-按实验设计拆分lane为文库)
    - [1.2. 按实验设计拆分文库为样品](#12-按实验设计拆分文库为样品)
    - [1.3. 样品双端合并、重命名、合并为单一文件](#13-样品双端合并重命名合并为单一文件)
    - [1.4. 切除引物与标签](#14-切除引物与标签)
    - [1.5. 质量控制](#15-质量控制)
    - [1.6. 序列去冗余](#16-序列去冗余)
    - [1.7. 挑选OTU](#17-挑选otu)
    - [1.8. 有参去嵌合体](#18-有参去嵌合体)
    - [1.9. 去除宿主](#19-去除宿主)
    - [1.10. 生成OTU表](#110-生成otu表)
    - [1.11. 过滤样本和OTUs](#111-过滤样本和otus)
    - [1.12. 物种注释](#112-物种注释)
    - [1.13. 物种统计](#113-物种统计)
    - [1.14. 多序列比对和进化树](#114-多序列比对和进化树)
    - [1.15. Alpha多样性指数计算](#115-alpha多样性指数计算)
    - [1.16. Beta多样性距离矩阵计算](#116-beta多样性距离矩阵计算)
    - [1.17. 有参考构建OTU表](#117-有参考构建otu表)
- [2. 统计绘图 Statistics and plot](#2-统计绘图-statistics-and-plot)
    - [2.1. Alpha多样性指数箱线图](#21-alpha多样性指数箱线图)
    - [2.2. Alpha丰富度稀释曲线](#22-alpha丰富度稀释曲线)
    - [2.3. 主坐标轴分析距离矩阵](#23-主坐标轴分析距离矩阵)
    - [2.4. 限制性主坐标轴分析](#24-限制性主坐标轴分析)
    - [2.5. 样品和组各级分类学堆叠柱状图](#25-样品和组各级分类学堆叠柱状图)
    - [2.6. 组间差异比较](#26-组间差异比较)
- [3. 高级分析 Advanced analysis](#3-高级分析-advanced-analysis)
    - [3.1. 添加可培养菌](#31-添加可培养菌)
    - [3.2. 查看Venn 4者共有菌，其中三者共有45，95，1，5](#32-查看venn-4者共有菌其中三者共有459515)
    - [3.3. 整理faprotax中菌的功能列表](#33-整理faprotax中菌的功能列表)
    - [3.4. 绘制网络](#34-绘制网络)
    - [3.5. /5/14 随机森林属水平区分籼粳稻](#35-514-随机森林属水平区分籼粳稻)
- [4. 个性化分析 Custom analysis](#4-个性化分析-custom-analysis)
    - [4.1. 低氮条件下挑选30个IND品种测宏基因组](#41-低氮条件下挑选30个ind品种测宏基因组)
    - [4.2. 挑选NRT样品测宏基因组](#42-挑选nrt样品测宏基因组)
    - [4.3. 新发现OTU_8的物种很像OTU_11](#43-新发现otu_8的物种很像otu_11)
    - [4.4. 筛选各品种最好的3个样品](#44-筛选各品种最好的3个样品)
    - [4.5. OTU或分类与PCoA轴的差异](#45-otu或分类与pcoa轴的差异)
    - [4.6. 检查GWAS是否可以发现Anaeromyxobacter与SNP的关联](#46-检查gwas是否可以发现anaeromyxobacter与snp的关联)
- [5. 图表整理 Figures and legends](#5-图表整理-figures-and-legends)
    - [籼粳稻分型](#籼粳稻分型)
        - [模式图](#模式图)
        - [alpha多样性](#alpha多样性)
        - [beta多样性](#beta多样性)
        - [门及纲水平差异](#门及纲水平差异)
    - [随机森林分类](#随机森林分类)
        - [纲水平建模](#纲水平建模)
        - [展示Top feature](#展示top-feature)
        - [在nrt和时间序列中验证](#在nrt和时间序列中验证)
        - [差异OTUs](#差异otus)
    - [亚种差异与氮相关](#亚种差异与氮相关)
    - [微生物与表型关联](#微生物与表型关联)

<!-- /TOC -->

# 1. 标准流程 Standard pipeline

	处理序列 Processing sequencing data

	# 1. 准备工作 Preparation

	## 1.1. 准备流程配置文件

	# Prepare config file of pipeline
	cd ~/other/guoxiaoxuan/180627_AQ/bac
	
	# 复制标准参数模板和操作指南至项目代码区：方便同步
	mkdir -p ~/github/Work/other/guoxiaoxuan/180627_AQ/bac
	cp ~/github/Work/rice/xianGeng/* ~/github/Work/other/guoxiaoxuan/180627_AQ/bac/
   
	# 链接代码至工作区
	ln -s ~/github/Work/other/guoxiaoxuan/180627_AQ/bac/
	ln -s ~/github/Work/other/guoxiaoxuan/180627_AQ/bac/manual.md manual.sh

	## 1.2. 初始化工作区

	# Initialize the working directory
	make init

	## 1.3. 准备原始数据

	# Prepare raw data
#	# 数据来自三个课题：miniCore + timecourse + nrt
#	
#	# 合并实验设计
#	# 合并三个项目的实验设计，检查样品是否有重名
#	cp ~/rice/miniCore/doc/design.txt doc/design_minicore.txt
#	cp ~/rice/timecourse/doc/design.txt doc/design_timecourse.txt
#	cp ~/rice/zjj.nitrogen/180116/doc/design.txt doc/design_nrt.txt
#	# 统计实验样品行和唯一行，确实样品名唯一
#	cat doc/design_* | grep -v 'SampleID'|wc -l
#	cat doc/design_* | grep -v 'SampleID'|cut -f 1| sort|uniq|wc -l
#	# 合并实验设计，前7列共有，只保留前7列
#	cat <(head -n1 doc/design_nrt.txt) <(cat doc/design_* | grep -v 'SampleID') | cut -f 1-7 > doc/design.txt
#
#	# 添加miniCore亚种
#	mv doc/design.txt doc/design0.txt
#	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$6]}' ~/rice/miniCore/180319/doc/design_group.txt doc/design0.txt|sed 's/\t$/\tsubspecies/'|less -S > doc/design1.txt # 添加亚种
#	awk 'BEGIN{FS=OFS="\t"} {print $0,$7$8}' doc/design1.txt|less -S>doc/design2.txt # 合并土壤类型和亚种
#	cp doc/design2.txt doc/design.txt
#	# 原始数据合并
#	cat ~/rice/miniCore/temp/seqs_usearch.fa ~/rice/timecourse/temp/seqs_usearch.fa ~/rice/zjj.nitrogen/180116/temp/seqs_usearch.fa | cut -f 1 -d ';' | sed 's/_/./g' > temp/filtered.fa
#	# 从2.7 fa_unqiue 开始

	# 原始数据
	cp /mnt/bai/xiaoning/xiaoxuan/180528/bac/clean_data/index??_?.fq seq/
	# 上传library.txt
	# 确定33/64
	determine_phred-score.pl seq/index11_1.fq
	# 如果33批量改名
	awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$3"_1.fq seq/"$1"_1.fq")}' <(grep -v '^$$' doc/library.txt|tail -n+2)
	awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$3"_2.fq seq/"$1"_2.fq")}' <(grep -v '^$$' doc/library.txt|tail -n+2)
	# 如果数据是64，转换为33并改名
	parallel --xapply -j 32 \
	"fastp -i seq/{1}_1.fq -I seq/{1}_2.fq -o seq/{2}_1.fq -O seq/{2}_2.fq -6 -A -G -Q -L -w 9" \
	::: `tail -n+2 doc/library.txt | cut -f 3` ::: `tail -n+2 doc/library.txt | cut -f 1`

## 1.1. 按实验设计拆分lane为文库

	# Split lane into libraries
	# lane文件一般为seq/lane_1/2.fq.gz
	# lane文库信息doc/library.txt：至少包括编号、Index和样品数量三列和标题
	head -n3 doc/library.txt
	#LibraryID	IndexRC	Samples
	#L1	CTCAGA	60
	
	# 按library.txt拆分lane为library # 如果无输入，应该不运行，而不是生成零字节
	make lane_split


## 1.2. 按实验设计拆分文库为样品

	# Prepare design of libraries
	
	# 情况1. 多文库实验设计拆分文库设计
	split_design.pl -i doc/design_raw.txt
 
	# 情况2. 从其它项目复制文库实验设计
	#cp ~/ath/jt.HuangAC/batch3/doc/L?.txt doc/
	#sed -i 's/ //g;s/\r/\n/' doc/*.txt # 删除多余空格


	# 拆分样品
	# 预览文库实验设计
	head -n3 doc/L1.txt
	# 按L1/2/3...txt拆分library为samples
	make library_split
 

## 1.3. 样品双端合并、重命名、合并为单一文件

	# Merge paired reads, renames and merge all samples
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt) <(cat doc/L* |grep -v -P '^SampleID') > doc/design.txt

	# 样品双端合并、重命名、合并为单一文件, 注意fastq为33格式，64位采用fastp转换
	make sample_merge


## 1.4. 切除引物与标签

	# Cut primers and lables
	# 切除左端标签和引物，右端 引物
	# Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
	make fq_trim


## 1.5. 质量控制

	# Quality control
	# 过滤序列中预期累计错误率>1%的序列
	make fq_qc


## 1.6. 序列去冗余
	
	# 从这里开始
	# ln ~/medicago/zjj170823/temp/seqs_usearch.fa temp/filtered.fa

	# Remove redundancy, get unique reads
	make fa_unqiue


## 1.7. 挑选OTU

	# Pick OTUs
	# unoise3速度比cluster_otus慢上百倍
	make otu_pick


## 1.8. 有参去嵌合体

	# Remove chimiras by silva database
	# 基于SILVA数据库去除
	make chimera_ref


## 1.9. 去除宿主

	# Remove host
	# 根据SILVA注释去除线粒体、叶绿体、真核生物18S和未知序列(非rRNA)
	make host_rm

## 1.10. 生成OTU表
	
	# Create OTUs table
	# 默认使用vsearch更快10倍，可选usearch10，线程不可超48
	make otutab_create


## 1.11. 过滤样本和OTUs

	# OTU table filter samples and OTU
	# 推荐过滤低测序量<5000的样本，筛选大于1RPM的OTU
	make otutab_filter 
	# 检查是否有样品缺失
	cat <(head -n1 result/otutab.txt|cut -f 2-|tr '\t' '\n') <(tail -n+2 doc/design.txt | cut -f 1) | sort | uniq -u > doc/missing_samples.txt 


## 1.12. 物种注释

	# Assign taxonomy
	# 默认使用RDP trainset快而准，GG太旧，Silva太慢
	# 推荐阈值为0.6保证注释更完整
	make tax_assign


## 1.13. 物种统计
	
	# Taxonomy summary
	# 必须所有物种有注释，否则有可能报错
	make tax_sum


## 1.14. 多序列比对和进化树
	
	# Multiply alignment and make_phylogeny
	# usearch10/culsterO结果不同可能影响多样性分析(usearch unifrac结果更可信)
	# 进化树，用于树图和多样性分析
	make tree_make

## 1.15. Alpha多样性指数计算
	
	# Calculate alpha diversity index
	# alpha指数计算结果为 result/alpha/index.txt
	# 稀释梯度结果位于 result/alpha/rare.txt
	make alpha_calc

## 1.16. Beta多样性距离矩阵计算
	
	# Beta diversity tree and distance matrix
	# 最好用usearch，结果unifrac分类更好；clustero+fastree结果PCoA较差
	make beta_calc
	# ---Fatal error--- ../calcdistmxu.cpp(32) assert failed: QueryUniqueWordCount > 0 致信作者; 改用qiime1

## 1.17. 有参考构建OTU表

	# Reference based OTU table
	# otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
	make otutab_gg



# 2. 统计绘图 Statistics and plot

## 2.1. Alpha多样性指数箱线图
	
	# Alpha index in boxplot
	make alpha_boxplot

## 2.2. Alpha丰富度稀释曲线
	
	# Alpha rarefracation curve
	make alpha_rare

## 2.3. 主坐标轴分析距离矩阵
	
	# PCoA of distance matrix
	make beta_pcoa

## 2.4. 限制性主坐标轴分析

	# Constrained PCoA / CCA of bray distance matrix
	# OTU表基于bray距离和CCA，至少3个组 
	make beta_cpcoa

## 2.5. 样品和组各级分类学堆叠柱状图

	# Stackplot showing taxonomy in each level
	make tax_stackplot

## 2.6. 组间差异比较 
	
	# Group compareing by edgeR or wilcox
	# 可选负二项分布，或wilcoxon秩和检验
	make DA_compare
	make DA_compare_tax

# 3. 高级分析 Advanced analysis

## 3.1. 添加可培养菌

	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $0,a[$1]}' temp/otutab.mean temp/culture_otu.blastn | cut -f 1-4,14 > result/41culture/otu.txt
	echo -ne "OTUs > 97% abundance :\t" >> result/41culture/summary.txt
	awk '$$3>=97 && $$4>=99' result/41culture/otu.txt | awk '{a=a+$$5} END {print a}' >> result/41culture/summary.txt
	cat result/41culture/summary.txt


## 3.2. 查看Venn 4者共有菌，其中三者共有45，95，1，5
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/41culture/otu.txt result/venn/3_45.txt |sort -k3,3nr|less -S
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/41culture/otu.txt result/venn/3_95.txt |sort -k3,3nr|less -S


## 3.3. 整理faprotax中菌的功能列表
	grep -v '^#' /mnt/bai/yongxin/software/FAPROTAX_1.1/FAPROTAX.txt|sed 's/\*/_/g' > culture/faprotax.tax


## 3.4. 绘制网络
	
	# 制作数据文件
	OTU按差异比较计算整理
	注释信息按network.R整理
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$1]}' result/taxonomy_8.txt <(cut -f 1,4 result/compare/database.txt) | cut -f 1-2,5- | sed '1 s/OTUID\tLTEJ/ID\ttotal\tphylum\tclass\torder\tfamily\tgenus\tspecies/' |less > network/database.txt

	# 绘制籼粳稻属水平模块和差异 co_network_genus_LIndTej.R
	# 以Burkholderia 相关菌 与 Anaeromyxobacter 的小网络 co_network_genus_LIndTej_core.R


## 3.5. /5/14 随机森林属水平区分籼粳稻

	数据集 | 丰度% | 数量 |  错误率%
	LTEJ/LIND | 0.1 | 11.7 | 11
	LTEJ/LIND | 0.5 | 28 | 11
	LTEJ/LIND | 1 | 12 | 13.4

	HTEJ/HIND | 0.1 | 110 | 18.2
	HTEJ/HIND | 0.5 | 28 |18.9
	HTEJ/HIND | 1 | 13 | 18.9

	TEJ/IND | 0.1 | 116 | 14.1
	TEJ/IND | 0.5 | 33 | 15 # 选择全局，features少还兼顾大多数
	TEJ/IND | 1 | 12 | 17.4
 
	# 查看Top3属的分布 Geobacter Tangfeifania Anaeromyxobacter Bacillus Burkholderia

	for i in `tail -n+2 result/randomForest/imp.txt|cut -f 1`; do; \
	alpha_boxplot.sh -i result/tax/sum_g.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND"' -m \"${i}\" -t TRUE -o result/randomForest/box_ -n TRUE	
	done

	tail -n+2 result/randomForest/imp_c.txt|cut -f 1|awk '{print "\""$1"\""}'|tr "\n" ","

	alpha_boxplot.sh -i result/tax/sum_c.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND"' -m '"Flavobacteriia","Bacilli","Gammaproteobacteria","Spirochaetia","Clostridia","Alphaproteobacteria","Actinobacteria","Bacteroidia","Caldilineae","Betaproteobacteria","Ignavibacteria","Holophagae","Acidobacteria_Gp1","Nitrospira","Deltaproteobacteria"' -t TRUE -o result/randomForest/c_ -n TRUE	

	# 发现Top feature贡献度不大。改为科、目、纲时差异较明显。尤其是丰度阈值0.3%时，Top1/2为两类氮相关菌，绘制15个features




# 4. 个性化分析 Custom analysis

## 4.1 2018/7/2 将海南、安徽两块地分开分析

	# 1. 修改海南 2. 统计绘图 部分
	# 筛选海南部分比较组，修改doc/compare.txt
	# 生成新的组列表并修改g1_list，和version为AQ1_Hn
	rm alpha_boxplot
	make rmd

	# 2. 安徽
	# 筛选安徽部分比较组，修改doc/compare.txt
	# 生成新的组列表并修改g1_list，和version为AQ1_Hn
	rm alpha_boxplot
	rm DA_compare
	make plot_manhattan
	make rmd


## 4.2 定量差异OTU比较

# 标准化脚本
/mnt/bai/xiaoning/xiaoxuan/180528/bac/script/wetD_0628.R

mkdir spikein

## 按正常分析比较，筛选共的OTU
# otu=/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/spikein/otutab.txt
# 2018/7/13 # 表有问题，更新如下
otu=/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/table/adjusted_absAbundance.txt 
cat <(head -n1 $otu) <(awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' $otu AQ1/result/compare/BacHnMH63dry-BacHnMH63wet_all.txt) | sed '/^$/d'| sed '1 s/^/OTU\t/' | less -S > spikein/otutab.txt
# 原始数据也制作同样筛选和标准化的RA表
cat <(head -n1 result/otutab_norm.txt) <(awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/otutab_norm.txt AQ1/result/compare/BacHnMH63dry-BacHnMH63wet_all.txt) | sed '/^$/d'| less -S > RA/otutab_norm.txt
## 筛选AA的OTU进行PCoA分析，与报告中的一致。报告中PCoA进行筛选了吗？没有呀

## 筛选AA的OTU进行PCoA分析
	biom convert -i spikein/otutab.txt -o spikein/otutab.biom --table-type="OTU table" --to-json
	# 计算4种距离矩阵 http://qiime.org/scripts/beta_diversity.html -s显示矩阵列表有34种距离可选
	beta_diversity.py -i spikein/otutab.biom -o spikein/beta/ -t result/otu.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
	# 删除文件名中多余字符，以方法.txt为文件名
	rename 's/_otutab//' spikein/beta/*.txt
	# Ah
beta_pcoa.sh -i `pwd`/spikein/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
        -d `pwd`/doc/design.txt -A groupID -B '"BacAhBulksoildry","BacAhBulksoilwet","BacAhMH63dry","BacAhMH63wet","BacAhMH63ZHdry","BacAhMH63ZHwet","BacAhWYJ7DEP1dry","BacAhWYJ7DEP1wet"' -E TRUE \
        -c `pwd`/doc/compare.txt \
        -o `pwd`/spikein/beta/Ah -h 5 -w 8
beta_pcoa.sh -i `pwd`/spikein/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
        -d `pwd`/doc/design.txt -A groupID -B '"BacHnBulksoildry","BacHnBulksoilwet","BacHnMH63dry","BacHnMH63wet","BacHnMH63ZHdry","BacHnMH63ZHwet","BacHnWYJ7DEP1dry","BacHnWYJ7DEP1wet"' -E TRUE \
        -c `pwd`/doc/compare.txt \
        -o `pwd`/spikein/beta/Hn -h 5 -w 8


# 比较RA
compare.sh -i result/otutab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
	-p 0.05 -q 0.05 -F 1.2 -t 0.1 \
	-d `pwd`/doc/design.txt -A groupID -B '"BacAhBulksoildry","BacAhBulksoilwet","BacAhMH63dry","BacAhMH63wet","BacAhMH63ZHdry","BacAhMH63ZHwet","BacAhWYJ7DEP1dry","BacAhWYJ7DEP1wet","BacHnBulksoildry","BacHnBulksoilwet","BacHnMH63dry","BacHnMH63wet","BacHnMH63ZHdry","BacHnMH63ZHwet","BacHnWYJ7DEP1dry","BacHnWYJ7DEP1wet"' \
	-o RA/compare/


# 比较AA(spikein)
compare.sh -i spikein/otutab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
	-p 0.05 -q 0.05 -F 1.2 -t 0 -N FALSE \
	-d `pwd`/doc/design.txt -A groupID -B '"BacAhBulksoildry","BacAhBulksoilwet","BacAhMH63dry","BacAhMH63wet","BacAhMH63ZHdry","BacAhMH63ZHwet","BacAhWYJ7DEP1dry","BacAhWYJ7DEP1wet","BacHnBulksoildry","BacHnBulksoilwet","BacHnMH63dry","BacHnMH63wet","BacHnMH63ZHdry","BacHnMH63ZHwet","BacHnWYJ7DEP1dry","BacHnWYJ7DEP1wet"' \
	-o spikein/compare/

## 4.2.2 定量下属比较

# 按差异OTU使用的560个OTUs进行合并，再差异比较
taxo_summary_RA.R
#
i="g"
for i in p c o f g; do
mkdir -p RA/compare_${i}; \
compare.sh -i RA/${i}.txt -c `pwd`/doc/compare.txt -m "wilcox" \
	-p 0.05 -q 0.05 -F 1.2 -t 0 -N FALSE -U 10000  \
	-d `pwd`/doc/design.txt -A groupID -B '"BacAhBulksoildry","BacAhBulksoilwet","BacAhMH63dry","BacAhMH63wet","BacAhMH63ZHdry","BacAhMH63ZHwet","BacAhWYJ7DEP1dry","BacAhWYJ7DEP1wet","BacHnBulksoildry","BacHnBulksoilwet","BacHnMH63dry","BacHnMH63wet","BacHnMH63ZHdry","BacHnMH63ZHwet","BacHnWYJ7DEP1dry","BacHnWYJ7DEP1wet"' \
	-o RA/compare_${i}/
done


Rscript script/tax_summary_spikein.R
# i="g"
for i in p c o f g; do
mkdir -p spikein/compare_${i}; \
compare.sh -i spikein/${i}.txt -c `pwd`/doc/compare.txt -m "wilcox" \
	-p 0.05 -q 0.05 -F 1.2 -t 0 -N FALSE \
	-d `pwd`/doc/design.txt -A groupID -B '"BacAhBulksoildry","BacAhBulksoilwet","BacAhMH63dry","BacAhMH63wet","BacAhMH63ZHdry","BacAhMH63ZHwet","BacAhWYJ7DEP1dry","BacAhWYJ7DEP1wet","BacHnBulksoildry","BacHnBulksoilwet","BacHnMH63dry","BacHnMH63wet","BacHnMH63ZHdry","BacHnMH63ZHwet","BacHnWYJ7DEP1dry","BacHnWYJ7DEP1wet"' \
	-o spikein/compare_${i}/
done



### Reclaculate RA

/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/doc/qPCR_new.txt


/mnt/bai/xiaoning/xiaoxuan/180528/180627_AQ/bac/spikein/AA.txt


师兄  我手动check了一下表  没有找到错误   需要师兄帮助 单独计算一下 用AA.txt的 OTU数 和 qPCR_new.txt 做相除  出一个 校正的表  如果咱两的一样  说明最开始给师兄的表可能是对的











