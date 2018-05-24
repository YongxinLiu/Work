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
	cd ~/github/Work/rice/xianGeng
	
	# 复制标准参数模板和操作指南至项目代码区：方便同步
	cp ~/github/Amplicon/16Sv2/parameter.md ./
	cp ~/github/Amplicon/16Sv2/manual.md ./
   
	# 链接代码至工作区
	ln -fs `pwd`/parameter.md ~/rice/xianGeng/makefile
	ln -fs `pwd`/manual.md ~/rice/xianGeng/manual.sh

	## 1.2. 初始化工作区

	# Initialize the working directory
	cd ~/rice/xianGeng
	make init

	## 1.3. 准备原始数据

	# Prepare raw data
	# 数据来自三个课题：miniCore + timecourse + nrt
	
	# 合并实验设计
	# 合并三个项目的实验设计，检查样品是否有重名
	cp ~/rice/miniCore/doc/design.txt doc/design_minicore.txt
	cp ~/rice/timecourse/doc/design.txt doc/design_timecourse.txt
	cp ~/rice/zjj.nitrogen/180116/doc/design.txt doc/design_nrt.txt
	# 统计实验样品行和唯一行，确实样品名唯一
	cat doc/design_* | grep -v 'SampleID'|wc -l
	cat doc/design_* | grep -v 'SampleID'|cut -f 1| sort|uniq|wc -l
	# 合并实验设计，前7列共有，只保留前7列
	cat <(head -n1 doc/design_nrt.txt) <(cat doc/design_* | grep -v 'SampleID') | cut -f 1-7 > doc/design.txt

	# 添加miniCore亚种
	mv doc/design.txt doc/design0.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$6]}' ~/rice/miniCore/180319/doc/design_group.txt doc/design0.txt|sed 's/\t$/\tsubspecies/'|less -S > doc/design1.txt # 添加亚种
	awk 'BEGIN{FS=OFS="\t"} {print $0,$7$8}' doc/design1.txt|less -S>doc/design2.txt # 合并土壤类型和亚种
	cp doc/design2.txt doc/design.txt


	# 原始数据合并
	cat ~/rice/miniCore/temp/seqs_usearch.fa ~/rice/timecourse/temp/seqs_usearch.fa ~/rice/zjj.nitrogen/180116/temp/seqs_usearch.fa | cut -f 1 -d ';' | sed 's/_/./g' > temp/filtered.fa
	# 从2.7 fa_unqiue 开始
	
	
## 1.1. 按实验设计拆分lane为文库

	# Split lane into libraries
	# lane文件一般为seq/lane_1/2.fq.gz
	# lane文库信息doc/library.txt：至少包括编号、Index和样品数量三列和标题
	head -n3 doc/library.txt
	#LibraryID	IndexRC	Samples
	#L1	CTCAGA	60
	
	# 按library.txt拆分lane为library
	make lane_split


## 1.2. 按实验设计拆分文库为样品

	# Prepare design of libraries
	
	# 情况1. 多文库实验设计拆分文库设计
	split_design.pl -i doc/design_raw.txt
 
	# 情况2. 从其它项目复制文库实验设计
	cp ~/ath/jt.HuangAC/batch3/doc/L?.txt doc/
	sed -i 's/ //g;s/\r/\n/' doc/*.txt # 删除多余空格

	# 拆分样品
	# 预览文库实验设计
	head -n3 doc/L1.txt
	# 按L1/2/3...txt拆分library为samples
	make library_split
 

## 1.3. 样品双端合并、重命名、合并为单一文件

	# Merge paired reads, renames and merge all samples
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt) <(cat doc/L* |grep -v '#') > doc/design.txt

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
	ln ~/medicago/zjj170823/temp/seqs_usearch.fa temp/filtered.fa

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
	# 此处结果统计只有3441个样品，而实验设计有3898个样品，少了哪些样品种？
	cat <(head -n1 result/otutab.txt|cut -f 2-|tr '\t' '\n') <(tail -n+2 doc/design.txt | cut -f 1) | sort | uniq -u > doc/missing_samples.txt 
	# 缺失 482样，主要是没有时间序列


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

## 4.1. 低氮条件下挑选30个IND品种测宏基因组

	# 目标：挑选基于LN条件下Weighted Unifrac结果中IND与TEJ差异明显的的30个代表品种
	# 方法：筛选LN下IND/TEJ样品，按组合并，添加标签后筛选。
	mkdir -p pick_variety
	cd pick_variety
	# 筛选LN下IND/TEJ样品
	Rscript /mnt/bai/yongxin/rice/miniCore/180319/scripts/otutab_sample_subset.r -i ../result/otutab.txt -d /mnt/bai/yongxin/rice/miniCore/180319/LN/design.txt -o LN_otutab0.txt
	# 按品种合并
	head -n1 LN_otutab0.txt|cut -f 2-|tr '\t' '\n'|awk '{print $1"\t"$1}'|cut -c1-13|less > LN_sample_variety.list
	usearch10 -otutab_group LN_otutab0.txt -labels LN_sample_variety.list -output LN_otutab1.txt
	# 计算PCoA
	usearch10 -cluster_agg ../result/otu.fa -treeout otu.tree
	usearch10 -beta_div LN_otutab1.txt -tree otu.tree -filename_prefix LN_ -metrics bray_curtis,unifrac
	# 提取实验设计
	grep -P 'soiltypesubspecies|LIND|LTEJ' ../doc/design.txt | cut -f 6- |uniq | less -S > design.txt
	beta_pcoa.sh -i LN_ -m '"bray_curtis","unifrac"' -d design.txt -A subspecies -B '"IND","TEJ"' -o LN_pcoa_ -w 8 -h 5

	Rscript /mnt/bai/yongxin/rice/miniCore/180319/scripts/beta_pcoa_group.r -i LN/bray_curtis.txt -d doc/design.txt -n GroupID -o LN/pcoa_bray

	# 绘制土壤PCoA
beta_pcoa.sh -i result/beta/ -m '"bray_curtis","weighted_unifrac"' \
-B '"HSoil1","HSoil2","LSoil1","LSoil2"' -E TRUE \
-o pick_variety/soil_ -h 5 -w 8



## 4.2. 挑选NRT样品测宏基因组

	# 按Anaeromyxobacter为新版OTU_11(1)，备选OTU_8(0.37)，必须 使用RDP注释，其它注释科faprotax无法识别
	# 查看OTU_11在CP 2016年中箱线图
	mkdir -p pick_nrt
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND","A50LnCp6","A56LnCp6","A50LnCp7","A56LnCp7","A50LnSz7","A56LnSz7","A50HnCp6","A56HnCp6","A50HnCp7","A56HnCp7","A50HnSz7","A56HnSz7","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","nrtHnCp7","ZH11HnCp7","ZH11LnCp7","nrtLnCp7","nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7","HSoil1","HSoil2","LSoil1","LSoil2","soilHnCp7","soilLnCp7","soilHnCp6","soilLnCp6","soilLnSz7","soilHnSz7"' -m '"OTU_8"' -t TRUE -o temp/all_ -n TRUE
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"A50HnCp6","A56HnCp6","A58HnCp6","IR24HnCp6","V3703HnCp6","ZH11HnCp6","A50LnCp6","A56LnCp6","A58LnCp6","IR24LnCp6","V3703LnCp6","ZH11LnCp6"' -m '"OTU_8"' -t TRUE -o pick_nrt/CP6_ -n true
	# OTU_8/11在CP7中变化
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND","A50LnCp7","A56LnCp7","A58LnCp7","IR24LnCp7","A50HnCp7","A56HnCp7","A58HnCp7","IR24HnCp7","nrtHnCp7","ZH11HnCp7","ZH11LnCp7","nrtLnCp7","nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7"' -m '"OTU_11"' -t TRUE -o pick_nrt/CP7_ -n TRUE

	# OTU_8在目标群体中准确度和丰度远小于OTU_11

	# OTU_11在时间序列中变化
	mkdir -p timecourse
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"A50Cp0","A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119","A50Sz0","A50Sz1","A50Sz2","A50Sz3","A50Sz5","A50Sz7","A50Sz10","A50Sz13","A50Sz27","A50Sz34","A50Sz41","A50Sz48","A50Sz56","A50Sz62","A50Sz69","A50Sz76","A50Sz83","A50Sz90","A50Sz97","A50Sz118"' -m '"OTU_11"' -t TRUE -o timecourse/A50_ -n TRUE
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"IR24Cp0","IR24Cp1","IR24Cp2","IR24Cp3","IR24Cp7","IR24Cp10","IR24Cp14","IR24Cp21","IR24Cp28","IR24Cp35","IR24Cp42","IR24Cp49","IR24Cp63","IR24Cp70","IR24Cp77","IR24Cp84","IR24Cp91","IR24Cp98","IR24Cp112","IR24Cp119","IR24Sz0","IR24Sz1","IR24Sz2","IR24Sz3","IR24Sz5","IR24Sz7","IR24Sz10","IR24Sz13","IR24Sz27","IR24Sz34","IR24Sz41","IR24Sz48","IR24Sz56","IR24Sz62","IR24Sz69","IR24Sz76","IR24Sz83","IR24Sz90","IR24Sz97","IR24Sz118"' -m '"OTU_11"' -t TRUE -o timecourse/IR24_ -n TRUE

	# 绘制土壤PCoA
	beta_pcoa.sh -i result/beta/ -m '"bray_curtis","weighted_unifrac"' \
	-B '"soilHnCp7","soilLnCp7","soilHnCp6","soilLnCp6","soilLnSz7","soilHnSz7"' -E TRUE \
	-o pick_nrt/soil_ -h 5 -w 8
	
	# OTU_11 在土壤中组间丰度
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"HSoil1","HSoil2","LSoil1","LSoil2","soilHnCp7","soilLnCp7","soilHnCp6","soilLnCp6","soilLnSz7","soilHnSz7"' -m '"OTU_11"' -t TRUE -o temp/soil_ -n TRUE


## 4.3. 新发现OTU_8的物种很像OTU_11
	less result/faprotax/report# 但按物种注释0.6注释到科无法识别；改为0.3阈值


## 4.4. 筛选各品种最好的3个样品

	# 结果备份，并筛选3个样品看结果
	cp -r result/ result180508/
	# 筛选实验设计为3个样品
	cp -r doc/design.txt doc/design.txt180508
	# 获得所有筛选样品
	cat ~/rice/miniCore/180319/HN/pcoa_bray_samples_top3.id <(cut -f 1 ~/rice/miniCore/180319/LN/pcoa_bray_samples_top3.group) > temp/temp1
	# 获得删除样品
	cat <(tail -n+2 doc/design_minicore.txt|cut -f 1) temp/temp1 | sort | uniq -u > doc/minicore_discard.id
	# 剔除点注释design
	for i in `cat doc/minicore_discard.id`; do;\
		sed -i "s/$i/#$i/" doc/design.txt;done


## 4.5. OTU或分类与PCoA轴的差异
	
	# 以LN的weighted unifrac PCoA1/2为例
	# 基于PCoA轴与OTU计算相关
	script/cor_pcao_otu.R 计算OTU与4unifrac前4轴spearman相关系数，保存为 result/cor_otu_pcoa_unifrac.txt
	# 添加至差异OTUs，和相关系数
	awk 'NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1]"\t"$0}' result/cor_otu_pcoa_unifrac.txt result/compare/LTEJ-LIND_sig.txt | sed '1 s/^/Cor_PC1\tCor_PC2/' > result/compare/LTEJ-LIND_sig_pcoa_unifrac.txt
	# 选择代表菌作为组的markers，我推荐按丰度选差异；可进一步结果与PCoA轴的相关性；或选3/4组共有

	# 2018/5/16 国家、经纬度、地区和亚种与PCoA间关系script/beta_pcoa_location.R 发现Latitudeg与PC1相关，0.4366
	
	# 绘制某个分类单元的箱线图，常用绘制按丰度取上下调的top3；也可按Pvalue选择，但高丰度并不占优
	wd=`pwd`
	Dct_tax=g
	mkdir -p result/otu_boxplot_${Dct_tax}
	Dct_input=${wd}/result/tax/sum_${Dct_tax}.txt
	plot_boxplot_ggpubr.sh -i ${wd}/result/tax/sum_${Dct_tax}.txt -d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND"' \
	-m '"Anaeromyxobacter","Curvibacter","Sideroxydans","Burkholderia","Bradyrhizobium","Rhizobium"' -t TRUE -o result/otu_boxplot_${Dct_tax}/ -n true
	


## 4.6. 检查GWAS是否可以发现Anaeromyxobacter与SNP的关联
	# 基于LN TEJ/IND中Anaeromyxobacter最低和最高的20个品种(保证数据可分开)，与nsSNP关联，最为仅为e-3，10m21759092,10,21759092,0.643641032317718,0.475,40,0.893753744914084,0.894396427377941,1。
	gapit_IndTej_Anaeromyxobacter_top20.R
	# 修改为30个品种有重合，肯定无法找到差异，再改为25个品种
	10m21759092,10,21759092,0.962618742911769,0.46,50,0.856299549567109,0.85630648832468,1
	# 可能从样本量、数据波动程度、SNP背景均无法满足


# 5. 图表整理 Figures and legends

##  图1. 籼粳稻分型

### 模式图
	
	地图，种植模式图——秦媛

### alpha多样性
	
	箱线图、稀释取线、样品稀释取线
#	# HN下TEJ和IND
#	alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
#	-d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND"' \
#	-o `pwd`/fig/1/HN -h 3 -w 5
#	# LN下TEJ和IND
#	alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
#	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND"' \
#	-o `pwd`/fig/1/LN -h 3 -w 5
	alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' \
	-o `pwd`/fig/1/alpha_ -h 3 -w 5
	# 稀释曲线
	alpha_rare.sh -i `pwd`/result/alpha/rare.txt \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' \
	-o `pwd`/fig/1/alpha_ -h 3 -w 5

	# 样品稀释曲线


### beta多样性
	
	# 不同加土，否则主要差异为土壤；不能两块地混合，否则主要差异为不同地块
	# HN下TEJ和IND
beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac"' \
	-d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND"' -E TRUE \
	-c `pwd`/doc/compare.txt \
	-o `pwd`/fig/1/beta_HN_ -h 3 -w 5
	# LN下TEJ和IND
beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac"' \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND"' -E TRUE \
	-c `pwd`/doc/compare.txt \
	-o `pwd`/fig/1/beta_LN_ -h 3 -w 5
	# Constrained LN + Soil
beta_cpcoa.sh -i `pwd`/result/otutab.txt -m '"bray","jaccard"' \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1"' -E TRUE \
	-o `pwd`/fig/1/beta_LN_ -h 3 -w 5
	# Constrained LN + HN
beta_cpcoa.sh -i `pwd`/result/otutab.txt -m '"bray","jaccard"' \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","HTEJ","HIND"' -E TRUE \
	-o `pwd`/fig/1/beta_ -h 3 -w 5

### Taxonomy 门+变形菌纲
	cut -f 3-4 result/taxonomy_8.txt|sort|uniq|grep 'Proteobacteria' # 为什么会有这么多结果，只选5类继续分析
	cat <(grep -v 'Proteobacteria' result/tax/sum_p.txt) <(grep 'proteobacteria' result/tax/sum_c.txt) > result/tax/sum_pc.txt

tax_stackplot.sh -i `pwd`/result/tax/sum_ -m '"pc"' -n 10 \
	-d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"' -O FALSE \
	-o `pwd`/fig/1/tax_pc_ -h 3 -w 5


### 门及纲水平差异


## 图2. 随机森林分类
	
	randomForest_class.R
	1. 纲水平建模，展示贡献度，和样品中热图
	使用高HN和HN下籼粳稻纲水平0.3%丰度的15个Feature机器学习；保存预测结果confusion.txt，整理16.4%错误率，TEJ 37.4%错误；
	2018/5/21 删除三个澳大亚利(纬度为负)粳稻, D4032, D4038, F4053; 标记A50/ZH11为TEJ，而IR24为IND，各分为Hn/Ln两种情况；
	错误率降低为15.3%，TEJ为36%;Top1 feature也变为了Nitrospira，Deltaproteobacteria

	2. Top feature：用各组柱状图/箱线图分类展示，再加梯度排序
	
	3. 在nrt和时间序列中验证
	# 都是IND比较准，TEJ不应

	4. 差异OTUs在时间序列中变化
	plot_boxplot_ggpubr.sh -i result/tax/sum_c.txt -d `pwd`/doc/design.txt -A groupID -B '"A50Cp0","A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119"' \
	-m '"Deltaproteobacteria","Actinobacteria","Alphaproteobacteria","Clostridia","Betaproteobacteria","Nitrospira"' -t TRUE -o result/randomForest/time_ -n TRUE # 42以后没有数据呢？改用alpha_boxplot.sh

	alpha_boxplot.sh -i result/tax/sum_c.txt -d `pwd`/doc/design.txt -A groupID -B '"A50Cp0","A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119"' \
	-m '"Deltaproteobacteria","Actinobacteria","Alphaproteobacteria","Clostridia","Betaproteobacteria","Nitrospira"' -t TRUE -o result/randomForest/time_ -n TRUE 

	# plot_heatmap_timecourse.R 绘制时间序列的图
	# 筛选HL/LN下TEJ-IND共同下调的菌
	cat result/compare/?TEJ-?IND_sig.txt | grep 'Depleted' | cut -f 1 | sort | uniq -d > fig/2/otu_IND_common_specific.txt


	
## 图3. 亚种差异与氮相关
	
	1. 差异菌存在氮代谢通路差异
	filter_otus_by_sample.sh -f result/faprotax/element_tab.txt -o result/faprotax/xiaogeng -d doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"'
	# 结果可用STAMP进一步探索
	
	# 筛选时间序列中上/下调两大类的OTU进行功能分析
	mkdir -p fig/3
	awk '$2>0' fig/2/otu_IND_common_specific_time_cor6.txt | cut -f 1 > fig/3/timecournse_increase.id
	awk '$2<0' fig/2/otu_IND_common_specific_time_cor6.txt | cut -f 1 > fig/3/timecournse_decrease.id
	# 以Incease为例, decrease
	type=decrease
	filter_otus_from_otu_table.py -i result/otutab_norm_tax.biom -o timecourse/${type}.biom --otu_ids_to_exclude_fp fig/3/timecournse_${type}.id --negate_ids_to_exclude
	/usr/bin/python2.7 /mnt/bai/yongxin/software/FAPROTAX_1.1/collapse_table.py -i timecourse/${type}.biom -o timecourse/${type}.faprotax -g /mnt/bai/yongxin/software/FAPROTAX_1.1/FAPROTAX.txt --collapse_by_metadata 'taxonomy' -v --force # --out_report result/faprotax/report 
	filter_otus_by_sample.sh -f timecourse/${type}.faprotax -o fig/3/timecournse_faprotax_${type} -d doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"'
	# increase的差异功能类型均为IND>TEJ>soil，且以芳香、氮 相关
	compare.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
	-p 0.01 -q 0.01 -F 1.2 -t 0.0005 \
	-d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"' \
	-o `pwd`/result/compare/
	# decrease的差异，stamp打开报错，但有时成功；下调无N循环相关，有
	
	# 可视化菌的功能有无
	## 筛选report为功能有无表
	grep 'OTU_' -B 1 result/faprotax/report | grep -v -P '^--$' > result/faprotax/report.clean
	faprotax_report_sum.pl -i result/faprotax/report.clean -o result/faprotax/report
	#OTU功能注释列表：result/faprotax/report.otu_func
	#功能包含OTU列表：result/faprotax/report.func_otu
	#OTU功能有无矩阵：result/faprotax/report.mat


	2. 差异菌与氮相关基因显著相关

	3. 关键氮高效基因不同形态、突变体可部分解析亚种差异


## 微生物与表型关联

	1. 微生物、多样性、PCoA主轴与表型相关分析