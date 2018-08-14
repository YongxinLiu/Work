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

	# 创建环境代码见~/github/Work/initial_project.sh
	# 本项目在xiangeng基础上继续
	cp ~/github/Work/rice/xianGeng/*md ~/github/Work/rice/integrate16s/

	## 1.2. 初始化工作区

	# Initialize the working directory
	cd ~/rice/miniCore/180718
	make init

	## 1.3. 准备原始数据

	# Prepare raw data
	# 数据来自三个课题：miniCore + timecourse + nrt1.1b + nrt1.1a + SL + epi + Gprotein + CTK + SD1

	
	# 合并实验设计
	# 合并三个项目的实验设计，检查样品是否有重名
	cp ~/rice/miniCore/doc/design.txt doc/design_minicore.txt
	cp ~/rice/timecourse/doc/design.txt doc/design_timecourse.txt
	cp ~/rice/zjj.nitrogen/180116/doc/design.txt doc/design_nrt.txt
	cp ~/rice/nrt1.1a/doc/design.txt doc/design_nrt1a.txt
	cp ~/rice/strigolactone.LiJY/doc/design.txt doc/design_SL.txt
	cp ~/rice/rice.epi/180409/doc/design.txt doc/design_epi.txt
	cp ~/rice/Gprotein/doc/design.txt doc/design_Gprotein.txt
	cp ~/rice/gxx_CTK/doc/design.txt doc/design_CTK.txt
	cp ~/rice/zn.sd1/doc/design.txt doc/design_sd1.txt

	# 统计实验样品行和唯一行，确实样品名唯一
	cat doc/design_* | grep -v 'SampleID'|wc -l
	cat doc/design_* | grep -v 'SampleID'|cut -f 1| sort|uniq|wc -l 
	# 9 projecct include 7253 samples

	# 合并实验设计，前7列共有，只保留前7列
	cat <(head -n1 doc/design_nrt.txt) <(cat doc/design_* | grep -v 'SampleID') | cut -f 1-7 > doc/design.txt

	# 添加miniCore亚种
	mv doc/design.txt doc/design0.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$6]}' ~/rice/miniCore/180319/doc/design_group.txt doc/design0.txt|sed 's/\t$/\tsubspecies/'|less -S > doc/design1.txt # 添加亚种
	awk 'BEGIN{FS=OFS="\t"} {print $0,$7$8}' doc/design1.txt|less -S>doc/design2.txt # 合并土壤类型和亚种
	cp doc/design2.txt doc/design.txt

	# 原始数据合并
	cat ~/rice/miniCore/temp/seqs_usearch.fa \
	~/rice/timecourse/temp/seqs_usearch.fa \
	~/rice/zjj.nitrogen/180116/temp/seqs_usearch.fa \
	~/rice/nrt1.1a/temp/filtered.fa \
	~/rice/strigolactone.LiJY/temp/seqs_usearch.fa \
	~/rice/rice.epi/180409/temp/seqs_usearch.fa \
	~/rice/Gprotein/temp/seqs_usearch.fa \
	~/rice/gxx_CTK/temp/seqs_usearch.fa \
	~/rice/zn.sd1/temp/seqs_usearch.fa | cut -f 1 -d ';' | sed 's/_/./g' > temp/filtered.fa
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
	# Remove redundancy, get unique reads
	ll temp/filtered.fa # 175 GB
	make fa_unqiue
	# 4亿条序列，100为164640条序列，97% OTU 12625；400为50167，97% OTU 6166，unoise 26756 太多


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


# 设置实验设计 ../xianGeng/doc/design.txt，计算soildtype, H、L下IND/TEJ/TRJ、AUS/
cp ../xianGeng/doc/design.txt doc/design_xiangeng.txt 
# 拥用solitype x subspecies列 soiltypesubspecies，更准确的手工修正为gorupID


# 2. 统计绘图 Statistics and plot

## 比较四大亚种，实验设计位于doc/xiangeng 目录中
	sub=xiangeng
	mkdir -p doc/${sub}
	cp doc/design_${sub}.txt doc/${sub}/design.txt
	cp doc/compare.txt doc/${sub}/compare.txt
	# 实验组改为"HIND","HTEJ","HAUS","HTRJ","LIND","LTEJ","LAUS","LTRJ"

## 比较SL，实验设计位于doc/SL目录中
	sub=SL
	mkdir -p doc/${sub}
	cp doc/design_${sub}.txt doc/${sub}/design.txt
	sed -i 's/genotypeID/groupID/g' doc/${sub}/design.txt
	cat doc/${sub}/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
	# 修改makefile中的sub和g1_list



## 2.1. Alpha多样性指数箱线图
	
	# Alpha index in boxplot
	make alpha_boxplot

# 绘制四大亚种，在H/L下的alpha多样性
alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
        -d `pwd`/doc/design_xiangeng.txt  -A groupID -B '"HIND","HTEJ","HAUS","HTRJ","HSoil1","LIND","LTEJ","LAUS","LTRJ","LSoil1"' \
        -o `pwd`/result/alpha/ -h 3 -w 5

## 2.2. Alpha丰富度稀释曲线
	
	# Alpha rarefracation curve
	make alpha_rare

## 2.3. 主坐标轴分析距离矩阵
	
	# PCoA of distance matrix
	make beta_pcoa

# 绘制四大亚种，在H/L下的beta多样性
beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
        -d `pwd`/doc/design_xiangeng.txt  -A groupID -B '"HIND","HTEJ","HTRJ","HAUS"' -E TRUE \
        -c `pwd`/doc/compare.txt \
        -o `pwd`/result/beta/H -h 3 -w 5
beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
        -d `pwd`/doc/design_xiangeng.txt  -A groupID -B '"LIND","LTEJ","LTRJ","LAUS"' -E TRUE \
        -c `pwd`/doc/compare.txt \
        -o `pwd`/result/beta/L -h 3 -w 5

## 2.4. 限制性主坐标轴分析

	# Constrained PCoA / CCA of bray distance matrix
	# OTU表基于bray距离和CCA，至少3个组 
	make beta_cpcoa

# 绘制四大亚种，在H/L下的限制性PCoA
beta_cpcoa.sh -i `pwd`/result/otutab.txt -m '"bray","jaccard"' \
        -d `pwd`/doc/design_xiangeng.txt  -A groupID -B '"HIND","HTEJ","HTRJ","HAUS"' -E TRUE \
        -o `pwd`/result/beta/H -h 3 -w 5
beta_cpcoa.sh -i `pwd`/result/otutab.txt -m '"bray","jaccard"' \
        -d `pwd`/doc/design_xiangeng.txt  -A groupID -B '"LIND","LTEJ","LTRJ","LAUS"' -E TRUE \
        -o `pwd`/result/beta/L -h 3 -w 5

## 2.5. 样品和组各级分类学堆叠柱状图

	# Stackplot showing taxonomy in each level
	make tax_stackplot

# 输出命令行
make -n -B tax_stackplot
# 修改为需要内容，绘制H/L下亚种和土壤的各层级相对丰度,sp代表亚种subspecies
tax_stackplot.sh -i `pwd`/result/tax/sum_ -m '"p","pc","c","o","f","g"' -n 10 \
        -d `pwd`/doc/design_xiangeng.txt  -A groupID -B '"HIND","HTEJ","HAUS","HTRJ","HSoil1","LIND","LTEJ","LAUS","LTRJ","LSoil1"' -O FALSE \
        -o `pwd`/result/tax/sp_ -h 3 -w 5

## 2.6. 组间差异比较 
	
	# Group compareing by edgeR or wilcox
	# 可选负二项分布，或wilcoxon秩和检验
	make DA_compare

make -B -n DA_compare
# 比对四大亚种间区别，先以LN下亚种间为例
compare.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare_sp.txt -m "wilcox" \
        -p 0.01 -q 0.05 -F 1.2 -t 0.1 \
        -d `pwd`/doc/design_xiangeng.txt  -A groupID -B '"HIND","HTEJ","HAUS","HTRJ","LIND","LTEJ","LAUS","LTRJ"' \
        -o `pwd`/result/compare/ -C "groupID2"



# 3. 高级分析 Advanced analysis

## 3.1. 添加可培养菌

	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $0,a[$1]}' temp/otutab.mean temp/culture_otu.blastn | cut -f 1-4,14 > result/41culture/otu.txt
	echo -ne "OTUs > 97% abundance :\t" >> result/41culture/summary.txt
	awk '$$3>=97 && $$4>=99' result/41culture/otu.txt | awk '{a=a+$$5} END {print a}' >> result/41culture/summary.txt
	cat result/41culture/summary.txt


## 3.2. 查看Venn 4者共有菌，其中三者共有45，95，1，5
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/41culture/otu.txt result/venn/3_45.txt |sort -k3,3nr|less -S
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/41culture/otu.txt result/venn/3_95.txt |sort -k3,3nr|less -S

	# 制作有平均丰度，和物种注释的表
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $1,$5,a[$1]}' result/taxonomy_2.txt result/41culture/otu.txt > result/41culture/otu_mean_tax.txt

## 3.3. 整理faprotax中菌的功能列表
	grep -v '^#' /mnt/bai/yongxin/software/FAPROTAX_1.1/FAPROTAX.txt|sed 's/\*/_/g' > culture/faprotax.tax


## 3.4. 绘制网络
	
	# 制作数据文件
	OTU按差异比较计算整理
	注释信息按network.R整理
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$1]}' result/taxonomy_8.txt <(cut -f 1,4 result/compare/database.txt) | cut -f 1-2,5- | sed '1 s/OTUID\tLTEJ/ID\ttotal\tphylum\tclass\torder\tfamily\tgenus\tspecies/' |less > network/database.txt

	# 绘制籼粳稻属水平模块和差异 co_network_genus_LIndTej.R
	# 以Burkholderia 相关菌 与 Anaeromyxobacter 的小网络 co_network_genus_LIndTej_core.R




# 4. 个性分析

## 4.1. 分蘖与菌相关性

	wd=/mnt/bai/yongxin/rice/integrate16s
	# 准备相关输入文件
	cd $wd
	# 硬链数据文件，保持可同步修改和可备份
	# miniCore分蘖数据整理
	ln ~/rice/xianGeng/doc/phenotype_sample_raw.txt doc/

	# 以miniCore 180319批次数据进行相关分析
	# LN otu表和实验设计
	mkdir -p data
	cp ~/rice/miniCore/180319/LN/otutab.txt data/LN_otutab.txt
	cp ~/rice/miniCore/180319/doc/design.txt doc/design_miniCore.txt
	mkdir -p data/cor/LN
	# 物种注释
	cp ~/rice/miniCore/180319/temp/otus_no_host.tax data/
	# 统计见script/cor_tiller_LN.Rmd
	# 相关系数，添加物种注释
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4} NR>FNR{print $0,a[$1]}' data/otus_no_host.tax data/cor/LN/otu_mean_pheno_cor.r.txt | less -S > data/cor/LN/otu_mean_pheno_cor.r.txt.tax
	# 再添加可培养相关菌

	# 以integrate16s批次数据进行相关分析
	# 统计见script/cor_tiller_LN.Rmd，输出top300 OTUs的相关；再输出all OTUs的相关结果
	# 相关系数，添加物种注释
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' result/taxonomy_2.txt result/cor/LN/otu_mean_pheno_cor.r.txt | less -S | sed '1 s/^/OTUID\t/' | sed '1 s/$/Taxonomy/' > result/cor/LN/otu_mean_pheno_cor.r.txt.tax
	# 再添加可培养相关菌
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result/39culture/otu.txt result/cor/LN/otu_mean_pheno_cor.r.txt.tax | less -S > result/cor/LN/otu_mean_pheno_cor.r.txt.tax.xls


## 4.2 样品选测宏基因组

	mkdir -p result/meta_select
	sed 's/#//g' doc/design_xiangeng.txt > doc/design_xiangeng_raw.txt
	# 添加亚种注释
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$8} NR>FNR{print $0,a[$1]}' doc/design_xiangeng_raw.txt ~/rice/miniCore/180319/LN/pcoa_bray_samples_all.txt | less -S > result/meta_select/LN_sample.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$8} NR>FNR{print $0,a[$1]}' doc/design_xiangeng_raw.txt ~/rice/miniCore/180319/HN/pcoa_bray_samples_all.txt | less -S > result/meta_select/HN_sample.txt
	

	# AUS、TRJ品种BC距离中心性排序，用于挑选样品
	cp ~/rice/miniCore/180319/scripts/beta_pcoa_group.r script/
	# 详见script/beta_pcoa_group.r 

## 4.3 鉴定核心OTU，100%，95%，90%，80%

	# 计算序列两两距离矩阵
	usearch10 -calc_distmx result/otu.fa -tabbedout temp/distmx.txt \
	  -sparsemx_minid 0.9 -termid 0.8
	# 计算核心OTUs
#	usearch11 -otutab_core result/otutab.txt -distmxin temp/distmx.txt -sintaxin temp/otu.fa.tax -tabbedout result/core.txt
#	usearch11 -otutab_core result/otutab.txt -distmxin temp/distmx.txt -tabbedout result/core.txt
	# 距离矩阵存在时报错 ../otutabcore.cpp(61) assert failed: SIZE(Fields) == 3
	# 去掉可计算频率，但缺少某些结果列，但总表受数据量影响，应该用抽平的更合理。
    # 以OTU总表为例
	usearch11 -otutab_core result/otutab.txt -sintaxin temp/otu.fa.tax -tabbedout result/core.txt
	# 采用抽平1万条的序列统计
	usearch11 -otutab_core result/otutab_norm.txt -sintaxin temp/otu.fa.tax -tabbedout result/core.txt
	# 样品数量，7173个
	cat result/otutab.log 
	# 添加统计各OTU频率
	awk '{print $0"\t"$2/7173*100}' result/core.txt > result/core_freq.txt
	# 90%以上42个，85%以上61个，统计各比例下的数量100 95 90 85 80
	mkdir result/core
    for i in `seq 1 100`; do
        echo -ne $i"\t"
        awk -v i="$i" '$14>=i' result/core_freq.txt|wc -l
        #awk -v i="$i" '$14>=i' result/core_freq.txt > result/core/$i.txt
    done

    # 筛选miniCore样品的OTU表，再进行计算核心OTUs
    # 筛选miniCore中根、精选3个样筛选后的结果，result/minicore
    Rscript script/filter_otutab_minicore.R
	usearch10 -otutab_stats result/minicore/otutab.txt -output result/minicore/otutab.stat
	cat result/minicore/otutab.stat
	# 1182 samples, 4996 OTUs, med 25362
	# 和qiime一样的抽样方式，可去年小于阈值的样品，还有抽样后为零的OTUs
	usearch11 -otutab_rare result/minicore/otutab.txt -sample_size 10000 -output result/minicore/otutab10k.txt
	usearch10 -otutab_stats result/minicore/otutab10k.txt -output result/minicore/otutab10k.stat
	cat result/minicore/otutab10k.stat
	# 抽平后，核心OTU数量下降明显，是否需要抽平？
	# 核心OTU，为是确定是否可检测，而不需样品间比较，无须抽平
	usearch11 -otutab_core result/minicore/otutab.txt -sintaxin temp/otu.fa.tax -tabbedout result/core/samples.txt
	awk '{print $0"\t"$2/1182*100}' result/core/samples.txt | sed '1 s/OTU/OTUID/;1 s/0/Core/' > result/core/freq.txt
    # 统计各比例100 95 90%下在miniCore 1182个样品中OTU数量，分别为54，276，370
    for i in 100 95 90; do
        echo -ne $i"\t"
        awk -v i="$i" '$14>=i' result/core/freq.txt|wc -l
        awk -v i="$i" '$14>=i' result/core/freq.txt > result/core/$i.txt
    done
    # 添加可培养和丰度注释


## 4.4 缺失品种缺失原因查找 N4129-412；R4159-448
	
	grep 'N4129' doc/design.txt
	# 有HN下只有2个样，LN下只有一个样，可能我选择3个样品时，只有1个样品的没有考虑而被丢弃
	head -n1 result/otutab.txt|tr '\t' '\n'|grep 'N4129'
	# 找到3个样
	grep 'N4129' result/otutab.biom.sum
	# 数据量为9488，9964，21182

	grep 'R4159' doc/design.txt
	# 有HN下只有1个样，LN下只有1个样
	head -n1 result/otutab.txt|tr '\t' '\n'|grep 'R4159'
	# 找到2	个样
	grep 'R4159' result/otutab.biom.sum
	# 数据量为43276



## 4.5 SL与分蘖相关菌共有关系

	# http://210.75.224.110/report/16Sv2/rice_SL_v1/result-otu.html#result-otu-sum，整理结果doc/SL/result.pptx
	# 合成途径D27, D17, D3突变体(分蘖增加)富集菌可能与分蘖正相关，Venny比较海南根D27-23/D10-34(D17-2太少不考虑)enriched OTUs与分蘖正相关>0.2(otu_top300_mean_pheno_cor.r.txt.tax.xls)的183个OTUs比较，有4个共有，分别为OTU_29，OTU_450，OTU_49，OTU_27，均可培养，3个为Burkholderiales、1个为Actinomycetales；可进一步提高相关性为>0.35的44个OTUs，有三个高丰度共有
	# 比较北京根的D24-/D10-, OTU_27, OTU_241, OTU_49, OTU_116, OTU_11, OTU_2405, OTU_106, OTU_29, OTU_105

    # 刚才分析只用了Top300丰度OTUs，丰度可以实验时进一步选择，现在用全部正相关>0.35的菌otu_all_mean_pheno_cor.r.txt.tax.xls共102个，与D27, D17, D10上调OTUs比较，中间结果保存于result\cor\LN\otu_all_mean_pheno_cor.r.txt.tax.xlsx


	# OTU_29，OTU_49，OTU_27，OTU_667，OTU_33，OTU_371三个高丰度共有菌的丰度分布，并使用OTU_7/9作为对照(NRT1.1b调OTU)
	# ggpubr绘制合成差异明显的2个基因型在两地间差别
    plot_boxplot_ggpubr.sh -i result/otutab.txt -d `pwd`/doc/"SL"/design.txt  -A groupID -B '"d10RtBj","d27RtBj","NpRtBj","d10RtHn","d27RtHn","NpRtHn"' -m '"OTU_49","OTU_29","OTU_27","OTU_7","OTU_9"' -t TRUE -o result/otu_boxplot/ -n true
    # ggplot2绘制SL所有基因型的丰度
	alpha_boxplot.sh -i `pwd`/result/otutab.txt -m '"OTU_667","OTU_33","OTU_371"' \
        -d `pwd`/doc/"SL"/design.txt -A groupID -B '"d27RtBj","d17RtBj","d10RtBj","d3AHLRtBj","d3NpRtBj","d3RtBj","d14AHLRtBj","d14RtBj","d53RtBj","NpRtBj","d27RtHn","d17RtHn","d10RtHn","d3RtHn","d14RtHn","d53RtHn","NpRtHn"' \
        -o `pwd`/result/otu_boxplot/ -h 3 -w 10 -t TRUE -n TRUE


## 4.6 可遗传OTUs

    mkdir -p heritability && cd heritability
    # 基因型数据bed文件，其中bim是SNP列表，fam是样品列表
    cp /mnt/zhou/chulab/miniCore/snp1.5x/T2.b* ./
    cp /mnt/zhou/chulab/miniCore/snp1.5x/T2.modi.fam ./T2.fam
    # 构建亲源关系
    gcta64 --bfile T2 --make-grm --out T2
    # 可遗传OTU分析表型：三列家族、个体ID(通常一样)、表型
    ## 获得LN下各品种OTU的均值LN_otu_mean.txt，并挑选 "OTU_49","OTU_29","OTU_27","OTU_7","OTU_9"
    Rscript ../script/filter_otutab_minicore.R
    # 批量计算遗传力，分别为0.19-0.93之间
    cut -f 1 pheno_test.txt.hsq | tr '\n' '\t' | sed 's/\t$/\n/' > heritable.txt
    for i in `seq 3 7`; do
        cut -f 1,2,${i} pheno_5.txt>pheno_test.txt
        gcta64 --reml --grm T2 --pheno pheno_test.txt --out pheno_test.txt
        echo -ne $i"\t" >>heritable.txt
        cut -f 2- pheno_test.txt.hsq|tail -n+2|sed 's/\t/+-/'|tr '\n' '\t' | sed 's/\t$/\n/' >>heritable.txt
    done
    # 计算所有OTU的遗传力作为属性
    #cut -f 1 pheno_test.txt.hsq | tr '\n' '\t' | sed 's/\t$/\n/' > heritable.txt 
    # 根据结果pheno_test.txt.hsq的行名修改为自定义的表头
    cp ../doc/heritable.header heritable.txt
    sed -i 's/Source/OTUID/' heritable.txt
    otu_num=`head -n1 LN_otu_gcta.txt|tr '\t' '\n'|awk 'END{print NR}'`
    for i in `seq 3 $otu_num`; do
        cut -f 1,2,$i LN_otu_gcta.txt | tail -n+2 > pheno_test.txt
        gcta64 --reml --grm T2 --pheno pheno_test.txt --out pheno_test.txt
        j=`cut -f $i LN_otu_gcta.txt | head -n1`
        echo -ne $j"\t" >>heritable.txt
        cut -f 2- pheno_test.txt.hsq|tail -n+2|sed 's/\t/+-/'|tr '\n' '\t' | sed 's/\t$/\n/' >>heritable.txt
    done
    sed 's/+-/\t/g' heritable.txt|less -S > heritable_otu.txt
    # 筛选显著可遗传的OTUs，基因型解析量V(G)2列，V(G)/Vp遗传力8列，Pval14列,显著性优先，再看解析率
    cut -f 1,2,8,14 heritable_otu.txt > heritable_otu.elite
    # 很多基因型贡献为0的值，筛选V(G)>0有1712个
    awk '$2>0' heritable_otu.elite | less -S | wc -l
    # 遗传力>0.15，1955个
    awk '$3>0.15' heritable_otu.elite | less -S | wc -l
    # 筛选Pvalue<0.05，有1328个显著
    awk '$4<0.001' heritable_otu.elite | less -S | wc -l 
    # 筛选以上三个条件，有943个
    awk '$2>0 && $3>0.15 && $4<0.001' heritable_otu.elite | less -S | wc -l
    # 按V(G)、V(G)/Vp和Pval排序查看
    awk '$2>0 && $3>0.15 && $4<0.001' heritable_otu.elite | sort -k2,2nr |less -S
    awk '$2>0 && $3>0.15 && $4<0.001' heritable_otu.elite | sort -k3,3nr |less -S
    awk '$2>0 && $3>0.15 && $4<0.001' heritable_otu.elite | sort -k4,4g |less -S


## 4.7 制作OTU注释表

    mkdir -p result/anno
    # 基于物种注释制作1_2列
    sed '1 i OTUID\tTaxonomy' result/taxonomy_2.txt | less -S > result/anno/otu.1_2tax
    # 添加所有样品平均丰度Mean，相似度pident，菌保ID和物种注释3——6列，无注释的OTU会为空
    #awk 'NR==FNR{a[$1]=$5"\t"$2"\t"$3"\t"$6} NR>FNR{print $0"\t"a[$1]}' result/39culture/otu.txt result/anno/otu.1_2tax > result/anno/otu.3_6culture
    awk 'NR==FNR{a[$1]=$5"\t"$2"\t"$3"\t"$6} NR>FNR {if(a[$1]==""){a[$1]="\t\t\t"};print $0"\t"a[$1]}' result/39culture/otu.txt result/anno/otu.1_2tax > result/anno/otu.3_6culture
    # 添加核心比例7列
    awk 'NR==FNR{a[$1]=$14} NR>FNR{print $0"\t"a[$1]}' result/core/freq.txt result/anno/otu.3_6culture > result/anno/otu.7core
    # 添加可遗传解析率、比例和p值，8-10列
    awk 'NR==FNR{a[$1]=$2"\t"$3"\t"$4} NR>FNR{print $0"\t"a[$1]}' heritability/heritable_otu.elite result/anno/otu.7core > result/anno/otu.8_10heritable
       

## 4.8 选菌验证
    
    date=wet/180813
    cluture_db=/mnt/bai/yongxin/culture/rice/result/culture_select.fasta
    format_seq2fasta.pl -i "${date}/*.seq" -o ${date}.fa 
    # 输出13列为coverage
    blastn -query ${date}.fa -db ${cluture_db} -out ${date}.xls -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 
    sed -i '1 i qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' ${date}.xls
    sed -i '1 s/ /\t/' ${date}.xls




# 5. 图表整理 Figures and legends

    # 位于doc目录中，方便同步

##  图1. 水稻miniCore品种和OTU描述 fig1.description.variety.sample.taxonomy.Rmd

### 1A. 203个水稻品种全球分布
	
	地图，种植模式图——秦媛

### 1B, 样品、品种整理多样性、OTU稀释取线 script/alpha_rare_sample.R

### 1C. 物种组成Taxonomy boxplot in genus and phylum

### 1D. 进化树

    # 选择丰度>0.1%的166个序列建树，并用iTOL注释物种门水平颜色、核心菌加粗和可遗传标为红色
    # select OTU abundance > 0.1% 
    awk '$3>0.1' result/anno/otu.8_10heritable | cut -f 1 > temp/otu_k1.id
    usearch10 -fastx_getseqs result/otu.fa -labels temp/otu_k1.id -fastaout script/fig1/1d.otu_k1.fa
    # check number
    grep -c '>' script/fig1/1d.otu_k1.fa
    # Multiply alignment
    clustalo -i script/fig1/1d.otu_k1.fa -o temp/1d.otu_k1_align.fa --seqtype=DNA --full --force --threads=9
    make_phylogeny.py -i temp/1d.otu_k1_align.fa -o script/fig1/1d.otu_k1.tree
    sed -i "s/'//g" script/fig1/1d.otu_k1.tree
    # 在iTOL中展示树
    # 注释物种四大门水平颜色
    # Phylum = c("Proteobacteria", "Actinobacteria", "Bacteroidetes", "Firmicutes", "Other"), Color = c("#85F29B", "#F58D8D", "#F7C875", "#91DBF6", "#AAAAAA" )
    Rscript ~/github/Amplicon/16Sv2/script/tree_color_iTOL.R
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3" "$2} NR>FNR{print $0" range "a[$1]}' result/taxonomy_8_color.txt temp/otu_k1.id
    # 标注核心OTUs $7>100% 蓝色加粗2号
    awk '$3>0.1 && $7>=100 ' result/anno/otu.8_10heritable | cut -f 1 | awk '{print $0" label #0000ff bold 2"}'
    # 标注可遗传OTUs枝标为红色
    awk '$3>0.1 && $9>=0.3 ' result/anno/otu.8_10heritable | tail -n+2 | cut -f 1 | awk '{print $0"|"$0" clade #ff0000 normal 1"}'





	# 不同加土，否则主要差异为土壤；不能两块地混合，否则主要差异为不同地块
	# HN下TEJ和IND
beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
	-d `pwd`/doc/design.txt -A groupID -B '"HIND","HTEJ"' -E TRUE \
	-c `pwd`/doc/compare.txt \
	-o `pwd`/fig1/1subspecies/beta_HN_ -h 3 -w 5
# 匹配非注释行，输出用于发表
grep -P '^\s*#' script/beta_pcoa.R | less
grep -P -v '^\s*#' script/beta_pcoa.R > fig1/script/beta_pcoa_fieldII.R
	# LN下TEJ和IND
beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
	-d `pwd`/doc/design.txt -A groupID -B '"LIND","LTEJ"' -E TRUE \
	-c `pwd`/doc/compare.txt \
	-o `pwd`/fig1/1subspecies/beta_LN_ -h 3 -w 5
grep -P -v '^\s*#' script/beta_pcoa.R > fig1/script/beta_pcoa_fieldI.R


### alpha多样性
	
	箱线图、稀释取线、样品稀释取线

	alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' \
	-o `pwd`/fig1/1/alpha_ -h 3 -w 5
	# 稀释曲线
	alpha_rare.sh -i `pwd`/result/alpha/rare.txt \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' \
	-o `pwd`/fig/1/alpha_ -h 3 -w 5

### Taxonomy 门+变形菌纲
	cut -f 3-4 result/taxonomy_8.txt|sort|uniq|grep 'Proteobacteria' # 为什么会有这么多结果，只选5类继续分析
	cat <(grep -v 'Proteobacteria' result/tax/sum_p.txt) <(grep 'proteobacteria' result/tax/sum_c.txt) > result/tax/sum_pc.txt
tax_stackplot.sh -i `pwd`/result/tax/sum_ -m '"pc"' -n 10 \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' -O FALSE \
	-o `pwd`/fig1/1/tax_pc_ -h 3 -w 5


### 门及纲水平差异


## 图2. 随机森林分类
	
	# 用family水平建模，用HN数据training，用LN验证。randomForest_family.R
	randomForest_class.R
	1. 纲水平建模，展示贡献度，和样品中热图
	使用高HN和HN下籼粳稻纲水平0.3%丰度的15个Feature机器学习；保存预测结果confusion.txt，整理16.4%错误率，TEJ 37.4%错误；
	2018/5/21 删除三个澳大亚利(纬度为负)粳稻, D4032, D4038, F4053; 标记A50/ZH11为TEJ，而IR24为IND，各分为Hn/Ln两种情况；
	错误率降低为15.3%，TEJ为36%;Top1 feature也变为了Nitrospira，Deltaproteobacteria
	2. Top feature：用各组柱状图/箱线图分类展示，再加梯度排序
	3. 在nrt和时间序列中验证

	
## 图3. 亚种差异与氮相关

	1. 差异OTUs曼哈顿图，维恩图
	由筛选组，改为筛选亚组(品种)中位数的OTUs: 原万5为343个OTUs，万一为942个；最终丰度为0.2%
	compare_sub.R # 修改丰度筛选group为groupID2，接下来 rm plot_volcano ; make plot_venn; make rmd

grep -P -v '^\s*#' script/compare_sub.R > fig1/script/compare.R

	差异OTUs在两块地曼哈顿、韦恩图;	曼哈顿图要写标颜色为门、纲,	plot_manhattan_pc.r
	plot_manhattan.sh -i result/compare/LTEJ-LIND_all.txt
	# 我们重点是突出IND，大多数是IND特异的，添加IND vs TEJ的组，重画曼哈顿图，让IND为向上实心三角

	# 曼哈顿图的代码和数据
	grep -P -v '^\s*#' script/plot_manhattan_pc.r > fig1/script/plot_manhattan_pc.r
	cp result/tax/sum_pc.txt fig1/data/
	cp result/compare/*IND-*TEJ_all.txt fig1/data/

	# 维恩图的代码和数据
	cp result/compare/diff.list* fig1/data/
	diff.list.vennHTEJ_HIND_DLTEJ_LIND_DCDE.r

	2. 差异OTUs在时间序列中变化
	alpha_boxplot.sh -i result/tax/sum_c.txt -d `pwd`/doc/design.txt -A groupID -B '"A50Cp0","A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119"' \
	-m '"Deltaproteobacteria","Actinobacteria","Alphaproteobacteria","Clostridia","Betaproteobacteria","Nitrospira"' -t TRUE -o result/randomForest/time_ -n TRUE # 42以后没有数据呢？改用alpha_boxplot.sh

	3. 绘制维恩图的共有饼图
	# fig1/3compare.rmd 绘制时间序列的图，最后添加共有的成份
	# 筛选HL/LN下TEJ-IND共同下调的菌
	cat result/compare/?TEJ-?IND_sig.txt | grep 'Depleted' | cut -f 1 | sort | uniq -d > fig1/3compare/otu_IND_common_specific.txt
	cat result/compare/?TEJ-?IND_sig.txt | grep 'Enriched' | cut -f 1 | sort | uniq -d > fig1/3compare/otu_TEJ_common_specific.txt
	# 两块地保守上调、下调的OTUs
	cat fig1/3compare/otu_IND_common_specific.txt fig1/3compare/otu_TEJ_common_specific.txt > fig1/3compare/otu_common.txt
	# 并用faprotax注释
	# 绘制nrt和A50差异与籼粳稻共有
	tail -n 70 ~/rice/xianGeng/fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.txt | cut -f 1 > fig1/4nrt/venn_nrt_indiaHL.txt
	tail -n 6 ~/rice/xianGeng/fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D.txt | cut -f 1 > fig1/4nrt/venn_NRTsnp_indiaHL.txt


	5. 差异菌功能有无热图 plot_heatmap_timecourse.R
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
	# plot_heatmap_timecourse.R 绘制时间序列的图，再添加相应菌的主要功能，
	# 同时对时间序列中不表达的也可视化功能:IND的功能绘制于 fig/2/otu_IND_common_specific_time_faprotax_noabundance.txt，TEJ单一条目录为 fig/2/otu_TEJ_common_specific_time_faprotax_noabundance.txt"
	# 再对时间序列中0点去掉重新计算，发现分为了4组，在原文件基础上添加-0标志


	# 菌种功能注释整体差异
	rm result/compare_far/diff.list
	compare.sh -i `pwd`/result/faprotax/element_tab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
	-p 0.01 -q 0.05 -F 1.2 -t 0.0005 \
	-d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1","V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50LnCp6","A56LnCp6"' \
	-o `pwd`/result/compare_far/ -N FALSE
	batch_venn.pl -i doc/venn.txt -d result/compare_far/diff.list
	# 注释比较结果
	rm result/compare_far/diff.list.venn*.xls.*
	batch2.pl -i 'result/compare_far/diff.list.venn*.xls' -d result/compare_far/database.txt -o result/compare_far/ -p vennNumAnno.pl
	# 绘制箱线图
	# 确定要展示的Features，两组共有
	tail -n 39 result/compare_far/diff.list.vennHTEJ_HIND_DLTEJ_LIND_D.xls.xls|cut -f 1 > result/compare_far/IND.list
	tail -n 9 result/compare_far/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls|cut -f 1 > result/compare_far/TEJ.list
	make plot_fa_barplot # 绘制单个功能的箱线图
	# 修改alpha_boxplot.R为alpha_boxplot_far.R

grep -P -v '^\s*#' script/alpha_boxplot_far.R > fig1/script/alpha_boxplot_far.R

# 图4. 菌群与nrt关系

	## 2. 差异菌/功能与氮相关基因显著相关
	# 来自胡斌整理的氮相关基因doc/N-related genes in rice.docx共9个基因，先在doc/rice_nitrogen_list.xlsx中惠惠相关ID，保存为doc/rice_nitrogen_list.txt, 其中第5列RAP_SNP的ID与SNP注释文件对应
	dos2unix doc/rice_nitrogen_list.txt # windows转换为linux
	# 整理出SNP数据中对应的基因型、提取相应的位点
	mkdir -p result/nitrogen_cor
	cut -f 5 doc/rice_nitrogen_list.txt | tail -n+2 | tr '\n' '|' # 提取ID并替换为|分隔
	grep -P 'OS10G0554200|OS08G0155400|OS02G0112100|OS02G0595900|OS01G0704100|OS01G0547600|OS03G0687000|OS04G0509600|OS06G0706400|OS06G0706500' /mnt/bai/yongxin/rice/miniCore/180319/gemma/snp.anno > result/nitrogen_cor/all_snp.list # 筛选到10个基因在miniCore中存在502个相关位点
	cut -f 3 result/nitrogen_cor/all_snp.list | sort | uniq -c # 8个MODERATE，467个MODIFIER和27个LOW
	grep -P 'HIGH|MODERATE' result/nitrogen_cor/all_snp.list > result/nitrogen_cor/good_snp.list # 其中重要SNP仅有8个，来自4个基因
	# 提SNP对应基因型
	cat /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_*.hmp.txt > /tmp/temp
	cut -f 1 result/nitrogen_cor/good_snp.list|tr '\n' '|' # 获取列表
	grep -P '2m655515\t|2m657013|6m29839102|6m29839240|8m3183208|10m21759092|10m21761740|10m21761997' /tmp/temp > /tmp/temp1
	cat <(head -n1 /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_1.hmp.txt) /tmp/temp1 > result/nitrogen_cor/good_snp.geno # 添加标题
	# 用excel转置 result/nitrogen_cor/good_snp.geno.t，添加注释
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$0} NR>FNR{print $0,a[$1]}' ../miniCore/doc/minicore_list.txt result/nitrogen_cor/good_snp.geno.t > result/nitrogen_cor/good_snp.geno.t.txt


	# 统计SNP基因型作为分组信息，来统计氮功能丰度组间P值
	nitrogen_cor.r # 保存实验设计+基因型，方便识别SNP不同基因型在籼粳稻中区别
	# 整理这4个重要SNP信息表，见SNP_list.xlsx
	# 统计基因型与亚种分布
	sed -i '1 s/^/SampleID\t/' result/nitrogen_cor/design.txt
	head -n1 result/nitrogen_cor/design.txt|tr '\t' '\n'|awk '{print NR,$1}'
	# NRT2.1 - 17; 1.1A - 21; 1.1B - 14，统计每个亚种内基因型的数量，可看到亚种内主要的SNP类型
	cut -f 2,8,20 result/nitrogen_cor/design.txt|sort|uniq|cut -f 2,3|sort|uniq -c
	grep -P '10m21759092|2m655515\t|8m3183208' /mnt/bai/yongxin/rice/miniCore/180319/gemma/T2.ann.vcf # 查询SNP变化位置、碱基和AA详细

grep -P -v '^\s*#' script/nitrogen_cor.r > fig1/script/nitrogen_cor.R


	3. 关键氮高效基因不同形态、突变体可部分解析亚种差异 2018/5/31
	"A50LnCp6","A56LnCp6","A50LnCp7","A56LnCp7","A50LnSz7","A56LnSz7","A50HnCp6","A56HnCp6","A50HnCp7","A56HnCp7","A50HnSz7","A56HnSz7","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","nrtHnCp7","ZH11HnCp7","ZH11LnCp7","nrtLnCp7","nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7"
	# nrt vs ZH11(TEJ): 主图："V3703HnCp6","ZH11HnCp6", 附图："V3703LnCp6","ZH11LnCp6",
	# 近等基因系：主图："A50LnCp7","A56LnCp7", 附图："A50LnCp6","A56LnCp6",
	# 主图/附图各分析一次：alpha, beta, 差异OTUs, venn: HTEJ_HIND_D LTEJ_LIND_D V3703HnCp6_ZH11HnCp6; HTEJ_HIND_D LTEJ_LIND_D A50HnCp7_A56HnCp7
	# 主图4组，xiangeng_wilcoxon_main
	# 附图4组，xiangeng_wilcoxon_supp


## 图4.e宏基因组KO注释

	从金桃处获得KO表和丰度，获取KO的功能描述；在kegg中没找到，picurst中没找到(输出结果有注释但不完整)，google搜索KO description download，找到biostar解答；https://www.genome.jp/kegg-bin/get_htext?ko00001.keg 中的Download htext/jason下载KO和描述；htext方便检索，而jason方便在线分析，如jason2table
	
	cd ~/rice/xianGeng/metagenome
	grep 'D      K' ko00001.keg | cut -c 8- | sed 's/  /\t/' > KO_description.txt #|  cat -A | less -S
	# 比较两组差异
compare.sh -i ko.txt -c compare.txt -m "wilcox"         -p 0.01 -q 0.05 -F 1.2 -t 0 -N FALSE         -d design.txt -A group -B '"HnNrt","HnZH11"'         -o compare/ -N FALSE -U 100 # Pvalue和FDR并不显著，可能秩和检验需要较多的样本数
	# 注释KO
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' KO_description.txt ko.txt > ko.anno

	# 绘制单个KO和箱线图
alpha_boxplot.sh -i ko.txt -m '"K02568","K02567","K00363","K10535","K00362"' \
-d design.txt -A group -B '"HnZH11","HnNrt"' \
-o ./ -h 2 -w 2.5 -t TRUE -n TRUE -U 1000000



## 图5. 微生物与表型关联

	1. 微生物、多样性、PCoA主轴与表型相关分析
	# 先使用胡斌整理表型数据+faprotax中氮通路相关
	# 方法1：script/phenotype_cor.R直接关联，spearman相关系数只有0-0.2，但能看到nitrogen_amonification正相关，而固
	
	# 方法2. 采用分组协变量关联，需要基因型的聚类/PCA信息
	# PCA信息 /mnt/bai/yongxin/rice/miniCore/180319/gemma/pca4.txt
	# 行名 head -n1 /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_1.hmp.txt
	head -n1 /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_1.hmp.txt|cut -f 12-|tr '\t' '\n' > temp.txt
	paste temp.txt /mnt/bai/yongxin/rice/miniCore/180319/gemma/pca4.txt | cut -f 1,3- | sed '1 i variety\tPC1\tPC2\tPC3\tPC4' | less > fig/4cor/genotype_pca4.txt
	
	# OTU/属水平(比OTU数量少至有描述)与张小宁整理表型关联(品种对应): HN/LN分开关联phenotype_cor.Rmd
	cp ~/rice/miniCore/mwas/phenotype/minicore?NPhenotype.txt doc/ # 准备miniCore中HN/LN表型
	
	# 表型-faprotax相关 phenotype_cor_faprotax.Rmd

	# 2018/6/5 相关热图+注释
	# LN条件下OTUs与表型关联，筛选>0.4相关的值进行注释物种，差异OTUs，和功能 phenotype_cor2.Rmd ；表型数据重新整理minicore低氮下为单株，excel整理
	sed -i '/\t0\t0\t0/d' doc/phenotype_sample_raw.txt # 删除缺失样品
	# 株的OTU对应株的表型，OTU_11与tiller相关仅为0.33(可能植株波动大，或不对应，会规律只有平均才能看出来)，而OTU_11与分均值还有0.44的相关。那OTUs的均值对应是否会更高呢？改为均值结果更好。

	# 2018/6/11 ## 筛选相关系数 > 0.4，# 添加faprotax有无 X OTU丰度，绘制泡泡图，line 310
	
	# 2018/6/22 用品种中全部OTUs计算相关系数，筛选相关系数>0.4，丰度大于均值0.1%的OTUs，物种名改为低级注释
	cp xiangeng4wilcox/result/compare/LTEJ-LIND_all.txt fig1/5pheno/ # 品种0.2%过滤OTU列表


附图.

1. 样本稀释曲线，平滑处理
alpha_rare_sample.R

2. 2块地所有样品IND/TEJ一起PCoA, CPCoA
3. 主图c/d的unifrac距离图
4. alpha diversity chao1 richness
5. 机器学习层级选择准确度柱状图
6. 共有OTUs饼图
7. 完整的OTUs时间序列和功能注释
8. 主图的功能差异分析完整版本
9. 图4a的，其它N相关功能与SNP关联
10. 图4b的beta其它距离， constraind
11. 附图专用4组多样性分析

附数据统计：

1. 品种数：68 IND
cd ~/rice/xianGeng/fig1/ST
     27 TEJ
sed -i 's/,/\t/g' 01.variety_geotable.csv # csv替换为tsv
# 数据筛选95个，并添加95个
awk '{FS=OFS="\t"} NR==FNR{a[$2]=$0} NR>FNR{print $0,a[$2]}' ~/rice/miniCore/doc/minicore_list.txt 01.variety_geotable.csv > 01.variety_geotable.txt 
cut -f 14 01.variety_geotable.txt|tail -n+2|sort|uniq|wc -l # 44个国家

2. 凌水海南样本数和测序量
# 需要使用compare_sub.R中来原代码，在3compare中汇总样品和测序量

3. 查看共有16个下调菌的丰度范围
cd ~/rice/xianGeng
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$14} NR>FNR{print a[$1]}' xiangeng4wilcox/result/compare/HTEJ-HIND_all.txt fig1/3compare/otu_TEJ_common_specific.txt|sort -nr
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$14} NR>FNR{print a[$1]}' xiangeng4wilcox/result/compare/LTEJ-LIND_all.txt fig1/3compare/otu_TEJ_common_specific.txt|sort -nr
awk '$2<0.05' fig1/ST/09.tiller_cor_p.txt|wc -l # 381个有214个显著P<0.05
cd fig1/ST
paste 09.tiller_cor.txt 09.tiller_cor_p.txt | cut -f 1,2,4> 09.tiller_cor_pr.txt



# 上传数据 PRJNA478068
# 整理的最终上传实验设计 fig1/metadata.txt
cut -f 16 fig1/metadata.txt|sort|uniq -c # 587个亚种数据来自miniCore，127+114=241来自nrt
wc -l fig1/metadata.txt # 检查是否有样品重名
cut -f 1 fig1/metadata.txt|sort|uniq|wc -l 
mkdir -p seq/submit # 创建提交数据目录

# 1. 籼粳稻样本，587个
for RPM in `grep 'minicore' fig1/ST/02.design.txt|cut -f 1`; do
	cp ../miniCore/clean_data/sample/${RPM}.fq.gz seq/submit/
done

# 2. 时间序列己上传，见 https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA435900
for RPM in `grep 'nrt' fig1/ST/02.design.txt|cut -f 1`; do
	cp /mnt/bai/yongxin/rice/zjj.nitrogen/180116/clean_data/sample/${RPM}.fq.gz seq/submit/
done



# 共享中间文件和分析流程代码 submit用于投稿，publish用于正式发表后共享

## 在fig1中创建index.Rmd并生成网页 http://210.75.224.110/submit/rice_microbiome, username: rice, password: microbiome
cd ~/rice/xianGeng/fig1 # 共享目录
ln -sf `pwd` /var/www/html/submit/rice_microbiome # 链接至外网
cp ~/github/Amplicon/16Sv2/rmd/.htaccess ./ # 加密
htpasswd /mnt/bai/yongxin/bin/config/users rice # 添加新用户和密码


## 获得分析流程
pipeline=fig1/pipeline.sh
make -n -B library_split_stat|grep -v '#' > $pipeline
make -n -B fq_qc|grep -v '#' >> $pipeline
make -n -B beta_calc|grep -v '#' >> $pipeline

## 根据最新实验设计和OTU表整理结果
cd fig1/
cp /mnt/bai/yongxin/github/Amplicon/16Sv2/script/stat_plot_functions.R script/
