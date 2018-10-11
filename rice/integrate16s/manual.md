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

    # 不同组的版本管理

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
	

    ## 比较SL_Hn
    sub=SL_Hn
	mkdir -p doc/${sub}
	grep 'Hn' doc/SL/compare.txt > doc/${sub}/compare.txt
    # 编号比较韦恩图venn.txt
	# 基于比较组的编号给alpha, beta小组绘图；cat doc/${sub}/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'，总体给g1_list保证不同组差异比较一致

    ## 比较SL_Bj
    sub=SL_Bj
	mkdir -p doc/${sub}
	grep 'Bj' doc/SL/compare.txt > doc/${sub}/compare.txt
    sed 's/Hn/Bj/g' doc/SL_Hn/venn.txt > doc/${sub}/venn.txt
	# 基于比较组的编号给alpha, beta小组绘图；cat doc/${sub}/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'，总体给g1_list保证不同组差异比较一致

    ## 比较nrt Bj 2016土
    cp -r doc/xiangeng doc/nrt
	sub=nrt
    # 手动编辑compare.txt, venn.txt
    rm alpha_boxplot
    make plot_venn # DA otu
    make DA_compare_tax2 # DA taxonomy
    make rmd # report
    mkdir -p result/nrt/2016/LN
    version=rice_nrt_wilcox_soil16LN2/result/compare
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$14} NR>FNR {print $0,a[$1]}' $version/ZH11LnCp6-soilLnCp6_all.txt $version/diff.list.vennZH11LnCp6_soilLnCp6_DIR24LnCp6_soilLnCp6_DA50LnCp6_soilLnCp6_D.xls.xls > result/nrt/2016/LN/soil_enriched1
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$14} NR>FNR {print $0,a[$1]}' $version/IR24LnCp6-soilLnCp6_all.txt result/nrt/2016/LN/soil_enriched1 > result/nrt/2016/LN/soil_enriched2
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$14"\t"$15} NR>FNR {print $0,a[$1]}' $version/ZH11LnCp6-soilLnCp6_all.txt result/nrt/2016/LN/soil_enriched2 > result/nrt/2016/LN/soil_enriched3.xls

    # 查看土壤里特异菌的
	alpha_boxplot.sh -i `pwd`/result/otutab.txt -m '"OTU_61","OTU_95","OTU_256","OTU_49","OTU_421","OTU_150","OTU_959","OTU_261","OTU_456","OTU_347","OTU_2881","OTU_991","OTU_3077"' \
        -d `pwd`/doc/"nrt"/design.txt -A groupID -B '"ZH11LnCp6","IR24LnCp6","A50LnCp6","soilLnCp6"' \
        -o `pwd`/result/otu_boxplot_soil/ -h 3 -w 5 -t TRUE -n TRUE

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

    # LN条件下OTU均值
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
    awk '$2>0 && $3>0.15 && $4<0.00001' heritable_otu.elite | less -S | wc -l # 615
    awk '$2>0 && $3>0.4 && $4<0.00001' heritable_otu.elite | less -S | wc -l # 428
    # 按V(G)、V(G)/Vp和Pval排序查看
    awk '$2>0 && $3>0.15 && $4<0.001' heritable_otu.elite | sort -k2,2nr |less -S
    awk '$2>0 && $3>0.15 && $4<0.001' heritable_otu.elite | sort -k3,3nr |less -S
    awk '$2>0 && $3>0.15 && $4<0.001' heritable_otu.elite | sort -k4,4g |less -S
    # 结果遗传力偏高，尝试用H/LN所有数据
    mv heritable.txt LN_heritable.txt # bak to specific name
    awk '{a=a+$3} END {print a}' heritable_otu.elite # total heritability is 1281.79

    # H/LN条件下OTU均值A，或A/L/H三种
    sub=H
    mkdir $sub
    # 根据结果pheno_test.txt.hsq的行名修改为自定义的表头
    cp ../doc/heritable.header $sub/heritable.txt
    # `head -n1 ${sub}N_otu_gcta.txt|tr '\t' '\n'|awk 'END{print NR}'`
    otu_num=4998
    for i in `seq 3 $otu_num`; do
        cut -f 1,2,$i ${sub}N_otu_gcta.txt | tail -n+2 > ${sub}/pheno_test.txt
        gcta64 --reml --grm T2 --pheno ${sub}/pheno_test.txt --out ${sub}/pheno_test.txt
        j=`cut -f $i ${sub}N_otu_gcta.txt | head -n1`
        echo -ne $j"\t" >> ${sub}/heritable.txt
        cut -f 2- ${sub}/pheno_test.txt.hsq|tail -n+2|sed 's/\t/+-/'|tr '\n' '\t' | sed 's/\t$/\n/' >>$sub/heritable.txt
    done
    sed 's/+-/\t/g' $sub/heritable.txt|less -S > $sub/heritable_otu.txt
    # 筛选显著可遗传的OTUs，基因型解析量V(G)2列，V(G)/Vp遗传力8列，Pval14列,显著性优先，再看解析率
    cut -f 1,2,8,14 $sub/heritable_otu.txt > $sub/heritable_otu.elite
    # 查看可遗传汇总，2263.95，平均之后更高了？？？因为没有初始化。实际为1118.84，较低氮下降了，14%；HN为982.707
    awk '{a=a+$3} END {print a}' $sub/heritable_otu.elite # total heritability is 1281.79
    # 查看3个分蘖正相关菌的遗传力
    grep -P 'OTU_27\t|OTU_29\t|OTU_49\t' $sub/heritable_otu.elite
    grep -P 'OTU_27\t|OTU_29\t|OTU_49\t' heritable_otu.elite
    # P值变显著了，但贡献在变小
    # 筛选以上三个条件：基因贡献>0，遗传率>0.15，P值<0.001，有820个
    awk '$2>0 && $3>0.15 && $4<0.001' $sub/heritable_otu.elite | less -S | wc -l
    awk '$2>0 && $3>0.15 && $4<0.00001' $sub/heritable_otu.elite | less -S | wc -l # 499
    awk '$2>0 && $3>0.4 && $4<0.00001' $sub/heritable_otu.elite | less -S | wc -l # 305
    # 按V(G)、V(G)/Vp和Pval排序查看
    awk '$2>0 && $3>0.15 && $4<0.001' $sub/heritable_otu.elite | sort -k2,2nr |less -S
    awk '$2>0 && $3>0.15 && $4<0.001' $sub/heritable_otu.elite | sort -k3,3nr |less -S
    awk '$2>0 && $3>0.15 && $4<0.001' $sub/heritable_otu.elite | sort -k4,4g |less -S

    # 遗传力分布：所有，和measureable OTUs
    # 筛选miniCore的measureable OTUs 0.1% in one of variety
    # script/filter_otutab_minicore.R己经筛选了minicore中挑选的植物，并去除土壤的样本
    # 此处仅筛选OTUs，按所有品种中位数1%筛选为513个
    filter_otus_by_group_median.sh -i result/minicore/otutab.txt \
        -d `pwd`/doc/"xiangeng"/design.txt  -A groupID -B '' -t 0.1 \
        -o result/minicore/otutab0.1.txt -C groupID2
    cut -f 1 result/minicore/otutab0.1.txt | tail -n+2 | less > result/minicore/otutab0.1.id
    # fig2. 2.D 遗传力分布，0.1以上最多，0.2起逐渐减少，但0.9以上又升高

    # 使用重复和添加协变量，不支持
    mkdir zhouyao
    cp T2.b* T2.f* LN_otu_gcta.txt zhouyao/
    tar zcvf zhouyao.tar.gz zhouyao/
    # 添加重复的数据
    cd zhouyao
    cp -r  /mnt/bai/qinyuan/test/rice/* ./
    # design.txt为doc/xiangeng/design.txt，all.txt为输出添加注释的OTU，otutab0.1.txt来自result/minicore中筛选measureable OTU，design2.txt为手动调整添加重复的样本
    # addNA.R 调整字母编号为1、2、3，没有的样本补为NA缺失phenotype.txt，
    # /mnt/bai/yongxin/rice/integrate16s/heritability/zhouyao 周姚计算为heritabity.txt，代码为Heritability.R第一个是广义遗传力，第二个和第三个分别使用blup和mean计算的遗传力，原文见邮件
    
    # 遗传力表整理，我公式调了多次，结果影响不大，以zhou为准继续，高遗传力的菌只在极少样本中出现，不可靠。




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
    # 180812/180815/180921(21拼接，30末拼接)
    date=wet/180930
    cluture_db=/mnt/bai/yongxin/culture/rice/result/culture_select.fasta
    # 条件1. 原始seq序列采用format_seq2fasta.pl合并为fa
    format_seq2fasta.pl -i "${date}/*.seq" -o ${date}.fa 
    # 条件2. 双端合并文件fa采用cat合并
    # cat ${date}/*.seq.txt | sed 's/>/\n>/' | sed '/^$/d' > ${date}.fa
    # 输出13列为coverage
    blastn -query ${date}.fa -db ${cluture_db} -out ${date}.xls -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 
    sed -i '1 i qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' ${date}.xls
    sed -i '1 s/ /\t/g' ${date}.xls
    blastn -query ${date}.fa -db ${cluture_db} -out ${date}.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1000 -evalue 1 -num_threads 9 


## 4.9 基因型与Beta多样性

    # PNAS-2018-多年多点5千样本鉴定玉米根际可遗传微生物
    # 在Hapmap2中查询NAM祖先并制作NAM亲源关系，过滤在至少10个祖先中出现的点。使用TASSEL产生遗传距离矩阵。使用beta距离矩阵和kinship矩阵计算R2。
    # gcta中计算的--make-grm(genetic relationship matrix)为二进制heritability/T2.grm.bin
    # gemma计算过kinship矩阵 /mnt/bai/yongxin/rice/miniCore/180319/gemma/kinship.txt
    # 计算kinship与beta_diversity相关性
    # 详见script/fig2.diversity.genetics.rmd ## 2.C
    # 结果相关性为-0.46-0.26，以负相关为居多，因为kinship为相关性，改为1-similarity=distance

    # 方法2. tassel计算kinship矩阵

## 4.10 筛选SL使用样品
    mkdir -p SL
    # 按beta距离，15个重重复保留12个
    Rscript ~/github/Amplicon/16Sv2/script/beta_pcoa_group.r -i temp/beta/bray_curtis.txt -N 12 -d doc/design_SL.txt -n groupID -o SL/pcoa_bray
	# 备份实验设计
	cp doc/design.txt doc/design.txt180821
	# 获得删除样品
    cut -f 1 SL/pcoa_bray* | sort | uniq -u > SL/discard.id
	# 剔除点注释design
    for i in `cat SL/discard.id`; do \
        sed -i "s/^$i/#$i/" doc/design.txt;done

## 4.11 SL 6个基因型热图展示均值或中位数 script/fig4.SL_mutant.rmd

    # 查看是否存在合成、信号通路一致或分类的结果
    mkdir -p script/fig4
    # 筛选差异OTU用于展示
    cat <(tail -n12 rice_SL_Bj_wilcox_v1/result/compare/diff.list.vennd27RtBj_NpRtBj_Dd17RtBj_NpRtBj_Dd10RtBj_NpRtBj_D.xls.xls) \
    <(tail -n71 rice_SL_Bj_wilcox_v1/result/compare/diff.list.vennd27RtBj_NpRtBj_Ed17RtBj_NpRtBj_Ed10RtBj_NpRtBj_E.xls.xls) \
    | cut -f 1 > script/fig4/Bj_d27_d17_d10_common.id
    # 所有Bj共有差异OTU
    cat rice_SL_Bj_wilcox_v3/result/compare/diff.list.vennd27RtBj_NpRtBj_Ed17RtBj_NpRtBj_Ed10RtBj_NpRtBj_E.xls.xls | grep 'OTU_' | cut -f 1 | sort| uniq | less -S > script/fig4/Bj_d27_d17_d10_common.id # rice_SL_Bj_wilcox_v3/result/compare/diff.list.vennd27RtBj_NpRtBj_Dd17RtBj_NpRtBj_Dd10RtBj_NpRtBj_D.xls.xls 
    # ## Heatmap of genotypes 中查看mat_mean_final变量，居然有的OTU_153|OTU_440|OTU_510在所有组中为零？差异基因不可能全为零的，检查原因
    grep -P 'OTU_153|OTU_440|OTU_510' rice_SL_Bj_wilcox_v3/result/compare/diff.list.vennd27RtBj_NpRtBj_Dd17RtBj_NpRtBj_Dd10RtBj_NpRtBj_D.xls.xls rice_SL_Bj_wilcox_v3/result/compare/diff.list.vennd27RtBj_NpRtBj_Ed17RtBj_NpRtBj_Ed10RtBj_NpRtBj_E.xls.xls 
    # all in depleted file，如OTU_153来自d10RtBj_NpRtBj_D，发现Np中有小部分不为零，但中位数为零

    # ggplot2绘制SL所有3组共有下调OTU的相对丰度"OTU_23","OTU_31","OTU_18","OTU_125","OTU_57","OTU_43","OTU_443","OTU_37","OTU_1436"
    tail -n12 rice_SL_Bj_wilcox_v1/result/compare/diff.list.vennd27RtBj_NpRtBj_Dd17RtBj_NpRtBj_Dd10RtBj_NpRtBj_D.xls.xls|cut -f 1|awk '{print "\""$$1"\""}'|tr "\n" ","|sed 's/,$$//'
	alpha_boxplot.sh -i `pwd`/result/otutab.txt -m '"OTU_23","OTU_31","OTU_18","OTU_125","OTU_57","OTU_43","OTU_443","OTU_37","OTU_1436"' \
        -d `pwd`/doc/"SL"/design.txt -A groupID -B '"d27RtBj","d17RtBj","d10RtBj","d3AHLRtBj","d3NpRtBj","d3RtBj","d14AHLRtBj","d14RtBj","d53RtBj","NpRtBj","d27RtHn","d17RtHn","d10RtHn","d3RtHn","d14RtHn","d53RtHn","NpRtHn"' \
        -o `pwd`/result/otu_boxplot/ -h 3 -w 10 -t TRUE -n TRUE

## 4.12 科水平差异展示台球图

    # 脚本来自timecourse的plot_pie_DA_Bplylum.r
    # 获得Top9 phylum+class，手工制作前10类Top6Phylum+Proteobacteria-4class
    cut -f 1 result/tax/sum_p.txt | tail -n+2 | head -n6 | sort | head -n5 > result/tax/tax_pc.top9
    cut -f 1 result/tax/sum_pc.txt | grep 'proteobacteria' | sort | grep -v 'Epsilon' >> result/tax/tax_pc.top9
    # 以rice_SL_Bj_wilcox_v3/result/compare/d27RtBj-NpRtBj_sig.txt rice_SL_Bj_wilcox_v3/result/compare_f/d27RtBj-NpRtBj_sig.txt为例
    awk '{if ($3!="Proteobacteria") {print $0"\t"$3} else {print $0"\t"$4}}' result/taxonomy_8.txt|less -S|sed '1 s/Phylum$/PC/' > result/taxonomy_9.txt
    cut -f 6,9 result/taxonomy_9.txt|tail -n+2|sort|uniq|sed '1 i Family\tPC'|less > result/taxonomy_familyPC.txt

    # family同名存在于多个phylum中，无法一一对应。需要使用完整门、属分类，改为多级合并式组合物种
    wc -l result/tax/sum_?.txt
    wc -l result/tax/count_?.txt
    # 改用完整物种时，各级都增加了2-3倍的数据
    mkdir -p SL/pie_family


## 4.13 GWAS关联可遗传菌

    # GAPIT LN heritable 0.15 —— gapit/gapit_LN_heratibleOTU.R 
    cd ~/rice/integrate16s && mkdir -p gwas && cd gwas
    # 4.6 A不分H/L氮下显著>0.001且可遗传的OTU(Vg>万一，H2>0.15)
    awk '$2>0.0001 && $3>0.15 && $4<0.001' ../heritability/A/heritable_otu.elite | sort -k2,2nr > heritability_H.15P.001.otu #  | wc -l # 245
    # 筛选LN下均值OTU表 ，D:\work\rice\miniCore\180319\result\gene_cor\SNPcorDiversity.xls有重点基因的SNP列表
    cd ~/rice/integrate16s/gwas/gapit/LN_heri/
    ls | grep 'OTU_' | cut -f 3 -d '.' | sort | uniq | wc -l # 一夜仅运行完了13个
    less -S GAPIT.MLM.OTU_11.GWAS.Results.csv # Proteobacteria  Betaproteobacteria      Gallionellales  Gallionellaceae Sideroxydans最高仅为1e-5查看nrt相关基因10m21759092，仅为0.36,FDR=1
    less -S GAPIT.MLM.OTU_29.GWAS.Results.csv #  Actinobacteria  Actinobacteria  Actinomycetales Kineosporiaceae Mobilicoccus    Mobilicoccus_pelagius最高，但D17/D27并不显著
    # OTU 29、331、34、450、99、存在>1e-8的SNP

    # GWAS正对照，显状位点的SNP、相关情况、分组和数据分析 —— fig0.test.Rmd ## GWAS
    cd ~/rice/integrate16s/gwas/gapit
    # 提SNP对应基因型
	cat /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_*.hmp.txt > ~/rice/integrate16s/gwas/minicroe.hmp.txt
	snp=`tail -n+2 LN_heri/GAPIT.MLM.OTU_29.GWAS.Results.csv|cut -f 1 -d ','|head -n 10|tr '\n' '|'|sed 's/|$//'`
	grep -P $snp ~/rice/integrate16s/gwas/minicroe.hmp.txt > /tmp/temp1
	cat <(head -n1 /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_1.hmp.txt) /tmp/temp1 > LN_heri/OTU_29_good_snp.geno # 添加标题
    # R代码统计，发现数据分布正常，结果有很多低频allele，以后要注释过滤掉。
    
    # GWAS分析所有的measureOTU, 0.1%中位数于H/L氮和品种
    # 由于gapit需要大量时间，1-2小时一个，需要并行，分十批，硬盘没空间的，重跑k1为k2
    mkdir -p k2
    for i in `seq 1 8`; do
        echo $i
        grep _$i ~/rice/integrate16s/result/minicore/otutab0.1.id > k2/otuid.txt
        Rscript ~/rice/integrate16s/script/gapit.R &
        sleep 300
    done

    # 检索基因SNP是否显著
    # 检索研究基因功能位点、所有功能基因SNP位点，再分析plink结果
    # 基因SNP列表文件 /mnt/bai/yongxin/rice/miniCore/180319/result/gene_cor/all_snp.list, good_snp.list
    list=/mnt/bai/yongxin/rice/miniCore/180319/result/gene_cor/good_snp.list
    cut -f 1 /mnt/bai/yongxin/rice/miniCore/180319/result/gene_cor/good_snp.list|tr '\n' '|'|sed 's/|/,|/g'|sed 's/|$//'
    # 匹配一个，且至少小于0.001
    # OTU替换为变量，并写为文件第一列
    mkdir -p snp1
    for i in `ls k1/GAPIT.MLM.*.GWAS.Results.csv|cut -f 3 -d '.'|sort|uniq`; do
    #i=OTU_9
    grep -P '1m39596400,|1m39596699,|1m39596713,|1m39596808,|1m39600096,|1m39602780,|1m39605678,|1m40659706,|1m40659949,|1m40663740,|2m28751006,|2m28751738,|2m28752575,|4m27567935,|4m27568586,|4m27568855,|5m5941257,|5m5941396,|5m5943888,|5m5944557,|5m5944851,|5m5944944,|8m3183208,|9m16411594,|9m16414341,|9m16414429,|9m16415032,|9m16415203,|9m16415254,|10m21759092,|10m21761740,|10m21761997,|11m22230369,|11m22230384,|11m22230417,|11m22230594,|11m22230660,|11m22230777,|11m22230789,' k1/GAPIT.MLM.${i}.GWAS.Results.csv | sed 's/,/\t/g' |awk '$4<0.001' | awk -v i=$i '{print i"\t"$0}' > snp1/${i}
    done
    # 删除零字节文件
    find ./snp1 -name "*" -type f -size 0c | xargs -n 1 rm -f
    cat ./snp1/OTU_* > ./snp1/all.txt
    # OTU_35有1e-6次方，而且SNP为D27对应Streptomyces，查看OTU_35 manhattan图，发现很多过e8的点，但没有成线的特征，e-6不是很显著
    # 物种注释和基因注释，只保留OTUID, SNPID, P值
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print $1,$2,$5,a[$2]}' /mnt/bai/yongxin/rice/miniCore/180319/result/gene_cor/alpha/pvalue.txt3 ./snp1/all.txt > ./snp1/all1.txt
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' ~/rice/integrate16s/result/taxonomy_2.txt ./snp1/all1.txt > ./snp1/all2.txt

    # 检索所有相关SNP
    list=/mnt/bai/yongxin/rice/miniCore/180319/result/gene_cor/all_snp.list
    cut -f 1 $list|tr '\n' '|'|sed 's/|/,|/g'|sed 's/|$//'
    # 匹配一个，且至少小于0.001
    # OTU替换为变量，并写为文件第一列
    mkdir -p snp1_a
    for i in `ls k1/GAPIT.MLM.*.GWAS.Results.csv|cut -f 3 -d '.'|sort|uniq`; do
    #i=OTU_9
    grep -P '1m31230247,|1m31230250,|1m31230333,|1m31230363,|1m31230368,|1m31230401,|1m31230445,|1m31230618,|1m31230624,|1m31230694,|1m31230761,|1m31230960,|1m31231213,|1m31231368,|1m31231377,|1m31231381,|1m31231417,|1m31231496,|1m31231524,|1m31232165,|1m31232397,|1m31232399,|1m31232412,|1m31232428,|1m31232527,|1m31232610,|1m38370381,|1m38370402,|1m38370498,|1m38370698,|1m38370715,|1m38370863,|1m38370966,|1m38371087,|1m38371126,|1m38371138,|1m38371194,|1m38371230,|1m38371375,|1m38371546,|1m38371599,|1m38371676,|1m38371748,|1m38371783,|1m38372432,|1m38372434,|1m38372441,|1m38372504,|1m38372568,|1m38372621,|1m38372803,|1m38372860,|1m38372948,|1m38372999,|1m38373007,|1m38373023,|1m38373102,|1m38373112,|1m38373113,|1m38373121,|1m38373136,|1m38373146,|1m38373455,|1m38373571,|1m38373622,|1m38373695,|1m38373719,|1m38373732,|1m38373743,|1m38373757,|1m38373758,|1m38373772,|1m38373792,|1m38373816,|1m38373878,|1m38373898,|1m38374214,|1m38374310,|1m38374330,|1m38374400,|1m38374589,|1m38374725,|1m38374790,|1m38374820,|1m38374921,|1m38375013,|1m38375048,|1m38376962,|1m38376963,|1m38377142,|1m38377155,|1m38377178,|1m38377227,|1m38377240,|1m38377288,|1m38377325,|1m38377337,|1m38377382,|1m38377428,|1m38377473,|1m38377718,|1m38377868,|1m38378231,|1m38378262,|1m38378302,|1m38378328,|1m38378332,|1m38378343,|1m38378368,|1m38378420,|1m38378470,|1m38378659,|1m38378668,|1m38378687,|1m38378750,|1m38378960,|1m38379068,|1m38379477,|1m38380261,|1m38380318,|1m38380369,|1m38380410,|1m38380563,|1m38380634,|1m38380836,|1m38380863,|1m38380873,|1m38380875,|1m38380886,|1m38380927,|1m38380929,|1m38380946,|1m38380981,|1m38380986,|1m38380997,|1m38381039,|1m38381175,|1m38381188,|1m38381358,|1m38381393,|1m38381483,|1m38381658,|1m38381676,|1m38381681,|1m38381706,|1m38381767,|1m38381785,|1m38381800,|1m38381934,|1m38381935,|1m38381991,|1m38382080,|1m38382084,|1m38382111,|1m38382137,|1m38382194,|1m38382292,|1m38383221,|1m38389314,|1m38389919,|1m39595895,|1m39595988,|1m39596071,|1m39596115,|1m39596269,|1m39596400,|1m39596699,|1m39596713,|1m39596808,|1m39596921,|1m39597779,|1m39598158,|1m39598992,|1m39600096,|1m39601505,|1m39601550,|1m39601906,|1m39602780,|1m39603303,|1m39603503,|1m39605678,|1m39606189,|1m39606807,|1m39606845,|1m39606849,|1m39606902,|1m39606917,|1m39606926,|1m39606946,|1m39606985,|1m39607030,|1m39607054,|1m39607082,|1m39607122,|1m39607126,|1m39607199,|1m39607264,|1m39607492,|1m39607495,|1m39607499,|1m39607524,|1m39607557,|1m39607565,|1m39607575,|1m39607608,|1m39607614,|1m39607655,|1m39607714,|1m39607719,|1m39607851,|1m39607941,|1m39607999,|1m39608034,|1m39608059,|1m39608078,|1m39608079,|1m39608086,|1m39608191,|1m39608214,|1m39608293,|1m39608372,|1m39608499,|1m39608565,|1m39609498,|1m39609782,|1m39609928,|1m39609967,|1m39609971,|1m39610271,|1m39610508,|1m39611111,|1m40649334,|1m40649417,|1m40649451,|1m40649606,|1m40649654,|1m40649657,|1m40649666,|1m40649668,|1m40649672,|1m40649704,|1m40649747,|1m40649799,|1m40649810,|1m40650110,|1m40650193,|1m40650195,|1m40650492,|1m40650547,|1m40650595,|1m40650671,|1m40650672,|1m40650738,|1m40651054,|1m40651194,|1m40651416,|1m40651511,|1m40651791,|1m40651873,|1m40651955,|1m40652366,|1m40652575,|1m40652711,|1m40652778,|1m40652784,|1m40652916,|1m40653029,|1m40653145,|1m40653187,|1m40653760,|1m40654191,|1m40654878,|1m40655838,|1m40657001,|1m40658900,|1m40659706,|1m40659949,|1m40663740,|1m40669137,|2m28749717,|2m28751006,|2m28751302,|2m28751738,|2m28752575,|2m28757766,|3m5422399,|3m5423417,|3m5423420,|3m5423574,|3m5426011,|3m5426371,|3m5426494,|3m5426662,|3m5426713,|3m5426823,|3m5426883,|3m5426909,|3m5426919,|3m5426954,|3m5427129,|3m5427130,|3m5427182,|3m5427226,|3m5427257,|3m5427287,|3m5427338,|3m5427346,|3m5427362,|3m5427376,|3m5427386,|3m5427434,|3m5427435,|3m5427488,|3m5427500,|3m5427553,|3m5427590,|3m5427614,|3m5427658,|3m5427662,|3m5427693,|3m5427739,|3m5427746,|3m5428618,|3m5428669,|3m5428670,|3m5428718,|4m25496921,|4m25499289,|4m25499381,|4m25499449,|4m25499615,|4m25507542,|4m25507833,|4m25507839,|4m25507877,|4m25507921,|4m25508042,|4m25508289,|4m25508297,|4m25508364,|4m25508409,|4m25508484,|4m25508500,|4m25508511,|4m25508556,|4m25508628,|4m25508694,|4m25508773,|4m25508834,|4m25508866,|4m25508886,|4m25508899,|4m25508949,|4m25509042,|4m25509068,|4m25509087,|4m25509114,|4m25509188,|4m25509189,|4m25509197,|4m25509209,|4m25509230,|4m25509585,|4m25510254,|4m27567935,|4m27568586,|4m27568602,|4m27568838,|4m27568855,|4m27569961,|4m27569991,|4m27570368,|5m5938960,|5m5939076,|5m5939344,|5m5939353,|5m5939420,|5m5939477,|5m5939644,|5m5939650,|5m5939739,|5m5939741,|5m5939770,|5m5939801,|5m5939905,|5m5939912,|5m5939916,|5m5939954,|5m5939965,|5m5940240,|5m5940458,|5m5940547,|5m5940684,|5m5941101,|5m5941257,|5m5941396,|5m5941964,|5m5943888,|5m5944557,|5m5944851,|5m5944944,|5m5946865,|5m5946869,|5m5946922,|5m5946947,|5m19871018,|5m19871081,|5m19884756,|5m19884773,|5m19884780,|5m19884790,|5m19884794,|5m19884798,|5m19884803,|5m19884846,|5m19885030,|5m19885032,|6m577626,|8m3182631,|8m3182703,|8m3182725,|8m3183208,|8m3184252,|8m3184413,|8m3184700,|8m3184729,|8m3184950,|8m3185444,|8m3185825,|8m3185896,|8m3186716,|8m3186889,|8m3187203,|8m3187354,|8m3187376,|8m3187406,|8m3187411,|8m3187421,|8m3187430,|8m3187491,|8m3187504,|8m3187516,|8m3187594,|8m3187602,|8m3187641,|8m3187669,|8m3187674,|8m3187685,|8m3187698,|8m3187699,|8m3187723,|8m3187740,|8m3187747,|8m3187757,|8m3187784,|8m3187789,|8m3187796,|8m3187827,|8m3187866,|8m3187867,|8m3187942,|8m3187976,|8m3188073,|8m3188097,|8m3188306,|8m3188396,|8m3188405,|8m3188419,|8m3188441,|8m3188473,|8m3188482,|8m3188493,|8m3188503,|8m3188523,|8m3188527,|8m3188541,|8m3188588,|8m3188602,|8m3188619,|8m3188679,|8m3189468,|8m3189495,|8m3189556,|8m3189724,|8m3189782,|8m3189860,|8m3189971,|8m3190010,|8m3190301,|8m3190328,|8m3190360,|8m3190550,|8m3190556,|8m3190559,|8m3190700,|8m3190708,|8m3190710,|8m3190731,|8m3190754,|8m3190757,|8m3190777,|8m3190782,|8m3190785,|8m3190799,|8m3190839,|8m3190847,|8m3190855,|8m3190858,|8m3190872,|8m3190949,|8m3190977,|8m3191013,|8m3191022,|8m3191028,|8m3191032,|8m3191040,|8m3191077,|8m3191083,|8m3191087,|8m3191121,|8m3191172,|8m3191179,|8m3191190,|8m3191216,|8m3191221,|8m3191240,|8m3191259,|8m3191290,|8m3191312,|8m3191349,|8m3191358,|8m3191465,|8m3191488,|8m3191549,|8m3191617,|8m3191625,|8m3191661,|8m3191669,|8m3191778,|8m3191874,|8m3191877,|8m3191977,|8m3192006,|8m3192009,|8m3192017,|8m3192357,|9m16406496,|9m16406528,|9m16406545,|9m16406555,|9m16406563,|9m16406577,|9m16406664,|9m16406690,|9m16406718,|9m16406753,|9m16406779,|9m16406821,|9m16406829,|9m16406839,|9m16406859,|9m16406883,|9m16406913,|9m16406922,|9m16406940,|9m16411594,|9m16413000,|9m16414341,|9m16414399,|9m16414429,|9m16414735,|9m16414739,|9m16415032,|9m16415104,|9m16415203,|9m16415254,|9m16415255,|9m16415391,|9m16415724,|9m16415831,|9m16415859,|9m16416012,|9m16416126,|9m16416275,|9m16416891,|9m16416892,|10m21757901,|10m21759092,|10m21761740,|10m21761946,|10m21761997,|10m21762382,|10m21762724,|10m21762787,|10m21762831,|10m21762917,|10m21763222,|10m21763261,|10m21763301,|10m21763315,|10m21763345,|10m21763361,|10m21763372,|10m21763380,|10m21763384,|10m21763396,|10m21763411,|10m21763413,|10m21763418,|10m21763450,|10m21763663,|10m21763687,|10m21763762,|10m21763782,|10m21763783,|10m21763799,|10m21763874,|10m21764084,|10m21764124,|10m21764165,|10m21765227,|10m21765231,|10m21765273,|10m21765292,|10m21765339,|10m21765387,|10m21765393,|10m21765413,|10m21765495,|10m21765591,|10m21765600,|10m21765611,|10m21765683,|10m21765691,|10m21766422,|10m21766437,|10m21766633,|10m21766668,|10m21766741,|10m21766780,|10m21766853,|10m21766868,|10m21766884,|10m21766893,|10m21766918,|10m21767257,|10m21767288,|10m21767291,|10m21767602,|10m21767725,|10m21767738,|10m21768009,|10m21768017,|10m21768028,|10m21768048,|10m21768062,|10m21768066,|10m21768091,|10m21768094,|10m21768131,|10m21768140,|10m21768247,|10m21768257,|10m21768452,|10m21768562,|10m21768660,|10m21768751,|10m21768768,|10m21768838,|10m21768905,|10m21768950,|10m21769062,|10m21769093,|10m21769101,|10m21769105,|10m21769119,|11m194222,|11m195792,|11m195980,|11m22230369,|11m22230384,|11m22230417,|11m22230594,|11m22230660,|11m22230777,|11m22230789,|11m22230893,|11m22230904,|11m22230908,' k1/GAPIT.MLM.${i}.GWAS.Results.csv | sed 's/,/\t/g' |awk '$4<0.001' | awk -v i=$i '{print i"\t"$0}' > snp1_a/${i} &
    done
    # 删除零字节文件
    find ./snp1_a -name "*" -type f -size 0c | xargs -n 1 rm -f
    cat ./snp1_a/OTU_* > ./snp1_a/all.txt
    # OTU_35有1e-6次方，而且SNP为D27对应Streptomyces，查看OTU_35 manhattan图，发现很多过e8的点，但没有成线的特征，e-6不是很显著
    # 物种注释和基因注释，只保留OTUID, SNPID, P值
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4} NR>FNR{print $1,$2,$5,a[$2]}' /mnt/bai/yongxin/rice/miniCore/180319/result/gene_cor/all_snp.list ./snp1_a/all.txt > ./snp1_a/all1.txt
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$3]=$1} NR>FNR{print $0,a[$4]}' ~/rice/miniCore/180319/doc/gwas_candidate_gene.txt ./snp1_a/all1.txt > ./snp1_a/all2.txt
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' ~/rice/integrate16s/result/taxonomy_2.txt ./snp1_a/all2.txt > ./snp1_a/all3.txt

    # 计算整理alpha, beta多样性，取显著SNP的基因，作GO分析
    # 计算LN下alpha, beta与SNP的关联 gapit_LN_diversity.R
    cd ~/rice/integrate16s/gwas/gapit/LN_diversity
    # 显示列表
    ls man/GAPIT.MLM.*|cut -f 3 -d '.'
    # richness为例, chao1, shannon_3
    i=chao1
    for i in `ls man/GAPIT.MLM.*|cut -f 3 -d '.'`; do
    # 筛选显著SNP
    head -n 10000 GAPIT.MLM.${i}.GWAS.Results.csv | sed 's/,/\t/g' | awk '$4<0.001' | cut -f 1,4 > ${i}_1_sig.txt
    wc -l ${i}_1_sig.txt
    # 基因注释
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' ~/rice/miniCore/180319/gemma/snp.anno ${i}_1_sig.txt > ${i}_1_anno.txt 
    # 提取第6列OS ID，基因间替换为回车，并去冗余
    # cut -f 6 ${i}_1_anno.txt|sed 's/-/\n/g' | sort|uniq | sed 's/OS/Os/g;s/G/g/g' >${i}_OS.txt # AgriGO的Oryza sative无法识别
    # 提取第8列LOC
    cut -f 8 ${i}_1_anno.txt|sort|uniq|grep 'LOC' >${i}_LOC.txt # AgriGO的Oryza sative无法识别
    done



    
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
    awk '$3>0.1 && $7>=100' result/anno/otu.8_10heritable | cut -f 1 | awk '{print $0" label #0000ff bold 2"}'
    # 标注可遗传OTUs枝标为红色
    awk '$3>0.1 && $9>=0.3 ' result/anno/otu.8_10heritable | tail -n+2 | cut -f 1 | awk '{print $0"|"$0" clade #ff0000 normal 1"}'

    # 选择核心菌95%序列，与greeengene比对建树，并用iTOL注释物种门水平颜色、核心菌加粗和可遗传标为红色
    # 核心菌95%序列
    awk '$7>=95' result/anno/otu.8_10heritable | cut -f 1 | grep -v 'OTUID'> temp/otu_k1.id
    wc -l temp/otu_k1.id # 276
    usearch10 -fastx_getseqs result/otu.fa -labels temp/otu_k1.id -fastaout script/fig1/1d.otu_k1.fa
    grep -c '>' script/fig1/1d.otu_k1.fa
    # Multiply alignment
    #clustalo -i script/fig1/1d.otu_k1.fa -o temp/1d.otu_k1_align.fa --seqtype=DNA --full --force --threads=9
    align_seqs.py -i script/fig1/1d.otu_k1.fa -t /mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta -o temp/aligned/
	filter_alignment.py -i temp/aligned/1d.otu_k1_aligned.fasta -o temp/aligned/
	make_phylogeny.py -i temp/aligned/1d.otu_k1_aligned_pfiltered.fasta -o script/fig1/1d.otu_core95.tree # generate tree by FastTree
    sed -i "s/'//g" script/fig1/1d.otu_core95.tree
    # 在iTOL中展示树
    # 注释物种四大门水平颜色
    # Phylum = c("Proteobacteria", "Actinobacteria", "Bacteroidetes", "Firmicutes", "Other"), Color = c("#85F29B", "#F58D8D", "#F7C875", "#91DBF6", "#AAAAAA" )
    anno=script/fig1/1d.otu_core95.txt
    cp ~/github/Amplicon/16Sv2/doc/iTOL.txt $anno
    Rscript ~/github/Amplicon/16Sv2/script/tree_color_iTOL.R
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3" "$2} NR>FNR{print $0" range "a[$1]}' result/taxonomy_8_color.txt temp/otu_k1.id >> $anno
    # 标注核心OTUs $7>100% 蓝色加粗2号
    awk '$3>0.1 && $7>=100 ' result/anno/otu.8_10heritable | cut -f 1 | awk '{print $0" label #0000ff bold 2"}' >> $anno
    # 标注可遗传OTUs枝标为红色
    awk '$3>0.1 && $9>=0.3 ' result/anno/otu.8_10heritable | tail -n+2 | cut -f 1 | awk '{print $0"|"$0" clade #ff0000 normal 1"}' >> $anno

    # 基于物种和丰度均值注释
mkdir -p tree
    # 筛选高丰度菌对应物种注释
    sed '1 i OTUID' temp/otu_k1.id > tree/otutab_high.id
    cp script/fig1/1d.otu_core95.tree tree/otu_core95.nwk
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/taxonomy_8.txt tree/otutab_high.id > tree/otutab_high.tax
    # 获得均值和最大值
    # Rscript otu_mean.R
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/taxonomy_8.txt tree/otutab_high.id > tree/otutab_high.tax
    awk 'NR==FNR{a[$1]=$3} NR>FNR{print $0"\t"a[$1]}' result/anno/otu.8_10heritable tree/otutab_high.tax > tree/annotation.txt
    # 注释
    Rscript ~/ehbio/train/04Amplicon/1809/26Evolution/table2itol/table2itol.R -a -c double -D tree/ -i OTUID -l Genus -t %s -w 0.5 tree/annotation.txt
    Rscript ~/ehbio/train/04Amplicon/1809/26Evolution/table2itol/table2itol.R -a -d -c none -D tree/ -b Phylum -i OTUID -l Genus -t %s -w 0.5 tree/annotation.txt
    # 用门背景色、属 目标签、丰度

##  图2. 微生物组与基因型 fig2.diversity.genetics.rmd

### 2A. Beta多样性

    # 按四大亚种绘制alpha, beta多样性
    beta_pcoa.sh -i `pwd`/temp/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
        -d `pwd`/doc/xiangeng/design.txt  -A groupID -B `cat doc/xiangeng/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'` -E TRUE \
        -c `pwd`/doc/xiangeng/compare.txt \
        -o `pwd`/result/beta/ -h 3 -w 5
    alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
        -d `pwd`/doc/"xiangeng"/design.txt  -A groupID -B '"LIND","LTEJ","LAUS","LTRJ"' \
        -o `pwd`/result/alpha/ -h 3 -w 5    
