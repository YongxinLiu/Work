<!-- TOC -->

- [1. 处理序列 Processing sequences](#1-处理序列-processing-sequences)
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
- [3. 高级分析](#3-高级分析)
- [4. 个性分析](#4-个性分析)
    - [4.1. 分蘖与菌相关性](#41-分蘖与菌相关性)

<!-- /TOC -->

# 1. 处理序列 Processing sequences

	# 0. 准备工作 Preparation

	## 0.1 准备流程配置文件

	# 设置工作目录
	wd=ath/integrate16s
	# 创建环境代码见~/github/Work/initial_project.sh
	

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	## 1.3. 准备原始数据

	# Prepare raw data
	# 数据来自三个课题：miniCore + timecourse + nrt1.1b + nrt1.1a + SL + epi + Gprotein + CTK + SD1
	
	# 合并实验设计
	# 合并三个项目的实验设计，检查样品是否有重名
	cp ~/ath/jt.terpene.16S/batch4_unoise/doc/design.txt doc/design_2.5.txt
	cp ~/ath/jt.HuangAC/batch3all/doc/design.txt doc/design_3.txt

	# 统计实验样品行和唯一行，确实样品名唯一
	cat doc/design_* | grep -v 'SampleID'|wc -l
	cat doc/design_* | grep -v 'SampleID'|cut -f 1| sort|uniq|wc -l 
	# 如果不一致，检查显示重名的方法
	cat doc/design_* | grep -v 'SampleID'|cut -f 1| sort|uniq -d
	# 检查土壤有异常，观察原因
	grep 'Soil' doc/design_2.5.txt
	grep 'Soil' doc/design_3.txt
	# 2.5中样本和数据都小，修改2.5中的土为Soil25r
	sed -i 's/^Soil/Soil25r/' doc/design_2.5.txt
	# Double check，序列也需要替换，再上面检查
	grep 'Soil' doc/design_2.5.txt

	# 统计各实验来源样本数量和总量
	wc -l doc/design*
	# 2 projecct include 915 samples (330, 585)

	# 合并实验设计，前7列共有，只保留前7列
	cat <(head -n1 doc/design_3.txt) <(cat doc/design_* | grep -v 'SampleID') | cut -f 1-11 > doc/design.txt

	# 原始数据合并
	cat <(sed 's/>Soil/>Soil25r/' ~/ath/jt.terpene.16S/batch4_unoise/temp/seqs_usearch.fa) \
	~/ath/jt.HuangAC/batch3all/temp/seqs_usearch.fa | cut -f 1 -d ';' | sed 's/_/./g' > temp/filtered.fa
	# 从1.6 fa_unqiue 开始
	

## 1.1. 按实验设计拆分lane为文库

	# Split lane into libraries
	# lane文件一般为seq/lane_1/2.fq.gz
	# lane文库信息doc/library.txt：至少包括编号、Index和样品数量三列和标题
	# head -n3 doc/library.txt
	#LibraryID	IndexRC	Samples
	#L1	CTCAGA	60
	
	# 按library.txt拆分lane为library
	# make lane_split


## 1.2. 按实验设计拆分文库为样品


	# 拆分样品
	head -n3 doc/L01.txt
	# 按L1/2/3...txt拆分library为samples
	make library_split
	make library_split_stat

## 1.3. 样品双端合并、重命名、合并为单一文件

	# Merge paired reads, renames and merge all samples
	# 样品双端合并、重命名、合并为单一文件, 注意fastq为33格式，64位采用fastp转换
	make sample_merge
	make sample_merge_stat


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

	# Remove redundancy, get unique reads
	make fa_unqiue
	# 143,053,957序列，>30为188654，改为143

## 1.7. 挑选OTU

	# Pick OTUs
	# unoise3速度比cluster_otus慢上百倍
	# 修改culster_otus为unoise3，143阈值下97%下为5000个OTUs，unoise3下1.3万
	cp -r result result_97
	rm otu_pick
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

	## 比较二半萜-第二批，实验设计位于doc/2.5 目录中
	sub=2.5
	mkdir -p doc/${sub}
	pwd=~/ath/jt.terpene.16S/batch4_unoise/doc
	# 编辑2、3批组名为groupID，添加b2/b3区分
	cp ${pwd}/design.txt doc/${sub}/design.txt
	sed -i 's/Soil/Soil25r/' doc/${sub}/design.txt
	# 第二批数据添加b2
	cp ${pwd}/group_compare.txt doc/${sub}/compare.txt
	sed -i 's/\t/b2\t/g;s/$/b2/g' doc/${sub}/compare.txt
	cat -A doc/${sub}/compare.txt
	cp ${pwd}/group_venn.txt doc/${sub}/venn.txt
	sed -i 's/vs/b2vs/g;s/_/b2_/g' doc/${sub}/venn.txt
	cat -A doc/${sub}/venn.txt

	## 比较二半萜-第三批，实验设计位于doc/2.5.3 目录中
	sub=2.5.3
	mkdir -p doc/${sub}
	pwd=~/ath/jt.terpene.16S/batch4_unoise/doc
	cp doc/2.5/design.txt doc/${sub}/design.txt
	# 获得比较和维恩文件，并添加b3与实验设计对应
	cp ${pwd}/b3/group_compare.txt doc/${sub}/compare.txt
	sed -i 's/\t/b3\t/g;s/$/b3/g' doc/${sub}/compare.txt
	cat -A doc/${sub}/compare.txt
	cp ${pwd}/b3/group_venn.txt doc/${sub}/venn.txt
	sed -i 's/vs/b3vs/g;s/_/b3_/g;s/enriched/E/g;s/depleted/D/g;s/vs/_/g' doc/${sub}/venn.txt
	cat -A doc/${sub}/venn.txt
	# 获得比较组中的组名
	cat doc/${sub}/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","

	## 2018/8/22 比较二半萜-第三批+Soil，实验设计位于doc/2.5.3soil 目录中
	sub=2.5.3soil
	mkdir -p doc/${sub}
	cp doc/2.5.3/* doc/${sub}
    # 添加WT4b3 vs Soil1b3
    head -n1 ~/ath/jt.terpene.16S/batch4_unoise/doc/b3soil/group_compare.txt|sed 's/\t/b3\t/g;s/$/b3/g' >> doc/${sub}/compare.txt
	cat doc/${sub}/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
    rm alpha_boxplot
    make plot_venn # DA otu
    make DA_compare_tax # DA taxonomy
    make rmd # report

	## 比较三萜-第三批，实验设计位于doc/3.1 目录中
	sub=3.1
	mkdir -p doc/${sub}
	pwd=~/ath/jt.HuangAC/batch3all/doc/b3_4
	# cp doc/design.txt doc/${sub}/design.txt
	cp ~/ath/jt.HuangAC/batch3all/doc/design.txt doc/${sub}/design.txt
    # 获得比较和维恩文件，并添加b3与实验设计对应
	cp ${pwd}/group_compare.txt doc/${sub}/compare.txt
	cp ${pwd}/group_venn.txt doc/${sub}/venn.txt
	sed -i 's/enriched/E/g;s/depleted/D/g;s/vs/_/g' doc/${sub}/venn.txt
	cat -A doc/${sub}/venn.txt
	# 获得比较组中的组名
	cat doc/${sub}/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
    # modify sub, then report
    rm alpha_boxplot
    make plot_venn # DA otu
    make DA_compare_tax # DA taxonomy
    make rmd # report



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
	make plot_volcano
	make plot_heatmap
	make plot_manhattan

# 2.11 plot_venn 维恩图
	
	make plot_venn



# 3. 高级分析


## 3.9 培养菌注释

    # 培养菌注释，采用ath root的菌库，COTU，目前只注释plot_veen的结果
    make culture



# 4. 个性化分析

    # 数据上传和可重复绘图

## 4.1 发表前样本合并和实验设计上传

# 根据结果PCoA/CPCoA/热图聚类筛选样品，删除异常点，保存于doc/3.1/outlier.txt
# Colr2热图异常，但CPCoA正常；ACT2KOr1热图和CPCoA离群，但与ACT2CR组拼接；
# 剔除点注释design
for i in `cut -f 1 doc/3.1/outlier.txt`; do sed -i "s/^${i}\t/#${i}\t/" doc/3.1/design.txt;done

## Fig3. 可重复计算~/ath/integrate16s/3T/fig3

    # 参考模板 D:\work\ath\jt.HuangAC\manuscript\Figure 3
    # 3A. CPCoA 五组()绘制CPCoA, 代重复；
    # 3B调整为附图；新增Thas/ACT2共有韦恩图4个
    # 3C 有重复的ACT2、Thas取交集,再做曼哈顿图；
    # 3D科水平维恩，人工制作，进附录？
    # 3E OTU水平4组维恩，与土壤相比富集的OTU

    # 7个基因型："b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO"
    # 精选5个基因型："b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO"
    # 精选5个基因型+Soil："b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO","b3BS"

### 3A. CPCoA 五组()绘制CPCoA, 代重复；
    
    # 生成图A脚本供修改
    beta_cpcoa.sh -i `pwd`/result/otutab.txt -m '"bray"' \
	-d `pwd`/doc/"3.1"/design.txt  -A groupID -B '"b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO"' -E TRUE \
	-o `pwd`/result/beta/ -h 3 -w 5

### 3B. 调整为附图；新增Thas/ACT2共有韦恩图4个

tax_stackplot.sh -i `pwd`/result/tax/sum_ -m '"p","pc"' -n 10 \
	-d `pwd`/doc/"3.1"/design.txt  -A groupID -B '"b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO","b3BS"' -O FALSE \
	-o `pwd`/result/tax/sum_ -h 3 -w 5

### 3C. 两基因型共有曼哈顿图
	
	plot_manhattan.sh -i result/compare/b3ACT2CR-b3Col_all.txt -Y 10
	plot_manhattan.sh -i result/compare/b3ACT2KO-b3Col_all.txt -Y 10

### 4D. Venny图比较，共有比特有

	# 追加两批共有OTU至result/compare/diff.list文件，再绘制维恩图
	# 制作新共有的两组的上调和下调OTU到新文件，不可追加老文件，防错误；运行完以下内容可马上make rmd保存至网页
	tail -n `grep -A1 'specific_to_others' result/compare/diff.list.vennb3ACT2CR_b3Col_Eb3ACT2KO_b3Col_E.xls.xls|tail -n1` result/compare/diff.list.vennb3ACT2CR_b3Col_Eb3ACT2KO_b3Col_E.xls.xls | awk '{print $1"\tACT2_E"}' > 3T/fig3/mutant_overlap.txt
	tail -n `grep -A1 'specific_to_others' result/compare/diff.list.vennb3ACT2CR_b3Col_Db3ACT2KO_b3Col_D.xls.xls|tail -n1` result/compare/diff.list.vennb3ACT2CR_b3Col_Db3ACT2KO_b3Col_D.xls.xls | awk '{print $1"\tACT2_D"}' >> 3T/fig3/mutant_overlap.txt
	tail -n `grep -A1 'specific_to_others' result/compare/diff.list.vennb3ThasKO1_b3Col_Eb3ThasKO2_b3Col_E.xls.xls|tail -n1` result/compare/diff.list.vennb3ThasKO1_b3Col_Eb3ThasKO2_b3Col_E.xls.xls | awk '{print $1"\tThas_E"}' >> 3T/fig3/mutant_overlap.txt
	tail -n `grep -A1 'specific_to_others' result/compare/diff.list.vennb3ThasKO1_b3Col_Db3ThasKO2_b3Col_D.xls.xls|tail -n1` result/compare/diff.list.vennb3ThasKO1_b3Col_Db3ThasKO2_b3Col_D.xls.xls | awk '{print $1"\tThas_D"}' >> 3T/fig3/mutant_overlap.txt
	# 备份旧文件再追加
	cp result/compare/diff.list result/compare/diff.list.bak
	cat result/compare/diff.list.bak 3T/fig3/mutant_overlap.txt > result/compare/diff.list

	# 修改Venn.txt文件添加比较组
	cp doc/3.1/venn.txt doc/3.1/venn.txt181016
	# batch_venn.pl -i `pwd`/doc/"3.1"/venn.txt -d result/compare/diff.list
	make plot_venn
	make rmd
	# Other downloaded from website: http://210.75.224.110/report/16Sv2/ath_3.1_edgeR_unoisev2

### 3D. 土壤中富集的菌 
	# 编辑doc/3.1/compare_soil.txt 带土，只比较土
	compare_OTU.sh -i `pwd`/result/otutab.txt -c doc/3.1/compare_soilonly.txt -m "edgeR" \
        -p 0.05 -q 0.2 -F 1.2 -t 0.01 \
        -d `pwd`/doc/"3.1"/design.txt  -A groupID -B `cat doc/"3.1"/compare_soil.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'` \
        -o `pwd`/result/compare/ -C groupID
	# 正常7个基因型才916个OTU，加土为2439个
	cp result/compare/b3Col-b3BS_all.txt 3T/fig3/3d.b3Col-b3BS_all.txt # 备份最终结果
	# 注释两类OTU
	cut -f 6 3T/fig3/3d.b3Col-b3BS_all.txt|sort|uniq -c
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4} NR>FNR{print $0,a[$1]}'
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print $0,a[$1]}' 3T/fig3/3d.b3Col-b3BS_all.txt result/compare/diff.list.vennThas_Eb3ThahKO_b3Col_Eb3ThadKO_b3Col_EACT2_E.xls.xls > 3T/fig3/3d.veen.enriched.soil.txt
	# http://210.75.224.110/report/16Sv2/ath_3.1_edgeR_unoisev2/result/compare/diff.list.vennThas_Db3ThahKO_b3Col_Db3ThadKO_b3Col_DACT2_D.xls.xls
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print $0,a[$1]}' 3T/fig3/3d.b3Col-b3BS_all.txt result/compare/diff.list.vennThas_Db3ThahKO_b3Col_Db3ThadKO_b3Col_DACT2_D.xls.xls > 3T/fig3/3d.veen.depleted.soil.txt
	# 提取Overlap ID
	tail -n `grep -A1 'specific_to_others' result/compare/diff.list.vennThas_Eb3ThahKO_b3Col_Eb3ThadKO_b3Col_EACT2_E.xls.xls|tail -n1` result/compare/diff.list.vennThas_Eb3ThahKO_b3Col_Eb3ThadKO_b3Col_EACT2_E.xls.xls | awk '{print $1"\tAll_E"}' > 3T/fig3/all_common.txt
	tail -n `grep -A1 'specific_to_others' result/compare/diff.list.vennThas_Db3ThahKO_b3Col_Db3ThadKO_b3Col_DACT2_D.xls.xls|tail -n1` result/compare/diff.list.vennThas_Db3ThahKO_b3Col_Db3ThadKO_b3Col_DACT2_D.xls.xls | awk '{print $1"\tAll_D"}' >> 3T/fig3/all_common.txt



## 附图15. Alpha多样性
	# 绘制8个基因型的alpha多样性
	alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
        -d `pwd`/doc/"3.1"/design.txt  -A groupID -B '"b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO","b3BS"' \
        -o `pwd`/result/alpha/ -h 5 -w 8
	# 排除土壤，否则基因型间差异过小
	alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
        -d `pwd`/doc/"3.1"/design.txt  -A groupID -B '"b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO"' \
        -o `pwd`/3T/fig3/S15. -h 5 -w 8


## 附图16. 热图绘制不聚类

	# 绘制单个比较组
	plot_heatmap.sh -i result/compare/b3ThahKO-b3Col_sig.txt -o result/compare/b3ThahKO-b3Col -w 7 -h 5 -C FALSE
	# 绘制多个比较组
	awk 'BEGIN{OFS=FS="\t"}{system("plot_heatmap.sh -i result/compare/"$1"-"$2"_sig.txt \
        -o result/compare/"$1"-"$2" -w 7 -h 5 -C FALSE");}' \
        <(grep -v '^$' `pwd`/doc/"3.1"/compare.txt)

## 附数据：上传三萜3批数据

	# 此目录三萜数据来自~/ath/jt.HuangAC/batch3all，它又来自~/ath/jt.HuangAC/batch2、3
	cd ~/ath/jt.HuangAC/batch3
	mkdir -p clean_data/submit
parallel --xapply -j 8 \
"split_libraries_fastq.py -i temp/{1}_barcode/reads.fastq \
 -b temp/{1}_barcode/barcodes.fastq \
 -o temp/{1}_barcode \
 -m doc/{1}.txt --store_demultiplexed_fastq --barcode_type 10
split_fastq_qiime.pl -i temp/{1}_barcode/seqs.fastq -o clean_data/submit/" \
::: `tail -n+2 doc/library.txt | cut -f 1`
# compress
cd clean_data/submit/
pigz *
# 检查ID是否唯一且一致
ls *.gz|cut -f 1 -d '.'|sort|uniq|wc -l
ls *.gz|cut -f 1 -d '.'|sort|uniq>sampleID1
cut -f 1 ../../doc/design.txt |tail -n+2|sort|uniq|wc -l
cut -f 1 ../../doc/design.txt |tail -n+2|sort|uniq>sampleID2
# 显示是否名称有不对应的样本
cat sampleID?|sort|uniq -u