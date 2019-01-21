# 2019/1/14 图片重出文件准备 fig2

wd=fig2
mkdir -p $wd
mkdir -p $wd/script $wd/data 

# 确定的组筛选样本
Rscript script/filter_otutab_by_sample.R # 按样本组筛选

# 图1subspecies Figure 1. Root microbiomes are distinct between indica and japonica varieties
    # 1a品种地理信息
    cp fig1/ST/01.variety_geotable.csv $wd/data/variety_geotable.txt
    # beta多样性矩阵
    cp result/beta/bray_curtis.txt $wd/data/
    cp result/beta/weighted_unifrac.txt $wd/data/
    cp result/beta/unweighted_unifrac.txt $wd/data/
    
    # alpha多样性更新usearch10的norm为11的rare
    usearch11 -otutab_rare result/otutab.txt -sample_size 5000 -output result/otutab_rare.txt 
    usearch10 -alpha_div result/otutab_rare.txt -output result/alpha/index.txt 
    ln result/alpha/index.txt $wd/data/alpha.txt

    # Taxonomy phylum + class of proteobacteria for stackplot
    cp result/tax/sum_pc.txt $wd/data/

# 图2model Figure 2. Random forest model detects bacterial taxa that accurately predict indica and japonica subspeciation.

    # family for randomForest
    for i in p pc c o f g; do
    cp result/tax/sum_${i}.txt $wd/data/
    done
    cp result/taxonomy_8.txt $wd/data/
    # 
    make DA_compare_tax
    cp result/compare_f/HTEJ-HIND_all.txt $wd/data/family_HTEJ-HIND_all.txt

# 图3compare Figure 3. Taxonomical and functional characteristics of differential bacteria between the indica and japonica root microbiome.

make -n -B DA_compare # 12个组，筛选每组中位数大于0.1%，共296个OTU；改为0.05%为521个，也包括更多的丰度
#  A50HnCp7   A50LnCp7   A56HnCp7   A56LnCp7       HIND       HTEJ       LIND 
#  62.46741   67.43844   64.26266   71.35957   41.52890   47.06853   47.05644 
#      LTEJ V3703HnCp6 V3703LnCp6  ZH11HnCp6  ZH11LnCp6 
#  58.25825   83.05020   83.27368   82.23098   82.49534 
#  A50HnCp7   A50LnCp7   A56HnCp7   A56LnCp7       HIND       HTEJ       LIND 
#  67.46607   71.38812   68.93914   74.77070   46.99615   51.61166   52.21140 
#      LTEJ V3703HnCp6 V3703LnCp6  ZH11HnCp6  ZH11LnCp6 
#  61.98529   84.73981   85.32541   84.55106   84.44299 

cp result/compare/database.txt temp/12geno_mean.abu

# 比对IND和TEJ相对于土壤中富集的情况，这296个OTU仅包括了TEJ       IND      Soil 14.343008 26.818806  5.542759%为什么变少了，应该占一半左右的？
compare_OTU_ref.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare_soil.txt -m "wilcox" \
        -p 0.05 -q 0.05 -F 1.2 -t 0.1 \
        -d `pwd`/doc/design.txt -A subspecies -B '"TEJ","IND","Soil"' \
        -o `pwd`/result/compare_ref/ -r result/compare/LIND-LTEJ_all.txt
cat result/compare/summary.txt # 几乎全部为植物富集，极少土壤富集，20%为末知
mkdir -p soil
cp result/compare/database.txt soil/mean_abundance.txt
# result/39culture/otu.txt 添加基因型丰度和土壤enrich与否信息
cp result/39culture/otu.txt result/39culture/otu.txt.bak
paste temp/12geno_mean.abu <(cut -f 6 result/compare/IND-Soil_all.txt) <(cut -f 6 result/compare/IND-Soil_all.txt) | sed '1 s/level\tlevel/IND-soilt\tTEJ-soil/' | less > soil/mean_type.txt
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' soil/mean_type.txt result/39culture/otu.txt.bak | less -S > result/39culture/otu.txt

# 确定选菌是否与新比对的ASV一致
cut -f 1 result/compare/database.txt | tail -n+2 > wet/otuid296.id
usearch10 -fastx_getseqs result/otu.fa -labels wet/otuid296.id -fastaout wet/otuid296.fa
# cp result/otu.fa wet/otuid296.fa
makeblastdb -in wet/otuid296.fa -dbtype nucl
blastn -query wet/Syncom25.fa -db wet/otuid296.fa -out temp/Syncom25.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 9 
# less temp/Syncom25.blastn # 18/24，6个相似并对，Co7无法比对；改为万5为3/24，万4为2/24，万1和3只差土全包括1/24
awk '$3>97' temp/Syncom25.blastn | wc -l
# 这里没有测试成功的Co_07替换为rice_230，相似度为96.53%；换成全部ASV所有菌为100%，土为97.35%

# 确定NRT与IND一至和菌存在于这些菌保
grep 'HTEJ_HIND_D_V3703HnCp6_ZH11HnCp6_D-specific_to_others' -A 1000 |
# 所有IND和NRT富集
cat xiangengASV_wilcox_0.03/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls.xls|cut -f 1|grep 'OTU' > wet/otuid296.id
# 所有IND富集
cat xiangengASV_wilcox_0.01/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_D.xls.xls|cut -f 1|grep 'OTU' > wet/otuid296.id
# 所有TEJ富集
cat xiangengASV_wilcox_0.01/result/compare/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls |cut -f 1|grep 'OTU' > wet/otuid296.id
wc -l wet/otuid296.id # 统计筛选部分数量109，重复运行上部确定一致性



# 图4









	# 快速分析 Quick Start(所需文件准备好)
	make library_split # 样本拆分、
	make fq_qc # 合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
    rm host_rm
    make host_rm # 序列去冗余、去噪、挑选代表序列、去嵌合、去宿主，获得OTU代表序列result/otu.fa
    rm otutab_create
    rm otutab_filter
	make beta_calc # 生成OTU表、过滤、物种注释、建库和多样性统计
	# 清除统计绘图标记(重分析时使用)
	rm -rf alpha_boxplot 
    rm DA_compare
    make DA_compare # 绘制alpha、beta、taxonomy和差异OTU比较
	rm -f plot_volcano # 删除OTU差异比较可化标记
	make plot_manhattan # 绘制差异比较的火山图、热图、曼哈顿图
	make plot_venn # 绘制OTU差异共有/特有维恩图
	make DA_compare_tax # 高分类级差异比较，维恩图绘制，2为reads count负二项分布统计
	make rmd # 生成网页报告，必须依赖的只有alpha, beta, taxonomy

	# 提取脚本
	submit=3T
	make -n -B fq_qc > pipeline.sh # 样本拆分、合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
	make -n -B host_rm >> pipeline.sh # 序列去冗余、去噪、挑选代表序列、去嵌合、去宿主，获得OTU代表序列result/otu.fa
	make -n -B beta_calc >> pipeline.sh # 生成OTU表、过滤、物种注释、建库和多样性统计
	grep -v '#' pipeline.sh > ${submit}/pipeline.sh

# 标准流程 Standard pipeline

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
	
	
## 按实验设计拆分lane为文库

	# Split lane into libraries
	# lane文件一般为seq/lane_1/2.fq.gz
	# lane文库信息doc/library.txt：至少包括编号、Index和样品数量三列和标题
	head -n3 doc/library.txt
	#LibraryID	IndexRC	Samples
	#L1	CTCAGA	60
	
	# 按library.txt拆分lane为library
	make lane_split


## 按实验设计拆分文库为样品

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
 

## 样品双端合并、重命名、合并为单一文件

	# Merge paired reads, renames and merge all samples
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt) <(cat doc/L* |grep -v '#') > doc/design.txt

	# 样品双端合并、重命名、合并为单一文件, 注意fastq为33格式，64位采用fastp转换
	make sample_merge


## 切除引物与标签

	# Cut primers and lables
	# 切除左端标签和引物，右端 引物
	# Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
	make fq_trim


## 质量控制

	# Quality control
	# 过滤序列中预期累计错误率>1%的序列
	make fq_qc


## 序列去冗余
	
	# 从这里开始
	ln ~/medicago/zjj170823/temp/seqs_usearch.fa temp/filtered.fa

	# Remove redundancy, get unique reads
	make fa_unqiue


## 挑选OTU

	# Pick OTUs
	# unoise3速度比cluster_otus慢上百倍
	make otu_pick


## 有参去嵌合体

	# Remove chimiras by silva database
	# 基于SILVA数据库去除
	make chimera_ref


## 去除宿主

	# Remove host
	# 根据SILVA注释去除线粒体、叶绿体、真核生物18S和未知序列(非rRNA)
	make host_rm


## 生成OTU表
	
	# Create OTUs table
	# 默认使用vsearch更快10倍，可选usearch10，线程不可超48
	make otutab_create


## 过滤样本和OTUs

	# OTU table filter samples and OTU
	# 推荐过滤低测序量<5000的样本，筛选大于1RPM的OTU
	make otutab_filter 
	# 此处结果统计只有3441个样品，而实验设计有3898个样品，少了哪些样品种？
	cat <(head -n1 result/otutab.txt|cut -f 2-|tr '\t' '\n') <(tail -n+2 doc/design.txt | cut -f 1) | sort | uniq -u > doc/missing_samples.txt 
	# 缺失 482样，主要是没有时间序列


## 物种注释

	# Assign taxonomy
	# 默认使用RDP trainset快而准，GG太旧，Silva太慢
	# 推荐阈值为0.6保证注释更完整
	make tax_assign


## 物种统计
	
	# Taxonomy summary
	# 必须所有物种有注释，否则有可能报错
	make tax_sum


## 多序列比对和进化树
	
	# Multiply alignment and make_phylogeny
	# usearch10/culsterO结果不同可能影响多样性分析(usearch unifrac结果更可信)
	# 进化树，用于树图和多样性分析
	make tree_make

## Alpha多样性指数计算
	
	# Calculate alpha diversity index
	# alpha指数计算结果为 result/alpha/index.txt
	# 稀释梯度结果位于 result/alpha/rare.txt
	make alpha_calc

## Beta多样性距离矩阵计算
	
	# Beta diversity tree and distance matrix
	# 最好用usearch，结果unifrac分类更好；clustero+fastree结果PCoA较差
	make beta_calc
	# ---Fatal error--- ../calcdistmxu.cpp(32) assert failed: QueryUniqueWordCount > 0 致信作者; 改用qiime1

## 有参考构建OTU表

	# Reference based OTU table
	# otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
	make otutab_gg



# 统计绘图 Statistics and plot

## Alpha多样性指数箱线图
	
	# Alpha index in boxplot
	make alpha_boxplot

## Alpha丰富度稀释曲线
	
	# Alpha rarefracation curve
	make alpha_rare

## 主坐标轴分析距离矩阵
	
	# PCoA of distance matrix
	make beta_pcoa

## 限制性主坐标轴分析

	# Constrained PCoA / CCA of bray distance matrix
	# OTU表基于bray距离和CCA，至少3个组 
	make beta_cpcoa

## 样品和组各级分类学堆叠柱状图

	# Stackplot showing taxonomy in each level
	make tax_stackplot

## 组间差异比较 
	
	# Group compareing by edgeR or wilcox
	# 可选负二项分布，或wilcoxon秩和检验
	make DA_compare


# 高级分析 Advanced analysis

## 元素循环预测

    make faprotax_calc
    ## 整理faprotax中菌的功能列表
	grep -v '^#' /mnt/bai/yongxin/software/FAPROTAX_1.1/FAPROTAX.txt|sed 's/\*/_/g' > culture/faprotax.tax


## 添加可培养菌

	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $0,a[$1]}' temp/otutab.mean temp/culture_otu.blastn | cut -f 1-4,14 > result/41culture/otu.txt
	echo -ne "OTUs > 97% abundance :\t" >> result/41culture/summary.txt
	awk '$$3>=97 && $$4>=99' result/41culture/otu.txt | awk '{a=a+$$5} END {print a}' >> result/41culture/summary.txt
	cat result/41culture/summary.txt


## 查看Venn 4者共有菌，其中三者共有45，95，1，5
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/41culture/otu.txt result/venn/3_45.txt |sort -k3,3nr|less -S
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/41culture/otu.txt result/venn/3_95.txt |sort -k3,3nr|less -S

	# 制作有平均丰度，和物种注释的表
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $1,$5,a[$1]}' result/taxonomy_2.txt result/41culture/otu.txt > result/41culture/otu_mean_tax.txt


## 绘制网络
	
	# 制作数据文件
	OTU按差异比较计算整理
	注释信息按network.R整理
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$1]}' result/taxonomy_8.txt <(cut -f 1,4 result/compare/database.txt) | cut -f 1-2,5- | sed '1 s/OTUID\tLTEJ/ID\ttotal\tphylum\tclass\torder\tfamily\tgenus\tspecies/' |less > network/database.txt

	# 绘制籼粳稻属水平模块和差异 co_network_genus_LIndTej.R
	# 以Burkholderia 相关菌 与 Anaeromyxobacter 的小网络 co_network_genus_LIndTej_core.R


## /5/14 随机森林属水平区分籼粳稻

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

# 个性化分析 Custom analysis

## 低氮条件下挑选30个IND品种测宏基因组

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



## 挑选NRT样品测宏基因组

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


## 筛选各品种最好的3个样品

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


## OTU或分类与PCoA轴的差异
	
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
	


## 检查GWAS是否可以发现Anaeromyxobacter与SNP的关联
	# 基于LN TEJ/IND中Anaeromyxobacter最低和最高的20个品种(保证数据可分开)，与nsSNP关联，最为仅为e-3，10m21759092,10,21759092,0.643641032317718,0.475,40,0.893753744914084,0.894396427377941,1。
	gapit_IndTej_Anaeromyxobacter_top20.R
	# 修改为30个品种有重合，肯定无法找到差异，再改为25个品种
	10m21759092,10,21759092,0.962618742911769,0.46,50,0.856299549567109,0.85630648832468,1
	# 可能从样本量、数据波动程度、SNP背景均无法满足

## 筛选HN/LN下可高丰度菌及是否可培养
	# core_microbiome_culture.R筛选中位数并排序 
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result/41culture/otu.txt culture/core.txt > culture/core_culture.txt

### 差异菌/功能与氮相关基因显著相关
	
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



# 图表整理 Figures and legends

http://210.75.224.110/submit/rice_microbiome/

## 图1. 籼粳稻分型

### a, b 模式图
	
	地图，种植模式图——秦媛

### beta多样性
	
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
#	alpha_rare.sh -i `pwd`/result/alpha/rare.txt \
#	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' \
#	-o `pwd`/fig/1/alpha_ -h 3 -w 5

### Taxonomy 门+变形菌纲
	cut -f 3-4 result/taxonomy_8.txt|sort|uniq|grep 'Proteobacteria' # 为什么会有这么多结果，只选5类继续分析
	cat <(grep -v 'Proteobacteria' result/tax/sum_p.txt) <(grep 'proteobacteria' result/tax/sum_c.txt) > result/tax/sum_pc.txt
tax_stackplot.sh -i `pwd`/result/tax/sum_ -m '"pc"' -n 10 \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' -O FALSE \
	-o `pwd`/fig1/1/tax_pc_ -h 3 -w 5


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
	由筛选组，改为筛选亚组(品种)中位数的OTUs: 原万5为343个OTUs，万一为942个；最终丰度为0.2%, 381个OTUs
	compare_sub.R # 修改丰度筛选group为groupID2，接下来 
    rm plot_volcano
    make plot_venn
    make rmd
    grep -P -v '^\s*#' script/compare_sub.R > fig1/script/compare.R
    # 差异见xiangengASV_wilcoxon_main/result/compare/
    source=result/compare/
    target=fig1/data/
    cp $source/* $target/
    cp result/tax/sum_pc.txt $target/
    # A50LnCp6-A56LnCp6原来仅有3个差异刚过0.05，现在没差异：和之前的对比，是因为OTU拆分为ASV差异变小，OTU数量增加导致FDR减少
    sort -k5,5n result/compare/A50LnCp6-A56LnCp6_all.txt|less -S
    sort -k5,5n ../xiangeng/xiangeng_wilcoxon_main/result/compare/A50LnCp6-A56LnCp6_all.txt | less -S
    # 尝试另两组比较
    # "HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","A50LnCp7","A56LnCp7","A50LnCp6","A56LnCp6"
    # 新增 "A50HnCp7","A56HnCp7","A50HnCp6","A56HnCp6","HSoil1","LSoil1"
    # 新增 "nrtHnCp7","ZH11HnCp7","nrtLnCp7","ZH11LnCp7",
    # 新增 "nrtHnCp7","ZH11HnCp7","nrtLnCp7","ZH11LnCp7",
    # 新增 "nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7","A50HnSz7","A56HnSz7","A50LnSz7","A56LnSz7", 结果差异极小，暂不考虑
    # cp doc/venn.txt doc/venn.txt.190109
    # cp doc/compare.txt doc/compare.txt.190109
    # 手动修改比较组和venn
compare_OTU_ref.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
        -p 0.01 -q 0.05 -F 1.2 -t 0.2 \
        -d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","nrtHnCp7","ZH11HnCp7","nrtLnCp7","ZH11LnCp7","A50LnCp7","A56LnCp7","A50LnCp6","A56LnCp6","A50HnCp7","A56HnCp7","A50HnCp6","A56HnCp6"' \
        -o `pwd`/result/compare/ -r xiangengASV_wilcoxon_main/result/compare/LIND-LTEJ_all.txt

# 新增四8组，ASV从361变为484；改为0.1增长至853
compare_OTU.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
        -p 0.05 -q 0.05 -F 1.2 -t 0.1 \
        -d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","nrtHnCp7","ZH11HnCp7","nrtLnCp7","ZH11LnCp7","A50LnCp7","A56LnCp7","A50HnCp7","A56HnCp7"' \
        -o `pwd`/result/compare/ -C "groupID2"
cat result/compare/summary.txt

# 2019/1/10 制作科水平维恩
    # 制作对应科和差异列表
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print a[$1],$2}' result/taxonomy_8.txt result/compare/diff.list \
    | tail -n+2 | sort | uniq > result/compare/Fdiff.list  #  less -S
    # 绘制维恩图
    batch_venn.pl -i `pwd`/doc/venn.txt -d result/compare/Fdiff.list
    # 批量展示
    out=xiangengASV_wilcoxon_main_k1F
report_16S_Family.pl -g groupID -b ${out} -m "wilcox" -d `pwd`/doc/design.txt -c `pwd`/doc/compare.txt -v `pwd`/doc/venn.txt -a 0.2 -F 1.2 -p 0.01 -q 0.05 # -D  -F  -l  -S  -s   -t 
ln -sf `pwd`/${out}/ /var/www/html/report/16Sv2/${out}
rm -f ${out}/${out}



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
    make -B faprotax_calc 
    # 筛选IND和TEJ的注释，比较找出功能差异类
#	filter_otus_by_sample.sh -f result/faprotax/element_tab.txt -o result/faprotax/xiangeng -d doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"'
#	# 结果可用STAMP进一步探索
#	
#	# 筛选时间序列中上/下调两大类的OTU进行功能分析
#	mkdir -p fig/3
#	awk '$2>0' fig/2/otu_IND_common_specific_time_cor6.txt | cut -f 1 > fig/3/timecournse_increase.id
#	awk '$2<0' fig/2/otu_IND_common_specific_time_cor6.txt | cut -f 1 > fig/3/timecournse_decrease.id
#	# 以Incease为例, decrease
#	type=decrease
#	filter_otus_from_otu_table.py -i result/otutab_norm_tax.biom -o timecourse/${type}.biom --otu_ids_to_exclude_fp fig/3/timecournse_${type}.id --negate_ids_to_exclude
#	/usr/bin/python2.7 /mnt/bai/yongxin/software/FAPROTAX_1.1/collapse_table.py -i timecourse/${type}.biom -o timecourse/${type}.faprotax -g /mnt/bai/yongxin/software/FAPROTAX_1.1/FAPROTAX.txt --collapse_by_metadata 'taxonomy' -v --force # --out_report result/faprotax/report 
#	filter_otus_by_sample.sh -f timecourse/${type}.faprotax -o fig/3/timecournse_faprotax_${type} -d doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"'
#	# increase的差异功能类型均为IND>TEJ>soil，且以芳香、氮 相关
#	compare.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
#	-p 0.01 -q 0.01 -F 1.2 -t 0.0005 \
#	-d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"' \
#	-o `pwd`/result/compare/
#	# decrease的差异，stamp打开报错，但有时成功；下调无N循环相关，有
#	
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
        -p 0.05 -q 0.05 -F 1.2 -t 0.0005 \
        -d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1","V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50LnCp6","A56LnCp6"' \
        -o `pwd`/result/compare_far/ -N FALSE
	batch_venn.pl -i doc/venn.txt -d result/compare_far/diff.list
	# 注释比较结果
	rm result/compare_far/diff.list.venn*.xls.*
	batch2.pl -i 'result/compare_far/diff.list.venn*.xls' -d result/compare_far/database.txt -o result/compare_far/ -p vennNumAnno.pl
	# 绘制箱线图
	# 确定要展示的Features，两组共有
	tail -n 36 result/compare_far/diff.list.vennHTEJ_HIND_DLTEJ_LIND_D.xls.xls|cut -f 1 > result/compare_far/IND.list
	tail -n 9 result/compare_far/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls|cut -f 1 > result/compare_far/TEJ.list
	make plot_fa_barplot # 绘制单个功能的箱线图
    grep -P -v '^\s*#' script/alpha_boxplot_far.R > fig1/script/alpha_boxplot_far.R
	# 修改alpha_boxplot.R为alpha_boxplot_far.R



## 菌群与nrt关系

### 关键氮高效基因NRT不同形态、突变体可部分解析亚种差异

	"A50LnCp6","A56LnCp6","A50LnCp7","A56LnCp7","A50LnSz7","A56LnSz7","A50HnCp6","A56HnCp6","A50HnCp7","A56HnCp7","A50HnSz7","A56HnSz7","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","nrtHnCp7","ZH11HnCp7","ZH11LnCp7","nrtLnCp7","nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7"
	# nrt vs ZH11(TEJ): 主图："V3703HnCp6","ZH11HnCp6", 附图："V3703LnCp6","ZH11LnCp6",
	# 近等基因系：主图："A50LnCp7","A56LnCp7", 附图："A50LnCp6","A56LnCp6",
	# 主图/附图各分析一次：alpha, beta, 差异OTUs, venn: HTEJ_HIND_D LTEJ_LIND_D V3703HnCp6_ZH11HnCp6; HTEJ_HIND_D LTEJ_LIND_D A50HnCp7_A56HnCp7

主图4组: http://bailab.genetics.ac.cn/report/16Sv2/xiangeng_wilcoxon_main/

附图4组: http://bailab.genetics.ac.cn/report/16Sv2xiangeng_wilcoxon_supp

主图3.10.2和3.10.3中共有菌数量但不，但丰度很高，计算其占全部菌的相对丰度。详见fig1/venn_overlap_abundance.R



### 宏基因组KO注释

	从金桃处获得KO表和丰度，获取KO的功能描述；在kegg中没找到，picurst中没找到(输出结果有注释但不完整)，google搜索KO description download，找到biostar解答；https://www.genome.jp/kegg-bin/get_htext?ko00001.keg 中的Download htext/jason下载KO和描述；htext方便检索，而jason方便在线分析，如jason2table
	
	cd ~/rice/xianGeng/meta
	grep 'D      K' ko00001.keg | cut -c 8- | sed 's/  /\t/' > KO_description.txt #|  cat -A | less -S
	# 比较两组差异
    compare.sh -i ko.txt -c compare.txt -m "wilcox"         -p 0.01 -q 0.05 -F 1.2 -t 0 -N FALSE         -d design.txt -A group -B '"HnNrt","HnZH11"'         -o compare/ -N FALSE -U 100 # Pvalue和FDR并不显著，可能秩和检验需要较多的样本数
	# 注释KO
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' KO_description.txt ko.txt > ko.anno

	# 绘制单个KO和箱线图
	alpha_boxplot.sh -i ko.txt -m '"K02568","K02567","K00363","K10535","K00362"' \
	-d design.txt -A group -B '"HnZH11","HnNrt"' \
	-o ./ -h 2 -w 2.5 -t TRUE -n TRUE -U 1000000


    # 2018/11/14 完整KO表
	cd ~/rice/xianGeng/meta
    # 制作kotab.txt, design.txt, compare.txt
    rm -r compare
    compare.sh -i kotab.txt -d design.txt -c compare.txt \
        -m "ttest" -p 0.05 -q 1 -F 1.2 -t 0 -N FALSE -A GroupID -B '"HnZH11","HnNrt","LnZH11","LnNrt"' \
        -o compare/ -N FALSE -U 100 
    cat compare/summary.txt
        # Pvalue和FDR并不显著，可能秩和检验需要较多的样本数，HN要比LN明显特别多

    ## https://www.kegg.jp 搜索Nitrogen metabolism ，打开PATHWAY map00910，点击Ortholog table打开KO与物种对应表，保存第一行KO和对应描述为nitrogene/ko.id
    # 设置分析本批关键字，同类分析替换即可 nitrogen/sulfur 有代谢通路, potassium 只有M00454中两个KO 	K07646 K07667, phosphorus 相关KO由张鹏帆贡献
    type=sulfur
    mkdir -p ${type}
    # 转换为KO ID对应基因名称，并添加差异比较对应表头
    sed ":a;N;s/\n/ /g;ta" ${type}/ko_raw.id |sed "s/\t/\n/g"| sed "s/ /\t/;s/(//;s/)//;" | cut -f 1 -d '['| sed '1 i HnZH11_HnNrt'|less > ${type}/ko.id
    # 筛选KO对应比较结果
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' compare/HnZH11-HnNrt_all.txt ${type}/ko.id | grep -v -P '^$' > ${type}/HnNrt_all.txt
    # 添加注释
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$7]=$8} NR>FNR{print $0,a[$1]}' ~/github/Metagenome/denovo1/kegg/ko00001.tsv  ${type}/HnNrt_all.txt > ${type}/HnNrt_all_anno.txt
    # 查看正对照K02568 K02567 K00363
    # 筛选显著差异并标颜色，用于在线展示
    # awk '$4<0.05' ${type}/HnNrt_all.txt | cut -f 1-4,7-8 | awk '{if ($2>0) {print $0"\tred"}else {print $0"\tgreen"}}' > ${type}/HnNrt_diff.txt
    # cat ${type}/HnNrt_diff.txt
    # KO添加描述
    # head -n1 ~/github/Metagenome/denovo1/kegg/ko00001.tsv|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$7]=$8} NR>FNR{print $0,a[$1]}' ~/github/Metagenome/denovo1/kegg/ko00001.tsv ${type}/HnNrt_diff.txt > ${type}/HnNrt_diff_anno.txt
    # N和S都有上升的基因，钾无差异，P有一个下调


## 图5. 可培养菌

    # 可培养比例与自然样品  ~/culture/jingmei171208.rice.med.soy/makefile.man # 2018/10/17 培养菌与自然样本比较
    # 统计科数量
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/rep_seqs_tax.txt culture_A50/otu_cultured.txt > culture_A50/otu_cultured.tax
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/rep_seqs_tax.txt culture_IR24/otu_cultured.txt > culture_IR24/otu_cultured.tax
    cat culture_A50/otu_cultured.tax culture_IR24/otu_cultured.tax | cut -f 6 | sort | uniq | wc -l
    # 分菌见 ~/culture/rice/makefile.man # 2018/10/19 水稻分菌整理
    # 总结：1098个单菌，去掉一个序列异常的1097，去冗余519，序列比对剩512条序列



## 图5. 微生物与表型关联——2018/11/19去掉

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

    # 2018/11/16 株高、鲜重与OTU相关性 fig1\5pheno\LN_otu_mean_pheno_cor.r.txt，用echart可视化为LN_otu_mean_pheno_cor.r.pdf


## syncom重组体系表型分析
    # 数据整理为 wet/Xu-NbSycomFunctionData.txt, 分析代码见 fig1/5Syncom.Rmd



## 附图

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

## 附数据统计

	# 1. 品种数：68 IND
	cd ~/rice/xianGeng/fig1/ST
		27 TEJ
	sed -i 's/,/\t/g' 01.variety_geotable.csv # csv替换为tsv
	# 数据筛选95个，并添加95个
	awk '{FS=OFS="\t"} NR==FNR{a[$2]=$0} NR>FNR{print $0,a[$2]}' ~/rice/miniCore/doc/minicore_list.txt 01.variety_geotable.csv > 01.variety_geotable.txt 
	cut -f 14 01.variety_geotable.txt|tail -n+2|sort|uniq|wc -l # 44个国家

	# 2. 凌水海南样本数和测序量
	# 需要使用compare_sub.R中来原代码，在3compare中汇总样品和测序量

	# 3. 查看共有16个下调菌的丰度范围
	cd ~/rice/xianGeng
	awk '{FS=OFS="\t"} NR==FNR{a[$1]=$14} NR>FNR{print a[$1]}' xiangeng4wilcox/result/compare/HTEJ-HIND_all.txt fig1/3compare/otu_TEJ_common_specific.txt|sort -nr
	awk '{FS=OFS="\t"} NR==FNR{a[$1]=$14} NR>FNR{print a[$1]}' xiangeng4wilcox/result/compare/LTEJ-LIND_all.txt fig1/3compare/otu_TEJ_common_specific.txt|sort -nr
	awk '$2<0.05' fig1/ST/09.tiller_cor_p.txt|wc -l # 381个有214个显著P<0.05
	cd fig1/ST
	paste 09.tiller_cor.txt 09.tiller_cor_p.txt | cut -f 1,2,4> 09.tiller_cor_pr.txt

## 上传数据 PRJNA478068

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


## 补充结果2018/10/29 

	目前文章需要补充的图和表

	主图

	1.分菌详细流程图：我基金中的那个流程图吗？
	原图参考2018面上项目基金D:\life\grant\NSFC2018\figure\figure.pptx中图3，复制到fig1/20180610-组合图.pptx中：和白老师定方案，秦媛绘制

	2.菌保图：树吗？
	~/culture/rice/makefile.man 中 	# 绘制物种树，不带丰度和可培养 ，绘制进化树，无论科62，目22，都不能确保位于树的同一分枝上。
	改为graphlan绘制物种树，并添加标签，保存为fig1/graphlan.ai/pdf


	附图

	1. 无机氮功能实验结果
	2. A50实验结果（同IR24，数据）Figure6S.pdf

	Table

	1.分离菌占16S结果比例及taxonomy组成：自然样品，先不整理

	2.菌保信息对应表：verified & COTU，详见clutre/rice/makefile.man中# 制作isolate(鉴定结果)物种注释、来源和序列表 181024/TableS10.xls	# 制作COUT注释和序列信息 181024/TableS9.xls

	3.功能实验设计及表型结果（包括有机氮/无机氮，培养基成分，浓度，sycoms菌组成，表型数据等）



## 补充结果2018/11/17

1. IND/TEJ差异共有OTU在表型相关系数；

共有差异OTU列表，见： http://bailab.genetics.ac.cn/report/16Sv2/xiangeng_wilcoxon_main/result-otu.html#htej_hind_dltej_lind_dcde-5 http://bailab.genetics.ac.cn/report/16Sv2/xiangeng_wilcoxon_main/result-otu.html#htej_hind_eltej_lind_ecde-5

    wd=result/subspecies_otu_cor
    mkdir -p $wd
    # 筛选IND/TEJ富集的OTU
    tail -n 141 xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_D.xls.xls | awk '{print $0"\tIND_enriched"}'> $wd/IND_enriched_otu.id
    tail -n 16 xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls | awk '{print $0"\tTEJ_enriched"}'> $wd/TEJ_enriched_otu.id
    cat $wd/IND_enriched_otu.id $wd/TEJ_enriched_otu.id | sed '1 i OTUID\tAbundance\tTaxonomy\tType' > $wd/IND_TEJ_otu.id
    # 添加OTU相关
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' ./fig1/5pheno/LN_otu_mean_pheno_cor.r.txt $wd/IND_TEJ_otu.id | less -S > $wd/IND_TEJ_otu.cor # 添加亚种
    # 在ehcart绘制热图 http://www.ehbio.com/ImageGP/index.php/Home/Index/PHeatmap.html

2. Venn中共有比例

共有和特有OTU列表见: http://bailab.genetics.ac.cn/report/16Sv2/xiangeng_wilcoxon_main/result-otu.html#htej_hind_dltej_lind_dv3703hncp6_zh11hncp6_dde-5
    
    wd=result/nrt_venn_abundance
    mkdir -p $wd
    # 脚本见 fig1/venn_overlap_abundance.R
    cp xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D.xls.xls $wd/HTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D.xls
    cp xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls.xls $wd/HTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls
    # 替换名称
    # 整理为附表 TableSX.nrt_veen_list.xlsx
    

3. Figure S10统计

    # 检索原代码
    /rice/xianGeng/fig1]$grep 'Tiller number' *.Rmd
    # 5phenotype_cor.Rmd 中有代码和图，带wilcoxon检验的P值

4. Figure 3.I/J 统计注释菌的总丰度，注释到N的总丰度；
    # ../result/faprotax/report.mat在菌和功能对应表
    # ../result/41culture/otu_mean_tax.txt 丰度和注释
    # 查看注释菌的丰度
    awk '{a=a+$2} END {print a}' ../result/41culture/otu_mean_tax.txt # 总量99.752%
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1]}' ../result/41culture/otu_mean_tax.txt ../result/faprotax/report.mat | awk '{a=a+$1} END {print a}' # 注释的菌有62.972
    # 用R脚本统计每类的丰度，和几类总和的丰度 sum_faprotax_anno.R，补充为3comapre_en.rmd中3E段落中末属


5. 机器学习方法写清楚：使用了那个包，调用了哪个函数。
    /mnt/bai/yongxin/rice/xianGeng/fig1/随机森林分类方法描述.docx

6. Nitrogene related KO STAMP: Welch's t-test, Storey FDR


# 实验绘图

    otu=`cat wet/181220functional.bac.id|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'`

	alpha_boxplot.sh -i result/otutab.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1","V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50LnCp6","A56LnCp6"' \
	-m $otu -t TRUE -o result/otu_boxplot/fuc_ -n TRUE 
