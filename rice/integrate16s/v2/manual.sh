
	# 快速分析 Quick Start(所需文件准备好)
	find . -name "*" -type f -size 0c | xargs -n 1 rm -f # 清理零字节文件，用于从头重新分析项目清空makefile点位文件
	make library_split # 样本拆分、
	make fq_qc # 合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
	make host_rm # 序列去冗余、去噪、挑选代表序列、去嵌合、去宿主，获得OTU代表序列result/otu.fa
	make beta_calc # 生成OTU表、过滤、物种注释、建库和多样性统计
	# 清除统计绘图标记(重分析时使用)
    # 实验设计来自上级目录 cp -r ../doc/SL_*/ doc/ # Hn Bj
	rm -rf alpha_boxplot 
    make DA_compare # 或下面两t地
    rm measurable_OTU
	make DA_compare2 # 绘制alpha、beta、taxonomy和差异OTU比较 或DA_compare
	rm -f plot_volcano # 删除OTU差异比较可化标记
	make plot_manhattan # 绘制差异比较的火山图、热图、曼哈顿图
	make plot_venn # 绘制OTU差异共有/特有维恩图
    make culture
    make culture_graphlan
	make DA_compare_tax # 高分类级差异比较，维恩图绘制，2为reads count负二项分布统计
	make rmd # 生成网页报告，必须依赖的只有alpha, beta, taxonomy


# 相关脚本

    GWAS —— 云协作 GWAS流程-emmax-rice


# 1. 标准流程 Standard pipeline

	处理序列 Processing sequencing data

	# 1. 准备工作 Preparation

	## 1.1. 准备流程配置文件

	# 创建环境代码见~/github/Work/initial_project.sh
	# 本项目在xiangeng基础上继续
	cd ~/rice/integrate16s/v2
	cp ../ma* ./

	## 1.2. 初始化工作区

	# Initialize the working directory
	make init

	## 1.3. 准备原始数据

	# Prepare raw data
	# 数据来自三个课题：miniCore + timecourse + nrt1.1b + nrt1.1a + SL + epi + Gprotein + CTK + SD1

	
	# 合并实验设计
	# 合并三个项目的实验设计，检查样品是否有重名
	cp ~/rice/zjj.nitrogen/180116/doc/design.txt doc/design_nrt.txt
	cp ~/rice/gxx_CTK/doc/design.txt doc/design_CTK.txt
	cp ~/rice/rice.epi/180409/doc/design.txt doc/design_epi.txt
	cp ~/rice/zn.sd1/v3/doc/design.txt doc/design_sd1.txt
	cp ~/rice/Gprotein/v2/doc/design.txt doc/design_Gprotein.txt
	cp ~/rice/miniCore/doc/design.txt doc/design_minicore.txt
	cp ~/rice/strigolactone.LiJY/doc/design.txt doc/design_SL.txt
	cp ~/rice/timecourse/doc/design.txt doc/design_timecourse.txt
	cp ~/rice/nrt1.1a/doc/design.txt doc/design_nrt1a.txt

	# 统计实验样品行和唯一行，确实样品名唯一
	cat doc/design_* | grep -v 'SampleID'|wc -l
	cat doc/design_* | grep -v 'SampleID'|cut -f 1| sort|uniq|wc -l 
	# 9 projecct include 9136 samples

	# 合并实验设计，前7列共有，只保留前7列
	cat <(head -n1 doc/design_nrt.txt) <(cat doc/design_* | grep -v 'SampleID') | cut -f 1-7 > doc/design.txt

	# 添加miniCore亚种
	mv doc/design.txt doc/design0.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$6]}' ~/rice/miniCore/180319/doc/design_group.txt doc/design0.txt|sed 's/\t$/\tsubspecies/'|less -S > doc/design1.txt # 添加亚种
	awk 'BEGIN{FS=OFS="\t"} {print $0,$7$8}' doc/design1.txt|less -S>doc/design2.txt # 合并土壤类型和亚种
	cp doc/design2.txt doc/design.txt

	# 原始数据合并
	cat ~/rice/zjj.nitrogen/180116/temp/seqs_usearch.fa \
        ~/rice/gxx_CTK/temp/seqs_usearch.fa \
        ~/rice/rice.epi/180409/temp/seqs_usearch.fa \
        ~/rice/zn.sd1/v3/temp/filtered.fa \
        ~/rice/Gprotein/v2/temp/filtered.fa \
        ~/rice/miniCore/temp/seqs_usearch.fa \
        ~/rice/strigolactone.LiJY/temp/seqs_usearch.fa \
        ~/rice/timecourse/v2/temp/filtered.fa \
        ~/rice/nrt1.1a/temp/filtered.fa \
        | cut -f 1 -d ';' | sed 's/_/./g' > temp/filtered.fa
	# 从1.7 fa_unqiue 开始
	
	
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
    #9518854 nt in 25283 seqs, min 90, max 429, avg 376
    #Matching query sequences: 451326033 of 610388499 (73.94%)
    #real    1019m11.764s
    #user    35924m56.004s


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

    ## 比较nrt 中突变体与野生型及近等基因系 (NBT文中8个组) 2020/2/25
    mkdir -p doc/nrt
    tail -n4 ~/rice/xianGeng/doc/compare.txt > doc/nrt/compare.txt
    # http://210.75.224.110/report/16Sv2/rice_OTU_nrt_wilcox_v1

    ## 比较nrt 中突变体与野生型及近等基因系 (全部) 2020/2/25
    mkdir -p doc/nrt2
    tail -n+3 ~/rice/xianGeng/doc/compare1.txt > doc/nrt2/compare.txt
    # http://210.75.224.110/report/16Sv2/rice_OTU_nrt_wilcox_v1


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
    # 绘制D27在北京和上海的丰度
	alpha_boxplot.sh -i `pwd`/result/otutab.txt -m '"OTU_39"' \
        -d `pwd`/doc/design_SL.txt -A genotypeID -B '"d27RtBj","NpRtBj","d27RtHn","NpRtHn"' \
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


## 4.14 抑制分蘖功能菌
    
    # 2019/7/11 来源 SL-Mutants-Depleted细菌定殖实验结果汇总 20190704.pptx，P6页有OTU物种注释
    grep -P 'OTU_148\t|OTU_18\t|OTU_21\t|OTU_24\t' result/taxonomy_2.txt  # 注释完全不对应，不此最新v2版
    cd ~/rice/integrate16s
    grep -P 'OTU_148\t|OTU_24\t|OTU_18\t|OTU_21\t' result/taxonomy_2.txt  # 注释基本对应，只有OTU_148的属种有差别
    # 观察来自148(4)21(5)两个肠杆菌科，18(8),21(11)两个金黄杆菌属；分别在OTU、属、科水平观察这些OTU在两块地点的分布
    # 参考报告 http://210.75.224.110/report/16Sv2/rice_SL_Bj_wilcox_v4 http://210.75.224.110/report/16Sv2/rice_SL_Hn_wilcox_v4/
    # 编写这4个OTU的信息描述 SL/otu.list
    OTU_list=`cut -f 1 SL/otu.list|tail -n+2|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'`
    group_list=`cat doc/SL_Bj/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'`
	alpha_boxplot.sh -i `pwd`/result/otutab.txt -m $OTU_list \
        -d `pwd`/doc/design.txt -A groupID -B $group_list \
        -o `pwd`/SL/boxplot/BJ_ -h 6 -w 8 -t TRUE -n TRUE
    # 查看科水平
	alpha_boxplot.sh -i `pwd`/result/tax/sum_f.txt -m '"Enterobacteriaceae","Flavobacteriaceae"' \
        -d `pwd`/doc/design.txt -A groupID -B $group_list \
        -o `pwd`/SL/boxplot/BJ_ -h 6 -w 8 -t TRUE -n TRUE

    # 获得序列后 SL/experiment1.fa，与v2版比较再展示OTU丰度
    makeblastdb -dbtype nucl -in result/otu.fa
    blastn -query SL/experiment1.fa -db result/otu.fa -out SL/experiment1_otu.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 100 -evalue 1 -num_threads 9 
    cat SL/experiment1_otu.txt
    # 筛选97%以上菌，并取前6列：StockID,OTUID,Similarity,Length,Mismacth,Gap，共249个候选
    awk '$3>97' SL/experiment1_otu.txt|cut -f 1-6>SL/experiment1_otu97.txt
    # 添加LN丰度大于0.03%的菌，共13个候选
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR {print $0,a[$2]}' ~/rice/integrate16s/v2OTU/LN/result/otutab_mean.txt SL/experiment1_otu97.txt > SL/experiment1_otu97_mean.txt
    awk '$7>0.03' SL/experiment1_otu97_mean.txt > SL/experiment1_otu97_mean0.03.txt
    # 添加注释
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR {print $0,a[$1]}' SL/experiment1.txt SL/experiment1_otu97_mean0.03.txt > SL/experiment1_otu97_mean_anno.txt
    less -S SL/experiment1_otu97_mean_anno.txt
    
    # 只有1个OTU被筛选出来，其它的肠杆菌呢？发现LN丰度文件错误，重复计算
    # 手动查找其被筛选掉的原因，相关性和丰度
    grep 'Enterobacteriaceae' result/taxonomy_8.txt > SL/taxonomy_entero.txt # 共有9个OTU注释为肠杆菌
    # 查看其丰度
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR {print $0,a[$1]}' LN/result/otutab_mean.txt SL/taxonomy_entero.txt  > SL/taxonomy_entero_mean.txt
    cat SL/taxonomy_entero_mean.txt
    # 查看其相似度
    cut -f 1 SL/taxonomy_entero.txt|tr '\n' '|'
    grep -P 'OTU_32\t|OTU_82\t|OTU_221\t|OTU_606\t|OTU_720\t' SL/experiment1_otu.txt
    # OTU_82没有，只有94%相似度；606和720仅为95,93%相似度


    OTU_list=`cut -f 2 SL/experiment1_otu.txt|tail -n+2|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'`
    group_list=`cat doc/SL_Bj/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'`
	alpha_boxplot.sh -i `pwd`/result/otutab.txt -m $OTU_list \
        -d `pwd`/doc/design.txt -A groupID -B $group_list \
        -o `pwd`/SL/boxplot/BJ_ -h 6 -w 8 -t TRUE -n TRUE
    # 查看科水平
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$2]}' result/taxonomy_2.txt SL/experiment1_otu.txt
	alpha_boxplot.sh -i `pwd`/result/tax/sum_f.txt -m '"Enterobacteriaceae","Flavobacteriaceae"' \
        -d `pwd`/doc/design.txt -A groupID -B $group_list \
        -o `pwd`/SL/boxplot/BJ_ -h 6 -w 8 -t TRUE -n TRUE
    # 查看对应.txt文件，黄杆菌NpRtBj与d14RtBj和d53RtBj显著差别，但肠杆菌不显著差异0.77, 0.88；是否为多重比较引起，仅有3个基因型尝试
	alpha_boxplot.sh -i `pwd`/result/tax/sum_f.txt -m '"Enterobacteriaceae","Flavobacteriaceae"' \
        -d `pwd`/doc/design.txt -A groupID -B '"NpRtBj","d53RtBj","d14RtBj"' \
        -o `pwd`/SL/boxplot/BJ_3 -h 6 -w 8 -t TRUE -n TRUE
    cat SL/boxplot/BJ_3Enterobacteriaceae.txt # 与野生型极显著-5次方，是多重比较校正所致，而突变体间无差异

    # 3株实验菌进行物种注释
    # 基于单端全长，4/5号菌f:Enterobacteriaceae属为Klebsiella、Enterobacter
    usearch10 -sintax SL/experiment1.fa \
        -db /mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb -sintax_cutoff 0.6 -strand both \
        -tabbedout SL/experiment1.fa.tax -threads 32
    # 基于799-1193区注释，匹配引物为正向
	# 正向引物如799F AACMGGATTAGATACCCKG ，检索GGATTAGATACCC位于第3行前方>250bp
	cutadapt -g AACMGGATTAGATACCCKG -e 0.2 SL/experiment1.fa -o temp/stock_rc1.P5.fa
	# 切除引物 http://themicrobiome.com/en/16s/16s-primers 
	# 反向时先切反向引物，1193R ACGTCATCCCCACCTTCC，F GGAAGGTGGGGATGACGT 正常为1492R GGTTACCTTGTTACGACTT找不到，改用1391R GACGGGCGGTGTGTRCA F TGYACACACCGCCCGTC
	cutadapt -a GGAAGGTGGGGATGACGT -e 0.2 temp/stock_rc1.P5.fa -o SL/experiment1V57.fa
    # 一个只能注释到Enterobacteriaceae科，属为0.4可信度的Enterobacter，另一个注释为Cedecea属
    usearch10 -sintax SL/experiment1V57.fa \
        -db /mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb -sintax_cutoff 0.6 -strand both \
        -tabbedout SL/experiment1V57.fa.tax -threads 32
    grep 'Enterobacteriales' SL/experiment1V57.fa.tax
    # R1452/R1483/R2129分别为4/5/7为肠杆菌
    cut -f 1,4 SL/experiment1V57.fa.tax|sort # 补充至实验菌信息

# 5. 图表整理 Figures and tables

    # 项目图表目录：fig，发表时在github改为一作姓-年代-杂志格式
    wd=/mnt/bai/yongxin/medicago/AMF/fig
    mkdir -p ${wd}
    cd ${wd}
    mkdir -p data script fig1 fig2 fig3 fig4


## 图0. 数据筛选 Filtering data and table

    # 筛选miniCore相关的实验设计和OTU表
    # 筛选OTU表与非种对应的HN/LN/平均值
    script/fig0.Rmd



##  图1. 水稻miniCore品种根系菌群物种和功能描述

    # fig1.description.variety.sample.taxonomy.Rmd
cp ../script/fig1.description.variety.sample.taxonomy.rmd fig/fig1/fig1.Rmd

### 1A. 203个水稻品种全球分布：地图，种植模式图——秦媛
    
    # 整理代码，绘制的原图和操作描述
cp /mnt/bai/qinyuan/rice/minicore/map203/worldmap_203.R script/
cp /mnt/bai/qinyuan/rice/minicore/map203/mapdata_mix.csv fig1/

### 1B. 样品、品种整理多样性、OTU稀释取线 fig/fig1/fig1.Rmd 2019/4/1

cp ~/rice/integrate16s/doc/xiangeng/design.txt data/design_xiangeng.txt
cp ~/rice/integrate16s/v2/result/otutab.txt data/
cp ~/rice/integrate16s/v2/result/otutab_norm.txt data/

### 1C. 物种组成Taxonomy boxplot in genus and phylum
cp ~/rice/integrate16s/v2/result/tax/sum_* data/
cp ~/rice/integrate16s/doc/minicore/design.txt data/design_miniCore.txt

### 1X. 功能组成 - 暂时只有 miniCore籼粳的两块点各60个样品

cd ~/rice/integrate16s/v2
mkdir meta
# 所有水稻宏基因组样本列表
scp yongxin@192.168.0.32:~/rice/miniCore/result/design.txt meta/
# 水稻miniCore 203个样品信息
cp /mnt/bai/yongxin/rice/miniCore/doc/minicore_list.txt doc/
# 修改表头EFD与品种编号Variety一致
sed -i '1 s/EFD/Variety/' doc/minicore_list.txt
# 添加亚种注释
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$13} NR>FNR{print $1,$2,$3,$4,a[$2]}' doc/minicore_list.txt meta/design.txt > meta/metadata.txt
cut -f 5 meta/metadata.txt|sort|uniq -c


### 1D. phylogenetics tree 进化树

#### 出现频率筛选核心OTU 2019/7/9

    # 筛选核心OTU：基于miniCore样品 1183，筛选90%存在的
    mkdir result/core
	usearch11 -otutab_core result/minicore/otutab.txt -sintaxin temp/otu.fa.tax -tabbedout result/core/samples.txt
    # 最后列添检查出现的频率
	awk '{print $0"\t"$2/1183*100}' result/core/samples.txt | sed '1 s/OTU/OTUID/;1 s/0/Core/' > result/core/freq.txt
    # 统计各比例100 95 90%下在miniCore 1183个样品中OTU数量，分别为5，56，85；较之前97%聚类少很多(54，276，370)
    # 统计usearch计算频率的单位，总频率为1
    awk '{a=a+$4}END{print a}' result/core/freq.txt 
    for i in `seq 50 5 100`; do
        echo -ne $i"\t"
        awk -v i="$i" '$14>=i' result/core/freq.txt|wc -l|tr '\n' '\t'
        awk -v i="$i" '$14>=i' result/core/freq.txt|awk '{a=a+$4}END{print a}'
    done
    # ASV方式非常琐碎，80%都不及OTU水平丰度高。暂用80%刚及格的149条序列
    # 绘制进化树
    awk '$14>=80' result/core/freq.txt > result/core/core.txt
    cut -f 1 result/core/core.txt > result/core/core.id
    usearch10 -fastx_getseqs result/otu.fa -labels result/core/core.id -fastaout result/core/otu.fa
    # check number
    grep '>' -c result/core/otu.fa
    # Multiply alignment
    clustalo -i result/core/otu.fa -o temp/1d.otu_k1_align.fa --seqtype=DNA --full --force --threads=9
    make_phylogeny.py -i temp/1d.otu_k1_align.fa -o result/core/otu.tree
    # 物种注释和丰度



#### 基于可遗传OTU计算
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

### 1N. 物种/功能与表型关联(分蘗、株高、鲜重)，见fig1.Rmd/fig1_meta.Rmd 结尾。转换到fig3.phenotype.Rmd开头
    # 相关图代码，来自 ~/rice/xianGeng/fig1/6phenotype_cor_en.rmd



##  图2. 微生物组多样性和遗传关系 2019/4/2

cd ~/rice/integrate16s/v2
cp ../script/fig2.diversity.genetics.rmd fig/fig2/fig2.Rmd
cp result/beta/* fig/data/
### 2A/B. 多样性
cp result/alpha/index.txt fig/data/alpha.txt
cp result/beta/* fig/data/ 
# beta多样性
beta_pcoa.sh -i `pwd`/temp/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
    -d `pwd`/doc/xiangeng/design.txt  -A groupID -B `cat doc/xiangeng/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'` -E TRUE \
    -c `pwd`/doc/xiangeng/compare.txt \
    -o `pwd`/result/beta/ -h 3 -w 5

# 按四大亚种绘制alpha多样性
cd fig2/
alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
    -d `pwd`/doc/"xiangeng"/design.txt  -A groupID -B '"LIND","LTEJ","LAUS","LTRJ"' \
    -o `pwd`/result/alpha/ -h 3 -w 5    


## 图3. SL相关基因

### 3.1 PCoA和CPCoA绘制SL所有基因型的丰度
    # "d27RtBj","d17RtBj","d10RtBj","d3AHLRtBj","d3NpRtBj","d3RtBj","d14AHLRtBj","d14RtBj","d53RtBj","NpRtBj","d27RtHn","d17RtHn","d10RtHn","d3RtHn","d14RtHn","d53RtHn","NpRtHn"
    # d3RtBj与D10混合，而d3NpRtBj两类有分开
    beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac"' \
        -d `pwd`/doc/SL/design.txt  -A groupID -B '"d27RtBj","d17RtBj","d10RtBj","d3NpRtBj","d14RtBj","d53RtBj","NpRtBj"' -E TRUE \
        -c `pwd`/doc/SL_Bj/compare.txt \
        -o `pwd`/result/beta/ -h 3 -w 5
    beta_cpcoa.sh -i `pwd`/result/otutab_norm.txt -m '"bray"' \
        -d `pwd`/doc/SL/design.txt  -A groupID -B '"d27RtBj","d17RtBj","d10RtBj","d3NpRtBj","d14RtBj","d53RtBj","NpRtBj"' -E TRUE \
        -o `pwd`/result/beta/ -h 3 -w 5

### 3.2 肠杆菌OTU在差异比较中是否显著下调

    # 北京地点出现9个全部显著下调
    grep -P 'OTU_32\t|OTU_221\t' rice_OTU_SL_Bj_wilcox_v1/result/compare/*_sig.txt | less -S | sed 's/rice_OTU_SL_Bj_wilcox_v1\/result\/compare\///;s/:/\t/;s/_sig.txt//' > SL/otu_Bj_enter.txt
    wc -l SL/otu_Bj_enter.txt
    # 海南差异小，基因型少，出现5次全为显著下调
    grep -P 'OTU_32\t|OTU_221\t' rice_OTU_SL_Hn_wilcox_v1/result/compare/*_sig.txt | less -S | sed 's/rice_OTU_SL_Hn_wilcox_v1\/result\/compare\///;s/:/\t/;s/_sig.txt//' > SL/otu_Hn_enter.txt
    wc -l SL/otu_Hn_enter.txt


## 菌与分蘖的关系 fig/fig3/fig3.phenotype.Rmd

### 以分蘖为属性进行constrained PCoA

    # 再将分蘖按5个一组，取连续分类，发现有非常好的规律。/mnt/bai/yongxin/rice/integrate16s/v2/LN/result/beta/cpcoa_tiller_number_5_bray.pdf

    ### 3. 科水平，分蘖划分连续型归类再求相关  /mnt/bai/yongxin/rice/integrate16s/v2/fig/fig3/pheno_LN_variety_cat5_Enterobacteriaceae_IND_tiller_cat.pdf 存在非常好的相关系数，改为pearson相关系数可以显著相关。
    # ASV和属水平的相关分析，查看肠杆菌中的属， Enterobacteriaceae 又分为Cedecea、Pantoea、Dickeya、Cronobacter和Citrobacter属
    grep 'Enterobacteriaceae' result/taxonomy_8.txt | cut -f 7|sort|uniq |awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'
    # 有Leminorella、Pantoea存在明显负相关，Unassigned显著正相关，即多样性越高，分蘖越多？Cedecea只有0.25负相关
    # 再筛选肠杆菌ASV层面的118个菌
    grep 'Enterobacteriaceae' result/taxonomy_8.txt | cut -f 1|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'
    # 采用相似度筛选，51个97%的ASV，找到14个相关>0.6，包括56 Unassigned/95 Citrobacter /100 Pantoea
    grep -v 'OTU_18_8' SL/experiment1_otu.txt | awk '$3>97' | cut -f 2 | sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'

### 2019/8/8 筛选与分蘖相关高丰度菌88个；2019/8/26 更新为99个

    # fig3.phenotype.LN.sample.Rmd和fig3.phenotype.HN.sample.Rmd，根据OTU与分蘖显著相关，并用FDR校正，找共有的上调和下调
    # 要先按HN/LN同时筛选丰度，的OTU再统计分析，见### 5. OTU水平，分蘖划分连续型归类再求相关
    # 筛选显著差异的注释可培养菌
    awk '$10<0.05' fig/fig3/HNsample/pheno_7_OTU_pearson_tiller_cat_.txt|cut -f 1,3,4,5,8,10 > temp/temp.txt
    awk '$10<0.05' fig/fig3/LNsample/pheno_5_OTU_pearson_tiller_cat_.txt|cut -f 1,3,4,5,8,10 >> temp/temp.txt
    cut -f 1 temp/temp.txt|sort|uniq > SL/tiller_sig_OTU.txt
    # 添加物种注释
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0"\t"a[$1]}' result/taxonomy_2.txt SL/tiller_sig_OTU.txt > SL/tiller_sig_OTU.tax
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/taxonomy_8.txt SL/tiller_sig_OTU.txt > SL/tiller_sig_OTU.tax8
    usearch10 -fastx_getseqs result/otu.fa -labels SL/tiller_sig_OTU.txt -fastaout SL/tiller_sig_OTU.fa
    # 添加可培养ID
    blastn -query SL/tiller_sig_OTU.fa -db ~/culture/rice/stock/sequence.fa -out SL/tiller_sig_OTU.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 5 -evalue 1 -num_threads 9 
    awk '$3>97' SL/tiller_sig_OTU.blastn|cut -f 1-6>SL/tiller_sig_OTU.blastn97
    # 添加培养菌的注释
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$2]}' ~/culture/rice/stock/rice_stock_190731_full.txt SL/tiller_sig_OTU.blastn97 > SL/tiller_sig_OTU.blastn97anno
    # 输出新菌ID对应的多行：读入ID，原文件检查是否存在，存在于ID中则输出
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$1} NR>FNR{if($1 in a){print $0}}' SL/newCorOTU.id SL/tiller_sig_OTU.blastn97anno > SL/tiller_sig_OTU.blastn97annoNew


    # 与原始分菌库比较
    blastn -query SL/tiller_sig_OTU.fa -db ~/culture/rice/result/culture_select.fa -out SL/tiller_sig_OTU.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 
    awk '$3>97' SL/tiller_sig_OTU.blastn|cut -f 1-6>SL/tiller_sig_OTU.blastn97
    sed 's/^/OTU_/' ~/culture/rice/result/culture_select.xls > ~/culture/rice/result/culture_select.xlsx
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$2]}' ~/culture/rice/result/culture_select.xlsx SL/tiller_sig_OTU.blastn97 > SL/tiller_sig_OTU.blastn97annobatch1
    # 输出新菌ID对应的多行：读入ID，原文件检查是否存在，存在于ID中则输出
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$1} NR>FNR{if($1 in a){print $0}}' SL/newCorOTU.id SL/tiller_sig_OTU.blastn97annobatch1 > SL/tiller_sig_OTU.blastn97annobatch1New

    # 与第二批菌库比较并追加
    blastn -query SL/tiller_sig_OTU.fa -db ~/culture/rice/190626/result/culture_select.fa -out SL/tiller_sig_OTU.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 
    awk '$3>97' SL/tiller_sig_OTU.blastn|cut -f 1-6>SL/tiller_sig_OTU.blastn97
    sed 's/^/COTU_/' ~/culture/rice/190626/result/culture_select.xls > ~/culture/rice/190626/result/culture_select.xlsx
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$2]}' ~/culture/rice/190626/result/culture_select.xlsx SL/tiller_sig_OTU.blastn97 > SL/tiller_sig_OTU.blastn97annobatch2
    # 输出新菌ID对应的多行：读入ID，原文件检查是否存在，存在于ID中则输出
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$1} NR>FNR{if($1 in a){print $0}}' SL/newCorOTU.id SL/tiller_sig_OTU.blastn97annobatch2 > SL/tiller_sig_OTU.blastn97annobatch2New


   
    # 结果差异OTU，D14，D53与相关均一致性比较，先以北京为例用Venny比较并统计于PPT doc/SL-新结果汇总190731.pptx
    http://210.75.224.110/report/16Sv2/rice_OTU_SL_Bj_wilcox_v1
    http://210.75.224.110/report/16Sv2/rice_OTU_SL_Hn_wilcox_v1

    ## 分析样本与OTU直接相关性 2019/8/12
    fig3.phenotype.LN.sample.raw.Rmd，与分蘖分组相关基本一致，只是更多正相关，而人工分组找到肠杆菌科的OTU
    ## 分析品种与OTU直接相关性


### 表型(分蘖)-微生物组-基因组-环境互作，解析率 (interaction/interaction.Rmd) 2019/8/13 余泓

    ## 分析表型与基因型、微生物组的解析模型
    # 需要1200个样本名，表型，微生物组PC前5，表型PC前5，以及环境因子
    # 分别制作HN/LN的PCoA表

    # 2019/10/10 获得分析代码继续分批分析
    cd LinearModel

    cp ../interaction/HN/HN_tiller_M_G.data ./
    sed -i "1 s/^/\t/" HN_tiller_M_G.data
    cut -f 1-6 HN_tiller_M_G.data > HN_tiller_M.data
    cut -f 1-2,7- HN_tiller_M_G.data > HN_tiller_G.data

    cp ../interaction/LN/LN_tiller_M_G.data ./
    sed -i "1 s/^/\t/" LN_tiller_M_G.data
    cut -f 1-6 LN_tiller_M_G.data > LN_tiller_M.data
    cut -f 1-2,7- LN_tiller_M_G.data > LN_tiller_G.data

    cut -f 3-4,6-9,16-20 ../interaction/data/tiller_microbiomePCo_genoPC_envNPKtxt | less > AN_tiller_M_G_E.data
    cut -f 1-6 AN_tiller_M_G_E.data > AN_tiller_M.data
    cut -f 1-2,7-10 AN_tiller_M_G_E.data > AN_tiller_G.data
    cut -f 1-2,11 AN_tiller_M_G_E.data > AN_tiller_E.data
    cut -f 1-10 AN_tiller_M_G_E.data > AN_tiller_M_G.data
    cut -f 1-2,7-11 AN_tiller_M_G_E.data > AN_tiller_G_E.data
    cut -f 1-6,11 AN_tiller_M_G_E.data > AN_tiller_M_E.data
    sed -i '1 s/rownames//' AN*.data
    wc -l *.data
    perl 3_calculate_regressino_lm.pl # 遍历*.data文件，计算1行与其它相关，输出脚本Run_LM.r
    rm result.summary 
    Rscript Run_LM.r # 输出结果为 result.summary 
    cat result.summary 



### 分蘖在LN/HN下显著相关菌可视化 fig3.phenotype.Rmd 2019/8/14
    # 结果见D:\work\rice\integrate16s\v2OTU\fig\fig3\LNsample\pheatmap_f_LN_HN_cor.pdf和pheatmap_OTU_LN_HN_cor.pdf

### 2019/8/28 network 网络分析 丰度、正负相关菌在网络中的位置和秩数量


# 宏基因组数据分析

    cd ~/rice/integrate16s/v2/meta
    # 实验设计
    scp 192.168.10.32:~/rice/miniCore/result/metadata.txt metadata1.txt 
    scp 192.168.10.32:~/rice/miniCore2/result/metadata.txt metadata2.txt 
    # 整理实验设计
    # 两表合并整理包括，品种，土壤类型，项目和批次为metadata12.txt
    # 修改品种信息表：ID名EFD为SampleID，与metadata.txt一致
    sed  '1 s/EFD/Variety/' /mnt/bai/yongxin/rice/miniCore/doc/minicore_list.txt | less > minicore_list.txt
    # 添加New_Structure，非miniCore材料补Unknown
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$2]=$13} NR>FNR{print $0"\t"a[$2]}' minicore_list.txt metadata12.txt | sed 's/\t$/\tUnknown/' > metadata13.txt
    awk '$5==0' metadata13.txt|wc -l # 第0批NRT有48个样，包括NRT(36)和NRTSoild(12)个
    awk '$5==1' metadata13.txt|wc -l # 第1批miniCore有132个样
    awk '$5==1' metadata13.txt|cut -f 3,6 | sort | uniq -c # 分H/L两类，分别为IND(30)、TEJ(30)、Soil(6);
    # metadata1.txt有180个样品，miniCore共132个样品，分别为IND(60)、TEJ(60)、Soil(12); 此外NRT有48个样，包括NRT(36)和NRTSoild(12)个
    awk '$5==2' metadata13.txt|wc -l # 第2批miniCore有66个样
    awk '$5==2' metadata13.txt|cut -f 3,6 | sort | uniq -c # L类，分别为AUS(30)、TRJ(30)、ARO(6)
    # metadata2.txt有66个样品，miniCore共66个样品，只有L类，分别为AUS(30)、TRJ(30)、ARO(6)
    ln -s metadata13.txt metadata.txt
    
## 功能通路表
    # scp 192.168.10.32:~/rice/miniCore/result/12humann2/uniref_relab_unstratified.tsv function1.tsv
    # scp 192.168.10.32:~/rice/miniCore2/result/12humann2/uniref_relab_unstratified.tsv function2.tsv
    # 两文件不对应，无法直接合并，从原始文件利用软件自带脚本合并
    mkdir -p humann2
    scp 192.168.10.32:~/rice/miniCore/temp/12humann2/*.tsv humann2/
    scp 192.168.10.32:~/rice/miniCore2/temp/12humann2/*.tsv humann2/
    humann2_join_tables --input humann2/ --file_name pathabundance --output pathwayabundance.tsv
    sed -i 's/_Abundance//g' pathwayabundance.tsv
    humann2_renorm_table --input pathwayabundance.tsv --units relab --output pathwayabundance_relab.tsv
    humann2_split_stratified_table --input pathwayabundance_relab.tsv --output ./
    sed -i 's/# Pathway/MetaCyc_pathway/' pathwayabundance_relab_*stratified.tsv
    wc -l pathwayabundance*
    # 共有545个通路

    # 功能组成绘图见 fig/fig1_meta.Rmd
    # uniref中89%为UNMAPPED，10%为UNINTEGRATED，其它只有1%为注释通路？是宿主没有去除干净吗？查看软件结果说明为末比对和末知通路，比例太高，应该去除再分析
    # 结果只有3个与nitr相关通路，均为硝酸盐还原，但并没有与tiller显著相关

    # 计算菌菌与LN下分蘖均值相关系数，来自 ~/rice/integrate16s/script/fig3.phenotype.Rmd
    # 菌的相关位于 ~/rice/integrate16s/script/fig3

## 基因家族表
    humann2_join_tables --input humann2/ --file_name genefamilies --output genefamilies.tsv
    sed -i '1 s/_Abundance-RPKs//g' genefamilies.tsv
    # humann2_renorm_table --input genefamilies.tsv --units relab --output genefamilies_relab.tsv
    humann2_split_stratified_table --input genefamilies.tsv  --output ./
    wc -l genefamilies* # 280万个基因家族
    # ID转换为KO
    # UniRef90 to UniPort https://www.uniprot.org/uploadlists/ UniRef90转换为UniPortKB，只有一小半可识别，转换后即为去掉UniRef90_的前缀
    # UniPort to KEGG https://www.kegg.jp/kegg/tool/conv_id.html 转换不成功


## prokka注释太慢了，手动跑后半部至新目录，再合并
    wc result/metadata.txt # 180个样
    ls temp/23prokka/ |wc # 109个运行完或在运行，跑土壤慢，留30个继续，后40个新建任务至新目录
    mkdir temp/23prokka2
    time parallel --xapply -j 8 \
        "prokka temp/22megahit/{1}/final.contigs.fa --outdir temp/23prokka2/{1} \
        -prefix mg --metagenome --force --cpus 5 \
        --kingdom Archaea,Bacteria,Mitochondria,Viruses" \
        ::: `tail -n+2 result/metadata.txt | cut -f 1 | tail -n 50`

# GWAS

## 基因型准备

    # 生成emmax的tped，原始bed文件见 /mnt/zhou/chulab/miniCore/snp1.5x/T2.*
    time plink --bfile /mnt/zhou/chulab/miniCore/snp1.5x/T2 \
        --recode 12 --output-missing-genotype 0 --transpose \
        --out emmax/snp 
    # 生成Kinship矩阵，2288569 SNPs，43s
    time emmax-kin -v emmax/snp 

    # 协变量(可选)
    # 玉米中协变量第4行全为1，而水稻中有正、有负，如-9为缺失 emmax/snp.tfam
    # awk 'NR==FNR{a[$1]=$0} NR>FNR {print $1,a[$2]}' snp/cubic_PopStructure.txt emmax/snp.tfam |sed 's/\t/ 1 /' | sed 's/\t/ /g' > emmax/snp.cov
    # sed替换多个空格为制表符，cut取前5列，排除数值
    paste <(sed 's/[ ][ ]*/\t/g' emmax/snp.tfam|cut -f 1-2) /mnt/zhou/chulab/miniCore/snp1.5x/sativa.pca > emmax/snp.cov
    cat -A emmax/snp.cov|less


## LN/HN子版本 2019/7/31 

	
### 以LN为例，基于beta bray_curtis挑选的3个代表样品列表 ~/rice/miniCore/180319/LN/beta_norm/sample_ids.txt
	# 参考 ~/rice/miniCore/180319/manual.sh L165继续
	wd=LN/result
	mkdir -p ${wd}
	# 按之前筛选的样品编号提取样品
	usearch10 -otutab_sample_subset result/otutab.txt -labels ~/rice/miniCore/180319/LN/beta_norm/sample_ids.txt -output ${wd}/otutab_raw.txt
	# 按组合并，usearch group直接求合
	usearch10 -otutab_group ${wd}/otutab_raw.txt -labels ~/rice/miniCore/180319/LN/beta_norm/sample_ids_group.txt -output ${wd}/otutab.txt
	# 统计，最小2万，最大24万
	usearch10 -otutab_stats ${wd}/otutab.txt -output ${wd}/otutab.txt.sum
    cat ${wd}/otutab.txt.sum

	# 采用标准流程生成alpha, beta多样性，按最小样本量抽平
    # 详见 /mnt/bai/yongxin/rice/integrate16s/v2/LN

    # 转换和关联详见子目录


### 以HN为例，基于beta bray_curtis挑选的3个代表样品列表 ~/rice/miniCore/180319/HN/beta_norm/sample_ids.txt
	# 参考 ~/rice/miniCore/180319/manual.sh L165继续
	wd=HN/result
	mkdir -p ${wd}
	# 按之前筛选的样品编号提取样品，删除 R4159Ha
	usearch10 -otutab_sample_subset result/otutab.txt -labels ~/rice/miniCore/180319/HN/beta_norm/sample_ids.txt -output ${wd}/otutab_raw.txt
	# 按组合并，usearch group直接求合
    paste ~/rice/miniCore/180319/HN/beta_norm/sample_ids.txt <(cut -c1-5 ~/rice/miniCore/180319/HN/beta_norm/sample_ids.txt) > ~/rice/miniCore/180319/HN/beta_norm/sample_ids_group.txt
	usearch10 -otutab_group ${wd}/otutab_raw.txt -labels ~/rice/miniCore/180319/HN/beta_norm/sample_ids_group.txt -output ${wd}/otutab.txt
    wc -l ${wd}/otutab.txt
	# 统计，最小2万，最大24万
	usearch10 -otutab_stats ${wd}/otutab.txt -output ${wd}/otutab.txt.sum
    cat ${wd}/otutab.txt.sum

	# 采用标准流程生成alpha, beta多样性，按最小样本量抽平
    # 详见 /mnt/bai/yongxin/rice/integrate16s/v2/HN
    cp LN/ma* HN/

    # 转换和关联详见子目录


# 2019/10/31 D:\work\rice\integrate16s\v2OTU\fig\191031\ 整理最新版图

### 群体结构SNP树

    # 参考 ~/ehbio/train/08ReSeq/reseq/population_genomics/pipeline.sh
    # vcf生成fasta，然后建树，itol可视化
    # 参考 ~/rice/miniCore/mwas/genotype
    # 需要先将vcf合并，再拆分为样品
    cd ~/rice/integrate16s/v2OTU
    mkdir -p genotype
    cd genotype
    
   # 按染色体合并vcf-concat -c 检查列名，-f指定列表
    #ln ~/rice/miniCore/mwas/genotype/*.vcf ./
    #vcf-concat Chr*.vcf > rice.vcf # 2GB 太大，无法建树
    # 提取单个样品，只提取基因区差异SNP
    # 参考~/rice/miniCore/ # 筛选基因型只留引起AA变化的SNP

    # 按指定顺序合并结果为空
#    ls *.vcf|sed 's/Chr//'|sort -n|sed 's/^/Chr/' > VcfChr.list
#    vcf-concat -c -f VcfChr.list > rice.vcf
    # 仅合并missense类型的SNP
    vcf-concat /mnt/bai/yongxin/rice/miniCore/mwas/genotype/missense/Chr*.vcf > rice.vcf # 80M
    # 拆分
    mkdir -p sample
    # 所有SNP要5min，只要60s，共5万多个位点
    time vcf-subset -c X4202_X4202 rice.vcf > sample/X4202.vcf
    # 以Chr9为例
    head -n 7 /mnt/bai/yongxin/rice/miniCore/mwas/genotype/Chr9.vcf | tail -n1 | cut -f 10- | sed 's/\t/\n/g'| cut -f 1 -d '_' > varieties.list # less > 
    mkdir -p sample
    for i in `cat varieties.list`; do
        time vcf-subset -c ${i}_${i} rice.vcf > sample/${i}.vcf
    done

    # 提取序列
    # 只提取了ref和Chr9，改用
    cd ~/rice/integrate16s/v2OTU/genotype/sample
    CreateSnpTreeSeq.py ../varieties.list .vcf > ../merge_seq.fsa
    cd ..
    grep '>' -c merge_seq.fsa # 207条序列，此时可以进一步筛选序列，或提取部分
    # 转换为phy格式
    fasta2phy.pl -i merge_seq.fsa -o merge_seq.phy
    # 方法1. fasttree
    
# 2019/10/31 D:\work\rice\integrate16s\v2OTU\fig\191031\ 整理最新版图

## 群体结构

    # 参考 ~/ehbio/train/08ReSeq/reseq/population_genomics/pipeline.sh
    # vcf生成fasta，然后建树，itol可视化
    # 参考 ~/rice/miniCore/mwas/genotype
    # 需要先将vcf合并，再拆分为样品
    cd ~/rice/integrate16s/v2OTU
    mkdir -p genotype
    cd genotype
    
   # 按染色体合并vcf-concat -c 检查列名，-f指定列表
    #ln ~/rice/miniCore/mwas/genotype/*.vcf ./
    #vcf-concat Chr*.vcf > rice.vcf # 2GB 太大，无法建树
    # 提取单个样品，只提取基因区差异SNP
    # 参考~/rice/miniCore/ # 筛选基因型只留引起AA变化的SNP

    # 按指定顺序合并结果为空
    #    ls *.vcf|sed 's/Chr//'|sort -n|sed 's/^/Chr/' > VcfChr.list
    #    vcf-concat -c -f VcfChr.list > rice.vcf
        # 仅合并missense类型的SNP
        vcf-concat /mnt/bai/yongxin/rice/miniCore/mwas/genotype/missense/Chr*.vcf > rice.vcf # 80M

    # vcf提取fasta
        # 方法1. 按ehbio/reseq示例制作单样本文件，vcf格式不同，此脚本不读取样本信息列
        # 拆分
    #    mkdir -p sample
    #    # 所有SNP要5min，只要60s，共5万多个位点
    #    time vcf-subset -c X4202_X4202 rice.vcf > sample/X4202.vcf
    #    # 以Chr9为例
    #    head -n 7 /mnt/bai/yongxin/rice/miniCore/mwas/genotype/Chr9.vcf | tail -n1 | cut -f 10- | sed 's/\t/\n/g'| cut -f 1 -d '_' > varieties.list # less > 
    #    mkdir -p sample
    #    for i in `cat varieties.list`; do
    #        time vcf-subset -c ${i}_${i} rice.vcf > sample/${i}.vcf
    #    done
    #    # 只提取了ref和Chr9，改用
    #    cd ~/rice/integrate16s/v2OTU/genotype/sample
    #    CreateSnpTreeSeq.py ../varieties.list .vcf > ../merge_seq.fsa # 脚本不参考0/1的值，结果全一样。
    #    cd ..
    #    grep '>' -c merge_seq.fsa # 207条序列，此时可以进一步筛选序列，或提取部分
    #    # 转换为phy格式
    #    fasta2phy.pl -i merge_seq.fsa -o merge_seq.phy


    # 方法2. vcf-consensus提取序列，tabix报错
    # vcf排序
    cat rice.vcf | vcf-sort > rice_sort.vcf
    tabix -p vcf rice_sort.vcf # tbx_index_build failed: rice_sort.vcf
    cat /mnt/bai/public/ref/rice/IRGSP1/IRGSP-1.0_genome.fasta|vcf-consensus rice_sort.vcf > merge_seq.fsa
    # vcf-consensus

    # 手动编写提取
    # perl脚本提取
    vcf2fasta.pl -i rice_sort.vcf -o merge_seq.fa -h 11
    # 方法1. fasttree,3s
    time fasttreeMP -nt -gtr merge_seq.fa > merge_seq.tree # MP是多线程版，nt核酸，-gtr适合核酸的广义时间重塑模型
    # fasttree -nt -gtr merge_seq.fa > merge_seq.tree
    # 方法2. iqtree，也很快，也准确
    # 只筛选四大亚种的序列
    mv merge_seq.fa merge_seq.fa.bak
    wc -l ~/rice/miniCore/doc/minicore_list.txt # 203个品种
    # R脚本筛选
    table_subset.sh -i ~/rice/miniCore/doc/minicore_list.txt \
        -A Subspecies -B '"IND","TEJ","TRJ","AUS"' \
        -o minicore4subspecies.txt # 203 to 178
    # csvtk筛选
    csvtk -t grep -f Subspecies -p AUS,IND,TEJ,TRJ ~/rice/miniCore/doc/minicore_list.txt > minicore4subspecies.txt
    tail -n+2 minicore4subspecies.txt|cut -f 2| awk '{print $1"_"$1}' > select.id
    # 与实验数据筛选
    tail -n+2 ../fig/data/miniCore_metadata.txt|cut -f 4|sort|uniq > miniCore165.id # 165个品种
    # cat miniCore165.id <(cut -f 1 -d '_' select.id) | sort | uniq -u
    # 找到G4064、X4201和X4202三个特异的编号
    # cat miniCore165.id <(cut -f 1 -d '_' select.id) | sort | uniq -d | wc  # 与miniCore165.id一致
    # 提取序列，找到173条序列
    awk '{print $1"_"$1}' miniCore165.id > select.id
    usearch10 -fastx_getseqs merge_seq.fa.bak -labels select.id -fastaout merge_seq.fa
    # 查看缺失的序列 L4102 / M4120
    # cat <(grep '>' merge_seq.fa|sed 's/>//') select.id | sort | uniq -u 
    # 提取最后使用的样本ID
    grep '>' merge_seq.fa | cut -f 1 -d '_' | sed 's/>//' > final_varieties.id
    # 统计使用的SNP数量
    head -n2 merge_seq1.fa|tail -n+2|awk '{print length($1)}' # 58450个SNP misense

    # 再用datafilter筛选品种
    # 建树 10h
    iqtree --version # 1.6.12
    time iqtree -s merge_seq.fa -st DNA -m TEST -bb 1000 -alrt 1000 -nt 20 -pre iqtree165 -quiet # 结果为genotype/iqtree.contree
    format_fasta_1line.pl -i merge_seq.fa -o merge_seq1.fa
    grep -v '>' merge_seq1.fa |sort|uniq |wc -l # 165非冗余序列
    
    # 注释树
    grep '>' merge_seq.fa|sed 's/>//' >merge_seq.id
    paste merge_seq.id <(cut -f 1 -d '_' merge_seq.id) | sed '1 i ID\tEFD' > merge_seq.id2
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$2]=$13} NR>FNR {print $0,a[$2]}' ~/rice/miniCore/doc/minicore_list.txt merge_seq.id2 | sed 's/\t$/\tNA/' > merge_seq.id3
    Rscript ~/ehbio/amplicon/02Software/script/table2itol.R -a -c double -D planA -i ID -l EFD -t %s -w 0.5 merge_seq.id3
    Rscript ~/ehbio/amplicon/02Software/script/table2itol.R -a -d -c none -D planB -b New_Structure -i ID -l EFD -t %s -w 0.5 merge_seq.id3


### 亚种建树

    cd ~/rice/integrate16s/v2OTU/genotype
    # 按亚种，把序列进行一致性序列计算，再构建进化树
    cut -f 1 -d '_' merge_seq.fa > merge_clean.fa
    consensus_fa.pl -i merge_seq.fa -o consensus_merge # 生成.fa一致序列和.stat统计
    cut -f 2,13 ../doc/minicore_list.txt > ../doc/minicore_subspecies.txt

    # 4大亚种无法确定根，需要有ARO在内当外类群
    for i in IND TEJ TRJ AUS ARO; do
    # i=AUS
    # 筛选亚种ID
    grep -P "\t${i}$" ../doc/minicore_subspecies.txt | cut -f 1 > subspecies_${i}.id
    # 提取序列
    usearch10 -fastx_getseqs merge_clean.fa -labels subspecies_${i}.id -fastaout subspecies_${i}.fa
    # 转换为单行fasta格式
    format_fasta_1line.pl -i subspecies_${i}.fa -o subspecies_${i}1.fa
    # 统计一致序列
    consensus_fa.pl -i subspecies_${i}1.fa -o subspecies_${i}_consensus
    done
    rm consensus*
    cat subspecies_*_consensus.fa > consensus.fa
    sed -i 's/subspecies_//;s/_consensus//' consensus.fa
    time iqtree -s consensus.fa -st DNA -m TEST -bb 1000 -alrt 1000 -nt 20 -pre subspecies5 -quiet -redo # 结果为 subspecies5.contree
    # 结果采用itol展示

    # 以LN下均值为例进行高低分蘖的差异比较，均值下分蘖和菌丰度更能代表品种的真实情况
    cd ~/rice/integrate16s/v2OTU/LN
    # 新实验设计，由D:\work\rice\integrate16s\v2OTU\fig\191031\fig1.Rmd ### LN分蘖分组 产生LNmean/metadata.txt
    http://210.75.224.110/report/16Sv2/rice_miniCore_LN_tiller_v1


### 实验菌与高丰度菌比较

    # 2019/11/20号邮件，实验用菌 wet/徐浩然-分蘖相关细菌两批实验结果.xlsx
    # 选择菌的序列为final列，保存ID和序列为experiment.fa，效果为experiment.txt
    cd ~/rice/integrate16s/v2OTU/wet
    sed -i 's/\t/\n/' experiment.fa
    sed -i 's/^/>/' experiment.fa
    grep -c '>' experiment.fa # 65个实验菌
    # 筛选高丰度菌3%
    cut -f 1 ../LN/rice_miniCore_LN_tiller_v0.3/result/compare/T6-T1_all.txt > otu.id
    usearch10 -fastx_getseqs ../result/otu.fa -labels otu.id -fastaout otu.fa # 46个高丰度菌
    makeblastdb -dbtype nucl -in otu.fa
    blastn -query experiment.fa -db otu.fa -out experiment_otu.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 
    # 比对中有63个可比对，36个大于97%
    awk '$3>97' experiment_otu.txt | cut -f 1-3 > experiment_otu97.txt
    # 添加差异和相关及作用类型
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$6} NR>FNR {print $0,a[$2]}' ../LN/rice_miniCore_LN_tiller_v0.3/result/compare/T6-T1_all.txt experiment_otu97.txt > temp
    dos2unix ../fig/191031/LNmean/pheno_variety_OTU_pearson_tiller.txt
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR {print $0,a[$2]}' ../fig/191031/LNmean/pheno_variety_OTU_pearson_tiller.txt temp > temp1
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$5} NR>FNR {print $0,a[$1]}' experiment.txt temp1 > experiment_otu97_anno.txt



## 与菌库对应示例代码

    # 提取差异OTUID和序列
    cd ~/rice/integrate16s/v2OTU/wet
    cut -f 1 ../fig/191031/otutab_mean.txt > otu.id
    usearch10 -fastx_getseqs ../result/otu.fa -labels otu.id -fastaout otu.fa # 65个高丰度菌
    # 与水稻菌库比对
    blastn -query otu.fa -db /mnt/bai/yongxin/culture/rice/stock/sequence.fa -out experiment_otu.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 
    # 筛选97%以上相似匹配结果
    awk '$3>97' experiment_otu.txt | cut -f 1-3 > experiment_otu97.txt
    sed -i 's/\r//g' ../fig/191031/K2DiffCor.txt
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR {print $0,a[$1]}' experiment_otu97.txt ../fig/191031/K2DiffCor.txt > ../fig/191031/K2DiffCor_culture.txt
    # 筛选更多的比对结果
    blastn -query otu.fa -db /mnt/bai/yongxin/culture/rice/stock/sequence.fa -out experiment_otu10.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 
    awk '$3>97' experiment_otu10.txt | cut -f 1-3 > experiment_otu97-10.txt
    # 添加来源
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$5"\t"$6"\t"$7} NR>FNR {print $0,a[$2]}' /mnt/bai/yongxin/culture/rice/stock/rice_stock_200103_full.txt experiment_otu97-10.txt  > experiment_otu97-10_anno.txt

    # K2DiffCor_culture.txt HN/LN的差异、相关共4个批量的结果整合，最后标注出了显著的次数和是否有菌保
    #experiment_otu97-10_anno.txt OTU对应菌库中多个菌株的信息

    # 与原始分菌库比较 2019/12/23
    mkdir -p hts
    blastn -query otu.fa -db ~/culture/rice/result/culture_select.fa -out hts/tiller_sig_OTU.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 
    awk '$3>97' hts/tiller_sig_OTU.blastn|cut -f 1-6>hts/tiller_sig_OTU.blastn97
    sed 's/^/OTU_/' ~/culture/rice/result/culture_select.xls > ~/culture/rice/result/culture_select.xlsx
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$2]}' ~/culture/rice/result/culture_select.xlsx hts/tiller_sig_OTU.blastn97 > hts/tiller_sig_OTU.blastn97annobatch1
    # 输出新菌ID对应的多行：读入ID，原文件检查是否存在，存在于ID中则输出
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' hts/tiller_sig_OTU.blastn97annobatch1 ../fig/191031/K2DiffCor_culture.txt> hts/tiller_sig_OTU.blastn97annobatch1New

    # 与第二批菌库比较并追加
    blastn -query otu.fa -db ~/culture/rice/190626/result/culture_select.fa -out hts/tiller_sig_OTU.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 
    awk '$3>97' hts/tiller_sig_OTU.blastn|cut -f 1-6>hts/tiller_sig_OTU.blastn97
    sed 's/^/COTU_/' ~/culture/rice/190626/result/culture_select.xls > ~/culture/rice/190626/result/culture_select.xlsx
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$2]}' ~/culture/rice/190626/result/culture_select.xlsx hts/tiller_sig_OTU.blastn97 > hts/tiller_sig_OTU.blastn97annobatch2
    # 输出新菌ID对应的多行：读入ID，原文件检查是否存在，存在于ID中则输出
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' hts/tiller_sig_OTU.blastn97annobatch2 hts/tiller_sig_OTU.blastn97annobatch1New > hts/tiller_sig_OTU.blastn97annobatch2New

    # 添加OTU物种注释
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]"\t"$0}' ../result/taxonomy_8.txt hts/tiller_sig_OTU.blastn97annobatch2New > hts/tiller_sig_OTU.blastn97annobatch2New.tax

    # Batch5分蘖实验 wet/batch5  2020/2/18
    # wet/batch5/R2413-B02.txt 为一株错误菌保但有正相关；保存序列为stock.fa
    cd ~/rice/integrate16s/v2OTU/wet
    sed -i 's/^/>/;s/\t/\n/' batch5/stock.fa
    blastn -query batch5/stock.fa  -db otu.fa -out batch5/stock.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 
    awk '$3>97' batch5/stock.blastn |cut -f 1-6 > batch5/stock.blastn97 # 与表比对基本一致，2413对应OTU_268，在tiller_sig_OTU.blastn97annobatch2New.xlsx中显著差异或相关
    less batch5/stock.blastn97 # 查看到原始表对应，重点看四个目标菌-正确
    # 查R2413对应OTU_268外其它序列和相关分类菌 2020/5/6
    blastn -query batch5/stock.fa  -db ../result/otu.fa -out batch5/stock.blastnA -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 

### 200509batch6verify 2020/5/9 

    # 根据表格A7/B7为R2413
    cd ~/rice/integrate16s/v2OTU/wet/200509batch6verify
    cat seq/*.txt > seq.fa
    blastn -query seq.fa  -db ../../result/otu.fa -out stock.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 


## GWAS结果与实验菌比较 2020/1/13
    # GWAS分析结果见HN 或 LN/emmax/target_symbol.txt 
    cd ~/rice/integrate16s/v2OTU
    cat ?N/emmax/target_symbol.txt | cut -f 1 | sed 's/log2.//' | grep -v 'PositiveControl' | sort | uniq | tr '\n' '|'
    grep -P 'Azospirillum|Flavobacterium|Hyphomicrobium|Methylobacterium|OTU_20|OTU_203|OTU_28|OTU_39|OTU_43|OTU_61|Paenibacillus|Spirochaeta' wet/hts/tiller_sig_OTU.blastn97annobatch2New.tax
    # 在tiller_sig_OTU.blastn97annobatch2New.xlsx标记匹配结果，同时查看wet/experiment_otu97_anno.txt文件
    # R2488抑制分蘖效果好，查找其与OTU_5相近， k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__Burkholderia_cenocepacia


## R2488抑制分蘖菌
    # 查看菌来源，与OTU对应关系
    grep 'R2488' wet/* # 为与OTU_5 98.39%相似，在 experiment_otu97_anno.txt 中一致
    # OTU的物种注释
    # OTU_5/OTU_39/OTU_4/OTU_268分别查看
    grep -P 'OTU_268\t' result/taxonomy_2.txt # c__Betaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__Burkholderia_cenocepacia
    # 查看在miniCore中根据分蘖分组的差异情况
    grep -P 'OTU_268\t' fig/191031/K2DiffCor_culture.txt # 它在HN/LN均显著差异和显著相关
    # 查看SL中差异情况
    grep -P 'OTU_268\t' rice_OTU_SL_Bj_wilcox_v1/result/compare/*sig.txt # 结果没有OTU_268，是在这两地中不表达吗？表达低，我看
    grep -P 'OTU_268\t' rice_OTU_SL_Hn_wilcox_v1/result/compare/*sig.txt
    grep -P 'OTU_4\t' rice_OTU_nrt_wilcox_v1/result/compare/*sig.txt
    grep -P 'OTU_4\t' rice_OTU_nrt2_wilcox_v1/result/compare/*sig.txt
    # 查看OTU_268在基因型中的分布，绝大多数无表达或极低；进一步看属也比较低
	alpha_boxplot.sh -i `pwd`/result/otutab_norm.txt -m '"OTU_268"' \
        -d `pwd`/doc/design_SL.txt -A genotypeID -B '"d27RtBj","d17RtBj","d10RtBj","d3AHLRtBj","d3NpRtBj","d3RtBj","d14AHLRtBj","d14RtBj","d53RtBj","NpRtBj","d27RtHn","d17RtHn","d10RtHn","d3RtHn","d14RtHn","d53RtHn","NpRtHn"' \
        -o `pwd`/result/otu_boxplot/ -h 3 -w 10 -t TRUE -n TRUE
    # 对应属为 Burkholderia/Pleomorphomonas/Exiguobacterium/Sphingomonas
    # 查看属水平差异
    grep -P 'Sphingomonas\t' rice_OTU_SL_Bj_wilcox_v1/result/compare_g/*sig.txt # 找到多个结果，都是下调
    grep -P 'Sphingomonas\t' rice_OTU_SL_Bj_wilcox_v1/result/compare_g/*sig.txt|cut -f 4 -d '/'|cut -f 1 -d '-'
    grep -P 'Sphingomonas\t' rice_OTU_SL_Hn_wilcox_v1/result/compare_g/*sig.txt # 在Hn无显著差异
    grep -P 'Pleomorphomonas' rice_OTU_nrt2_wilcox_v1/result/compare_g/*sig.txt # 在Hn无显著差异
	alpha_boxplot.sh -i `pwd`/result/tax/sum_g.txt -m '"Sphingomonas"' \
        -d `pwd`/doc/design_SL.txt -A genotypeID -B '"d27RtBj","d17RtBj","d10RtBj","d3AHLRtBj","d3NpRtBj","d3RtBj","d14AHLRtBj","d14RtBj","d53RtBj","NpRtBj","d27RtHn","d17RtHn","d10RtHn","d3RtHn","d14RtHn","d53RtHn","NpRtHn"' \
        -o `pwd`/result/otu_boxplot/ -h 3 -w 10 -t TRUE -n TRUE
    # 只看北京在该属在SL中变化，有显著下调的结果
	alpha_boxplot.sh -i `pwd`/result/tax/sum_g.txt -m '"Sphingomonas"' \
        -d `pwd`/doc/design.txt -A groupID -B '"A50HnCp6","A50HnCp7","A50HnSz7","A50LnCp6","A50LnCp7","A50LnSz7","A56HnCp6","A56HnCp7","A56HnSz7","A56LnCp6","A56LnCp7","A56LnSz7","nrtHnCp7","nrtHnSz7","nrtLnCp7","nrtLnSz7","V3703HnCp6","V3703LnCp6","ZH11HnCp6","ZH11HnCp7","ZH11HnSz7","ZH11LnCp6","ZH11LnCp7","ZH11LnSz7"' \
        -o `pwd`/result/otu_boxplot/ -h 3 -w 10 -t TRUE -n TRUE

## 优秀表型菌 2020/02/24
    # 代码参考上面 ## R2488抑制分蘖菌

    # 查看GWAS关联
    # /mnt/bai/yongxin/rice/integrate16s/v2OTU/LN/emmax 中 OTU_5 Burkholderia，并在HN/LN中查找统计结果
    grep -P 'OTU_5|Burkholderia' emmax/target_symbol.txt
    # LN下OTU_5位于emmax目录，log2.OTU_5位于emmax_log2中; HN下全在emmax目录中
    
    # 查看t-test结果
    代码参考云笔记：GWAS流程-emmax-rice，表型与基因SNP分组统计
    grep -P 'HIGH|MODERATE' ${out}/${i}.anno | cut -f 2-4,7 | uniq | less -S > ${out}/${i}.anno.good


## miniCore核心微生物组 2019/11/24
    
    mkdir -p miniCore
    # 203个品种在两块地中的分布，详见 fig/CoreMicrobiome.Rmd

### 绘制cladogram
    # 筛选核心菌ID及对应物种注释
    awk '$3>=0.9' otutab_metadata.txt | cut -f 4 > otu0.9.id # 289个，288OTU，去除末注释的乘95个  taxonomy.txt
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR {print a[$1]}' ../result/taxonomy_8.txt otu0.9.id | awk '$6!="Unclassified"' > taxonomy.txt # grep -v 'Unclassified' 
    
    # 生成骨架和背景色、标签和方向 fig/CoreMicrobiome.Rmd ### 制作graphlan文件
    cat /mnt/bai/yongxin/culture/rice/graphlan/global.cfg 2_annotation_Family.txt > track0
    # 结点ID
    cut -f 5 -d '.' 1_tree_plain.txt > genus.id

    # 1. 添加属合并丰度均值 ###属水平均计算，OTU按属合并，再筛选即可
    study=miniCore
    track=1
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR {print $1,a[$1]}' RAmean_genus.txt genus.id > taxonomy.txt.RA${track}${study}
    cat <(cat ~/github/Note/R/format2graphlan/cfg/heat3.cfg|sed "s/3/${track}/;s/green/yellow/") <(sed "s/\t/\tring_alpha\t${track}\t/g" taxonomy.txt.RA${track}${study} | tail -n+2 ) > track${track}${study}
    # 柱状用log2，热图用zscore效果最好？默认也挺好，好像已经内部进行了转换
    awk '{a=a+$2} END {print a}' taxonomy.txt.RA${track}${study} # 80.3%

    # 2. 添加NBT北京野生型数据
    mkdir -p pubNBT2019
    ## 中国北京2019 NBT ~/rice/integrate16s/v2OTU/fig/CoreMicrobiome.Rmd # pubNBT2019
    study=pubNBT2019
    track=2
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR {print $1,a[$1]}' RAmean_genus_${study}.txt genus.id | sed 's/\t$/\t0/' > taxonomy.txt.RA${study} 
    awk '$2>0' taxonomy.txt.RA${study}|wc -l # 115个属，有115个检测到值
    awk '{a=a+$2} END {print a}' taxonomy.txt.RA${study} # 73.4%
    # 绘制丰度热度
    # sed "s/\t/\tring_alpha\t${track}\t/g" taxonomy.txt.RA${study} | tail -n+2 > taxonomy.txt.RA${study}${track} # 柱状用log2，热图用zscore
    # cat  <(cat /mnt/bai/yongxin/culture/rice/graphlan/abundance_heat.cfg|sed "s/3/${track}/") taxonomy.txt.RA${study}${track} > track${study}
    # 绘制有无，track为位置整数，m为颜色,0控制丰度筛选，R为形状
    cat <(cat /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg|sed "s/1/${track}/;s/m/b/") <(awk '$2>0' taxonomy.txt.RA${study} | cut -f 1 | sed "s/$/\tring_shape\t${track}\tR/") > track${track}${study}

    # 3. 添加日本公共数据热图，需要继续研究实验设计筛选
    study=pubRice2016
    track=3
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR {print $1,a[$1]}' RAmean_genus_${study}.txt genus.id | sed 's/\t$/\t0/' > taxonomy.txt.RA${study}${track} 
    awk '$2>0' taxonomy.txt.RA${study}${track}|wc -l # 115个属，有104个检测到值
    awk '{a=a+$2} END {print a}' taxonomy.txt.RA${study}${track} # 18.6%
    #sed "s/\t/\tring_alpha\t${track}\t/g" taxonomy.txt.RA${study} | tail -n+2 > taxonomy.txt.RA${study}${track} # 柱状用log2，热图用zscore
    #cat  <(cat /mnt/bai/yongxin/culture/rice/graphlan/abundance_heat.cfg|sed "s/3/${track}/") taxonomy.txt.RA${study}${track} > track${track}${study}
    cat <(cat /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg|sed "s/1/${track}/;s/m/g/") <(awk '$2>0' taxonomy.txt.RA${study}${track} | cut -f 1 | sed "s/$/\tring_shape\t${track}\tR/") > track${track}${study}
    
    # 4. 添加美国公共数据热图
    study=pubGB2019
    track=4
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR {print $1,a[$1]}' RAmean_genus_pubGB2019.txt genus.id | sed 's/\t$/\t0/' > taxonomy.txt.RA${study}${track}
    # sed 's/\t/\tring_alpha\t4\t/g' taxonomy.txt.RA2 | tail -n+2 > 4_annotation_match.txt # 柱状用log2，热图用zscore
    awk '$2>0' taxonomy.txt.RA2|wc -l # 115个属，有107个检测到值
    awk '{a=a+$2} END {print a}' taxonomy.txt.RA2 # 29.6%
    #cat  <(cat /mnt/bai/yongxin/culture/rice/graphlan/abundance_heat.cfg|sed 's/3/4/') 4_annotation_match.txt > track2
    # 添加公共数据有无
    #awk '$2>0' taxonomy.txt.RA2 | cut -f 1 | sed 's/$/\tring_shape\t11\tR/' > 11taxonomy.txt.RA2_R
    #cat <(cat /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg|sed 's/1/11/') 11taxonomy.txt.RA2_R > track${track}${study}
    cat <(cat /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg|sed "s/1/${track}/;s/m/m/") <(awk '$2>0' taxonomy.txt.RA${study}${track} | cut -f 1 | sed "s/$/\tring_shape\t${track}\tR/") > track${track}${study}

    # 5. 添加意大利属有无
    study=pubEMR2016
    track=5
    cat <(cat /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg|sed "s/1/${track}/;s/m/c/") <(cat RAmean_genus_${study}.txt | sed "s/$/\tring_shape\t${track}\tR/") > track${track}${study}
    wc -l track${track}${study} # 74/115
    # 添加列表，有填1
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=1} NR>FNR {print $1,a[$1]}' RAmean_genus_${study}.txt genus.id | sed 's/\t$/\t0/' > taxonomy.txt.RA${study}${track}


    # 6. 添加埃及属有无(太差，删除)
    study=pubME2019
    track=6
    cat <(cat /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg|sed "s/1/${track}/;s/m/y/") <(cat RAmean_genus_${study}.txt | sed "s/$/\tring_shape\t${track}\tR/") > track${track}${study}
    wc -l track${track}${study} # 46/115

    # 汇总
    # 生成注释文件
    cat track* > graphlan_annotate.txt
    # 注释文件修改树
    graphlan_annotate.py --annot graphlan_annotate.txt 1_tree_plain.txt graphlan.xml
    # 绘图，size是图片大小，越大，字相对越小
    graphlan.py graphlan.xml graphlan5.pdf --size 5

### 检查目标13个属是否表达和丰度
    
    cd ~/rice/integrate16s/v2OTU/miniCore
    # 查看属
    tail -n+2 ../fig/data/genus_list_1sig_simCor_2dot.txt|cut -f1
    # 按fig3中保存列表cor_genus.id
    # 查看数据格式
    head -n3 RAmean_genus_*
    # 以一个为例 pubNBT2019 pubRice2016 pubGB2019 pubEMR2016 pubME2019
    pub=pubME2019
    sed -i 's/\r//' cor_genus.id
    for i in `cat cor_genus.id`; do
      grep $i RAmean_genus_${pub}.txt
    done
    # 整理至minicore/cor_genus.xlsx


## 附录1. 公共数据整理比对核心菌

### 美国Genome Biology 2019 (GB2019) 2019/11/20

    # 下载github https://github.com/bulksoil/SoilDomestication
    cd ~/github/SoilDomestication
    git clone https://github.com/bulksoil/SoilDomestication
    # 数据全完没有，只有一个目录中有RDS数据，读取后有119个样品名，也无物种信息
    # Zenodo数据
    wget https://zenodo.org/record/3372822#.XdSzUpozYUE
    # SoilDomesticationData/OTU_Files/soil_domestication_sample_info.tsv # 442个样品
    #NCBI under project no. PRJNA548898 [56]. All processed datasets have been deposited in a Zenodo repository https://doi.org/10.5281/zenodo.3372822 [57]. R notebooks for the full analyses are freely available under the GNU General Public License v3.0 in the GitHub r
    # PRJNA548898在NCBI有样本214个，附表15个中没有实验设计，https://www.ncbi.nlm.nih.gov/bioproject/PRJNA548898
    # 点击样本数量 214进行列表，再点击Send results to Run selector，下载样本列表 Table
    # 筛选水稻且非土壤，找到95个水稻的根内和根际土样品
    # 下载数据
    prefetch SRR9617953
    # -A按编号下载数据无效；可以转换
    fastq-dump SRR9617953/SRR9617953.sra --split-spot -O seq/ --gzip --split-3 # 2.10.0

    # 批量下载
    cd ~/rice/integrate16s/v2OTU/pubGB2019/
    mkdir -p seq/sample
    for i in `tail -n+2 doc/metadata.txt|cut -f 1`; do
        # prefetch ${i} --output-directory seq/
        fastq-dump seq/${i}/${i}.sra --split-spot -O seq/ --split-3
    done
    rename 's/fastq/fq/' seq/sample/* # 改名至符合流程，从merge开始
    # GTGCCAGCMGCCGCGGTAA GGACTACHVGGGTWTCTAAT # 正反向引物均没匹配成功，是否已经去除了吗？
    # 接下来采用流程进行 pubGB2019/

### 日本2016 Rice Greengene
    mkdir -p pubRice2016
    cd pubRice2016    
    ln ../doc/list/07* ./
    biom convert -i 07JapanRice.biom -o otutab_gg85.txt --to-tsv
    sed -i '/# Construc/d;s/#OTU ID/OTUID/' otutab_gg85.txt
    wc -l otutab_gg85.txt # 有152300行
    # 只筛选greengeneID的，去除New后10797，也绝不可能是85_otus，而且基本无法匹配
    grep -v 'New' otutab_gg85.txt > otutab.txt # |wc -l
    grep '>' -c ~/ref/greengenes/gg_13_8_otus/rep_set/85_otus.fasta # 才有5088条序列，只有97和99比此表大，分别20万和40万条
    # 比较见CoreMicrobiome.Rmd ## pubRice2016 

### 意大利 2014 Rice 物种注释，与我们的注释文字比对
    mkdir -p pubEMR2016
    cd pubEMR2016
    cut -f 1 taxonomy.txt|tail -n+2|awk '{print "OTU_"NR"\t"$0}'|less>taxonomy2.txt
    sed 's/;/\t/g' taxonomy2.txt | sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' > otus.tax
    # 与树比对，见pubEMR2016

### 埃及 2019 ME 
    # PRJNA526033 但检索不到，可能末公开
    # 查看20个附件，没有目录；在文中有ESM中描述，需要检索；找到ESM3保存为taxonomy
    mkdir -p pubME2019/
    cd pubME2019/
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["d"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
        split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
        print $1,a["d"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
        taxonomy.txt > otus.tax
    sed -i '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' otus.tax


### 与核心OTU比对：gg提取v5-v7区域，sintax物种注释；按科和属合并丰度 2019/11/28
    # 提取采用QIIME2导入，提取和导出
    cd ~/ref/greengenes/gg_13_8_otus/rep_set
    conda activate qiime2-2019.7
    time qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path 97_otus.fasta \
      --output-path 97_otus.qza
    time qiime feature-classifier extract-reads \
      --i-sequences 97_otus.qza \
      --p-f-primer AACMGGATTAGATACCCKG \
      --p-r-primer ACGTCATCCCCACCTTCC \
      --o-reads 97_otusV5-7.qza # 7m
    qiime tools export \
      --input-path 97_otusV5-7.qza \
      --output-path 97_otusV5-7
    # 物种注释
    usearch10 -sintax 97_otusV5-7/dna-sequences.fasta \
        -db /mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb -sintax_cutoff 0.8 -strand both \
        -tabbedout 97_otusV5-7/dna-sequences.fasta.tax -threads 32
    cut -f 1,4 97_otusV5-7/dna-sequences.fasta.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > 97_otusV5-7/taxonomy_2.txt
    # 生成物种表格：注意OTU中会有末知为空白，补齐分类未知新物种为Unclassified
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="Unclassified";a["p"]="Unclassified";a["c"]="Unclassified";a["o"]="Unclassified";a["f"]="Unclassified";a["g"]="Unclassified";a["s"]="Unclassified"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' 97_otusV5-7/taxonomy_2.txt | sed '1 i #OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > 97_otusV5-7/taxonomy_8.txt
    # 去除#号和空格，会引起读取表格分列错误
    sed -i 's/#//g;s/ //g' 97_otusV5-7/taxonomy_8.txt
    # 在fig/CoreMicrobiome.Rmd中按科、属汇总求均值



## 选菌合成群落方案SynComm 2020/6/3

    # 筛选4个基因中千分之一共有上调47个，下调14个，不变12个。
    ll result/compare/*.xls.xls
    tail -n 47 result/compare/*E.xls.xls | cut -f 1,2,3,5 | sed 's/$/\tEnriched/' > wet/SL4_conserved_OTU.txt
    tail -n 14 result/compare/*D.xls.xls | cut -f 1,2,3,5 | sed 's/$/\tDepleted/' >> wet/SL4_conserved_OTU.txt
    tail -n 12 result/compare/*N.xls.xls | cut -f 1,2,3,5 | sed 's/$/\tNotSig/' >> wet/SL4_conserved_OTU.txt
    sed -i '1 i OTUID\tStockID\tSimilarity\tRA\tLevel' wet/SL4_conserved_OTU.txt
    # 筛选可培养菌，25个
    awk '$3>=97' wet/SL4_conserved_OTU.txt > wet/SL4_conserved_OTU_cultured.txt

    # 建树可视化选菌
    # 筛选高丰度，conserved建树，可培养，并树上注释门、可培养
    # 标注列的序列
    #### 方案1. 千分之1保守的菌
    mkdir -p wet/syncom
    # 提取序列
    cut -f 1 wet/SL4_conserved_OTU_cultured.txt > temp/syncom_OTU.id
    usearch10 -fastx_getseqs result/otu.fa -labels temp/syncom_OTU.id -fastaout temp/syncom_OTU.fa
    # 建树
    threshold=0.1
    mkdir -p wet/syncom/p${threshold}
    time muscle -in temp/syncom_OTU.fa -out temp/syncom_OTU_aligned.fas
    time iqtree -s temp/syncom_OTU_aligned.fas -st DNA -m TEST -bb 1000 -alrt 1000 -nt 20 -pre wet/syncom/p${threshold}/tree -quiet -redo
    # 制作注释文件
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result/taxonomy_8.txt wet/SL4_conserved_OTU_cultured.txt |cut -f 1-5,7-> wet/syncom/p${threshold}/annotation.txt
    Rscript ~/bin/table2itol/table2itol.R -a -c double -D  wet/syncom/p${threshold} -i OTUID -l Genus -t %s -w 0.5 wet/syncom/p${threshold}/annotation.txt
    # 输出更多候选
    blastn -query temp/syncom_OTU.fa -db /mnt/bai/yongxin/culture/rice/stock/sequence.fa -out temp/culture_otu.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 999 -evalue 1 -num_threads 9
    awk '$3>97' temp/culture_otu.blastn | cut -f 1-3,6,7 | sed '1 i OTUID\tStockID\tSimilarity\tMismatch\tGaps'> wet/syncom/p${threshold}/culture_otu.blastn.txt


    #### 人工挑选菌去冗余
    cd ~/rice/integrate16s/v2OTU/wet
    i=200617
    cat ${i}/plan1.id | sort | uniq -d # 无重复
    cat ${i}/plan1.id | sort | uniq > ${i}/plan1.id.nr
    wc -l ${i}/plan1.id.nr # 67个非冗余
    usearch10 -fastx_getseqs /mnt/bai/yongxin/culture/rice/stock/sequenceV5-7.fa -labels ${i}/plan1.id.nr -fastaout ${i}/plan1.fa
    makeblastdb -in ${i}/plan1.fa -dbtype nucl
    blastn -query ${i}/plan1.fa -db ${i}/plan1.fa -out ${i}/plan1.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 100 -evalue 1 -num_threads 9
    grep 100.00 ${i}/plan1.blastn|less|wc -l # 153重复
    # 结果用wet/blast_self_heatmap.Rmd可视化