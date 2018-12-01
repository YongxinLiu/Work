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

	# 快速分析 Quick Start(所需文件准备好)
	make fq_qc # 样本拆分、合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
	make host_rm # 序列去冗余、去噪、挑选代表序列、去嵌合、去宿主，获得OTU代表序列result/otu.fa
	make beta_calc # 生成OTU表、过滤、物种注释、建库和多样性统计
	# 清除统计绘图标记(重分析时使用)
	rm -rf alpha_boxplot 
	make DA_compare # 绘制alpha、beta、taxonomy和差异OTU比较
	#rm -f plot_volcano # 删除OTU差异比较可化标记
	make plot_manhattan # 绘制差异比较的火山图、热图、曼哈顿图
	make plot_venn # 绘制OTU差异共有/特有维恩图
	make DA_compare_tax # 高分类级差异比较，维恩图绘制，2为reads count负二项分布统计
	make rmd # 生成网页报告，必须依赖的只有alpha, beta, taxonomy



# 1. 处理序列 Processing sequences

	# 0. 准备工作 Preparation

	## 0.1 准备流程配置文件

	# 设置工作目录
	wd=maize/magic
	# 创建环境代码见~/github/Work/initial_project.sh

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt(含表头)
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/181024WangChao/"$4"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/181024WangChao/"$4"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz

	# 标准多文库实验设计拆分，保存模板中design页为doc/design_raw.txt
	split_design.pl -i doc/design_raw.txt
	# 从其它处复制实验设计
	# cp ~/ath/jt.HuangAC/batch3/doc/L*.txt doc/
	# 删除多余空格，windows换行符等
	sed -i 's/ //g;s/\r/\n/' doc/*.txt 
	head -n3 doc/L1.txt
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt | sed 's/#//g') <(cat doc/L* |grep -v '#'|grep -v -P '^SampleID\t') > doc/design.txt
	# 检查是否相等
	wc -l doc/design.txt
	cut -f 1 doc/design.txt|sort|uniq|wc -l

	## 准备原始数据

	# 拆lane和质量转换归为原始seq目录中处理
	# Prepare raw data
	#ln ~/seq/180210.lane9.ath3T/Clean/CWHPEPI00001683/lane_* ./
	#cp ~/ath/jt.HuangAC/batch3/doc/library.txt doc/
	
	# 检查数据质量，转换为33
	#determine_phred-score.pl seq/lane_1.fq.gz
	# 如果为64，改原始数据为33
	# rename 's/lane/lane_33/' seq/lane_*
	# 关闭质量控制，主要目的是格式转换64至33，不然usearch无法合并
	#time fastp -i seq/lane_64_1.fq.gz -I seq/lane_64_2.fq.gz \
	#	-o seq/lane_1.fq.gz -O seq/lane_2.fq.gz -6 -A -G -Q -L -w 9
	# 1lane 80GB, 2 threads, 102min

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
	# 按L1/2/3...txt拆分library为samples
	# 输入为seq/L*.fq，输出为seq/sample/*.fq
	make library_split
	make library_split_stat
	# 统计结果见result/split有txt/pdf/png，推荐看png方便快速查看每张位图
	# 查看样本量排序
	sort -k2,2n result/sample_split.log|less # 有11个小于5000，15个小于10000

## 1.3. 样品双端合并、重命名、合并为单一文件

	# Merge paired reads, renames and merge all samples
	# 样品双端合并、重命名、合并为单一文件, 注意fastq为33格式，64位采用fastp转换
	# 输入为seq/sample/*.fq，输出为seq/all.fq
	make sample_merge
	make sample_merge_stat
	# result/sample_merge.log中有每个样本合并后的序列数量


## 1.4. 切除引物与标签

	# Cut primers and lables
	# 切除左端标签和引物，右端 引物
	# Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
	# 输入为seq/all.fq，输出为temp/stripped.fq
	make fq_trim


## 1.5. 质量控制

	# Quality control
	# 过滤序列中预期累计错误率>1%的序列
	# 输入为temp/stripped.fq，输出为temp/filtered.fa
	make fq_qc



    # (第一阶段结束，获得纯净扩增子序列temp/filtered.fa，可提供此文件从下面开始)


## 1.6. 序列去冗余

	# Remove redundancy, get unique reads
	# 输入为temp/filtered.fa，输出为temp/uniques.fa
	make fa_unqiue


## 1.7. 挑选OTU

	# Pick OTUs
	# unoise3速度比cluster_otus慢上百倍，更精细但结果也更多
	# 输入为temp/uniques.fa，输出为temp/Zotus.fa
	make otu_pick


## 1.8. 有参去嵌合体

	# Remove chimiras by silva database
	# 基于SILVA数据库去除
	make chimera_ref


## 1.9. 去除宿主

	# Remove host
	# 根据SILVA注释去除线粒体、叶绿体、真核生物18S和未知序列(非rRNA)
	make host_rm


    # (第二阶段结束，获得OTU代表序列result/otu.fa，可提供此文件和测序数据temp/filtered.fa从下方起始)


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

# 3. 高级分析

## 3.9 培养菌注释

	# 默认为水稻，包括相似度、覆盖度、丰度和物种注释





# 4. 个性分析

## 4.1 基因型数据 

    mkdir -p snp/
    # snp/cubic_*共5个文件，hmp2plink代表hmp格式转换为plink格式，maf0.02为min allel frequency 0.02，bed/bim/ped为一组，分别为二进制的样本与SNP表，snp列表，样本对应的表型列表；_1404_Kinship.txt为kinship相似度矩阵，最相近即本身为2；_PopStructure.txt为群体结构数据，即前10个主成分
    # snp和样本数量，样式分别为chr1.s_717和MG_49
    wc -l snp/cubic_1404_hmp2plink_maf0.02.bim snp/cubic_1404_hmp2plink_maf0.02.fam # 1404个样品，11,825,030个SNP
    # 1404个样本的列表
    cut -f 2 snp/cubic_1404_hmp2plink_maf0.02.fam -d ' ' > snp/cubic_1404.id

    # 准备gemma分析要求格式
    # 链接至根目录并简化名称
    ln `pwd`/snp/cubic_1404_hmp2plink_maf0.02.b* ./
    rename 's/cubic_1404_hmp2plink_maf0.02/snp/' cubic_1404_hmp2plink_maf0.02.b*
    # 群体结构PC1-10与id对应
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' snp/cubic_PopStructure.txt snp/cubic_1404.id > snp.pca
    # kinship矩阵
    # 原始kinship去除行列名与gemma示例一致，下游报错矩阵是单数
    tail -n+2 snp/cubic_1404_Kinship.txt|cut -f 2-|less -S > snp.kin # gsl: lu.c:262: ERROR: matrix is singular
    # gemma生成kinship矩阵，第一次-o snp结果没有结果，后来发现输出至.output目录中，snp.cXx/log.txt，1h
    cd snp
    # 用原始数据生成kinship，error! number of analyzed individuals equals 0. ERROR: matrix dimension n1 must be positive integer，尝试修改空格为制作符仍报错，-9为1可运行
    sed -i 's/ /\t/g;s/-9/1/g' cubic_1404_hmp2plink_maf0.02.fam
    gemma -bfile cubic_1404_hmp2plink_maf0.02 -gk 1 -o snp 
    ## number of total individuals = 1404
    ## number of analyzed individuals = 1404
    ## number of covariates = 1
    ## number of phenotypes = 1
    ## number of total SNPs = 11825030
    ## number of analyzed SNPs = 8584966
    cd ..
    gemma -bfile snp -gk 1 -o snp
    cp output/snp.cXX.txt snp.kin
    mv output/snp.log.txt snp/
    rmdir output/

    # 准备gapit输入hmp格式
    cd snp
    # 转换bed为ped，30min, 还包括map, nosex, log文件
    time plink --bfile cubic_1404_hmp2plink_maf0.02 --out cubic_1404_hmp2plink_maf0.02 --recode
    # 转换ped为hmp，直接使用了60个线程，55m
    time ~/bin/TASSEL5/run_pipeline.pl -plink -ped  cubic_1404_hmp2plink_maf0.02.ped -map cubic_1404_hmp2plink_maf0.02.map  -export cubic_1404_hmp_maf0.02 -exportType Hapmap -Xmx50g
    cd ..
    


## 4.2 微生物组数据

    mkdir -p 16S/
    # 对应OTU表的编号与snp/cubic_1404_hmp2plink_maf0.02.fam第二列样本名一致
    # 对应实验设计中的 ID_MG 列，为26列
    head -n1 doc/design.txt | sed 's/\t/\n/g' | awk '{print NR"\t"$0}'
    # 查看编号唯一情况
    tail -n+2 doc/design.txt|cut -f 26|sort|uniq -u |wc -l # 唯一编号仅有1278
    tail -n+2 doc/design.txt|cut -f 26|sort|uniq -c|awk '{print $2"\t"$1}'|awk '$2>1'|sort -k2,2nr # 查看冗余情况
    tail -n+2 doc/design.txt|cut -f 26|sort|uniq -u > 16S/MG_unique.id # 保存非冗余编号
    cat 16S/MG_unique.id snp/cubic_1404.id|sort|uniq|wc -l|1415 # 仍有不一致的11个
    # 链接表型数据和SNP文件至根目录
    ln `pwd`/snp/cubic_1404_hmp2plink_maf0.02.b* ./
    rename 's/cubic_1404_hmp2plink_maf0.02/snp/' cubic_1404_hmp2plink_maf0.02.b*

    # 筛选otutab与fam列表对应样本
    # 提取样本名、组名和基因型对应ID
    cut -f 1,3,26 doc/design.txt |less > doc/design.txt3
    # 备份原始OTU表
    cp result/otutab.txt result/otutab.txt.181105
    # 对应OTU样本名至基因型
    otutab2genotype.Rmd
    # 重新计算物种和多样性
    usearch10 -otutab_stats result/otutab.txt -output result/otutab.stat
    usearch10 -otutab_norm result/otutab.txt -sample_size 10000 -output result/otutab_norm.txt 
    make tax_sum # 物种
    make beta_calc # alpha, beta

    # Alpha, beta, taxonomy，生成gemma输入fam格式性状
    Rscript ~/github/Amplicon/16Sv2/script/gemma_fam.R -i result/alpha/index.txt -d snp/cubic_1404_hmp2plink_maf0.02.fam
    # result/alpha/index.txt.fam 和 result/alpha/index.txt.fam.header
    



## 4.3 基因型与微生物组关联

    # gemma
    mkdir -p gemma
    cp result/alpha/index.txt.fam snp.fam
    # 样本顺序为字母顺序，不是fam中顺序
    i=alpha
	tail -n+6 result/alpha/index.txt.fam.header|head -n99|sed 's/[\(\)\"]//g'|awk '{print NR"\t"NR$0}'| sed "s/\t/\t${i}/" > gemma/${i}.list
    # 计算alpha多样性中第9列的richness关联
    gemma -bfile snp -k snp.kin -c snp.pca -lmm 4 -n 9 -o gemma/alpha9richness # le

    parallel -j 33 "gemma -bfile T2 -k gemma/kinship.txt -c snp.pca -lmm 4 -n {1} -o gemma/{1}" ::: `cut -f 1 gemma/${i}.list`

	parallel -j 33 "gemma -bfile snp/cubic_1404_hmp2plink_maf0.02 -k snp/cubic_1404_Kinship.txt -c snp/cubic_PopStructure.txt -lmm 4 -n {1} -o {1}" ::: `cut -f 1 ${dir}/HN/${i}.list`

    # gapit
    script/gapit.r # 使用内存超800G，手动终止

    # tassel
    # 转换ped格式为hmp
    time ~/bin/TASSEL5/run_pipeline.pl -plink -ped  cubic_1404_hmp2plink_maf0.02.ped -map cubic_1404_hmp2plink_maf0.02.map  -export cubic_1404_hmp_maf0.02 -exportType Hapmap -Xmx50g
    # 
    # 生成tassel格式
    mkdir -p tassel/16s
    cut -f 1,4,10,13 result/alpha/index.txt|sed '1 s/Sample/<Trait>/' > tassel/16s/alpha.txt
    cut -f 1,3  tassel/16s/alpha.txt >  tassel/16s/alpha_richness.txt
    # 运行tassel于alpha多样性，10G/30G内存溢出，改为200/600G
    rm -r tassel/alpha/
    out=tassel/alpha
    mkdir -p $out
    ~/software/tassel-5-standalone/run_pipeline.pl \
        -Xms200g -Xmx600g \
        -fork1 -h snp/cubic_1404_hmp_maf0.02.hmp.txt \
        -fork2 -r tassel/16s/alpha_richness.txt \
        -fork3 -q snp/cubic_PopStructure.txt \
        -fork4 -k snp/cubic_1404_Kinship.txt \
        -combine5 -input1 -input2 -input3 -intersect \
        -combine6 -input5 -input4 \
        -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
        -mlmOutputFile $out/mlmOut \
        -export $out/mlmExport \
        -runfork1 -runfork2 -runfork3 -runfork4 \
        > $out/out


# meta-gwas 在meta上分析GWAS
    
    [meta@meta:]scp -r yongxin@210.75.224.110:~/software/tassel-5-standalone/ ~/soft/tassel5
    screen -R gwas
    ssh yongxin@210.75.224.32
    cd ~/maize/magic/
    scp -r yongxin@210.75.224.110:~/maize/magic/snp ./
    scp -r yongxin@210.75.224.110:~/maize/magic/tassel ./
    scp -r yongxin@210.75.224.110:~/maize/magic/result ./


    cut -f 1,4,10,13 result/alpha/index.txt|sed '1 s/Sample/<Trait>/' > tassel/16s/alpha.txt
    cut -f 1,3  tassel/16s/alpha.txt >  tassel/16s/alpha_richness.txt
    # 运行tassel于alpha多样性，10G/30G内存溢出，改为200/600G
    rm -r tassel/alpha/
    out=tassel/alpha
    mkdir -p $out
    # 样本量和SNP太多，2天也运行不完
    /home/meta/soft/tassel5/run_pipeline.pl \
        -Xms200g -Xmx600g \
        -fork1 -h snp/cubic_1404_hmp_maf0.02.hmp.txt \
        -fork2 -r tassel/16s/alpha_richness.txt \
        -fork3 -q snp/cubic_PopStructure.txt \
        -fork4 -k snp/cubic_1404_Kinship.txt \
        -combine5 -input1 -input2 -input3 -intersect \
        -combine6 -input5 -input4 \
        -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
        -mlmOutputFile $out/mlmOut \
        -export $out/mlmExport \
        -runfork1 -runfork2 -runfork3 -runfork4 \
        > $out/out


### 拆分snp批量运行tassel

    cd ~/maize/magic/snp
    # 11M行，输出1位数字，每个1.2M SNP拆为10份，一夜tassel只运行20%，改为100份
    split -a 2 -d -l 120000 cubic_1404_hmp_maf0.02.hmp.txt split_hmp
    # 提取ID
    ls split_hmp??|cut -c10- > split_list.txt
    # 添加表头
    header=`head -n1 split_hmp00`
    for i in `tail -n+2 split_list.txt`;do sed -i "1 i $header" split_hmp$i; done
    cd ..

    out=tassel/alpha
    mkdir -p $out
    parallel --xapply -j 30 \
    "/home/meta/soft/tassel5/run_pipeline.pl \
        -Xms10g -Xmx30g \
        -fork1 -h snp/split_hmp{1} \
        -fork2 -r tassel/16s/alpha_richness.txt \
        -fork3 -q snp/cubic_PopStructure.txt \
        -fork4 -k snp/cubic_1404_Kinship.txt \
        -combine5 -input1 -input2 -input3 -intersect \
        -combine6 -input5 -input4 \
        -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
        -mlmOutputFile $out/mlmOut{1} \
        -export $out/mlmExport{1} \
        -runfork1 -runfork2 -runfork3 -runfork4 \
        > $out/out{1}" \
        ::: `cat snp/split_list.txt`
    # 类似的结果文件 tassel/alpha/mlmOut00_split_hmp00_+_alpha_richness_+_cubic_PopStructure_stats.txt


## Reference

### Tassel参考流程

TRAIT=$1
ID=$2

perl ~/software/tassel3-standalone/run_pipeline.pl \
	-Xms1g -Xmx3g \
	-fork1 -h ~/hjliu/cubic/geno/cubic_1404_maf0.02_ID$ID.hmp.gz \
	-fork2 -r ~/hjliu/cubic/AllTrait/${TRAIT}.txt \
	-fork3 -q ~/hjliu/cubic/GWAS/cubic_PopStructure.txt \
	-fork4 -k ~/hjliu/cubic/GWAS/cubic_1404_Kinship.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
	-mlmOutputFile mlmOut_cubic_${TRAIT}_ID${ID} \
	-export mlmExport_${TRAIT}_ID${ID} \
	-runfork1 -runfork2 -runfork3 -runfork4 \
	> out_${TRAIT}_ID${ID}

### emma参考流程

#!/bin/sh

i=$1

~/software/emmax-beta-07Mar2010/emmax -v -d 10 -t genotype/Plink/cubic_391_maf0.05.plk -p phenotype/${i}.rld.txt -k genotype/cubic_391_Kinship_emmax.txt -c genotype/cubic_PopStructure_391_peer_emmax.txt.bak -o GWAS/${i}.rld.peer
gzip GWAS/${i}.rld.peer.ps
