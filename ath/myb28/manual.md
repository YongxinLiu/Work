
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
	#rm -f plot_volcano # 删除OTU差异比较可化标记
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

# 1. 处理序列 Processing sequences

	# 0. 准备工作 Preparation

	## 0.1 准备流程配置文件

	# 设置工作目录
	wd=ath/myb28
	# 创建环境代码见~/github/Work/initial_project.sh

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/L181122/L181122_"$2"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/L181122/L181122_"$2"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz

	# 标准多文库实验设计拆分，保存模板中design页为doc/design_raw.txt
	split_design.pl -i doc/design_raw.txt
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
	rename 's/lane/lane_33/' seq/lane_*
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
	head -n3 doc/L1.txt
	# 按L1/2/3...txt拆分library为samples
	# 输入为seq/L*.fq，输出为seq/sample/*.fq
	make library_split
	make library_split_stat
	# 统计结果见result/split有txt/pdf/png，推荐看png方便快速查看每张位图
	# 查看样本量排序
	sort -k2,2n result/sample_split.log|less

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

	# 默认为水稻，包括相似度、覆盖度、丰度和物种注释，请修改参数处菌库位置和注释文件
	make culture

# 4. 个性分析

## 4.1

结果列表：http://210.75.224.110/report/16Sv2/ath_myb__v1 # 关闭了叶、线叶体过滤，所有样本
http://210.75.224.110/report/16Sv2/ath_myb__v1 # 过滤OTU_1叶绿体，OTU表过滤>2000 reads

    # 缺少样本，改为1000条两种方法分析
http://210.75.224.110/report/16Sv2/ath_myb_edgeR_v2
http://210.75.224.110/report/16Sv2/ath_myb_wilcox_v2

    # 过滤所有宿主，1000条抽样
http://210.75.224.110/report/16Sv2/ath_myb_wilcox_v3
http://210.75.224.110/report/16Sv2/ath_myb_edgeR_v3


## 文中统计

grep -P 'ColTot1|MybTot1|ColPhyl2|MybPhyl2' doc/design.txt | cut -f 3 | uniq -c # 每个组12个样
grep -P 'ColTot1|MybTot1|ColPhyl2|MybPhyl2' doc/design.txt | grep -v '#'| cut -f 3 | uniq -c # 每个组样品数
grep -P 'ColTot1|MybTot1|ColPhyl2|MybPhyl2' doc/design.txt | grep -v '#' > doc/filter_metadata.tsv
wc -l doc/filter_metadata.tsv # 45个样品
dos2unix doc/*
# 样品测序数据量
awk 'BEGIN{FS=OFS="\t"} ARGIND==1{a[$1]=$2} ARGIND==2{print $1,a[$1]}' result/sample_split.log doc/filter_metadata.tsv | awk '{a=a+$2} END {print a}' # 4119595
awk 'BEGIN{FS=OFS="\t"} ARGIND==1{a[$1]=$2} ARGIND==2{print $1,a[$1]}' result/sample_split.log doc/filter_metadata.tsv | sort -k2,2nr
# OTU表中数据总量
sed -i 's/: /\t/' result/otutab.biom.sum
awk 'BEGIN{FS=OFS="\t"} ARGIND==1{a[$1]=$2} ARGIND==2{print $1,a[$1]}' result/otutab.biom.sum doc/filter_metadata.tsv | awk '{a=a+$2} END {print a}'


# 数据上传 2019/7/5

    # 扩增子 根据筛选后的整理数据并上传
    mkdir -p ~/ath/myb28/submit/16s
    awk '{system("ln seq/sample/"$1"_1.fq.gz submit/16s/")}' doc/filter_metadata.tsv
    awk '{system("ln seq/sample/"$1"_2.fq.gz submit/16s/")}' doc/filter_metadata.tsv

    cd submit/16s
    md5sum *_1.fq.gz > /tmp/md5sum1.txt
    md5sum *_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt

    # 细菌转录组
    cd ~/ath/myb28/submit/RNA-Seq/
    md5sum *_1.fq.gz > /tmp/md5sum1.txt
    md5sum *_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt



# 结果中要输出筛选后的实验设计和OTU表，用于统计和发表文章
# 把数据量、alpha多样性、beta多样性等属性追加到metadata.txt


# 问题汇总

## otutab_filter过滤小于5000后样本减少，但样本数据量大于5万

    # 可能原因为宿主叶绿体、线粒体含量过高，关掉过滤叶绿体
    grep -P 'Mitochondria|Chloroplast|Eukaryota' temp/otus_no_chimeras.tax > temp/host.tax # OTU1/9都是叶、线粒体
    # 关掉过滤，即修改 gg_silva为none，比对率也上升为90.1%

    # 样本过滤前后数据量
    less temp/otutab.biom.sum
    less result/otutab.biom.sum

## 删除OTU_1继续分析
    cp temp/otutab.txt temp/otutab.txt.bak
    sed '/OTU_1\t/d' temp/otutab.txt
    rm otutab_filter