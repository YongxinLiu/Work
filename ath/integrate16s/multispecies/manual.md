
	# 快速分析 Quick Start(所需文件准备好)
	make library_split # 样本拆分、
	make fq_qc # 合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
	make host_rm # 序列去冗余、去噪、挑选代表序列、去嵌合、去宿主，获得OTU代表序列result/otu.fa
	make beta_calc # 生成OTU表、过滤、物种注释、建库和多样性统计
	# 清除统计绘图标记(重分析时使用)
	rm -rf alpha_boxplot 
    rm DA_compare
	make DA_compare # 绘制alpha、beta、taxonomy和差异OTU比较
    cat result/compare/summary.txt
	#rm -f plot_volcano # 删除OTU差异比较可化标记
	make plot_manhattan # 绘制差异比较的火山图、热图、曼哈顿图
	rm plot_venn
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
	wd=ath/integrate16s/multispecies
	# 创建环境代码见~/github/Work/initial_project.sh

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/181009.lane14/"$4"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/181009.lane14/"$4"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz

	# 标准多文库实验设计拆分，保存模板中design页为doc/design_raw.txt
	split_design.pl -i doc/design_raw.txt
	# 从其它处复制实验设计
	cp ~/ath/jt.HuangAC/batch3/doc/L*.txt doc/
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
	# 数据准备CP多物种数据： ath3T + Ler + rice + wheat
    # 1. 拟南芥三萜最新版 ~/ath/integrate16s/ 包括2.5和3T，引用原为~/ath/jt.HuangAC/batch3all/
    # 1.1 三萜和二半萜
    id=athT
    wd=ath/integrate16s
    cp ~/${wd}/doc/design.txt doc/design_${id}.txt
    ln ~/${wd}/temp/filtered.fa temp/seq_${id}.fa
    # 1.2 傅向东Gibberellin突变体gai
    id=athGA
    wd=ath/jt.FuXD
    sed 's/\.//g' ~/${wd}/doc/design.txt > doc/design_${id}.txt
    sed 's/\.//g' ~/${wd}/temp/seqs_usearch.fa | cut -f 1 -d ';'|sed 's/_/\./g'|less > temp/seq_${id}.fa
    # 1.3 SA 
    id=athSA
    wd=ath/SA
    cp ~/${wd}/doc/design.txt doc/design_${id}.txt
    ln ~/${wd}/temp/filtered.fa temp/seq_${id}.fa

    # 2. 水稻
    # 2.1 水稻时间序列CP ~/rice/timecourse/v2
    id=riceTC
    wd=rice/timecourse/v2
    cp ~/${wd}/doc/design.txt doc/design_${id}.txt
    ln ~/${wd}/temp/filtered.fa temp/seq_${id}.fa
    cut -f 5 doc/design_${id}.txt |grep A50Cp|uniq # 统计组类别
    # 3. 小麦
    # 3.1 小麦2016 CP 时间序列
    id=wheatTC
    wd=wheat/profile
    cp ~/${wd}/doc/design.txt doc/design_${id}.txt
    ln ~/${wd}/temp/seqs_usearch.fa temp/seq_${id}.fa
    cut -f 5 doc/design_${id}.txt |grep -v "soil"|grep 'XY54'|grep 'L1' | uniq # 统计组类别


    # 合并
    cat doc/design_* > doc/design_raw.txt # 手动修改
    cut -f 1 doc/design.txt |sort|uniq -d
    cat temp/seq_* > temp/filtered.fa


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

# 4.1 6个基因型共用，且野生型>0.01%的OTU
# 筛选2.5T中6个基因型中共有
cat result/compare/TPS*-WT4_sig.txt result/compare/DM*-WT4_sig.txt > temp/ath2t_sig.txt
grep 'Enriched' temp/ath2t_sig.txt | cut -f 1 |sort|uniq| sed 's/$/\tAll2T_E/' >>result/compare/diff.list
grep 'Depleted' temp/ath2t_sig.txt | cut -f 1 |sort|uniq| sed 's/$/\tAll2T_D/' >>result/compare/diff.list
sed -i '/^\t/d' result/compare/diff.list
grep 'All2T' result/compare/diff.list|cut -f 1 |sort|uniq -d # 有10个OTU在上调和下调中有都？

# 筛选3T中4个基因型中共有
cat result/compare/b3*_sig.txt > temp/ath3t_sig.txt
grep 'Enriched' temp/ath3t_sig.txt | cut -f 1 |sort|uniq| sed 's/$/\tAll3T_E/' >>result/compare/diff.list
grep 'Depleted' temp/ath3t_sig.txt | cut -f 1 |sort|uniq| sed 's/$/\tAll3T_D/' >>result/compare/diff.list
sed -i '/^\t/d' result/compare/diff.list
grep 'All3T' result/compare/diff.list|cut -f 1 |sort|uniq -d 

# 筛选野生型WT4中高丰度
wc -l result/compare/database.txt
# awk '$14>0.01' result/compare/database.txt > result/compare/athCol.txt # |wc -l # 1834中的798个
# 筛选WT4和b3Col中共有的
wc -l result/compare/database.txt # 2038个共有,1017在两个野生型在0.01%
awk '$3>0.01 || $19>0.01' result/compare/database.txt > result/compare/athCol.txt
awk '{a=a+$3} END{print a}' result/compare/athCol.txt # $3-b3Col 73.6481
awk '{a=a+$19} END{print a}' result/compare/athCol.txt # $19-WT4 71.6073
cp result/compare/diff.list result/compare/diff.list.bak
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=1} NR>FNR{print $0,a[$1]}' result/compare/athCol.txt result/compare/diff.list | awk '$3==1' | less | grep -P 'RiceCp35_WT4|WheatD35L1_WT4|All2T|All3T' > result/compare/diff.list.col.txt
cp result/compare/diff.list.col.txt result/compare/diff.list
cp result/compare/database.txt result/compare/database.txt.bak
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' result/taxonomy_2.txt result/compare/database.txt.bak | cut -f 1-7,10,13 > result/compare/database.txt


# 手动维恩分析筛选后的组，指定了新的-d差异列表，仍然还指向旧文件
batch_venn.pl -i `pwd`/doc/""/venn.txt -d result/compare/diff.list
rm -f result/compare/*.xls.xls
batch2.pl -i 'result/compare/diff.list.venn*.xls' -d result/compare/database.txt -o result/venn/ -p vennNumAnno.pl
make rmd
# 2.5T和3T比较
http://210.75.224.110/report/16Sv2/multispecies_evolve_2.5T_v2


# 2019/2/8 土壤三点比较，35d rice ath and wheat
# 三萜中的土壤为b3BS，但分为；水稻的土壤0，42，49，35天没取样，暂用42天；小麦取样封为soilD35L1："b3BS","soilCp42","soilD35L1"
# 放在soil/compare.txt中按标准流程分析，修改sub, version
rm DA_compare
make DA_compare
# 筛选各组中存在的的ASV做维恩 
awk '$2>0.0' result/compare/database.txt|grep 'OTU_'|cut -f 1 |sort|uniq| sed 's/$/\tAth3TSoil/' >>result/compare/diff.list
awk '$3>0.0' result/compare/database.txt|grep 'OTU_'|cut -f 1 |sort|uniq| sed 's/$/\tRiceSoil/' >>result/compare/diff.list
awk '$4>0.0' result/compare/database.txt|grep 'OTU_'|cut -f 1 |sort|uniq| sed 's/$/\tWheatSoil/' >>result/compare/diff.list
sed -i '/^\t/d' result/compare/diff.list



# 2019/2/11 三物种比较，采用Cp35天比较，来自~/ath/jt.HuangAC/coevolve
mkdir -p doc/coevolve
cp ~/ath/jt.HuangAC/coevolve/doc/design.txt doc/coevolve/
cp ~/ath/jt.HuangAC/coevolve/doc/compare.txt doc/coevolve/
cp ~/ath/jt.HuangAC/coevolve/doc/venn.txt doc/coevolve/
# 2019/2/22 cp doc/coevolve/compare.txt doc/coevolve/compare.txt.bak
# 修改sub版本为 coevolve, 实验组为groupID2
make DA_compare
make plot_venn
# 添加4个基因型比较共有结果
grep 'OTU' result/compare/diff.list.vennb3ACT2KO_b3Col_Db3ThadKO_b3Col_Db3ThahKO_b3Col_Db3ThasKO_b3Col_D.xls.xls | cut -f 1 | sed 's/$/\tAll3T_D/' >>result/compare/diff.list
grep 'OTU' result/compare/diff.list.vennb3ACT2KO_b3Col_Eb3ThadKO_b3Col_Eb3ThahKO_b3Col_Eb3ThasKO_b3Col_E.xls.xls | cut -f 1 | sed 's/$/\tAll3T_E/' >>result/compare/diff.list
sed -i '/^\t/d' result/compare/diff.list
## 筛选拟南芥中高丰度的OTU
#cp ~/ath/jt.HuangAC/coevolve/script/compare_athCol.R script/compare_athCol.R
## script/compare_athCol.R 输出 result/compare/athCol.txt，共768个OTU，新版本为773个
#awk '{a=a+$3} END{print a}' result/compare/athCol.txt # 72.31%
## 筛选diff.list只在拟南芥中高表达的
#awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=1} NR>FNR{print $0,a[$1]}' result/compare/athCol.txt result/compare/diff.list | awk '$3==1' | less | grep -P 'RiceCp35_b3Col|WheatD35L1_b3Col|All3T' > result/compare/diff.list.col.txt
#mv result/compare/diff.list result/compare/diff.list190211
#cp result/compare/diff.list.col.txt result/compare/diff.list
## 在线Venn比较，看到了显著的趋势，本地绘制并注释OTU
#cp result/compare/database.txt result/compare/database.txt.bak
#awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' result/taxonomy_2.txt result/compare/database.txt.bak | cut -f 1-7,10,13 > result/compare/database.txt
# 重新绘图
rm plot_venn
make plot_venn
# 可视化
make rmd 
http://210.75.224.110/report/16Sv2/multispecies_evolve_coevolve_v2


## Fig5. C/E 单菌对应OTU的箱线图
    # 姜婷发我的实验的菌测序结果 /mnt/bai/yongxin/ath/integrate16s/multispecies/3T/culture.txt
    grep -P '^A' 3T/culture.txt|sed 's/^A/>A/'|sed 's/\t/\n/'| less -S> 3T/culture.fa
    makeblastdb -in 3T/culture.fa -dbtype nucl

    # OTU太多，只筛选用于差异比较的高丰度OTU
    cut -f 1 multispecies_evolve_3.1_v1/result/compare/b3ThasKO1-b3Col_all.txt|tail -n+2 > 3T/otu.id
    usearch10 -fastx_getseqs result/otu.fa -labels 3T/otu.id -fastaout 3T/otu.fa
    # 比对otu至单菌序列
    blastn -query 3T/otu.fa -db 3T/culture.fa -out temp/3T.otu.blast \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' \
        -num_alignments 1 -evalue 1 -num_threads 9
    # 筛选>97%相似的OTU
    rm -r 3T/otu_box
    mkdir -p 3T/otu_box
    awk '$3>97' temp/3T.otu.blast|sort -k2,2 -k3,3nr > 3T/otu_box/3T.otu.blast.txt
    # 绘制每个OTU的箱线图，需要标准化+转置 -t -n TRUE; 修改compare为compare_soil加土，调整顺序为："b3BS","b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO"
    alpha_boxplot.sh -i result/otutab.txt -m `cut -f 1 3T/otu_box/3T.otu.blast.txt|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'
` \
        -d `pwd`/doc/"3.1"/design.txt  -A groupID -B '"b3BS","b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO"' \
        -o `pwd`/3T/otu_box/ -h 5 -w 8 -t TRUE -n TRUE
    # 查看A388和A224的100%，为OTU_309 OTU_231

## 附图22. Alpha多样性
	# 绘制8个基因型的alpha多样性
    mkdir -p fig/S22
	alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
        -d `pwd`/doc/"3.1"/design.txt  -A groupID -B '"b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR","b3ACT2KO","b3BS"' \
        -o `pwd`/fig/S22/ -h 5 -w 8

## 附图23. 土壤中富集的菌 
	# 编辑doc/3.1/compare_soil.txt 带土，只比较土
    mkdir -p fig/S23
	compare_OTU.sh -i `pwd`/result/otutab.txt -c doc/3.1/compare_soilonly.txt -m "edgeR" \
        -p 0.05 -q 0.2 -F 1.2 -t 0.01 \
        -d `pwd`/doc/"3.1"/design.txt  -A groupID -B `cat doc/"3.1"/compare_soil.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'` \
        -o `pwd`/result/compare/ -C groupID
	# 正常7个基因型才916个OTU，加土为2656个
	cp result/compare/b3Col-b3BS_all.txt fig/S23/b3Col-b3BS_all.txt # 备份最终结果
	# 注释两类OTU
	cut -f 6 fig/S23/b3Col-b3BS_all.txt|sort|uniq -c
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print $0,a[$1]}' fig/S23/b3Col-b3BS_all.txt multispecies_evolve_3.1_v1/result/compare/diff.list.vennThas_Eb3ThahKO_b3Col_Eb3ThadKO_b3Col_EACT2_E.xls.xls > fig/S23/veen.enriched.soil.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print $0,a[$1]}' fig/S23/b3Col-b3BS_all.txt multispecies_evolve_3.1_v1/result/compare/diff.list.vennThas_Db3ThahKO_b3Col_Db3ThadKO_b3Col_DACT2_D.xls.xls > fig/S23/veen.depleted.soil.txt
	# 提取Overlap ID
	tail -n `grep -A1 'specific_to_others' multispecies_evolve_3.1_v1/result/compare/diff.list.vennThas_Eb3ThahKO_b3Col_Eb3ThadKO_b3Col_EACT2_E.xls.xls|tail -n1` multispecies_evolve_3.1_v1/result/compare/diff.list.vennThas_Eb3ThahKO_b3Col_Eb3ThadKO_b3Col_EACT2_E.xls.xls | awk '{print $1"\tAll_E"}' > fig/S23/all_common.txt
	tail -n `grep -A1 'specific_to_others' multispecies_evolve_3.1_v1/result/compare/diff.list.vennThas_Db3ThahKO_b3Col_Db3ThadKO_b3Col_DACT2_D.xls.xls|tail -n1` multispecies_evolve_3.1_v1/result/compare/diff.list.vennThas_Db3ThahKO_b3Col_Db3ThadKO_b3Col_DACT2_D.xls.xls | awk '{print $1"\tAll_D"}' >> fig/S23/all_common.txt


## 附图24. 热图绘制不聚类
	# 绘制多个比较组，但图例不一致，保证着色一致尝试注释修改level也无效
    mkdir -p fig/S24
	awk 'BEGIN{OFS=FS="\t"}{system("plot_heatmap.sh -i multispecies_evolve_3.1_v1/result/compare/"$1"-"$2"_sig.txt \
        -o fig/S24/"$1"-"$2" -w 7 -h 5 -C FALSE");}' \
        <(grep -v '^$' `pwd`/doc/"3.1"/compare.txt)

## 2019/2/21 重新提取 design, otutab, otu.fa和taxonomy
# fig/index.Rmd ### Metadata, OTU table and taxonomy

# 2019/2/26 数据上传基因组所，拟南芥，水稻，小麦
# 样品列表见fig/table/metadata.txt
mkdir -p GSA
wc -l fig/table/metadata.txt
grep 'ath3t' fig/table/metadata.txt|wc -l
grep 'riceTC' fig/table/metadata.txt|wc -l
grep 'wheat' fig/table/metadata.txt|wc -l
# 共140个样品，分别为108，12和20
	# 1. 拟南芥3t
    # 找原始数据没找到，找上次提交NCBI目录也没找到，找到原始数据来自lane9
    # 从头目录中分析
	for i in `grep 'minicore' fig1/ST/02.design.txt|cut -f 1`; do
        ln /mnt/bai/yongxin/rice/miniCore/clean_data/sample/${i}.fq.gz seq/submitGSA/ ;done

	# 2. 时间序列己上传，见 https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA435900
	for i in `grep 'nrt' fig1/ST/02.design.txt|cut -f 1`; do
        ln /mnt/bai/yongxin/rice/zjj.nitrogen/180116/clean_data/sample/${i}.fq.gz seq/submitGSA/ ; done
    md5sum seq/submitGSA/* > seq/submitGSA_md5.txt
    sed -i 's/seq\/submitGSA\///' seq/submitGSA_md5.txt
    sed 's/.fq.gz//' seq/submitGSA_md5.txt > seq/temp.txt
    paste seq/temp.txt seq/submitGSA_md5.txt |sed 's/  /\t/g'| cut -f 2-4| less > seq/temp1.txt
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' seq/temp1.txt fig1/metadata.txt > fig1/metadata_md5.txt 

