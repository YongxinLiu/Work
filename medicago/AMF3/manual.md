
	# 快速分析 Quick Start(所需文件准备好)
	find . -name "*" -type f -size 0c | xargs -n 1 rm -f # 清理零字节文件，用于从头重新分析项目清空makefile点位文件
	make library_split # 样本拆分、
	make fq_qc # 合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
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
	wd=medicago/AMF3
	# 创建环境代码见~/github/Work/initial_project.sh

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt
	#sed -i 's/\t/\tL171121_/' doc/library.txt # time check SeqLibraryList.xlsx
    cp ../AMF2/doc/library.txt doc/
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz

	# 标准多文库实验设计拆分，保存模板中design页为doc/design_raw.txt

    dos2unix doc/*.txt
	split_design.pl -i doc/design_raw.txt
	# 从其它处复制实验设计
	# cp ~/ath/jt.HuangAC/batch3/doc/L*.txt doc/
	# 删除多余空格，windows换行符等(苹果用户勿用)
	sed -i 's/ //g;s/\r//' doc/*.txt 
	head -n3 doc/L1.txt
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt | sed 's/#//g') <(cat doc/L* |grep -v '#'|grep -v -P '^SampleID\t') > doc/design.txt
	# 检查是否相等
	wc -l doc/design.txt
	cut -f 1 doc/design.txt|sort|uniq|wc -l
	cut -f 1 doc/design.txt|sort|uniq -d


## 1.2. 按实验设计拆分文库为样品

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
	# Cut barcode 10bp + ITS1F 22bp in left， and ITS2 20bp in right
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
	# 输入为temp/filtered.fa，输出为temp/uniques.fa，根据上一步结果选择整数，此处为97
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

## 1.17. 有参考构建OTU表-非16S不可用

	# Reference based OTU table
	# otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
	make otutab_gg



# 2. 统计绘图 Statistics and plot

## 2.1. Alpha多样性指数箱线图
	
	# Alpha index in boxplot
    rm -rf alpha_boxplot
    make alpha_boxplot

## 2.2. Alpha丰富度稀释曲线
	
	# Alpha rarefracation curve
    rm -rf alpha_rare
	make alpha_rare

## 2.3. 主坐标轴分析距离矩阵
	
	# PCoA of distance matrix
    rm -rf beta_pcoa
	make beta_pcoa

## 2.4. 限制性主坐标轴分析

	# Constrained PCoA / CCA of bray distance matrix
	# OTU表基于bray距离和CCA，至少3个组 
    rm -rf beta_cpcoa
	make beta_cpcoa

## 2.5. 样品和组各级分类学堆叠柱状图

	# Stackplot showing taxonomy in each level
    rm -rf tax_stackplot
	make tax_stackplot

## 2.6. 组间差异比较 
	
	# Group compareing by edgeR, wilcox or ttest
	# 可选负二项分布GLM、Wilcoxon秩和检验或T检验
    rm -rf DA_compare
	make DA_compare
	make DA_compare_tax
	make plot_volcano
	make plot_heatmap
	make plot_manhattan

## 2.7 绘制维恩图和生成报告
	make plot_venn # 绘制OTU差异共有/特有维恩图
	make rmd # 生成网页报告，必须依赖的只有alpha, beta, taxonomy





# 3. 高级分析

## 3.9 培养菌注释

	# 默认为水稻，包括相似度、覆盖度、丰度和物种注释，请修改参数处菌库位置和注释文件
	make culture

# 4. 个性分析

## 4.1. 分蘖与菌相关性

	# 准备相关输入文件
	cd ~/rice/miniCore/180718
	# 硬链数据文件，保持可同步修改和可备份
	# miniCore分蘖数据整理
	ln ~/rice/xianGeng/doc/phenotype_sample_raw.txt doc/
	# LN otu表和实验设计
	mkdir -p data
	cp ~/rice/miniCore/180319/LN/otutab.txt data/LN_otutab.txt
	cp ~/rice/miniCore/180319/doc/design.txt doc/design_miniCore.txt
	mkdir -p data/cor/LN
	# 物种注释
	cp ~/rice/miniCore/180319/temp/otus_no_host.tax data/

	# 统计见script/cor_tiller_LN.Rmd
	# 相关系数，添加物种注释
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4} NR>FNR{print $0,a[$1]}' result/otus_no_host.tax data/cor/LN/otu_mean_pheno_cor.r.txt | less -S > result/cor/LN/otu_mean_pheno_cor.r.txt.tax
	# 再添加可培养相关菌
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result/39culture/otu.txt data/cor/LN/otu_mean_pheno_cor.r.txt.tax | less -S > data/cor/LN/otu_mean_pheno_cor.r.txt.tax



# 5. 发表图版

    cd ~/medicago/AMF3
    mkdir -p fig && cd fig
    mkdir -p fig1 fig2 fig3 fig4 data script && cd ..
    # 数据筛选 script/DataFilter.Rmd ## 筛选发表使用样本，otu表，alpha，beta, taxonomy
    # OTU和物种注释不变
    cp result/otu.fa fig/data/
    cp result/taxonomy_* fig/data/

## 1. 样本描述description fig/fig1/fig.Rmd

    cp ~/medicago/AMF2/fig/fig1/fig1.Rmd fig/fig1/

## 2. 差异比较compare fig2.Rmd

    cp /mnt/bai/yongxin/medicago/AMF2/fig/fig2/fig2.Rmd fig/fig2/
    ### 差异比较曼哈顿图(2/3批混合筛选点后6个基因型4组wilcox比较) http://210.75.224.110/report/16Sv2/med_AMF3_b23r6_v1
    cp med_AMF3_b23r6_v1/result/compare/*_all.txt fig/data/

### 绘制单菌的丰度图 Pseudomonadaceae 和 Bacillus
    # 前10个菌有5个重点关注， Bacillus 2,5; Pseudomonas 3
    # 绘制OTU2,5 Bacillus
    alpha_boxplot.sh -i fig/data/otutab.txt -d fig/data/design.txt -A genocomp -B '"A17R","AnfpR","dmi3R","R108R","lyk9R","lyr4R","A17Rs","AnfpRs","dmi3Rs","R108Rs","lyk9Rs","lyr4Rs","soilS"' -m '"OTU_2"' -t TRUE -o fig/fig2/boxplot_ -n TRUE -U 100
    # 绘制Bacillus属
    alpha_boxplot.sh -i fig/data/sum_g.txt -d fig/data/design.txt -A genocomp -B '"A17R","AnfpR","dmi3R","R108R","lyk9R","lyr4R","A17Rs","AnfpRs","dmi3Rs","R108Rs","lyk9Rs","lyr4Rs","soilS"' -m '"Bacillus"' -t TRUE -o fig/fig2/boxplot_ -n TRUE -U 100

## 3. 分菌

### 3.1 自然样品与分菌结果比较

    # 与项目中A17/R180野生型根样本与总体菌库比较
    # 修改makefile中3.9部分，参考~/medicago/culture_start/makefile
    # 注意与doc/design.txt中列和组对应，注意与分菌中文件名和品种名对应
    conda deactivate
    rm culture
    make culture
    rm culture_graphlan
    make culture_graphlan
    # A17 0.001只有58个ASV太少，改为0.0005有88个


    # 3.9 culture和culture_graphlan，只需修改A17r或R108r，基于实验中大量样本的graphlan
    # 2019/3/25 与使用culture_start自然样品 ~/medicago/culture_start # culture_graphlan
    # 2019/4/9 更新 pipeline.md 中分菌部分详细注释和代码优化
    # 统计分菌代码见文末“## Sanger测序验证菌”段落，sanger测序绘制代码参考~/culture/rice/makefile.man "# 总结：1098个单菌" 段落
    cd ~/medicago/AMF3/fig/fig3
    cp ~/culture/rice/verify/graphlan_culture.R ./
    cp ~/medicago/AMF3/wet/stock_sanger/taxonomy_9.txt ./
 
    Rscript graphlan_culture.R # 生成1树, 2科注释，和来原环
    cat /mnt/bai/yongxin/culture/rice/graphlan/global.cfg 2_annotation_family.txt /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg 3_annotation_match.txt > 5_annotation.txt
    graphlan_annotate.py --annot 5_annotation.txt 1_tree_plain.txt graphlan.xml
    graphlan.py graphlan.xml culture_graphlan.pdf --size 5
    # 提取相应序列
    tail -n+2 taxonomy_9.txt|wc -l # 325 seqs
    cut -f 1 taxonomy_9.txt|tail -n+2 > seq325.id
    usearch10 -fastx_getseqs ~/medicago/AMF3/wet/stock_sanger/16s_full_length_list1.fa -labels seq325.id -fastaout seq325.fa

    # 2019/7/2 与最新10个菌库比较
    cwd=culture10
    mkdir -p ${cwd}
    blastn -query ~/medicago/AMF3/result/otu.fa -db ~/culture/medicago/190626/result/culture_select.fa -out ${cwd}/otu_culture.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 # 输出13列为coverage
    awk '$3*$13>=9700' ${cwd}/otu_culture.blastn |cut -f 1 > ${cwd}/otu_cultured.txt # 197 OTU cultured, 197/535=36.8%
    # 我们关注的菌OTU_2 k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_1;g__Bacillus，在新库中${cwd}/otu_culture.blastn 对应 COTU_61
    # L8-10为新库，精选信息可查culture_select.xls中编号61，L8,L10中为最优解；更多孔详见 culture_bacteria.xls
    # 添加孔的信息
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$8"\t"$9"\t"$10} NR>FNR {print $0,a[$1]}' doc/design.txt result/culture_bacteria.xls > result/culture_bacteria_anno.xls

    # awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]="culture"} NR>FNR {print $0,a[$1]}' temp/otu_cultured.txt ${nature}/otu_table.txt.mean |grep 'culture'|awk '{a=a+$2} END {print a}' # 可培养的丰度 0.847 




# 附录
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


# 常见问题

## alpha_rare 出错，第一行不完整列，可删除
	sed -i '/-/d' result/alpha/rare.txt
