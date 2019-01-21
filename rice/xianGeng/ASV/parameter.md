SHELL:=/bin/bash

	# 16S扩增子第二版流程配置文件，请根据项目具体情况修改 v2.0 2018/4/3
	# Config file of 16S Amplicon pipeline version 2, please modify according to experiment design. v2.0 2018/4/3


# 1. 标准流程参数 Parameters of standard pipeline

	## 工作目录 Working directory
	# 修改wd为当前工作目录pwd
	wd=`pwd`
	# make init # 建立分析所需子目录

	# 设置任务最大运行任务/线程数，超过CPU数量效率反而会降低
	p=64
	
	# 数据库
	# Greengene 13 May database, fa for usearch format, udb for usearch index
	usearch_gg=/mnt/bai/public/ref/gg_13_5_otus/97_otus_usearch.udb
	## Silva 132 database, fa for usearch format, udb for usearch index
	usearch_silva=/mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.udb
	usearch_rdp=/mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb

## 1.1. lane_split 拆分下机数据为文库

	# Split lane into library
	# 需要将lane文件放入seq目录，对应的index和文库放在doc/library.txt
	# 注意：务必查看文库文件中Index具体格式，默认为#Index，其它情况需修改主流程源代码main_pipeine.sh
	# 文库名
	lane=lane

## 1.2. library_split 拆分文库为样品

	# Split library into sample
	# 默认只拆分单左端barcode类型的样品:先匹配左端，再提取序列ID，再提取右端，最后改名，注意实验设计要严格规范无空格

## 1.3. sample_merge 双端序列合并

	# Merge pair-end reads
	# 如果是pair-end reads是phred64，需要先使用fastp转换为33且关闭质控(质控影响序列长度)，再使用usearch10 mergepair

## 1.4. fq_trim 切除引物和标签

	# Cut barcode 10bp + primer V5 19bp in left, and primer V7 18bp in right
	stripleft=29
	stripright=18

## 1.5. fq_qc 质量控制
	
	# fastq filter
	# 默认错误率<0.01 keep reads error rates less than 1%
	fastq_maxee_rate=0.01

## 1.6. fa_unqiue 序列去冗余

	# Remove redundancy
	# 最小序列频率miniuniqusize默认为8，去除低丰度，增加计算速度，整lane的序列可更改为30，甚至100
	minuniquesize=211
	# 100条件下，211147447 seqs, 25377669 uniques, 21184572 singletons (83.5%)，78706个，而OTU只有6000；
	# 改为211，1/M，有23446个唯一序列，14619个ASV

## 1.7. otu_pick 挑选OTU

	# Pick OTUs
	# 可选97% cluster_otus 和 unoise3 ，默认unoise3
	otu_method=unoise3
	# OTU日志文件，记录从挑选至过滤为最终的过程
	otu_log=result/otu.log

## 1.8. chimera_ref 参考去嵌合

	# Remove chimeras
	# 此处推荐使用大数据，如SILVA132，其它数据库如rdp_gold是错误的
	# /mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.fasta # 99%非冗余1.1G, 内存使用5.8G, 8min, 16.2% chimeras
	# /mnt/bai/public/ref/silva/SILVA_132_SSURef_tax_silva.fasta # 全部3.3G, 内存使用15.9G, 30min, 16.2% chimeras
	# 使用非冗余99%的省内存和时间，但结果差不多
	chimera_ref=${usearch_silva}
	# 模式，嵌合体比例由大到小: high_confidence specific balanced sensitive sensitive
	chimera_mode=balanced

## 1.9. host_rm 去宿主

	# Remove host original sequences
	# 去宿主方法选择 blast / sintax_gg / sintax_silva，推荐：sintax_silva
	host_method=sintax_silva
	# 方法1. blast宿主基因组(含叶绿体/线粒体)去除同源序列，如水稻微生物，需要提供水稻基因组；可调相似度和覆盖度的阈值(百分数)
	host=/mnt/bai/public/ref/rice/msu7/all.con
	host_similarity=90
	host_coverage=90
	# 方法2. 基于gg注释结果筛选, 去除叶绿体Chloroplast、线粒体mitochondria，默认为usearch_gg数据库
	# 方法3. silva注释可识线粒体Mitochondria、叶绿体Chloroplast和真核生物Eukaryota(包括宿主、真菌、原生动物等)，默认为usearch_silva数据库


## 1.10. otutab_create 生成OTU表

	# Creat OTUs table
	# 有 usearch10 和 vsearch 两个软件可选，默认usearch10，vsearch多线程会更快些
	map_method=vsearch
	map_identify=0.97

## 1.11. otutab_filter/norm OTU表筛选

	# Filter OTU table
	# OTU表筛选日志文件
	log_otutable=result/otutab.log
	# 按样本量筛选，默认5000，根据otu_stats结果调整
	min_sample_size=5000
	# 按矩阵中每个点count, freq筛选，低于阈值变为0
	# 按OTU丰度和频率筛选，如OTU测序量至少8次，相对丰度百万分之一(建议挑选序列去冗余部分调高阈值更合理)
	min_otu_size=8
	# 按频率筛选，推荐十万分之一0.00001，范围千一至百分一0.001 - 0.000001之间
	min_otu_freq=0.000001
	# 抽样标准化的值，推荐最小10000，根据统计结果选择筛选后最小值或可保留大部分样品的值
	sample_size=5000

## 1.12. tax_assign 物种注释

	# Assign taxonomy
	# 物种注释推荐使用小而准的数据库，如rdp trainset 16(由Robert整理)
	# 可选gg, rdp, silva，分别从官网下载并shell调整格式
	sintax_db=${usearch_rdp}
	# 分类准确度阈值，默认0.8，注释太少最小可改0.5，推荐0.6注释较完整，发现有明显错误可最高上升为0.95
	sintax_cutoff=0

## 1.13. tax_sum 物种注释统计

	# Taxonomy summary
	# 按门、纲、目、科、属水平分类汇总，结果位于：result/tax/sum_*.txt
	# 输出可读的taxonomy与OTU对应文件，有2列和8列版，result/taxonomy_2/8.txt

## 1.14. tree_make 多序列比对和进化树

	# Multiply alignment and make_phylogeny

## 1.15. alpha_calc Alpha多样性指数

	# Calculate alpha diversity index
	# alpha指数计算结果为 result/alpha/index.txt
	# 稀释梯度抽样方法 richness (observed OTUs)-method fast / with_replacement / without_replacement , 结果位于 result/alpha/rare.txt
	rare_method=without_replacement

## 1.16. beta_calc Beta多样性距离矩阵

	# Beta diversity tree and distance matrix
	# 距离矩阵计算方法，34种可选： abund_jaccard, binary_chisq, binary_chord, binary_euclidean, binary_hamming, binary_jaccard, binary_lennon, binary_ochiai, binary_otu_gain, binary_pearson, binary_sorensen_dice, bray_curtis, bray_curtis_faith, bray_curtis_magurran, canberra, chisq, chord, euclidean, gower, hellinger, kulczynski, manhattan, morisita_horn, pearson, soergel, spearman_approx, specprof, unifrac, unifrac_g, unifrac_g_full_tree, unweighted_unifrac, unweighted_unifrac_full_tree, weighted_normalized_unifrac, weighted_unifrac
	# 默认使用4种
	dis_method=bray_curtis,weighted_unifrac,unweighted_unifrac
	tree_method=qiime

## 1.17. otutab_ref 有参比对生成OTU表

	# 如Greengenes，可用于picurst, bugbase分析
	# 比对方法和相似度同1.10 mapping
	otutab_gg=/mnt/bai/public/ref/gg_13_5_otus/rep_set/97_otus.fasta



# 2. 统计绘图

	# 绘图通用参数
	# 实验设计文件位置，全局，其它图默认调此变量，也可单独修改；并选择表中的组列和具体分组 head -n2 doc/design.txt
	design=${wd}/doc/design.txt
	g1=groupID
	# tail -n+2 doc/design.txt|cut -f 9|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
	# 绘图使用的实验组，顺序即图中显示顺序；为空时使用所有组和默认顺序"HIND","HTEJ","LIND","LTEJ" 
	# 根据实验比对较获得实验组列表: cat doc/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" "," # ,"A50LnCp6","A56LnCp6","A50LnCp7","A56LnCp7","A50LnSz7","A56LnSz7","A50HnCp6","A56HnCp6","A50HnCp7","A56HnCp7","A50HnSz7","A56HnSz7","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","nrtHnCp7","ZH11HnCp7","ZH11LnCp7","nrtLnCp7","nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7"
	# 本研究中使用的组名
	# 籼粳稻："HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"
	# RF验证："IR24HnCp7", "IR24LnCp7", "IR24HnSz7", "IR24LnSz7","ZH11HnCp6","ZH11LnCp6", "ZH11HnSz7", "ZH11LnSz7"
	# NRT1.1b："V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50HnCp7","A56HnCp7"
	## 主图："V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7"
	## 附图："V3703LnCp6","ZH11LnCp6","A50HnCp7","A56HnCp7"
	# 时间序列："A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119"
	g1_list='"HTEJ","HIND","LTEJ","LIND","V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50HnCp7","A56HnCp7"'
	# 去除土壤作为筛选条件，只关注根系高丰度,"HSoil1","LSoil1"

	# 组间比较列表
	compare=${wd}/doc/compare.txt
	# 组间共有、特有比较列表
	venn=${wd}/doc/venn.txt
	# 图片长宽，按nature全面版、半版页面设置
	# 图片长宽和字体大小，7组以下用默认，7组以上改为8x5或更大； figure size, recommend 4x2.5, 5x3(default), 8x5, 16x10, text_size 6, 7(default), 8
	width=5
	height=3
	text_size=7

	# 图中显示legend, 如taxonomy的数量，5，8(default)，10
	legend_number=10
	# 差异统计按丰度过滤 abundance filter，单位为百分比，如丰度按每组最位数大于万分之5过滤，即0.05，减少计算量，降低OTU的FDR值，可选千一，或万一
	# 0.01 - 1711, 0.02 - 1017, 0.03 - 757, 0.04 - 609, 0.05 - 561, 0.1 - 296
	abundance_thre=0.01
	# 差异比较方法，默认是 edgeR ，可选 wilcox 秩和检验
	compare_method="wilcox"
	# 显著性P值过滤 threshold of P-value，可选0.05, 0.01, 0.001。采用FDR校正，此参数意义不大，即使0.001也没有FDR < 0.2过于严格
	pvalue=0.05
	# 统计检验方式fdr
	FDR=0.05
	# 差异变化倍数常用1.5, 2, 4倍，对应logFC为0.585, 1, 2；菌丰度变化倍数不明显，还可用1.3和1.7倍对应0.379和0.766
	FC=1.2

	# 统计绘图和网页报告版本控制
	sub=""
	doc=doc/${sub}
	# 报告输出目录
	version=xiangengASV_wilcox_${abundance_thre}


## 2.1 alpha_boxplot Alpha多样性指数箱线图 Alpha index in boxplot

	# alpha箱线图绘制参数
	ab_input=${wd}/result/alpha/index.txt
	# alpha指数种类14种 head -n1 result/alpha/index.txt | tr '\t' '\n'|tail -n+2|awk '{print "\""$1"\""}'|tr "\n" ","
	# "berger_parker","buzas_gibson","chao1","dominance","equitability","jost","jost1","reads","richness","robbins","simpson","shannon_e","shannon_2","shannon_10"
	# 默认使用chao1, richness和shannon_e三种，可修改ab_type添加或减少，shannon2/10与e相似，仅y值略有差异
	ab_method='"chao1","richness","shannon_e"'
	ab_design=${design}
	ab_group_name=${g1}
	# 默认 ab_group_list=${g1_list}
	# 主图 ab_group_list='"V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7"'
	ab_group_list=${g1_list}
	ab_output=${wd}/result/alpha/
	ab_width=${width}
	ab_height=${height}


## 2.2 alpha_rare Alpha丰富度稀释曲线 Alpha rarefracation curve

	ar_input=${wd}/result/alpha/rare.txt
	ar_design=${design}
	ar_group_name=${g1}
	ar_group_list=${ab_group_list}
	ar_output=${wd}/result/alpha/
	ar_width=${width}
	ar_height=${height}


## 2.3 beta_pcoa 主坐标轴分析距离矩阵 PCoA of distance matrix "bray_curtis","binary_jaccard","weighted_unifrac","unweighted_unifrac"

	bp_input=${wd}/result/beta/
	bp_method='"bray_curtis","weighted_unifrac","unweighted_unifrac"'
	bp_design=${design}
	bp_group_name=${g1}
	bp_group_list=${ab_group_list}
	bp_output=${wd}/result/beta/
	bp_width=${width}
	bp_height=${height}
	# 散点图是否按组添加置信椭圆，TRUE添加，FALSE不添加，默认T
	bp_ellipse=TRUE
	# 实验比较组，可用默认，也可设置单独文件，没有则不计算
	bp_compare=${compare}


## 2.4 beta_cpcoa 限制性主坐标轴分析: OTU表基于bray距离和CCA  CCA of bray distance matrix
	# 输入的OTU表，可原始count，也可以标准化的结果，无差异
	bc_input=${wd}/result/otutab.txt
	# Method from vegdist() of vegan: "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis"
	bc_method='"bray","jaccard"'
	bc_design=${design}
	bc_group_name=${g1}
	bc_group_list=${ab_group_list}
	bc_output=${wd}/result/beta/
	bc_width=${width}
	bc_height=${height}
	# 散点图是否按组添加置信椭圆，TRUE添加，FALSE不添加，默认T
	bc_ellipse=TRUE


## 2.5 tax_stackplot 样品和组分类学各级别的堆叠柱状图 Stackplot showing taxonomy in each level
	ts_input=${wd}/result/tax/sum_
	ts_level='"p","c","o","f","g"'
	ts_design=${design}
	ts_group_name=${g1}
	ts_group_list=${ab_group_list}
	ts_output=${ts_input}
	ts_width=${width}
	ts_height=${height}
	# 显列图例的数量，推荐6，8，10，默认10
	ts_number=${legend_number}
	# 设置图例的顺序，默认FALSE按分类单元字母顺序排列，TRUE则按丰度由到大小排列
	ts_order=FALSE


## 2.6 DA_compare 组间差异比较
	Dc_input=${wd}/result/otutab.txt
	# 差异比较方法edgeR or wilcox，默认edgeR
	Dc_compare=${compare}
	Dc_method=${compare_method}
	Dc_pvalue=${pvalue}
	Dc_FDR=${FDR}
	Dc_FC=${FC}
	Dc_thre=${abundance_thre}
	Dc_design=${design}
	Dc_group_name=${g1}
	Dc_group_list=${g1_list}
	# 按丰度筛选列，默认为${g1}
	Dc_group_name2="groupID2"
	Dc_output=${wd}/result/compare/


## 2.6.1 DA_compare_tax 组间差异比较,phylum/order/genus
	# 物种分类级，p/c/o/f/g五级可选
	#Dct_tax=g
	#Dct_input=${wd}/result/tax/sum_${Dct_tax}.txt
	# 差异比较方法edgeR or wilcox，默认edgeR
	Dct_compare=${compare}
	Dct_method=${compare_method}
	Dct_pvalue=${pvalue}
	Dct_FDR=${FDR}
	Dct_FC=${FC}
	Dct_thre=0.001
	Dct_design=${design}
	Dct_group_name=${g1}
	Dct_group_list=${g1_list}
	#Dct_output=${wd}/result/compare_${Dct_tax}/


## 2.7 plot_volcano 基于差异OTU表绘制火山图
	pv_input=${wd}/result/tax/sum_
	pv_design=${design}
	pv_output=${pv_input}
	pv_width=5
	pv_height=7
	# 显列图例的数量，推荐6，8，10，默认10
	pv_number=${legend_number}
	# 设置图例的顺序，默认FALSE按分类单元字母顺序排列，TRUE则按丰度由到大小排列
	pv_order=FALSE


## 2.8 plot_heatmap 基于差异OTU表绘制火山图
	ph_design=${design}
	ph_output=${ph_input}
	ph_width=7
	ph_height=5
	# 显列图例的数量，推荐6，8，10，默认10
	ph_number=${legend_number}
	# 设置图例的顺序，默认FALSE按分类单元字母顺序排列，TRUE则按丰度由到大小排列
	ph_order=FALSE
	# 绘制不同分类级热图，p,c,o,f,g
	ph_tax=g
	# 列聚类，默认TRUE
	cluster_cols=TRUE

## 2.9 plot_manhattan 绘制OTU按门着色曼哈顿图
	pm_yax=20

## 2.10 plot_boxplot 基于差异OTU表绘制火山图
	pb_input=result/otutab.txt
	pb_list='"OTU_4","OTU_3","OTU_6","OTU_11","OTU_7","OTU_9"'
	pb_design=${design}
	pb_group_name=$g1
	pb_group_list=$g1_list
	pb_output=result/otu_boxplot/
	pb_width=7
	pb_height=5
	# 显列图例的数量，推荐6，8，10，默认10
	pb_number=${legend_number}
	# 设置图例的顺序，默认FALSE按分类单元字母顺序排列，TRUE则按丰度由到大小排列
	pb_order=FALSE
	# 设置标准化
	pb_norm=TRUE
	# 转置，OTU表默认需要转置
	pb_trans=TRUE


# 2.10 plot_venn 维恩图

	# venn OTU注释数据库，如差异比较result/compare/database.txt、菌库result/39culture/otu.txt
	venn_anno=result/39culture/otu.txt
# 3 高级分析

## 3.3 faprotax 元素循环预测

	fapro_list='"nitrate_ammonification","nitrogen_fixation"'



## 3.9 culture_graphlan 可培养菌
	 
	# 可培养菌库类型，如组织root / rhizosphere / leaf, 品种A50 / IR24
	type=""
	# 指定可培养菌库位置，fa为OTU，fasta为物种如rice
	culture_db=/mnt/bai/yongxin/culture/rice/result/${type}culture_select.fasta
	# 可培养菌结果输出文件
	# 绘制Graphlan图的筛选阈值
	graph_thre=0.001
	filter=culture_${type}
	# 过滤方法，默认median，可选max, mean, median, min，数据依次减少
	filter_method=max
	otu_table=`pwd`/result/otutab.txt
	# 指定具体的实验设计、列、组筛选
	cg_design=${design}
	cg_group_name=${g1}
	cg_group_list=${ab_group_list}



# 9. 其它

	# 整合主流程
	include /mnt/bai/yongxin/github/Amplicon/16Sv2/pipeline.md
