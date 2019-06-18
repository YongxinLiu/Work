
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
	wd=medicago/AMF2
	# 创建环境代码见~/github/Work/initial_project.sh

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz

	# 从其它处复制实验设计
	cp ~/medicago/AMF/doc/L* doc/
	#cp ~/medicago/culture_start/doc/L1.txt doc/L12.txt
    # L12与之前列数不一致，调为一致
	# 删除多余空格，windows换行符等
	sed -i 's/ //g;s/\r//' doc/*.txt 
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt | sed 's/#//g') <(cat doc/L* |grep -v '#'|grep -v -P '^SampleID\t') > doc/design.txt
	# 检查是否相等
	wc -l doc/design.txt
	cut -f 1 doc/design.txt|sort|uniq|wc -l


## 1.2. 按实验设计拆分文库为样品

	# 按L1/2/3...txt拆分library为samples
	# 输入为seq/L*.fq，输出为seq/sample/*.fq
	make library_split
	make library_split_stat
	# 统计结果见result/split有txt/pdf/png，推荐看png方便快速查看每张位图
	# 查看样本量排序
	sort -k2,2n result/sample_split.log|less
    # 只有2个低于5k，4个低于5万

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
    # 100阈值去冗余 1/1M，4万多，但unoise读只有1.9d万，去噪后为1.1万

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
	# 默认过滤低测序量<5000的样本，可，筛选大于1RPM的OTU
	make otutab_filter 
    # 抽平sample_size=50000，根据 less result/otutab.biom.sum 591个样品只有12个小于5万
	make otutab_norm


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


# 复制之前筛选过样本的、二三批合并的实验设计和比较组
cp -r ~/medicago/AMF/doc/b23r10/ doc/b23r10
# 选择的8个基因型，去除dmi2和lyk9nfp
cp -r doc/b23r10 doc/b23r8
# 最后选择的6个基因型，去除lyk3, dmi2, Rnfp和lyk9nfp
cp -r doc/b23r8 doc/b23r6



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
	
	# Group compareing by edgeR, wilcox or ttest
	# 可选负二项分布GLM、Wilcoxon秩和检验或T检验
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



# 5. 发表图版

    cd ~/medicago/AMF2/fig/
    mkdir -p fig && cd fig
    mkdir -p fig1 fig2 fig3 fig4 data script
    cp /mnt/bai/yongxin/medicago/AMF/doc/design.txt data/
    cp /mnt/bai/yongxin/medicago/AMF2/result/alpha/index.txt data/
    cp /mnt/bai/yongxin/medicago/AMF2/result/beta/bray_curtis.txt data/
    cp /mnt/bai/yongxin/medicago/AMF2/result/beta/unweighted_unifrac.txt data/
    cp /mnt/bai/yongxin/medicago/AMF2/result/beta/weighted_unifrac.txt data/
    cp /mnt/bai/yongxin/medicago/AMF2/result/otutab.txt data/
    cp /mnt/bai/yongxin/medicago/AMF2/result/otutab_norm.txt data/
    cp /mnt/bai/yongxin/medicago/AMF2/result/tax/sum_* data/

## 数据筛选 
    
    # 共586个样品，包括3批(180,189,190)10个基因型根、根际土的重复，和4批重测(27)个验证Bacillus的真实性
    tail -n+2 data/design.txt|wc -l
    tail -n+2 data/design.txt|cut -f 5|uniq -c

    mkdir -p table && cd table
    cp /mnt/bai/yongxin/medicago/AMF/fig/table/table.Rmd ./
    # 详细见table/table.Rmd
    # 228个样品，9023个ASV，筛选对应的序列和物种注释
    cut -f 1 table.txt | tail -n+2 > ASV.id
    usearch10 -fastx_getseqs ~/medicago/AMF2/result/otu.fa -labels ASV.id -fastaout rep_seqs.fa
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ~/medicago/AMF2/result/taxonomy_8.txt table.txt > taxonomy.txt 
    tail -n+2 taxonomy.txt | wc -l 
    # 统计筛选后的metadata，各种样本量
    tail -n+2 metadata.txt|wc -l # 样本数量
    tail -n+2 metadata.txt|cut -f 5|uniq -c # 按批量统计
    tail -n+2 metadata.txt|cut -f 4|sort|uniq -c|awk '{print $2,$1}' # 按基因型统计
    # 统计筛选后OTU表
    usearch10 -otutab_stats table.txt -output table.stat
    cat table.stat


## 1. 样本描述description fig/fig1/fig.Rmd

    cd ~/medicago/AMF2/fig/
    cp /mnt/bai/yongxin/medicago/AMF/fig/fig1/fig1.Rmd ./fig1/
    # 2018/12/11 使用7个基因型2，3批进行分析 "A17","Anfp","lyk3","dmi3","R108","lyk9","lyr4" 去掉"dmi2","Rnfp","lyk9nfp"，备份为fig-7
    # 2019/1/3 删除lyk3共6个基因型分析 "A17","Anfp","dmi3","R108","lyk9","lyr4" 
    # 2019/5/27 删除lyk3共6个基因型分析 "A17","Anfp","dmi3","R108","lyk9","lyr4" 

### PCoA compartment shape, genotype color
    # 代码 fig/beta_pcoa_all.R 绘制根、根际、土间差异
    # fig/beta_pcoa_root.R 表现基因型可变
    # 组间距离箱线图的 beta_boxplot.R, 结果 箱线图太单一，用echart绘制beta_boxplot.txt
    

### 物种组成
    # tax_stackplot_all.R

### 附图
    # Constrained PCoA: fig/beta_cpcoa_all.R 绘制根、根际、土间差异
    # Constrained PCoA: fig/beta_cpcoa_root.R 表现基因型可变
    # Alpha diversity: fig/alpha_boxplot.R 7个基因型两批次


## 2. 差异比较compare fig2.Rmd

    cd ~/medicago/AMF2/fig/
    cp /mnt/bai/yongxin/medicago/AMF/fig/fig2/fig2.Rmd ./fig2/
    # 绘制饼形图 radius 设置中空无效，可能是版本升级原因吗？
### 差异比较曼哈顿图(2/3批混合筛选点后6个基因型4组edgeR比较) http://210.75.224.110/report/16Sv2/med_AMF_b23r6_edgeR_v1/
    cp ../script/plot_manhattan.r script/plot_manhattan.R

### 饼形图，见 fig/plot_bar_pie.R
    # 饼形图数据，edgeR明显有大量错误，改用wilcox
    cp ../med_AMF_b23r6_wilcox_v1/result/compare/*_all.txt data/

### 绘制单菌的丰度图 Pseudomonadaceae 和 Bacillus
    # 前10个菌有5个重点关注，Bacillus 2,5 Pseudomonadaceae 3
    # 绘制菌在6个基因型中变化
    cut -f 4 table/metadata.txt|tr '\t' '\n'|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'
    cd ~/medicago/AMF2/fig
    # 绘制OTU2,5 Bacillus
    alpha_boxplot.sh -i table/table.txt -d table/metadata.txt -A genocomp -B '"A17r","nfpr","dmi3r","R108r","lyk9r","lyr4r","A17rs","nfprs","dmi3rs","R108rs","lyk9rs","lyr4rs","soils"' -m '"OTU_2"' -t TRUE -o fig2/boxplot_ -n TRUE -U 100
    # 绘制Bacillus属
    alpha_boxplot.sh -i data/sum_g.txt -d table/metadata.txt -A genocomp -B '"A17r","nfpr","dmi3r","R108r","lyk9r","lyr4r","A17rs","nfprs","dmi3rs","R108rs","lyk9rs","lyr4rs","soils"' -m '"Bacillus"' -t TRUE -o fig2/boxplot_ -n TRUE -U 100
    # 绘制OTU2,5 Bacillus
    alpha_boxplot.sh -i table/table.txt -d table/metadata.txt -A genocomp -B '"A17r","nfpr","dmi3r","R108r","lyk9r","lyr4r","A17rs","nfprs","dmi3rs","R108rs","lyk9rs","lyr4rs","soils"' -m '"OTU_3"' -t TRUE -o fig2/boxplot_ -n TRUE -U 100
    # 绘制Bacillus属
    alpha_boxplot.sh -i data/sum_g.txt -d table/metadata.txt -A genocomp -B '"A17r","nfpr","dmi3r","R108r","lyk9r","lyr4r","A17rs","nfprs","dmi3rs","R108rs","lyk9rs","lyr4rs","soils"' -m '"Pseudomonas"' -t TRUE -o fig2/boxplot_ -n TRUE -U 100


## 3. 分菌
    

### 3.1 自然样品与分菌结果比较

    # 与项目中A17/R180野生型根样本与总体菌库比较
    # 修改makefile中3.9部分，参考~/medicago/culture_start/makefile
    # 注意与doc/design.txt中列和组对应，注意与分菌中文件名和品种名对应
    make culture
    make culture_graphlan
    # A17 0.001只有58个ASV太少，改为0.0005有88个


    # 3.9 culture和culture_graphlan，只需修改A17r或R108r，基于实验中大量样本的graphlan
    # 2019/3/25 与使用culture_start自然样品 ~/medicago/culture_start # culture_graphlan
    # 2019/4/9 更新 pipeline.md 中分菌部分详细注释和代码优化
    # 统计分菌代码见文末“## Sanger测序验证菌”段落，sanger测序绘制代码参考~/culture/rice/makefile.man "# 总结：1098个单菌" 段落
    cd ~/medicago/AMF2/fig/fig3
    cp ~/culture/rice/verify/graphlan_culture.R ./
    cp ~/medicago/AMF2/wet/stock_sanger/taxonomy_9.txt ./
 
    Rscript graphlan_culture.R # 生成1树, 2科注释，和来原环
    cat /mnt/bai/yongxin/culture/rice/graphlan/global.cfg 2_annotation_family.txt /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg 3_annotation_match.txt > 5_annotation.txt
    graphlan_annotate.py --annot 5_annotation.txt 1_tree_plain.txt graphlan.xml
    graphlan.py graphlan.xml culture_graphlan.pdf --size 5
    # 提取相应序列
    tail -n+2 taxonomy_9.txt|wc -l # 325 seqs
    cut -f 1 taxonomy_9.txt|tail -n+2 > seq325.id
    usearch10 -fastx_getseqs ~/medicago/AMF2/wet/stock_sanger/16s_full_length_list1.fa -labels seq325.id -fastaout seq325.fa



# 参数讨论

## 1. minisize讨论
    # http://www.drive5.com/usearch/manual/cmd_fastx_uniques.html
    # result/otu.log， 我设置 1/M 为 100 时，40991个Unique，11039个ASV
    # 1, 默认1为不过滤，6275201 seqs, 21500654 uniques, 18638709 singletons (86.7%)
    # 2， 2861945 uniques written, 289290 clusters size < 2 discarded (1.3%)
    # 绘制1-100的Unique序列曲线，观察规律
    for i in `seq 1 10`; do
    #i=1
    usearch11 -fastx_uniques temp/filtered.fa \
        -minuniquesize ${i} -sizeout \
        -fastaout temp/uniques.fa -threads 64
    echo -ne "Unique reads > ${i}\t" >> result/otu.log
    grep -c '>' temp/uniques.fa >> result/otu.log
    done
    cat result/otu.log
    # 统计10次以上各梯度的序列数量
    grep '>' temp/uniques.fa | cut -f 2 -d '=' | sed 's/;/\t/' > temp/temp.txt
    # awk筛选时，必须有\t分隔才有列
    for i in `seq 10 10 300`; do
    awk '$1>='$i temp/temp.txt | wc -l 
    done
    # 统计图表见 AMF2/doc/discussion.xlsx
    # 对于100M的reads有unique reads，21M，非Singlton 286K，> 8 518K，>10 410K，>100 40K

    # unoise聚类时间：基于>=10的41万unique Reads
    time usearch10 -unoise3 temp/uniques.fa -zotus temp/Zotus.fa
    # 139334 amplicons, 17874350 bad (size >= 10), 35393 good, 103942 chimeras 非冗余序列中有1/3为扩增子，挑选出1/4非嵌合；11小时
    # 新版本在unoise3在结果、计算时间上没有区别
    time usearch11 -unoise3 temp/uniques.fa -zotus temp/Zotus.fa -minsize 100
    # 40991 挑选出 100.0% 18986 amplicons, 11714428 bad (size >= 100), 11039 good, 7947 chimeras 非冗余序列中有1/2为扩增子，挑选出3/5非嵌合；14分钟
    # 数据量下降7倍，速度提高了47倍，大约为数据量平方倍


## 2. rhizosphere根际比较
    cp -r doc/b23r6 doc/b23rs6
    sed -i 's/b2r/b2rs/;s/b3r/b3rs/' doc/b23rs6/design.txt