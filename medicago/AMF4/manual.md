
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
	wd=medicago/AMF4
	# 创建环境代码见~/github/Work/initial_project.sh

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt
    cp ../AMF3/doc/library.txt doc/
    cat doc/library.txt
    # 共12个文库，来自17个8，11，12三批
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz

	# 标准多文库实验设计拆分，保存模板中design页为doc/design_raw.txt
    cp ../AMF3/doc/design_raw.txt doc/
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
	# 查看样本量排序，仅有3上个小于4.8万
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

    # 链接之前的OTU表
    ln ../AMF3/temp/otutab.txt temp/
	# OTU table filter samples and OTU
	# 推荐过滤低测序量<5000的样本，筛选大于1RPM的OTU
	make otutab_filter
    # 其余最小值为23907，仅有2个2万，3个3万，7个4万，3个5万，10个6万。
    # 筛选非第一批的结果，仅有2个2万，3个2万，1个4万，2个5万，10个6万。
    tail -n+16 result/otutab.biom.sum |  grep -v 'B0' | less
    # 采用6万抽平
    make otutab_norm


## 1.12. 物种注释

    # 链接之前的OTU序列
    ln ../AMF3/result/otu.fa result/
	# Assign taxonomy
	# 默认使用RDP trainset快而准，GG太旧，Silva太慢
	# 推荐阈值为0.6保证注释更完整
	make tax_assign


## 1.13. 物种统计
	
	# Taxonomy summary(ASV多，时间长)
	# 必须所有物种有注释，否则有可能报错
	make tax_sum


## 1.14. 多序列比对和进化树
	
    touch tree_make
    ln ../AMF3/result/otu.tree temp/
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
    
    # 筛选OTU
    make measurable_OTU
    # 差异比较
    rm -rf DA_compare2
	make DA_compare2
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

## 4.2 2019/10/14 根际土与根差异共有情况

    # http://210.75.224.110/report/16Sv2/med_AMF3_b23r6_v1 根差异；
    # http://210.75.224.110/report/16Sv2/med_AMF3_b23rs6_v1 根际差异；有趋势，但没有根强
    # 大于0.01%的OTU数量：根 545，根际1385，多了两倍多；共有452个；
    tail -n+2 med_AMF3_b23r6_v1/result/compare/Anfp-A17_all.txt|wc -l # 根 545 OTU大于万一
    tail -n+2 med_AMF3_b23rs6_v1/result/compare/Anfp-A17_all.txt|wc -l # 根际 1385 大于万一
    # 找两批共有的ASV，共452个
    mkdir -p comare_r_rs
    cat <(tail -n+2 med_AMF3_b23r6_v1/result/compare/Anfp-A17_all.txt) <(tail -n+2 med_AMF3_b23rs6_v1/result/compare/Anfp-A17_all.txt) | cut -f 1 |sort|uniq -d > comare_r_rs/common_OTU.id #|wc -l
    # v1筛选共有OTU的维恩列表，只选这452个中共有的(但这不太好和别人解析)
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]="TRUE"} NR>FNR{print $0,a[$1]}' \
        comare_r_rs/common_OTU.id med_AMF3_b23rs6_v1/result/compare/diff.list | \
        grep 'TRUE'| cut -f 1-2 | sed 's/$/_RS/' | less > comare_r_rs/diff.list.RS # | wc -l
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]="TRUE"} NR>FNR{print $0,a[$1]}' \
        comare_r_rs/common_OTU.id med_AMF3_b23r6_v1/result/compare/diff.list | \
        grep 'TRUE'| cut -f 1-2 | sed 's/$/_R/' | less > comare_r_rs/diff.list.R # | wc -l
    cat comare_r_rs/diff.list.RS comare_r_rs/diff.list.R | sed "1 i ID\ttype" > comare_r_rs/diff.list
    # Venn比较根，根际上、下变化OTU，编写venn.txt
	batch_venn.pl -i comare_r_rs/venn.txt -d comare_r_rs/diff.list
    # 结果重合非常好
    # 另一种筛选数据的方法，判断有无，则输出
    awk 'NR==FNR{a[$1]; next} ($1 in a){print}' comare_r_rs/common_OTU.id med_AMF3_b23r6_v1/result/compare/diff.list > comare_r_rs/diff.list.R2

    # v2筛选根的OTU在根际中的变化，需要依赖根的OTU来比较根际
    # 查找compare_OTU_ref.sh的用法
    for i in `find ./ -name manual.md`; do
        grep -A 3 'compare_OTU_ref.sh' $i
    done
    # 指定参考ID差异分析根际土
    rm -fr `pwd`/result/compare/
    mkdir -p `pwd`/result/compare/
    compare_OTU_ref.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/"b23rs6"/compare.txt -m "wilcox" \
        -p 0.05 -q 0.05 -F 1.1 -t 0.01 \
        -d `pwd`/doc/"b23rs6"/design.txt  -A genotype -B '"A17","Anfp","dmi3","R108","lyk9","lyr4"' \
        -o `pwd`/result/compare/ -r med_AMF3_b23r6_v1/result/compare/Anfp-A17_all.txt
    # 根和根际混合维恩
    tail -n+2 med_AMF3_b23r6_v1/result/compare/diff.list | sed 's/$/_R/' > comare_r_rs/diff.list.R # | wc -l
    tail -n+2 result/compare/diff.list | sed 's/$/_RS/' > comare_r_rs/diff.list.RS # | wc -l
    cat comare_r_rs/diff.list.RS comare_r_rs/diff.list.R | sed "1 i ID\ttype" > comare_r_rs/diff.list
    # Venn比较根，根际上、下变化OTU，编写venn.txt
	batch_venn.pl -i comare_r_rs/venn.txt -d comare_r_rs/diff.list
    # 结果重合非常好



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

    # 与项目中A17/R180野生型根样本与总体菌库比较 2019/10/15
    # 修改makefile中3.9部分，参考~/medicago/culture_start/makefile
    # 注意与doc/design.txt中列和组对应，注意与分菌中文件名和品种名对应
    conda deactivate
    rm culture
    make culture
    rm culture_graphlan
    make culture_graphlan
    # A17 0.001只有58个ASV太少，改为0.0005有88个；这是两个样R1研磨和R2根最大值；
    # 当只有A17 R2根时，改为万三 95
    # R108： 万一 213程序自动选Top100, 万三 67；感沉图形与实际样品不一致，查看OTU_2，在A17中有0.01，在R108中为0
    # 而该OTU_2中样品根中大多数非常高，高达60%



    # 3.9 culture和culture_graphlan，只需修改A17r或R108r，基于实验中大量样本的graphlan
    # 2019/3/25 与使用culture_start自然样品 ~/medicago/culture_start # culture_graphlan
    # 2019/4/9 更新 pipeline.md 中分菌部分详细注释和代码优化
    # 统计分菌代码见文末“## Sanger测序验证菌”段落，sanger测序绘制代码参考~/culture/rice/makefile.man "# 总结：1098个单菌" 段落

    # 2019/10/22 分菌两批总结，数据源~/culture/medicago/190626/：详见 ~/culture/medicago/190626/manual.sh ## Sanger鉴定
    # 绘图
    cd ~/medicago/AMF3/fig/fig3
    cp ~/culture/rice/verify/graphlan_culture.R ./
    mkdir -p result
    cp ~/culture/medicago/190626/wet/result/taxonomy_9.txt result/
    # 脚本需要根据实验的文件修改，本次生成3，4两个文件
    Rscript graphlan_culture.R 
	cat /mnt/bai/yongxin/culture/rice/graphlan/global.cfg 2_annotation_family.txt /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg 3_annotation_match.txt /mnt/bai/yongxin/culture/rice/graphlan/ring2.cfg 4_annotation_match.txt > 5_annotation.txt 
    graphlan_annotate.py --annot 5_annotation.txt 1_tree_plain.txt graphlan.xml
    graphlan.py graphlan.xml culture_graphlan.pdf --size 5 # 243
    # 提取相应序列
#    tail -n+2 taxonomy_9.txt|wc -l # 325 seqs
#    cut -f 1 taxonomy_9.txt|tail -n+2 > seq325.id
#    usearch10 -fastx_getseqs ~/medicago/AMF3/wet/stock_sanger/16s_full_length_list1.fa -labels seq325.id -fastaout seq325.fa


    ### 与分菌桑格比较 2020/4/29仅使用了191022的第二批，2020/5/14更新为三批
    # 提取序列
    cut -f 1 wet/R108_OTU.txt | tail -n +2 > wet/R108_OTU.id
    usearch10 -fastx_getseqs result/otu.fa -labels wet/R108_OTU.id -fastaout wet/R108_OTU.fa
    blastn -query wet/R108_OTU.fa -db ~/culture/medicago/190626/wet/merge.fa -out wet/R108_OTU.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 100 -evalue 1 -num_threads 9 # 输出13列为coverage
    # 筛选可培养
    awk '$3>=97' wet/R108_OTU.blastn|cut -f 1-3|uniq > wet/R108_OTU.cultured.txt
    # 追加至表格
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]="cultured"} NR>FNR {print $0,a[$1]}' wet/R108_OTU.cultured.txt wet/R108_OTU.txt > wet/R108_OTU_cultured.txt
    # 添加OTU孔信息
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]="cultured"} NR>FNR {print $0,a[$1]}' culture10/otu_cultured.id wet/R108_OTU_cultured.txt | sed '1 s/\tcultured\t$/\tCulturedStock\tCulturedWell/' > wet/R108_OTU_cultured2.txt
    # 添加菌保信息
    sed -i '1 i OTUID\tStockID\tSimilarity' wet/R108_OTU.cultured.txt
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print $0,a[$2]}' \
        ~/culture/medicago/190626/wet/IDinformation.txt \
        wet/R108_OTU.cultured.txt \
        > wet/R108_OTU.cultured_anno.txt 

    ### R108与土壤差异比 2020/5/10
    cut -f 1 med_AMF3_b23r6S_v1/result/compare/soil-R108_all.txt | tail -n +2 > wet/soil.id
    ### 土壤与R108和突变体比较 2020/5/12
    cut -f 1 med_AMF3_b23r6S_v2/result/compare/diff.list.vennsoil_lyk9_Esoil_lyr4_Esoil_R108_E.xls.xls | tail -n 144 | sed '1 i OTUID' > wet/soil.id
    usearch10 -fastx_getseqs result/otu.fa -labels wet/soil.id -fastaout wet/soil.fa
    blastn -query wet/soil.fa -db ~/culture/medicago/190626/wet/191022V5-7.fa -out wet/soil.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 # 输出13列为coverage
    # 筛选可培养
    awk '$3>=97' wet/soil.blastn|cut -f 1-3|uniq > wet/soil.cultured.txt
    # 土壤差异比较排序，
    # sed -i 's/\r//' med_AMF3_b23r6S_v2/result/compare/soil-R108_all_sort.txt 
    # 土壤中富集的ID，添加物种丰度和物种注释
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/taxonomy_8.txt wet/soil.id > wet/soil.tax
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/compare/database.txt wet/soil.id> wet/soil.mean
    paste wet/soil.tax wet/soil.mean | cut -f 1-8,10- > wet/soil.tax_mean
    # 添加是否可培养
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]="cultured"} NR>FNR {print $0,a[$1]}' wet/soil.cultured.txt wet/soil.tax_mean | grep -P 'OTUID|cultured'> wet/soil_cultured2.txt # 144个共有，只有15个有Stocks
    # 添加OTU孔信息
    # awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]="cultured"} NR>FNR {print $0,a[$1]}' culture10/otu_cultured.id wet/soil_cultured.txt > wet/soil_cultured2.txt

    # 整理结果文件
    wd=wet/200510
    mkdir -p $wd
    cp wet/R108_OTU_cultured2.txt $wd/1.R108_OTU_cultured.txt # 最后两列cultured分别代表菌保，菌高通量测序well中存在
    cp wet/R108_OTU.cultured.txt $wd/2.R108_OTU_cultured_stocks.txt # OTU对应菌保列表，1对多
    cp culture10/otu_cultured.txt $wd/3.R108_OTU_COTU.txt # OTU对应COTU，1对多
    sed 's/^/COTU_/' /mnt/bai/yongxin/culture/medicago/190626/result/culture_select2.xls > $wd/4.R108_COTU_wells.txt # COTU对应孔，1对多
    cp wet/soil_cultured2.txt $wd/5.Soil_OTU_cultured.txt # 最后两列cultured分别代表菌保，菌高通量测序well中存在，stock列表见5，well信息同3，4
    cp wet/soil.cultured.txt $wd/6.Soil_OTU_cultured_stocks.txt # 土壤OTU对应菌保列表，1对多

    ### 聚类可培养菌
    # 提取高丰度序列
    usearch10 -fastx_getseqs result/otu.fa -labels wet/R108_OTU.id -fastaout wet/R108_OTU.fa # 546个高丰度菌
    # 聚类结果，545个序列，获得213个簇，最多18条序列聚一起
    usearch10 -cluster_fast wet/R108_OTU.fa -id 0.97 -centroids wet/R108_OTU.nr.fa -uc wet/R108_OTU.clusters.uc -sizeout 
    # 提取序列，聚类和相似度
    cut -f 4,9,10 wet/R108_OTU.clusters.uc | awk 'BEGIN{FS=OFS="\t"}{if($3=="*"){print $2,$2,100} else{print $2,$3,$1}}' > wet/R108_OTU.clusters
    # 追加聚类和相似度
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3} NR>FNR {print $0,a[$1]}' wet/R108_OTU.clusters wet/R108_OTU_cultured2.txt | sed '1 s/$/Cluster\tSimilarity/' > wet/R108_OTU_cultured3.txt
 
    # 建树可视化选菌
    # 筛选高丰度，conserved建树，可培养，并树上注释门、可培养
    # 标注列的序列
    #### 方案1. 千分之3的菌
    sed -i 's/\r//' wet/R108_OTU_cultured3.txt
    head -n1 wet/R108_OTU_cultured3.txt|tr '\t' '\n'|awk '{print NR,$1}'
    # 筛选conserved，再筛选可培养，再筛选丰度，1%为21个，0.5为34个，0.4为40，0.3为56个，0.2为76个，0.1为117个；0.3且有Stocked的
    grep 'Conserved' wet/R108_OTU_cultured3.txt|grep 'cultured'|awk '$9>0.3 || $10>0.3 || $14>0.3'|awk '$19=="cultured"'|less -S| wc -l 
    grep 'Conserved' wet/R108_OTU_cultured3.txt|grep 'cultured'|awk '$9>0.3 || $10>0.3 || $14>0.3'|less -S > temp/syncom_OTU.txt # wc -l 
    # 提取序列
    cut -f 1 temp/syncom_OTU.txt > temp/syncom_OTU.id
    usearch10 -fastx_getseqs result/otu.fa -labels temp/syncom_OTU.id -fastaout temp/syncom_OTU.fa
    # 建树
    threshold=0.3
    mkdir -p wet/syncom/p${threshold}
    time muscle -in temp/syncom_OTU.fa -out temp/syncom_OTU_aligned.fas
    time iqtree -s temp/syncom_OTU_aligned.fas -st DNA -m TEST -bb 1000 -alrt 1000 -nt 20 -pre wet/syncom/p${threshold}/tree -quiet -redo
    # 制作注释文件
    cat <(head -n1 wet/R108_OTU_cultured3.txt) temp/syncom_OTU.txt | cut -f 1,3-7,9,10,14,17-19,21> wet/syncom/p${threshold}/annotation.txt #  sed '1 s/\t$/CulturedStock\tCluster/'| less 
    Rscript ~/bin/table2itol/table2itol.R -a -c double -D  wet/syncom/p${threshold} -i OTUID -l Genus -t %s -w 0.5 wet/syncom/p${threshold}/annotation.txt

    #### 方案2. 组间差异菌Top30的可培养部分，需要对wet/R108_OTU_cultured3.txt求三组均值再排序
    # R108-9, lyk9-10, lyr4-14，按Enriched、Depleted和NotSig三类筛选
    type=NotSig
    # 共140
    grep 'Conserved' wet/R108_OTU_cultured3.txt|awk '$9>0.1 || $10>0.1 || $14>0.1'|less -S| wc -l 
    # 13个Depleted, 99个Enriched, 28个NotSig；可培养分别为7，18，27
    grep 'Conserved' wet/R108_OTU_cultured3.txt|awk '$9>0.1 || $10>0.1 || $14>0.1'|less -S|cut -f 17|sort|uniq -c
    # 筛选各类Top30
    grep 'Conserved' wet/R108_OTU_cultured3.txt|awk '$9>0.1 || $10>0.1 || $14>0.1'|grep $type |head -n 30| less -S > temp/syncom_OTU_${type}.txt
    # 提取序列
    cut -f 1 temp/syncom_OTU_${type}.txt > temp/syncom_OTU_${type}.id
    usearch10 -fastx_getseqs result/otu.fa -labels temp/syncom_OTU_${type}.id -fastaout temp/syncom_OTU_${type}.fa
    # 建树
    mkdir -p wet/syncom/${type}
    time muscle -in temp/syncom_OTU_${type}.fa -out temp/syncom_OTU_${type}_aligned.fas
    rm -rf wet/syncom/${type}/${type}*
    time iqtree -s temp/syncom_OTU_${type}_aligned.fas -st DNA -m TEST -bb 1000 -alrt 1000 -nt 20 -pre wet/syncom/${type}/${type} -quiet -redo
    # 制作注释文件
    cat <(head -n1 wet/R108_OTU_cultured3.txt) temp/syncom_OTU_${type}.txt | cut -f 1,3-7,9,10,14,17-19,21| less > wet/syncom/${type}/annotation.txt # | sed '1 s/\t$/CulturedStock\tCluster/'
    Rscript ~/bin/table2itol/table2itol.R -a -c double -D  wet/syncom/${type} -i OTUID -l Genus -t %s -w 0.5 wet/syncom/${type}/annotation.txt
    grep -c 'cultured' wet/syncom/${type}/iTOL_colorstrip-CulturedStock.txt
    

    #### 人工挑选菌去冗余
    cd ~/medicago/AMF3/wet

    # Plan1
    cat 200518/plan1.id | sort | uniq -d # M260和M6重复
    cat 200518/plan1.id | sort | uniq > 200518/plan1.id.nr
    wc -l 200518/plan1.id.nr # 35个非冗余
    usearch10 -fastx_getseqs /mnt/bai/yongxin/culture/medicago/190626/wet/merge.fa -labels 200518/plan1.id.nr -fastaout 200518/plan1.fa
    # 建索引比对
    makeblastdb -in 200518/plan1.fa -dbtype nucl
    blastn -query 200518/plan1.fa -db 200518/plan1.fa -out 200518/plan1.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 100 -evalue 1 -num_threads 9 # 输出13列为coverage
    grep 100.000 200518/plan1.blastn|less|wc -l # 只有35个，没有重复的

    # Plan2
    cat 200518/plan2.id | sort | uniq -d # 无
    cat 200518/plan2.id | sort | uniq > 200518/plan2.id.nr
    wc -l 200518/plan2.id.nr # 26个非冗余
    usearch10 -fastx_getseqs /mnt/bai/yongxin/culture/medicago/190626/wet/merge.fa -labels 200518/plan2.id.nr -fastaout 200518/plan2.fa
    # 建索引比对
    makeblastdb -in 200518/plan2.fa -dbtype nucl
    blastn -query 200518/plan2.fa -db 200518/plan2.fa -out 200518/plan2.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 100 -evalue 1 -num_threads 9 # 输出13列为coverage
    grep 100.000 200518/plan2.blastn|less|wc -l # 28个百分百，有2个重复
    grep 100.000 200518/plan2.blastn | cut -f 1 | sort | uniq -d # MN39和52序列相同

    # Plan3，Plan1+8个新菌
    cat 200518/plan3.id | sort | uniq -d # M260和M6重复
    cat 200518/plan3.id | sort | uniq > 200518/plan3.id.nr
    wc -l 200518/plan3.id.nr # 43个非冗余
    usearch10 -fastx_getseqs /mnt/bai/yongxin/culture/medicago/190626/wet/merge.fa -labels 200518/plan3.id.nr -fastaout 200518/plan3.fa
    # 建索引比对
    makeblastdb -in 200518/plan3.fa -dbtype nucl
    blastn -query 200518/plan3.fa -db 200518/plan3.fa -out 200518/plan3.blastn -num_alignments 10 -evalue 1 -num_threads 9 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' # 输出13列为coverage
    grep 100.00 200518/plan3.blastn|less|wc -l # 51个，(51-43)/2=4个重复
    grep 100.00 200518/plan3.blastn | cut -f 1 | sort | uniq -d > 200518/plan3.dup # 显示序列相同的ID

    # plan4，Plan2+6个新菌
    cat 200518/plan4.id | sort | uniq -d # 无
    cat 200518/plan4.id | sort | uniq > 200518/plan4.id.nr
    wc -l 200518/plan4.id.nr # 32个非冗余
    usearch10 -fastx_getseqs /mnt/bai/yongxin/culture/medicago/190626/wet/merge.fa -labels 200518/plan4.id.nr -fastaout 200518/plan4.fa
    # 建索引比对
    makeblastdb -in 200518/plan4.fa -dbtype nucl
    blastn -query 200518/plan4.fa -db 200518/plan4.fa -out 200518/plan4.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 100 -evalue 1 -num_threads 9 # 输出13列为coverage
    grep 100.00 200518/plan4.blastn|less|wc -l # 40个，(40-32)/2=4个重复
    grep 100.00 200518/plan4.blastn | cut -f 1 | sort | uniq -d > 200518/plan4.dup # 显示序列相同的ID

    # plan5，Plan3+3个新菌
    cat 200518/plan5.id | sort | uniq -d # M260和M6重复
    cat 200518/plan5.id | sort | uniq > 200518/plan5.id.nr
    wc -l 200518/plan5.id.nr # 45个非冗余
    usearch10 -fastx_getseqs /mnt/bai/yongxin/culture/medicago/190626/wet/merge.fa -labels 200518/plan5.id.nr -fastaout 200518/plan5.fa
    # 建索引比对
    makeblastdb -in 200518/plan5.fa -dbtype nucl
    blastn -query 200518/plan5.fa -db 200518/plan5.fa -out 200518/plan5.blastn -num_alignments 10 -evalue 1 -num_threads 9 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' # 输出13列为coverage
    grep 100.00 200518/plan5.blastn|less|wc -l # 55个，(55-45)/2=5个重复
    grep 100.00 200518/plan5.blastn | cut -f 1 | sort | uniq -d > 200518/plan5.dup # 显示序列相同的ID
    # 见脚本 blast_self_heatmap.Rmd

    # 200525/plan1
    cat 200525/plan1.id | sort | uniq -d # 无
    cat 200525/plan1.id | sort | uniq > 200525/plan1.id.nr
    wc -l 200525/plan1.id.nr # 39个非冗余
    usearch10 -fastx_getseqs /mnt/bai/yongxin/culture/medicago/190626/wet/merge.fa -labels 200525/plan1.id.nr -fastaout 200525/plan1.fa
    # 建索引比对
    makeblastdb -in 200525/plan1.fa -dbtype nucl
    blastn -query 200525/plan1.fa -db 200525/plan1.fa -out 200525/plan1.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 100 -evalue 1 -num_threads 9 # 输出13列为coverage
    grep 100.00 200525/plan1.blastn|less|wc -l # 45个，(45-39)/2=3个重复
    grep 100.00 200525/plan1.blastn | cut -f 1 | sort | uniq -d > 200525/plan1.dup # 显示序列相同的ID

    # 200525/plan2
    cat 200525/plan2.id | sort | uniq -d # 无
    cat 200525/plan2.id | sort | uniq > 200525/plan2.id.nr
    wc -l 200525/plan2.id.nr # 39个非冗余
    usearch10 -fastx_getseqs /mnt/bai/yongxin/culture/medicago/190626/wet/merge.fa -labels 200525/plan2.id.nr -fastaout 200525/plan2.fa
    # 建索引比对
    makeblastdb -in 200525/plan2.fa -dbtype nucl
    blastn -query 200525/plan2.fa -db 200525/plan2.fa -out 200525/plan2.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 100 -evalue 1 -num_threads 9 # 输出13列为coverage
    grep 100.00 200525/plan2.blastn|less|wc -l # 28个，(32-28)/2=2个重复
    grep 100.00 200525/plan2.blastn | cut -f 1 | sort | uniq -d > 200525/plan2.dup # 显示序列相同的ID

    # 合并三组,48个(暂不用)
    type=root
    cat temp/syncom_OTU_Enriched.id temp/syncom_OTU_Depleted.id temp/syncom_OTU_NotSig.id > temp/syncom_OTU_${type}.id
    usearch10 -fastx_getseqs result/otu.fa -labels temp/syncom_OTU_${type}.id -fastaout temp/syncom_OTU_${type}.fa
    # 建树
    mkdir -p wet/syncom/${type}
    time muscle -in temp/syncom_OTU_${type}.fa -out temp/syncom_OTU_${type}_aligned.fas
    time iqtree -s temp/syncom_OTU_${type}_aligned.fas -st DNA -m TEST -bb 1000 -alrt 1000 -nt 20 -pre wet/syncom/${type}/${type} -quiet -redo
    # 制作注释文件
    cat temp/syncom_OTU_Enriched.txt temp/syncom_OTU_Depleted.txt temp/syncom_OTU_NotSig.txt > temp/syncom_OTU_${type}.txt
    cat <(head -n1 wet/R108_OTU_cultured3.txt) temp/syncom_OTU_${type}.txt | cut -f 1,3-7,9,10,14,17-19,21| less > wet/syncom/${type}/annotation.txt # | sed '1 s/\t$/CulturedStock\tCluster/'
    Rscript ~/bin/table2itol/table2itol.R -a -c double -D  wet/syncom/${type} -i OTUID -l Genus -t %s -w 0.5 wet/syncom/${type}/annotation.txt

    # 土壤中富集(暂不用)
    type=soil
    sed -i 's/\r//' wet/soil_cultured3.txt
    head -n1 wet/soil_cultured3.txt|tr '\t' '\n'|awk '{print NR,$1}'
    # 筛选土壤中特异的OTUID，276个，可培养149，丰度>0.1共42个，富集的38个
    grep 'cultured' wet/soil_cultured3.txt |awk '$14>0.1' | grep 'Enriched'| cut -f 1|less -S > temp/soil_cultured.id
    # 与root的比较找特有
    cat temp/soil_cultured.id temp/syncom_OTU_root.id | sort | uniq -d > temp/soil_cultured_nonUniq.id
    # 提取土壤Uniq 
    cat temp/soil_cultured.id temp/soil_cultured_nonUniq.id | sort | uniq -u > temp/soil_cultured_Uniq.id
    # 标记并筛选
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]="noUniq"} NR>FNR {print $0,a[$1]}' temp/soil_cultured_nonUniq.id wet/soil_cultured3.txt > wet/soil_cultured3_marker.txt
    grep 'cultured' wet/soil_cultured3_marker.txt |awk '$14>0.1' | grep -v 'noUniq' | head -n20 > temp/syncom_OTU_${type}.txt 
    # 提取序列
    cut -f 1 temp/syncom_OTU_${type}.txt > temp/syncom_OTU_${type}.id
    usearch10 -fastx_getseqs result/otu.fa -labels temp/syncom_OTU_${type}.id -fastaout temp/syncom_OTU_${type}.fa
    # 建树
    mkdir -p wet/syncom/${type}
    time muscle -in temp/syncom_OTU_${type}.fa -out temp/syncom_OTU_${type}_aligned.fas
    time iqtree -s temp/syncom_OTU_${type}_aligned.fas -st DNA -m TEST -bb 1000 -alrt 1000 -nt 20 -pre wet/syncom/${type}/${type} -quiet -redo
    # 制作注释文件
    cat <(head -n1 wet/soil_cultured3.txt) temp/syncom_OTU_${type}.txt | cut -f 1,6-18|less > wet/syncom/${type}/annotation.txt
    Rscript ~/bin/table2itol/table2itol.R -a -c double -D  wet/syncom/${type} -i OTUID -l Genus -t %s -w 0.5 wet/syncom/${type}/annotation.txt


    # 比较两种方法共有
    # 0.3%下的56个，分上下不变挑选的48个。共有41个
    cat temp/syncom_OTU.id temp/syncom_OTU_root.id | sort | uniq -d |wc -l




    ### 与最新10个菌库比较 2019/7/2
    cwd=culture10
    mkdir -p ${cwd}
    blastn -query ~/medicago/AMF3/result/otu.fa -db ~/culture/medicago/190626/result/culture_select.fa -out ${cwd}/otu_culture.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 # 输出13列为coverage
    awk '$3>=97' ${cwd}/otu_culture.blastn |cut -f 1-3 > ${cwd}/otu_cultured.txt # 197 OTU cultured, 197/535=36.8%
    awk '$3>=97' ${cwd}/otu_culture.blastn |cut -f 1 | sort |uniq > ${cwd}/otu_cultured.id
    # 我们关注的菌OTU_2 k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_1;g__Bacillus，在新库中${cwd}/otu_culture.blastn 对应 COTU_61
    # L8-10为新库，精选信息可查culture_select.xls中编号61，L8,L10中为最优解；更多孔详见 culture_bacteria.xls
    # 添加孔的信息
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$8"\t"$9"\t"$10} NR>FNR {print $0,a[$1]}' doc/design.txt result/culture_bacteria.xls > result/culture_bacteria_anno.xls




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

