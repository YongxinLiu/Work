	# 快速分析 Quick Start(所需文件准备好)
	make library_split # 样本拆分、
	make fq_qc # 合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
	make host_rm # 序列去冗余、去噪、挑选代表序列、去嵌合、去宿主，获得OTU代表序列result/otu.fa
    make beta_calc # 生成OTU表、过滤、物种注释、建库和多样性统计
	# 清除统计绘图标记(重分析时使用)
	rm -rf alpha_boxplot 
	make DA_compare # 绘制alpha、beta、taxonomy和差异OTU比较
	#rm -f plot_volcano # 删除OTU差异比较可化标记
	make plot_manhattan # 绘制差异比较的火山图、热图、曼哈顿图
	make DA_compare_tax # 高分类级差异比较，维恩图绘制，2为reads count负二项分布统计
	make plot_venn # 绘制OTU差异共有/特有维恩图
	make rmd # 生成网页报告，必须依赖的只有alpha, beta, taxonomy

	# 提取脚本
	submit=3T
	make -n -B fq_qc > pipeline.sh # 样本拆分、合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
	make -n -B host_rm >> pipeline.sh # 序列去冗余、去噪、挑选代表序列、去嵌合、去宿主，获得OTU代表序列result/otu.fa
	make -n -B beta_calc >> pipeline.sh # 生成OTU表、过滤、物种注释、建库和多样性统计
	grep -v '#' pipeline.sh > ${submit}/pipeline.sh

# 1. 准备工作 Preparation

## 1.1. 准备流程配置文件

    # Prepare config file of pipeline
	mkdir ~/medicago/AMF
    cd ~/medicago/AMF
    
    # 复制标准参数模板和操作指南至项目代码区：方便同步
	mkdir ~/github/Work/medicago/AMF
    cp ~/github/Amplicon/16Sv2/parameter.md ~/github/Work/medicago/AMF/
    cp ~/github/Amplicon/16Sv2/manual.md ~/github/Work/medicago/AMF/
   
    # 链接代码至工作区
    ln -sf ~/github/Work/medicago/AMF/parameter.md makefile
    ln -sf ~/github/Work/medicago/AMF/manual.md manual.sh

## 1.2. 初始化工作区

    # Initialize the working directory
    make init

	# 标准多文库实验设计拆分
	# split_design.pl -i doc/design_raw.txt
	# 从其它处复制实验设计 10个库
	cp ~/medicago/zjj170823/doc/L*.txt doc/
	# 删除多余空格，windows换行符等
	sed -i 's/ //g;s/\r/\n/' doc/*.txt 
	head -n3 doc/L1.txt



## 1.3. 准备原始数据

    # Prepare raw data
    ln ~/medicago/zjj170823/clean_data/*.gz seq/
    cp ~/medicago/zjj170823/doc/library.txt doc/
     # 补测27个样命名为L11
    ln ~/medicago/zjj170823/171225test/clean_data/L1_1.fq.gz seq/L11_1.fq.gz
    ln ~/medicago/zjj170823/171225test/clean_data/L1_2.fq.gz seq/L11_2.fq.gz 
    cp ~/medicago/zjj170823/171225test/doc/L1.txt doc/L11.txt
    # 补测的样本，本体为/mnt/bai/yongxin/rice/strigolactone.LiJY/clean_data/L30_1.fq.gz或L171220_CACCGG；是小璇在水稻中给苜蓿带的样品
    # 手动修改L11符合前面的格式，不重名，且可比较   
    
    # Merge paired reads, renames and merge all samples
    # 依据各文库L*.txt文件生成实验设计
    cat <(head -n1 doc/L1.txt | sed 's/#//g') <(cat doc/L* |grep -v '#') > doc/design.txt
    # 检查文件名是否唯一
    cut -f 1 doc/design.txt|sort|uniq|wc -l 
    wc -l doc/design.txt
    gunzip -f seq/*.gz

    # 检查数据质量，如果64则转换为33
    # determine_phred-score.pl seq/lane_1.fq.gz
    determine_phred-score.pl seq/L1_1.fq.gz
    # 如果为64，改原始数据为33
    rename 's/lane/lane_33/' seq/lane_*
    # 关闭质量控制，主要目的是格式转换64至33，不然usearch无法合并
	time fastp -i seq/lane_64_1.fq.gz -I seq/lane_64_2.fq.gz \
        -o seq/lane_1.fq.gz -O seq/lane_2.fq.gz -6 -A -G -Q -L -w 9
    # 1lane 80GB, 2 threads, 102min


# 2. 处理序列 Processing sequencing data

## 2.1. 按实验设计拆分lane为文库

    # Split lane into libraries
    # lane文件一般为seq/lane_1/2.fq.gz
    # lane文库信息doc/library.txt：至少包括编号、Index和样品数量三列和标题
    head -n3 doc/library.txt
    #LibraryID	IndexRC	Samples
    #L1	CTCAGA	60
    
    # 按library.txt拆分lane为library
    # make lane_split


## 2.2. 按实验设计拆分文库为样品

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
 

## 2.3. 样品双端合并、重命名、合并为单一文件

    # 样品双端合并、重命名、合并为单一文件, 注意fastq为33格式，64位采用fastp转换
    make sample_merge


## 2.4. 切除引物与标签

    # Cut primers and lables
    # 切除左端标签和引物，右端 引物
    # Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
    make fq_trim


## 2.5. 质量控制

    # Quality control
    # 过滤序列中预期累计错误率>1%的序列
    make fq_qc


## 2.6. 序列去冗余

    # Remove redundancy, get unique reads
    make fa_unqiue


## 2.7. 挑选OTU

    # Pick OTUs
    # unoise3速度比cluster_otus慢上百倍
    make otu_pick


## 2.8. 有参去嵌合体

    # Remove chimiras by silva database
    # 基于SILVA数据库去除
    make chimera_ref


## 2.9. 去除宿主

    # Remove host
    # 根据SILVA注释去除线粒体、叶绿体、真核生物18S和未知序列(非rRNA)
    make host_rm


## 2.10. 生成OTU表
    
    # Create OTUs table
    # 默认使用vsearch更快10倍，可选usearch10，线程不可超48
    make otutab_create


## 2.11. 过滤样本和OTUs

    # OTU table filter samples and OTU
    # 推荐过滤低测序量<5000的样本，筛选大于1RPM的OTU
    make otutab_filter 


## 2.12. 物种注释

    # Assign taxonomy
    # 默认使用RDP trainset快而准，GG太旧，Silva太慢
    # 推荐阈值为0.6保证注释更完整
    make tax_assign


## 2.13. 物种统计
    
    # Taxonomy summary
    # 必须所有物种有注释，否则有可能报错
    make tax_sum


## 2.14. 多序列比对和进化树
    
    # Multiply alignment and make_phylogeny
    # usearch10/culsterO结果不同可能影响多样性分析(usearch unifrac结果更可信)
    # 进化树，用于树图和多样性分析
    make tree_make

## 2.15. Alpha多样性指数计算
    
    # Calculate alpha diversity index
    # alpha指数计算结果为 result/alpha/index.txt
    # 稀释梯度结果位于 result/alpha/rare.txt
    make alpha_calc

## 2.16. Beta多样性距离矩阵计算
    
    # Beta diversity tree and distance matrix
    # 最好用usearch，结果unifrac分类更好；clustero+fastree结果PCoA较差
    make beta_calc
    # ---Fatal error--- ../calcdistmxu.cpp(32) assert failed: QueryUniqueWordCount > 0 致信作者; 改用qiime1

## 2.17. 有参考构建OTU表

    # Reference based OTU table
    # otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
    make otutab_gg



# 3. 统计绘图 Statistics and plot

    # 批量处理6个批次的分析 b1r b2r b3r b1rs b2rs b3rs
    for sub in ; do
    sub="b1r"
    sed -i "s/sub=.*/sub=$sub/" makefile
    doc=doc/${sub}
    mkdir -p $doc
    head -n8 doc/compare.txt | sed "s/b1r/$sub/g" > $doc/compare.txt
    sed "s/b1r/$sub/g" doc/venn.txt > $doc/venn.txt
    # modify sub, then report
    rm alpha_boxplot
    make plot_venn # DA otu
    make DA_compare_tax # DA taxonomy
    make rmd # report
    done

    # 制作A17/R108 ：Rhizosphere, Root比土的6个表根富集Enriched/depleted
    sub="compartment"
    doc=doc/${sub}
    mkdir -p $doc
    cp doc/compare.txt $doc/
    cp doc/venn.txt $doc/
    # 手动编辑
    rm alpha_boxplot
    make plot_venn # DA otu
    make DA_compare_tax # DA taxonomy
    make rmd # report

    # 恢复makefile至某版本，如退到b2r版本v3版本，更新为RDP注释为v4
    cp med_b2r_edgeR_v3/makefile /mnt/bai/yongxin/github/Work/medicago/AMF/parameter.md


## 3.1. Alpha多样性指数箱线图
    
    # Alpha index in boxplot
    make alpha_boxplot

## 3.2. Alpha丰富度稀释曲线
    
    # Alpha rarefracation curve
    make alpha_rare

## 3.3. 主坐标轴分析距离矩阵
    
    # PCoA of distance matrix
    make beta_pcoa

## 3.4. 限制性主坐标轴分析

    # Constrained PCoA / CCA of bray distance matrix
    # OTU表基于bray距离和CCA，至少3个组 
    make beta_cpcoa

## 3.5. 样品和组各级分类学堆叠柱状图

    # Stackplot showing taxonomy in each level
    make tax_stackplot

## 3.6. 组间差异比较 
    
    # Group compareing by edgeR or wilcox
    # 可选负二项分布，或wilcoxon秩和检验
    make DA_compare

## 2.7 绘制维恩图和生成报告
	make plot_venn # 绘制OTU差异共有/特有维恩图
	make rmd # 生成网页报告，必须依赖的只有alpha, beta, taxonomy

# 手工比较
wd=/mnt/bai/yongxin/medicago/AMF
/mnt/bai/yongxin/github/Amplicon/16Sv2/script/compare.sh -i $wd/result/otutab.txt -c doc/compare.txt -m "edgeR" \
	-p 0.01 -q 0.05 -F 1.3 -t 0.0001 \
	-d $wd/doc/design.txt -A groupID -B '"A17b1rs","Anfpb1rs","dmi2b1rs","dmi3b1rs","lyk3b1rs","lyk9b1rs","lyk9nfpb1rs","lyr4b1rs","R108b1rs","Rnfpb1rs","A17b1r","Anfpb1r","dmi2b1r","dmi3b1r","lyk3b1r","lyk9b1r","lyk9nfpb1r","lyr4b1r","R108b1r","Rnfpb1r","soilB1S","A17b2rs","Anfpb2rs","dmi2b2rs","dmi3b2rs","lyk3b2rs","lyk9b2rs","lyk9nfpb2rs","lyr4b2rs","R108b2rs","Rnfpb2rs","A17b2r","Anfpb2r","dmi2b2r","dmi3b2r","lyk3b2r","lyk9b2r","lyk9nfpb2r","lyr4b2r","R108b2r","Rnfpb2r","soilB2S","A17b3rs","Anfpb3rs","dmi2b3rs","dmi3b3rs","lyk3b3rs","lyk9b3rs","lyk9nfpb3rs","lyr4b3rs","R108b3rs","Rnfpb3rs","A17b3r","Anfpb3r","dmi2b3r","dmi3b3r","lyk3b3r","lyk9b3r","lyk9nfpb3r","lyr4b3r","R108b3r","Rnfpb3r","soilB3S"' \
	-o compare/

# 分三批分析，以第一批为例，添加KEGG及比较
修改 version 和 ab_group_list 处为第一组

# 样本筛选：根据PCoA和热图排除异常点，保存exclude.txt
cp doc/design.txt doc/design.txt.181127
for i in `grep -v '#' doc/exclude.txt`; do sed -i "s/^$i/#$i/" doc/design.txt;done
for i in `grep -v '#' doc/exclude.txt`; do sed -i "s/^$i/#$i/" doc/b23r/design.txt;done
grep '#' doc/design.txt

# 清除#注册样本
sed -i 's/#//' doc/design.txt

## 继续删除 A17b2r4 点，重新绘制PCoA beta_pcoa_compartment_batch.R


# 3 高级分析

## 3.2 picrust_compare KO组间比较

    make picrust_compare
    # 查询N nitrogen, P phosphorus 的KO结果，以 lyk9-R108_all.txt 为例
    sed '1 s/HnZH11_HnNrt/OTUID/' ~/rice/xianGeng/meta/nitrogen/ko.id > PICRUSt/nitrogen.id
    # sed '1 s/HnZH11-HnNrt/OTUID/' ~/rice/xianGeng/meta/phosphorus/ko.id > PICRUSt/phosphorus.id 手动注释和添加了KO ID
    # 追加Phosphate的功能注释
    # 批量提取
    t=phosphorus
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' /db/kegg/ko_description.txt PICRUSt/${t}.id | sed 's/\t\t/\t/' > PICRUSt/${t}.id.txt
    for i in `ll PICRUSt/ko/*sig.txt|cut -f 3 -d '/'|cut -f 1 -d '_'`; do
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $2,$3,a[$1]}' PICRUSt/ko/${i}_all.txt PICRUSt/${t}.id.txt|grep -P -v '$^' |cut -f 1-8,16- > PICRUSt/ko/${t}_${i}_all.txt        
    done
    # 结果没有可与实验功能相对应的结果


## 3.9 培养菌注释

    make culture 
    make culture_graphlan 

## 实验数据绘图Wet
    
    # 自然土
    cut -f 6 wet/1NatureSoil.txt|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'


## 与PNAS图5比较突变体下调目
    Gp10,Acidimicrobidae,Flavobacteriales,Elusimicrobiales,Caulobacterales,Rhizobiales,Sphingomonadales,Burkholderiales,Neisseriales,Rhodocyclales,Myxococcales,Syntrophobacterales,Chromatiales,Pseudomonadales,Spirochaetales,
    # 与差异分析的结果比较
    grep 'Depleted' med_b3r_edgeR_v4/result/compare/Anfpb3r-A17b3r_sig.txt | cut -f 10 | sort| uniq
    grep 'Depleted' med_b3r_edgeR_v4/result/compare/Rnfpb3r-R108b3r_sig.txt | cut -f 10 | sort| uniq
    grep 'Depleted' med_b2r_edgeR_v4/result/compare/Anfpb2r-A17b2r_sig.txt | cut -f 10 | sort| uniq
    grep 'Depleted' med_b2r_edgeR_v4/result/compare/Rnfpb2r-R108b2r_sig.txt | cut -f 10 | sort| uniq
    # 查看各组中在丰度
    grep -P 'OTUID|Flavobacteriales|Caulobacterales|Rhizobiales|Burkholderiales|Pseudomonadales' result/compare_o/database.txt
    # 查看来自这些目OTU的相对丰度
    head -n1 med_b3r_edgeR_v4/result/compare/Anfpb3r-A17b3r_sig.txt|tr '\t' '\n' | awk '{print NR,$0} # 查看列编号
    grep 'Depleted' med_b3r_edgeR_v4/result/compare/Anfpb3r-A17b3r_sig.txt | grep -P 'Flavobacteriales|Caulobacterales|Rhizobiales|Burkholderiales|Pseudomonadales' | cut -f 15 | awk '{a=a+$1} END{print a}' # 14为A，15为B
    grep 'Depleted' med_b3r_edgeR_v4/result/compare/Rnfpb3r-R108b3r_sig.txt |  grep -P 'Flavobacteriales|Caulobacterales|Rhizobiales|Burkholderiales|Pseudomonadales' | cut -f 15 | awk '{a=a+$1} END{print a}' 
    grep 'Depleted' med_b2r_edgeR_v4/result/compare/Anfpb2r-A17b2r_sig.txt |  grep -P 'Flavobacteriales|Caulobacterales|Rhizobiales|Burkholderiales|Pseudomonadales' | cut -f 15 | awk '{a=a+$1} END{print a}' 
    grep 'Depleted' med_b2r_edgeR_v4/result/compare/Rnfpb2r-R108b2r_sig.txt |  grep -P 'Flavobacteriales|Caulobacterales|Rhizobiales|Burkholderiales|Pseudomonadales' | cut -f 15 | awk '{a=a+$1} END{print a}' 

# 4. 高级分析


# 5. 发表图版

    mkdir -p fig && cd fig
    mkdir -p fig1 fig2 fig3 fig4 data script
    cp /mnt/bai/yongxin/medicago/AMF/doc/design.txt data/
    cp /mnt/bai/yongxin/medicago/AMF/result/alpha/index.txt data/
    cp /mnt/bai/yongxin/medicago/AMF/result/beta/bray_curtis.txt data/
    cp /mnt/bai/yongxin/medicago/AMF/result/beta/unweighted_unifrac.txt data/
    cp /mnt/bai/yongxin/medicago/AMF/result/beta/weighted_unifrac.txt data/
    cp /mnt/bai/yongxin/medicago/AMF/result/otutab.txt data/
    cp /mnt/bai/yongxin/medicago/AMF/result/otutab_norm.txt data/
    cp /mnt/bai/yongxin/medicago/AMF/result/tax/sum_* data/

## 数据筛选 
    
    # 共586个样品，包括3批(180,189,190)10个基因型根、根际土的重复，和4批重测(27)个验证Bacillus的真实性
    tail -n+2 data/design.txt|wc -l
    tail -n+2 data/design.txt|cut -f 5|uniq -c

    mkdir -p table && cd table
    # 详细见table/table.Rmd
    # 228个样品，9598个ASV，筛选对应的序列和物种注释
    cut -f 1 table.txt | tail -n+2 > ASV.id
    usearch10 -fastx_getseqs ~/medicago/AMF/result/otu.fa -labels ASV.id -fastaout rep_seqs.fa
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ~/medicago/AMF/result/taxonomy_8.txt table.txt > taxonomy.txt 
    wc -l taxonomy.txt 
    # 统计筛选后的metadata
    tail -n+2 table/metadata.txt|wc -l # 样本数量
    tail -n+2 table/metadata.txt|cut -f 5|uniq -c # 按批量统计
    tail -n+2 table/metadata.txt|cut -f 4|sort|uniq -c|awk '{print $2,$1}'|tr ' ' '('|tr '\n' '),' # 按基因型统计
    # 统计筛选后OTU表
    usearch10 -otutab_stats table/table.txt -output table/table.stat
    cat table/table.stat


## 1. 样本描述description fig/fig1/fig.Rmd

    # 2018/12/11 使用7个基因型2，3批进行分析 "A17","Anfp","lyk3","dmi3","R108","lyk9","lyr4" 去掉"dmi2","Rnfp","lyk9nfp"，备份为fig-7
    # 2019/1/3 删除lyk3共6个基因型分析 "A17","Anfp","dmi3","R108","lyk9","lyr4" 

### PCoA compartment shape, genotype color
    # 代码 fig/beta_pcoa_all.R 绘制根、根际、土间差异
    cp /mnt/bai/yongxin/medicago/AMF/doc/design.txt
    # fig/beta_pcoa_root.R 表现基因型可变
    # 组间距离箱线图的 beta_boxplot.R, 结果 箱线图太单一，用echart绘制beta_boxplot.txt
    

### 物种组成
    # tax_stackplot_all.R

### 附图
    # Constrained PCoA: fig/beta_cpcoa_all.R 绘制根、根际、土间差异
    # Constrained PCoA: fig/beta_cpcoa_root.R 表现基因型可变
    # Alpha diversity: fig/alpha_boxplot.R 7个基因型两批次


## 2. 差异比较compare fig2.Rmd

### 差异比较曼哈顿图(2/3批混合筛选点后10个基因型8组edgeR比较) http://210.75.224.110/report/16Sv2/med_b23r10_edgeR_v1/
    cp ../script/plot_manhattan.r script/plot_manhattan.R

### 饼形图，见 fig/plot_bar_pie.R
    cp ../med_b23r10_edgeR_v1/result/compare/*_all.txt data/

### 绘制单菌的丰度图
    # 前10个菌有5个重点关注，Pseudomonadaceae 1(相反但不显著), **3(目标菌)**, 7(不显著); **Bacillus 2**, 5(趋势一致)
    # 绘制菌在7个基因型中变化
    cut -f 4 table/metadata.txt|tr '\t' '\n'|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'
    cd ~/medicago/AMF/fig
    # 绘制OTU2 Bacillus
    alpha_boxplot.sh -i table/table.txt -d table/metadata.txt -A genocomp -B '"A17r","nfpr","dmi3r","R108r","lyk9r","lyr4r","A17rs","nfprs","dmi3rs","R108rs","lyk9rs","lyr4rs","soils"' -m '"OTU_2"' -t TRUE -o fig2/boxplot_ -n TRUE -U 100
    # 绘制Bacillus属
    alpha_boxplot.sh -i data/sum_g.txt -d table/metadata.txt -A genocomp -B '"A17r","nfpr","dmi3r","R108r","lyk9r","lyr4r","A17rs","nfprs","dmi3rs","R108rs","lyk9rs","lyr4rs","soils"' -m '"Bacillus"' -t TRUE -o fig2/boxplot_ -n TRUE -U 100
    # 绘制单个菌在不同条件下的箱线图(景美发来wet/序列与OTU比较，100%一致)
    alpha_boxplot.sh -i result/otutab.txt -d doc/design.txt -A groupID -B '"soilB3S","R108b3rs","R108b3r","lyk9b3rs","lyk9b3r","soilB2S","R108b2rs","R108b2r","lyk9b2rs","lyk9b2r"' -m '"OTU_2","OTU_3"' -t TRUE -o fig/OTU_box -n TRUE -U 100


## 3. 分菌
    
    # 工作量和稀释曲线 /mnt/bai/yongxin/culture/medicago/result/sample_rarefracation_boxplot.pdf
    # 代码详见：/mnt/bai/yongxin/culture/medicago/makefile.man

    # 与项目中A17/R180野生型根样本与总体菌库比较
    # 3.9 culture和culture_graphlan，只需修改A17r或R108r，基于实验中大量样本的graphlan
    # 2019/3/25 与使用culture_start自然样品 ~/medicago/culture_start # culture_graphlan
    # 2019/4/9 更新 pipeline.md 中分菌部分详细注释和代码优化
    # 统计分菌代码见文末“## Sanger测序验证菌”段落，sanger测序绘制代码参考~/culture/rice/makefile.man "# 总结：1098个单菌" 段落
    cd ~/medicago/AMF/fig/fig3
    cp ~/culture/rice/verify/graphlan_culture.R ./
    cp ~/medicago/AMF/wet/stock_sanger/taxonomy_9.txt ./
 
    Rscript graphlan_culture.R # 生成1树, 2科注释，和来原环
    cat /mnt/bai/yongxin/culture/rice/graphlan/global.cfg 2_annotation_family.txt /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg 3_annotation_match.txt > 5_annotation.txt
    graphlan_annotate.py --annot 5_annotation.txt 1_tree_plain.txt graphlan.xml
    graphlan.py graphlan.xml culture_graphlan.pdf --size 5
    # 提取相应序列
    tail -n+2 taxonomy_9.txt|wc -l # 325 seqs
    cut -f 1 taxonomy_9.txt|tail -n+2 > seq325.id
    usearch10 -fastx_getseqs ~/medicago/AMF/wet/stock_sanger/16s_full_length_list1.fa -labels seq325.id -fastaout seq325.fa




## 4. 实验数据分析

    # 2018/12/19 绘制定殖实验批次1-5的箱线+统计，文件来源云：苜宿-2.colonization_data_20181218
    # 保存 wet/colonization181218 ，代码见 phenotype.Rmd
    # 先分析第4批，保存为文本 b4.txt，

## 其它


### 文中数据统计
    grep 'b[23]' doc/design.txt|wc -l # 360个样，有一个数据量太少，359
    grep 'b[23]' result/split/L*txt|awk '{a=a+$2} END {print a}' # 75245220
    cat result/otu.log  # 查看OTU数据 9622
    grep 'b[23]' result/alpha/index.txt|cut -f 9|awk '{a=a+$1} END {print a/359}' # 9531，这些也基本都有了



### 数据上传
    # GSA注册项目
    http://bigd.big.ac.cn/gsub/submit/bioproject/new

    ## Sanger测序验证菌 2019/1/30
    cd ~/medicago/AMF/wet/stock_sanger
    ls contig/*.txt|wc # 307条拼接好的序列，查看为fasta格式
    ls raw_27F_515R_1492R/*|cut -f 2 -d '/'|cut -f 1 -d '.'|uniq > list.txt
    wc -l list.txt # 测了有376条序列
    # cap3拼接
    mkdir -p cap3
    for file in `cat list.txt`; do
        echo $file
        format_seq2fasta.pl -i "raw_27F_515R_1492R/${file}.*.seq" -o ${file}.fa
        cap3 ${file}.fa > /tmp/temp.txt
        # 改名
        sed -i "1 s/Contig1/${file}/" ${file}.fa.cap.contigs
        # 移至目录
        mv ${file}.fa.cap.contigs cap3/
        # 删除其它临时文件
        rm ${file}.*
    done
    ls cap3/*|wc # 329个文件
    # 序列合并
    cat cap3/* > 16s_full_length_list.fa
    format_fasta_1line.pl -i 16s_full_length_list.fa -o 16s_full_length_list1.fa # 还生成16s_full_length_list1.fa.tsv
    sed -i 's/>P03A10/>MP03A10/;s/>P03A1$/>MP03A1/' 16s_full_length_list1.fa 
    grep '>' -c 16s_full_length_list1.fa # 326条序列
    # 2019/4/9 序列来源信息，物种注释
    file=16s_full_length_list1
    usearch10 -sintax ${file}.fa -db /mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb -sintax_cutoff 0.6 -strand both -tabbedout ${file}.tax -threads 9
    cut -f 1,4 ${file}.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > taxonomy_2.txt
    # 生成物种表格：注意OTU中会有末知为空白，补齐分类未知新物种为Unassigned
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' taxonomy_2.txt | sed '1 i OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > taxonomy_8.txt
    # 添加来源，参考 /mnt/bai/yongxin/culture/rice/makefile.man # 整理分离来源的标签先分四类，再分两类
    # 景美提供 list_detail.txt 是对应来源信息，4，5，6列分别为基因型，compartment和培养基，主要为A17，Root和TSB
    # list中有354条注释，而taxonomy中有326条序列；
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' list_detail.txt taxonomy_8.txt|cut -f 4|sort|uniq -c # 有312个A17,11个R108和4个空；只划了A17为主
    # 筛选空行ID，有3个没注释，找景美查原因，
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $1,a[$1]}' list_detail.txt taxonomy_8.txt|cut -f 1,5|grep -P -v 'A17|R108'
    # 两个P03A1*名称不完整，修改补全；另一MP05D6缺失删除
    sed -i 's/P03A10/MP03A10/;s/P03A1\t/MP03A1\t/' taxonomy_8.txt
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4} NR>FNR{print $0,a[$1]}' list_detail.txt taxonomy_8.txt|grep -P 'A17|R108|Source' > taxonomy_9.txt
