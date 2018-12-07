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

    # 培养菌注释，采用ath root的菌库，COTU，目前只注释plot_venn的结果
    make culture # 总体均值为 74.092, 修改为R108根(b23r)为97.69，A17根为94.598

    # 默认计算所有OTU的平均丰度，这里只想用WT的平均丰度
    make -n -B culture

    # 绘图
    make -n -B culture_graphlan 

    # 筛选根际土、根的k1 OTU,并在相应库中匹配培养比例；
mkdir -p culture_Root
filter_otus_from_otu_table.sh -t 0.001 -o culture_Root -f `pwd`/result/otutab.txt -d `pwd`/doc/b23r/design.txt -F 'TRUE' -A genotype -B '"R108"' -F mean
# 筛选WT组中0.1%丰度OTU 153条
#awk '$2>0.1' temp/otutab.mean | awk '{print $1"\t"$2/100}' |tail -n+2 > culture_"Root"/otu_table_ha.mean
#awk '$2>0.1' temp/otutab.mean | cut -f 1 |tail -n+2 > culture_"Root"/otu_table_ha.id
filter_fasta.py -f result/otu.fa -o culture_"Root"/rep_seqs.fa.top -s culture_"Root"/otu_table_ha.id
echo -ne "Nature_HA_OTUs:\t" > culture_"Root"/culture.sum
grep -c '>' culture_"Root"/rep_seqs.fa.top >> culture_"Root"/culture.sum
# 分析这些OTU中可培养的比例
blastn -query culture_"Root"/rep_seqs.fa.top -db /mnt/bai/yongxin/culture/ath/result/"Root"culture_select.fa -out culture_"Root"/rep_seqs.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 # 输出13列为coverage
awk '$3*$13>=9700' culture_"Root"/rep_seqs.blastn|cut -f 1 > culture_"Root"/otu_cultured.txt
echo -ne "Stocked_OTUs:\t" >> culture_"Root"/culture.sum
grep -c 'OTU' culture_"Root"/otu_cultured.txt >> culture_"Root"/culture.sum
echo -ne "Nature_HA_abundance:\t" >> culture_"Root"/culture.sum
awk '{a=a+$2} END {print a}' culture_"Root"/otu_table_ha.mean >> culture_"Root"/culture.sum # total is 0.835
echo -ne "Stocked_abundance:\t" >> culture_"Root"/culture.sum
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]="culture"} NR>FNR {print $0,a[$1]}' culture_"Root"/otu_cultured.txt culture_"Root"/otu_table_ha.mean |grep 'culture'|awk '{a=a+$2} END {print a}' >> culture_"Root"/culture.sum 
# 绘制graphlan
sed 's/\t/\;/g' result/taxonomy_8.txt|sed 's/\;/\t/' > temp/taxonomy_2.txt
graphlan_culture.pl -i culture_"Root"/otu_table_ha.id -d culture_"Root"/otu_cultured.txt -t temp/taxonomy_2.txt -o 0_ha_otu_culture.txt
Rscript /mnt/bai/yongxin/bin/graphlan_culture.R # 生成1树, 2科注释, 3培养注释文件
sed 's/\t/\tring_alpha\t3\t/g' culture_"Root"/otu_table_ha.zscore > culture_"Root"/abundance_heat.txt # 柱状用log2，热图用zscore
cat /mnt/bai/yongxin/culture/rice/graphlan/global.cfg 2_annotation_family.txt /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg 3_annotation_match.txt /mnt/bai/yongxin/culture/rice/graphlan/abundance_heat.cfg culture_"Root"/abundance_heat.txt > culture_"Root"/5_annotation.txt
graphlan_annotate.py --annot culture_"Root"/5_annotation.txt 1_tree_plain.txt culture_"Root"/graphlan.xml
graphlan.py culture_"Root"/graphlan.xml culture_"Root"/graphlan.pdf --size 5
cat culture_"Root"/culture.sum



## 实验数据绘图Wet
    
    # 自然土
    cut -f 6 wet/1NatureSoil.txt|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'


## 与PNAS图5比较突变体下调目
    Gp10
    Acidimicrobidae
    Flavobacteriales
    Elusimicrobiales
    Caulobacterales
    Rhizobiales
    Sphingomonadales
    Burkholderiales
    Neisseriales
    Rhodocyclales
    Myxococcales
    Syntrophobacterales
    Chromatiales
    Pseudomonadales
    Spirochaetales

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


# 5. 发版图表

## 样本描述

    # 使用7个基因型2，3批进行分析 "R108", "Rnfp","lyk9","lyr4","lyk9nfp","A17","dmi3"

### PCoA compartment shape, genotype color

    # beta_pcoa_all.R 绘制根、根际、土间差异，表现基因型可变


### 绘制单个菌在不同条件下的箱线图，景美发来wet/序列与OTU比较，100%一致
    alpha_boxplot.sh -i result/otutab.txt -d doc/design.txt -A groupID -B '"soilB3S","R108b3rs","R108b3r","lyk9b3rs","lyk9b3r","soilB2S","R108b2rs","R108b2r","lyk9b2rs","lyk9b2r"' -m '"OTU_2","OTU_3"' -t TRUE -o fig/OTU_box -n TRUE -U 100

## 差异比较

## 分菌
    
    # 工作量和稀释曲线
    详细：/mnt/bai/yongxin/culture/medicago/makefile.man

    # 与原生样本比较

    # 与自然样本比较