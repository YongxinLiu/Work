	# 快速分析 Quick Start(所需文件准备好)
	make library_split # 样本拆分、
	make fq_qc # 合并、去接头和引物、质控，获得纯净扩增子序列temp/filtered.fa
	make host_rm # 序列去冗余、去噪、挑选代表序列、去嵌合、去宿主，获得OTU代表序列result/otu.fa
	make beta_calc # 生成OTU表、过滤、物种注释、建库和多样性统计
	# 清除统计绘图标记(重分析时使用)
	rm -rf alpha_boxplot 
	make DA_compare # 绘制alpha、beta、taxonomy和差异OTU比较
	rm -f plot_volcano # 删除OTU差异比较可化标记
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

# 标准流程 Standard pipeline

# 2019/1/15 添加土壤中富集菌

# 比对IND和TEJ相对于土壤中富集的情况，这296个OTU仅包括了TEJ       IND      Soil 29.14521 47.47679 58.75728，这384个菌，土壤更富集？
compare_OTU_ref.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare_soil.txt -m "wilcox" \
        -p 0.05 -q 0.05 -F 1.2 -t 0.1 \
        -d `pwd`/doc/design.txt -A subspecies -B '"TEJ","IND","Soil"' \
        -o `pwd`/result/compare_ref/ -r xiangeng_wilcoxon_main/result/compare/A50LnCp6-A56LnCp6_all.txt
cat result/compare_ref/summary.txt # 多数为植物，小部分为土，极少为不变
cp result/compare_ref/IND-Soil_all.txt fig1/data/
cp result/compare_ref/TEJ-Soil_all.txt fig1/data/


# 注释更换组，维恩图添加丰度注释，参考 http://210.75.224.110/report/16Sv2/xiangeng_wilcoxon_main

# 比对IND和TEJ相对于土壤中富集的情况，这296个OTU仅包括了TEJ       IND      Soil 29.14521 47.47679 58.75728，这384个菌，土壤更富集？
# 修正：A50LnCp6-A56LnCp6共30个样本只有3个差异，添加 A50HnCp7 - A56HnCp7 有差异
# 只有有参，p0.01, q0.05, F1.2, t0.1与网页中结果一致
compare_OTU_ref.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
        -p 0.01 -q 0.05 -F 1.2 -t 0.1 \
        -d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","A50LnCp7","A56LnCp7","A50HnCp7","A56HnCp7"' \
        -o `pwd`/result/compare/ -r xiangeng_wilcoxon_main/result/compare/A50LnCp6-A56LnCp6_all.txt
cat result/compare/summary.txt # 多数为植物，小部分为土，极少为不变

# 制作带培养菌丰度和土壤富集信息的注释
cp result/compare_ref/database.txt temp/soil_mean_abundance.txt
paste result/compare/database.txt <(cut -f 6 result/compare_ref/IND-Soil_all.txt) <(cut -f 6 result/compare_ref/TEJ-Soil_all.txt) | sed '1 s/level\tlevel/IND-soil\tTEJ-soil/' | less > temp/mean_type.txt
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' temp/mean_type.txt result/39culture/otu.txt | less -S > result/39culture/otu.txt.bak

# 统计IND中对应数量，在soil中汇总对应的丰度
tail -n 141 xiangeng_wilcoxon_main1/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_D.xls.xls|cut -f 22|sort|uniq -c
tail -n 141 xiangeng_wilcoxon_main1/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_D.xls.xls|cut -f 23|sort|uniq -c
tail -n 16 xiangeng_wilcoxon_main1/result/compare/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls|cut -f 22|sort|uniq -c
tail -n 16 xiangeng_wilcoxon_main1/result/compare/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls|cut -f 23|sort|uniq -c


# 2019/1/15 测试nrt overlap IND enriched与选菌序列相似性

# 所有IND富集
grep 'HTEJ_HIND_D_V3703HnCp6_ZH11HnCp6_D-specific_to_others' -A 1000 xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls.xls |cut -f 1|grep 'OTU' > wet/otuid296.id
wc -l wet/otuid296.id # 统计筛选部分数量119，重复运行上部确定一致性

# 所有TEJ富集
grep 'HTEJ_HIND_D_V3703HnCp6_ZH11HnCp6_D-specific_to_others' -A 1000 xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls.xls |cut -f 1|grep 'OTU' > wet/otuid296.id
wc -l wet/otuid296.id # 统计筛选部分数量109，重复运行上部确定一致性


# 确定选菌是否与OTU一致
usearch10 -fastx_getseqs result/otu.fa -labels wet/otuid296.id -fastaout wet/otuid296.fa
# cp result/otu.fa wet/otuid296.fa
makeblastdb -in wet/otuid296.fa -dbtype nucl
blastn -query ASV/wet/Syncom25.fa -db wet/otuid296.fa -out temp/Syncom25.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 9 
# less temp/Syncom25.blastn # 18/24，6个相似并对，Co7无法比对；改为万5为3/24，万4为2/24，万1和3只差土全包括1/24
awk '$3>97' temp/Syncom25.blastn | wc -l

# TS6. 添加土壤对应的菌
# 整理好的两个品种比较的见 result/39culture/otu.txt.bak
head -n 1 result/39culture/otu.txt.bak|sed 's/\t/\n/g'|awk '{print NR,$0}' # 22/23列
cp xiangeng_wilcoxon_main/result/compare/LIND-LTEJ_all.txt fig1/ST/05.LIND-LTEJ_all.txt
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$22"\t"$23} NR>FNR{print $1,a[$1]}' result/39culture/otu.txt.bak fig1/ST/05.LIND-LTEJ_all.txt > fig1/ST/05.LIND-LTEJ_all_soil.txt
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$22"\t"$23} NR>FNR{print $1,a[$1]}' result/39culture/otu.txt.bak fig1/ST/05.HIND-HTEJ_all.txt > fig1/ST/05.HIND-HTEJ_all_soil.txt

# 图4

# 图4a. 基因型与Nrt1.1b对应关系
cut -f 8,10,14 result/nitrogen_cor/design.txt | awk '{print $2,$1,$3}' | less
# IND为A，TEJ为G；找相反的如IND G有L4105,N4126；TEJ为A的研究中没有

# 图4对应的附图8 replication beta diversity重新绘图
mkdir -p fig1/3S
beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
        -d `pwd`/doc/design.txt -A groupID -B '"V3703LnCp6","ZH11LnCp6","A50HnCp7","A56HnCp7"' -E TRUE \
        -c `pwd`/doc/compare.txt \
        -o `pwd`/fig1/3S/ -h 3 -w 5
beta_cpcoa.sh -i `pwd`/result/otutab.txt -m '"bray","jaccard"' \
        -d `pwd`/doc/design.txt -A groupID -B '"V3703LnCp6","ZH11LnCp6","A50HnCp7","A56HnCp7"' -E TRUE \
        -o `pwd`/fig1/3S/  -h 3 -w 5

# 用新菌库注释 IND enriched overlap with NRT1.1b enriched，查看可培养比例
wget http://bailab.genetics.ac.cn/culture_collection/data/16S_rice_culture_collection.fasta -O culture_/16S_rice_culture_collection.fasta 
dos2unix culture_/16S_rice_culture_collection.fasta 
makeblastdb -in culture_/16S_rice_culture_collection.fasta -dbtype nucl

blastn -query result/otu.fa -db culture_/16S_rice_culture_collection.fasta -out temp/culture_otu.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 36
# 添加blastn结果表头，最主要前4列：OTUID，培养菌ID，相似度，覆盖度
sed -i '1 i OTUID\tsseqid\tpident\tqcovs\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' temp/culture_otu.blastn
# 统计可培养菌所占种类和丰度比例
echo -ne "Total OTUs\t" > culture_/summary.txt
grep '>' -c result/otu.fa >> culture_/summary.txt
echo -ne "Cultured OTUs\t" >> culture_/summary.txt
awk '$3>=97 && $4>=99' temp/culture_otu.blastn|wc -l >> culture_/summary.txt
# IND overlap NRT原来存在48个可培养的菌
tail -n+203 xiangeng_wilcoxon_main1/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls.xls|awk '$3>=97 && $4>=99' |wc -l
# 与培养菌只有30个可培养的菌
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' temp/culture_otu.blastn <(tail -n+203 xiangeng_wilcoxon_main1/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls.xls) | awk '$3>=97 && $4>=99' |wc -l
# TEJ 原来存在10个可培养的菌 diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls
tail -n 16 xiangeng_wilcoxon_main1/result/compare/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls|awk '$3>=97 && $4>=99' |wc -l
# 与培养菌只有9个可培养的菌
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' temp/culture_otu.blastn <(tail -n 16 xiangeng_wilcoxon_main1/result/compare/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls) | awk '$3>=97 && $4>=99' |wc -l



## 2019/1/27 IND和TEJ选菌与OTU关系详细，同 # 2019/1/15 测试nrt overlap IND enriched与选菌序列相似性

# 比对stock有多少可培养, stock_1492r.fa 函数truncate：取反向，并单行，截取V5-V7为 stock_1492r_v57.fa
cd ~/rice/xianGeng/select_bac
makeblastdb -in stock_1492r_v57.fa -dbtype nucl
# indica(IR24) isolate and japonica(Nipponbare) isolate
makeblastdb -in stock_ind.fa -dbtype nucl
makeblastdb -in stock_tej.fa -dbtype nucl

# IND富集与nrt下调仍意两者共有的OTU分别为3，46，70，共119，见附表8。
grep 'HTEJ_HIND_D_V3703HnCp6_ZH11HnCp6_D-specific_to_others' -A 1000 ../xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls.xls |cut -f 1|grep 'OTU' > ../wet/otuid119.id
wc -l ../wet/otuid119.id # 统计筛选部分数量119，确定一致性
usearch10 -fastx_getseqs ../result/otu.fa -labels ../wet/otuid119.id -fastaout ../wet/otuid119.fa
blastn -query ../wet/otuid119.fa -db stock_ind.fa -out indica_nrt_otu119.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 9 
awk '$3>=97' indica_nrt_otu119.blastn | wc -l # 30个stocks in All， IND(IR24) 28个

# TEJ富集与nrt下调仍意两者共有的OTU分别为3，46，70，共119，见附表8。
grep 'HTEJ_HIND_E_V3703HnCp6_ZH11HnCp6_D-specific_to_others' -A 1000 ../result/compare/diff.list.vennHTEJ_HIND_ELTEJ_LIND_EV3703HnCp6_ZH11HnCp6_D.xls.xls |cut -f 1|grep 'OTU' > ../wet/otuid16.id
wc -l ../wet/otuid16.id # 统计筛选部分数量16，确定一致性
usearch10 -fastx_getseqs ../result/otu.fa -labels ../wet/otuid16.id -fastaout ../wet/otuid16.fa
blastn -query ../wet/otuid16.fa -db stock_tej.fa -out indica_nrt_otu16.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 9 
awk '$3>=97' indica_nrt_otu16.blastn|wc -l # 11个stock, TEJ(Nipponbare)10个stocks



# 菌保与实验菌交叉验证

cd ~/rice/xianGeng/select_bac

# 1. Symcom使用菌
# select_bac/水稻syncom第一批测序结果.xlsx 中有20个菌的编号，对应关系和正反向, 保存1492R的序列为 syncom_1492r.fa
# 建立重组20个菌数据库
makeblastdb -in syncom_1492r.fa -dbtype nucl
# 函数truncate：取反向，并单行(存在冗余序列)，截取V5-V7为 syncom_1492r_v57.fa
makeblastdb -in syncom_1492r_v57.fa -dbtype nucl


# 2. Stock对应菌：提取对应菌库序列数据库
# 提取对应菌保编号 syncom_stock.id，原始菌保序列 stock_1492r.fa
usearch10 -fastx_getseqs stock_1492r.fa -labels syncom_stock.id -fastaout syncom_stock_1492r.fa
# 函数truncate：取反向，并单行，截取V5-V7为 syncom_stock_1492r_v57.fa
makeblastdb -in 2019/4/7 -dbtype nucl

# 3. OTU序列
usearch10 -fastx_getseqs ../result/otu.fa -labels otu.id -fastaout otu.fa
makeblastdb -in otu.fa -dbtype nucl
makeblastdb -in ../result/otu.fa -dbtype nucl

# 4. Rice COTU序列，仅查看rice_700 /mnt/bai/yongxin/culture/rice/result/culture_select.fasta


# 2vs1 比对序列菌保列表与实验是否一致（1492R）
blastn -query syncom_stock_1492r.fa -db syncom_1492r.fa -out syncom_stock_1492r.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 20 -evalue 1 -num_threads 9 
# 查看结果，发现I1/6/7/10，S20不是最佳匹配没有在结果中；放在上面num_alignments为20，又发现了I1/10正确；但I6和S20相差较多，且I7没有发现
sort -k2,2 -k3,3nr syncom_stock_1492r.blastn|less -S 
awk '$3>97' syncom_stock_1492r.blastn | wc -l

# 2vs1 比对序列: Symcom使用菌与Stock对应菌是否一致（V5-V7）
blastn -query syncom_stock_1492r_v57.fa -db syncom_1492r_v57.fa -out syncom_stock_1492r_v57.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 20 -evalue 1 -num_threads 9 
# 查看结果，发现I1/6/7/10，S20不是最佳匹配没有在结果中；放在上面num_alignments为20，又发现了I1/10正确；但I6和S20相差较多，且I7没有发现
sort -k2,2 -k3,3nr syncom_stock_1492r_v57.blastn|less -S 

# 3vs2 OTU对菌库v5-v7
blastn -query otu.fa -db syncom_stock_1492r_v57.fa -out otu.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 20 -evalue 1 -num_threads 9 
# 查看结果，只有OTU_95与R2217a相似不大于97%
sort -k2,2 -k3,3nr otu.blastn|less -S 

# 3vs4 OTU对应COTU
blastn -query otu.fa -db /mnt/bai/yongxin/culture/rice/result/culture_select.fasta -out otu_cotu.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1000 -evalue 1 -num_threads 9 
# 查看结果OTU95 100%对应rice_262，97.84%对应635，与rice700仅有77%相似
less otu_cotu.blastn

# 1vs3 实验菌vsOTU，鉴定错误的实验菌来自那些OTU?，找原因
blastn -query syncom_1492r_v57.fa -db ../result/otu.fa -out syncom_otu.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 100 -evalue 1 -num_threads 9 
# 查看结果OTU95 100%对应rice_262，97.84%对应635，与rice700仅有77%相似
less syncom_otu.blastn # OTU_377与I12最好，OTU_371

# 2vs3 实验菌vsOTU，鉴定错误的实验菌来自那些OTU?，找原因
blastn -query syncom_stock_1492r_v57.fa -db ../result/otu.fa -out syncom_stock_otu.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 100 -evalue 1 -num_threads 9 
# 查看结果OTU95 100%对应rice_262，97.84%对应635，与rice700仅有77%相似
less syncom_stock_otu.blastn


# 函数truncate，快速截取V7-V7, 799-1192
	# 设置起始文件名
	i=stock_1492r
	# 取反向互补
	revseq -sequence ${i}.fa -outseq /tmp/seq_rc.fa -tag N
	# 转换为单行
	format_fasta_1line.pl -i /tmp/seq_rc.fa -o /tmp/seq_rc1.fa
	# 检查序列冗余度，最好没有，有就存在重复序列
	grep -v '>' /tmp/seq_rc1.fa|sort|uniq -c|awk '$1>1'
	# 长度分布，正常为800-1200
	grep -v '>' /tmp/seq_rc1.fa|awk '{print length($0)}'|sort -n|uniq -c
	# 正向引物如799F AACMGGATTAGATACCCKG ，检索GGATTAGATACCC位于第3行前方>250bp
	cutadapt -g AACMGGATTAGATACCCKG -e 0.2 /tmp/seq_rc1.fa -o /tmp/seq_rc1.P5.fa
	# 反向时先切反向引物，1193R ACGTCATCCCCACCTTCC RC GGAAGGTGGGGATGACGT，正常为1492R GGTTACCTTGTTACGACTT找不到，改用1391R GACGGGCGGTGTGTRCA F TGYACACACCGCCCGTC
	cutadapt -a GGAAGGTGGGGATGACGT -e 0.2 /tmp/seq_rc1.P5.fa -o ${i}_v57.fa


# 图5 菌保

## Table S11统计OTU数量和Order数量
sed -i 's/^/R/' fig1/181118/TableS11-old.txt
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$9} NR>FNR{print $2,a[$2]}' fig1/181118/TableS11-old.txt fig1/190117/TableS11.txt | tail -n +3 | cut -f 2 | sort|uniq -c |wc -l # 438个OTU
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$9} NR>FNR{print $2,a[$2]}' fig1/181118/TableS11-old.txt fig1/190117/TableS11.txt | tail -n +3 | cut -f 3 | sort|uniq -c |wc -l # 22个目
tail -n+3 fig1/190117/TableS11.txt|cut -f 10|sort|uniq|wc -l # 非冗余1077个吗？
tail -n+3 fig1/181118/TableS11-old.txt|cut -f 13|sort|uniq|wc -l # 非冗余1096个吗？
# 原来文中截取1098个V5-V7非冗余序列519条，现在1079条是？？？，参考~/culture/rice/makefile.man 519部分
tail -n+3 fig1/190117/TableS11.txt|cut -f 1,10|sed 's/^/>/'|sed 's/\t/\n/' > temp/stock.fa
revseq -sequence temp/stock.fa -outseq temp/stock_rc.fa -tag N
        # 转换为单行
        format_fasta_1line.pl -i temp/stock_rc.fa -o temp/stock_rc1.fa
cutadapt -g AACMGGATTAGATACCCKG -e 0.2 temp/stock_rc1.fa -o temp/stock_rc1.P5.fa
cutadapt -a TGYACACACCGCCCGTC -e 0.2 temp/stock_rc1.P5.fa -o temp/stock_rc1.P3.fa
grep -v '>' temp/stock_rc1.P3.fa | sort|uniq -c|wc -l # 515条非冗余序列V5-V7，少了18个菌，少了4个非冗余序列

## 附网站

These culture collections also stock in the following National Culture Collection Center:

China National GeneBank (Shenzhen) http://www.accc.org.cn/

Agricultural Culture Collection of China (Beijing) http://www.accc.org.cn/


# 附表

## 附表2. 品种信息添加
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$7} NR>FNR{print $0,a[$2]}' ~/rice/miniCore/doc/minicore_list.txt fig1/190117/modified/TableS2.txt | sed '1 s/$/Variety/' > fig1/190117/modified/TableS2.csv
## 附表4. 更新完整的门纲数据









# 处理序列 Processing sequencing data

	# 1. 准备工作 Preparation

	## 1.1. 准备流程配置文件

	# Prepare config file of pipeline
	cd ~/github/Work/rice/xianGeng
	
	# 复制标准参数模板和操作指南至项目代码区：方便同步
	cp ~/github/Amplicon/16Sv2/parameter.md ./
	cp ~/github/Amplicon/16Sv2/manual.md ./
   
	# 链接代码至工作区
	ln -fs `pwd`/parameter.md ~/rice/xianGeng/makefile
	ln -fs `pwd`/manual.md ~/rice/xianGeng/manual.sh

	## 1.2. 初始化工作区

	# Initialize the working directory
	cd ~/rice/xianGeng
	make init

	## 1.3. 准备原始数据

	# Prepare raw data
	# 数据来自三个课题：miniCore + timecourse + nrt
	
	# 合并实验设计
	# 合并三个项目的实验设计，检查样品是否有重名
	cp ~/rice/miniCore/doc/design.txt doc/design_minicore.txt
	cp ~/rice/timecourse/doc/design.txt doc/design_timecourse.txt
	cp ~/rice/zjj.nitrogen/180116/doc/design.txt doc/design_nrt.txt
	# 统计实验样品行和唯一行，确实样品名唯一
	cat doc/design_* | grep -v 'SampleID'|wc -l
	cat doc/design_* | grep -v 'SampleID'|cut -f 1| sort|uniq|wc -l
	# 合并实验设计，前7列共有，只保留前7列
	cat <(head -n1 doc/design_nrt.txt) <(cat doc/design_* | grep -v 'SampleID') | cut -f 1-7 > doc/design.txt

	# 添加miniCore亚种
	mv doc/design.txt doc/design0.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$6]}' ~/rice/miniCore/180319/doc/design_group.txt doc/design0.txt|sed 's/\t$/\tsubspecies/'|less -S > doc/design1.txt # 添加亚种
	awk 'BEGIN{FS=OFS="\t"} {print $0,$7$8}' doc/design1.txt|less -S>doc/design2.txt # 合并土壤类型和亚种
	cp doc/design2.txt doc/design.txt


	# 原始数据合并
	cat ~/rice/miniCore/temp/seqs_usearch.fa ~/rice/timecourse/temp/seqs_usearch.fa ~/rice/zjj.nitrogen/180116/temp/seqs_usearch.fa | cut -f 1 -d ';' | sed 's/_/./g' > temp/filtered.fa
	# 从2.7 fa_unqiue 开始
	
	
## 按实验设计拆分lane为文库

	# Split lane into libraries
	# lane文件一般为seq/lane_1/2.fq.gz
	# lane文库信息doc/library.txt：至少包括编号、Index和样品数量三列和标题
	head -n3 doc/library.txt
	#LibraryID	IndexRC	Samples
	#L1	CTCAGA	60
	
	# 按library.txt拆分lane为library
	make lane_split


## 按实验设计拆分文库为样品

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
 

## 样品双端合并、重命名、合并为单一文件

	# Merge paired reads, renames and merge all samples
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt) <(cat doc/L* |grep -v '#') > doc/design.txt

	# 样品双端合并、重命名、合并为单一文件, 注意fastq为33格式，64位采用fastp转换
	make sample_merge


## 切除引物与标签

	# Cut primers and lables
	# 切除左端标签和引物，右端 引物
	# Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
	make fq_trim


## 质量控制

	# Quality control
	# 过滤序列中预期累计错误率>1%的序列
	make fq_qc


## 序列去冗余
	
	# 从这里开始
	ln ~/medicago/zjj170823/temp/seqs_usearch.fa temp/filtered.fa

	# Remove redundancy, get unique reads
	make fa_unqiue


## 挑选OTU

	# Pick OTUs
	# unoise3速度比cluster_otus慢上百倍
	make otu_pick


## 有参去嵌合体

	# Remove chimiras by silva database
	# 基于SILVA数据库去除
	make chimera_ref


## 去除宿主

	# Remove host
	# 根据SILVA注释去除线粒体、叶绿体、真核生物18S和未知序列(非rRNA)
	make host_rm


## 生成OTU表
	
	# Create OTUs table
	# 默认使用vsearch更快10倍，可选usearch10，线程不可超48
	make otutab_create


## 过滤样本和OTUs

	# OTU table filter samples and OTU
	# 推荐过滤低测序量<5000的样本，筛选大于1RPM的OTU
	make otutab_filter 
	# 此处结果统计只有3441个样品，而实验设计有3898个样品，少了哪些样品种？
	cat <(head -n1 result/otutab.txt|cut -f 2-|tr '\t' '\n') <(tail -n+2 doc/design.txt | cut -f 1) | sort | uniq -u > doc/missing_samples.txt 
	# 缺失 482样，主要是没有时间序列


## 物种注释

	# Assign taxonomy
	# 默认使用RDP trainset快而准，GG太旧，Silva太慢
	# 推荐阈值为0.6保证注释更完整
	make tax_assign


## 物种统计
	
	# Taxonomy summary
	# 必须所有物种有注释，否则有可能报错
	make tax_sum


## 多序列比对和进化树
	
	# Multiply alignment and make_phylogeny
	# usearch10/culsterO结果不同可能影响多样性分析(usearch unifrac结果更可信)
	# 进化树，用于树图和多样性分析
	make tree_make

## Alpha多样性指数计算
	
	# Calculate alpha diversity index
	# alpha指数计算结果为 result/alpha/index.txt
	# 稀释梯度结果位于 result/alpha/rare.txt
	make alpha_calc

## Beta多样性距离矩阵计算
	
	# Beta diversity tree and distance matrix
	# 最好用usearch，结果unifrac分类更好；clustero+fastree结果PCoA较差
	make beta_calc
	# ---Fatal error--- ../calcdistmxu.cpp(32) assert failed: QueryUniqueWordCount > 0 致信作者; 改用qiime1

## 有参考构建OTU表

	# Reference based OTU table
	# otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
	make otutab_gg



# 统计绘图 Statistics and plot

## Alpha多样性指数箱线图
	
	# Alpha index in boxplot
	make alpha_boxplot

## Alpha丰富度稀释曲线
	
	# Alpha rarefracation curve
	make alpha_rare

## 主坐标轴分析距离矩阵
	
	# PCoA of distance matrix
	make beta_pcoa

## 限制性主坐标轴分析

	# Constrained PCoA / CCA of bray distance matrix
	# OTU表基于bray距离和CCA，至少3个组 
	make beta_cpcoa

## 样品和组各级分类学堆叠柱状图

	# Stackplot showing taxonomy in each level
	make tax_stackplot

## 组间差异比较 
	
	# Group compareing by edgeR or wilcox
	# 可选负二项分布，或wilcoxon秩和检验
	make DA_compare


# 高级分析 Advanced analysis

## 添加可培养菌

	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $0,a[$1]}' temp/otutab.mean temp/culture_otu.blastn | cut -f 1-4,14 > result/41culture/otu.txt
	echo -ne "OTUs > 97% abundance :\t" >> result/41culture/summary.txt
	awk '$$3>=97 && $$4>=99' result/41culture/otu.txt | awk '{a=a+$$5} END {print a}' >> result/41culture/summary.txt
	cat result/41culture/summary.txt


## 查看Venn 4者共有菌，其中三者共有45，95，1，5
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/41culture/otu.txt result/venn/3_45.txt |sort -k3,3nr|less -S
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/41culture/otu.txt result/venn/3_95.txt |sort -k3,3nr|less -S

	# 制作有平均丰度，和物种注释的表
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $1,$5,a[$1]}' result/taxonomy_2.txt result/41culture/otu.txt > result/41culture/otu_mean_tax.txt

## 整理faprotax中菌的功能列表
	grep -v '^#' /mnt/bai/yongxin/software/FAPROTAX_1.1/FAPROTAX.txt|sed 's/\*/_/g' > culture/faprotax.tax


## 绘制网络
	
	# 制作数据文件
	OTU按差异比较计算整理
	注释信息按network.R整理
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print $0"\t"a[$1]}' result/taxonomy_8.txt <(cut -f 1,4 result/compare/database.txt) | cut -f 1-2,5- | sed '1 s/OTUID\tLTEJ/ID\ttotal\tphylum\tclass\torder\tfamily\tgenus\tspecies/' |less > network/database.txt

	# 绘制籼粳稻属水平模块和差异 co_network_genus_LIndTej.R
	# 以Burkholderia 相关菌 与 Anaeromyxobacter 的小网络 co_network_genus_LIndTej_core.R


## /5/14 随机森林属水平区分籼粳稻

	数据集 | 丰度% | 数量 |  错误率%
	LTEJ/LIND | 0.1 | 11.7 | 11
	LTEJ/LIND | 0.5 | 28 | 11
	LTEJ/LIND | 1 | 12 | 13.4

	HTEJ/HIND | 0.1 | 110 | 18.2
	HTEJ/HIND | 0.5 | 28 |18.9
	HTEJ/HIND | 1 | 13 | 18.9

	TEJ/IND | 0.1 | 116 | 14.1
	TEJ/IND | 0.5 | 33 | 15 # 选择全局，features少还兼顾大多数
	TEJ/IND | 1 | 12 | 17.4
 
	# 查看Top3属的分布 Geobacter Tangfeifania Anaeromyxobacter Bacillus Burkholderia

	for i in `tail -n+2 result/randomForest/imp.txt|cut -f 1`; do; \
	alpha_boxplot.sh -i result/tax/sum_g.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND"' -m \"${i}\" -t TRUE -o result/randomForest/box_ -n TRUE	
	done

	tail -n+2 result/randomForest/imp_c.txt|cut -f 1|awk '{print "\""$1"\""}'|tr "\n" ","

	alpha_boxplot.sh -i result/tax/sum_c.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND"' -m '"Flavobacteriia","Bacilli","Gammaproteobacteria","Spirochaetia","Clostridia","Alphaproteobacteria","Actinobacteria","Bacteroidia","Caldilineae","Betaproteobacteria","Ignavibacteria","Holophagae","Acidobacteria_Gp1","Nitrospira","Deltaproteobacteria"' -t TRUE -o result/randomForest/c_ -n TRUE	

	# 发现Top feature贡献度不大。改为科、目、纲时差异较明显。尤其是丰度阈值0.3%时，Top1/2为两类氮相关菌，绘制15个features

# 个性化分析 Custom analysis

## 低氮条件下挑选30个IND品种测宏基因组

	# 目标：挑选基于LN条件下Weighted Unifrac结果中IND与TEJ差异明显的的30个代表品种
	# 方法：筛选LN下IND/TEJ样品，按组合并，添加标签后筛选。
	mkdir -p pick_variety
	cd pick_variety
	# 筛选LN下IND/TEJ样品
	Rscript /mnt/bai/yongxin/rice/miniCore/180319/scripts/otutab_sample_subset.r -i ../result/otutab.txt -d /mnt/bai/yongxin/rice/miniCore/180319/LN/design.txt -o LN_otutab0.txt
	# 按品种合并
	head -n1 LN_otutab0.txt|cut -f 2-|tr '\t' '\n'|awk '{print $1"\t"$1}'|cut -c1-13|less > LN_sample_variety.list
	usearch10 -otutab_group LN_otutab0.txt -labels LN_sample_variety.list -output LN_otutab1.txt
	# 计算PCoA
	usearch10 -cluster_agg ../result/otu.fa -treeout otu.tree
	usearch10 -beta_div LN_otutab1.txt -tree otu.tree -filename_prefix LN_ -metrics bray_curtis,unifrac
	# 提取实验设计
	grep -P 'soiltypesubspecies|LIND|LTEJ' ../doc/design.txt | cut -f 6- |uniq | less -S > design.txt
	beta_pcoa.sh -i LN_ -m '"bray_curtis","unifrac"' -d design.txt -A subspecies -B '"IND","TEJ"' -o LN_pcoa_ -w 8 -h 5

	Rscript /mnt/bai/yongxin/rice/miniCore/180319/scripts/beta_pcoa_group.r -i LN/bray_curtis.txt -d doc/design.txt -n GroupID -o LN/pcoa_bray

	# 绘制土壤PCoA
beta_pcoa.sh -i result/beta/ -m '"bray_curtis","weighted_unifrac"' \
-B '"HSoil1","HSoil2","LSoil1","LSoil2"' -E TRUE \
-o pick_variety/soil_ -h 5 -w 8



## 挑选NRT样品测宏基因组

	# 按Anaeromyxobacter为新版OTU_11(1)，备选OTU_8(0.37)，必须 使用RDP注释，其它注释科faprotax无法识别
	# 查看OTU_11在CP 2016年中箱线图
	mkdir -p pick_nrt
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND","A50LnCp6","A56LnCp6","A50LnCp7","A56LnCp7","A50LnSz7","A56LnSz7","A50HnCp6","A56HnCp6","A50HnCp7","A56HnCp7","A50HnSz7","A56HnSz7","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","nrtHnCp7","ZH11HnCp7","ZH11LnCp7","nrtLnCp7","nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7","HSoil1","HSoil2","LSoil1","LSoil2","soilHnCp7","soilLnCp7","soilHnCp6","soilLnCp6","soilLnSz7","soilHnSz7"' -m '"OTU_8"' -t TRUE -o temp/all_ -n TRUE
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"A50HnCp6","A56HnCp6","A58HnCp6","IR24HnCp6","V3703HnCp6","ZH11HnCp6","A50LnCp6","A56LnCp6","A58LnCp6","IR24LnCp6","V3703LnCp6","ZH11LnCp6"' -m '"OTU_8"' -t TRUE -o pick_nrt/CP6_ -n true
	# OTU_8/11在CP7中变化
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","LTEJ","LIND","A50LnCp7","A56LnCp7","A58LnCp7","IR24LnCp7","A50HnCp7","A56HnCp7","A58HnCp7","IR24HnCp7","nrtHnCp7","ZH11HnCp7","ZH11LnCp7","nrtLnCp7","nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7"' -m '"OTU_11"' -t TRUE -o pick_nrt/CP7_ -n TRUE

	# OTU_8在目标群体中准确度和丰度远小于OTU_11

	# OTU_11在时间序列中变化
	mkdir -p timecourse
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"A50Cp0","A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119","A50Sz0","A50Sz1","A50Sz2","A50Sz3","A50Sz5","A50Sz7","A50Sz10","A50Sz13","A50Sz27","A50Sz34","A50Sz41","A50Sz48","A50Sz56","A50Sz62","A50Sz69","A50Sz76","A50Sz83","A50Sz90","A50Sz97","A50Sz118"' -m '"OTU_11"' -t TRUE -o timecourse/A50_ -n TRUE
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"IR24Cp0","IR24Cp1","IR24Cp2","IR24Cp3","IR24Cp7","IR24Cp10","IR24Cp14","IR24Cp21","IR24Cp28","IR24Cp35","IR24Cp42","IR24Cp49","IR24Cp63","IR24Cp70","IR24Cp77","IR24Cp84","IR24Cp91","IR24Cp98","IR24Cp112","IR24Cp119","IR24Sz0","IR24Sz1","IR24Sz2","IR24Sz3","IR24Sz5","IR24Sz7","IR24Sz10","IR24Sz13","IR24Sz27","IR24Sz34","IR24Sz41","IR24Sz48","IR24Sz56","IR24Sz62","IR24Sz69","IR24Sz76","IR24Sz83","IR24Sz90","IR24Sz97","IR24Sz118"' -m '"OTU_11"' -t TRUE -o timecourse/IR24_ -n TRUE

	# 绘制土壤PCoA
	beta_pcoa.sh -i result/beta/ -m '"bray_curtis","weighted_unifrac"' \
	-B '"soilHnCp7","soilLnCp7","soilHnCp6","soilLnCp6","soilLnSz7","soilHnSz7"' -E TRUE \
	-o pick_nrt/soil_ -h 5 -w 8
	
	# OTU_11 在土壤中组间丰度
	alpha_boxplot.sh -i result/otutab_norm.txt -d doc/design.txt -A groupID -B '"HSoil1","HSoil2","LSoil1","LSoil2","soilHnCp7","soilLnCp7","soilHnCp6","soilLnCp6","soilLnSz7","soilHnSz7"' -m '"OTU_11"' -t TRUE -o temp/soil_ -n TRUE


## 新发现OTU_8的物种很像OTU_11
	less result/faprotax/report# 但按物种注释0.6注释到科无法识别；改为0.3阈值


## 筛选各品种最好的3个样品

	# 结果备份，并筛选3个样品看结果
	cp -r result/ result180508/
	# 筛选实验设计为3个样品
	cp -r doc/design.txt doc/design.txt180508
	# 获得所有筛选样品
	cat ~/rice/miniCore/180319/HN/pcoa_bray_samples_top3.id <(cut -f 1 ~/rice/miniCore/180319/LN/pcoa_bray_samples_top3.group) > temp/temp1
	# 获得删除样品
	cat <(tail -n+2 doc/design_minicore.txt|cut -f 1) temp/temp1 | sort | uniq -u > doc/minicore_discard.id
	# 剔除点注释design
	for i in `cat doc/minicore_discard.id`; do;\
		sed -i "s/$i/#$i/" doc/design.txt;done


## OTU或分类与PCoA轴的差异
	
	# 以LN的weighted unifrac PCoA1/2为例
	# 基于PCoA轴与OTU计算相关
	script/cor_pcao_otu.R 计算OTU与4unifrac前4轴spearman相关系数，保存为 result/cor_otu_pcoa_unifrac.txt
	# 添加至差异OTUs，和相关系数
	awk 'NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1]"\t"$0}' result/cor_otu_pcoa_unifrac.txt result/compare/LTEJ-LIND_sig.txt | sed '1 s/^/Cor_PC1\tCor_PC2/' > result/compare/LTEJ-LIND_sig_pcoa_unifrac.txt
	# 选择代表菌作为组的markers，我推荐按丰度选差异；可进一步结果与PCoA轴的相关性；或选3/4组共有

	# 2018/5/16 国家、经纬度、地区和亚种与PCoA间关系script/beta_pcoa_location.R 发现Latitudeg与PC1相关，0.4366
	
	# 绘制某个分类单元的箱线图，常用绘制按丰度取上下调的top3；也可按Pvalue选择，但高丰度并不占优
	wd=`pwd`
	Dct_tax=g
	mkdir -p result/otu_boxplot_${Dct_tax}
	Dct_input=${wd}/result/tax/sum_${Dct_tax}.txt
	plot_boxplot_ggpubr.sh -i ${wd}/result/tax/sum_${Dct_tax}.txt -d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND"' \
	-m '"Anaeromyxobacter","Curvibacter","Sideroxydans","Burkholderia","Bradyrhizobium","Rhizobium"' -t TRUE -o result/otu_boxplot_${Dct_tax}/ -n true
	


## 检查GWAS是否可以发现Anaeromyxobacter与SNP的关联
	# 基于LN TEJ/IND中Anaeromyxobacter最低和最高的20个品种(保证数据可分开)，与nsSNP关联，最为仅为e-3，10m21759092,10,21759092,0.643641032317718,0.475,40,0.893753744914084,0.894396427377941,1。
	gapit_IndTej_Anaeromyxobacter_top20.R
	# 修改为30个品种有重合，肯定无法找到差异，再改为25个品种
	10m21759092,10,21759092,0.962618742911769,0.46,50,0.856299549567109,0.85630648832468,1
	# 可能从样本量、数据波动程度、SNP背景均无法满足

## /5/28 筛选HN/LN下可高丰度菌及是否可培养
	# core_microbiome_culture.R筛选中位数并排序 
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result/41culture/otu.txt culture/core.txt > culture/core_culture.txt



# 图表整理 Figures and legends

## 图1. 籼粳稻分型

### a, b 模式图
	
	地图，种植模式图——秦媛


### beta多样性
	
	# 不同加土，否则主要差异为土壤；不能两块地混合，否则主要差异为不同地块
	# HN下TEJ和IND
	beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
	-d `pwd`/doc/design.txt -A groupID -B '"HIND","HTEJ"' -E TRUE \
	-c `pwd`/doc/compare.txt \
	-o `pwd`/fig1/1subspecies/beta_HN_ -h 3 -w 5
	# 匹配非注释行，输出用于发表
	grep -P '^\s*#' script/beta_pcoa.R | less
	grep -P -v '^\s*#' script/beta_pcoa.R > fig1/script/beta_pcoa_fieldII.R
		# LN下TEJ和IND
	beta_pcoa.sh -i `pwd`/result/beta/ -m '"bray_curtis","weighted_unifrac","unweighted_unifrac"' \
		-d `pwd`/doc/design.txt -A groupID -B '"LIND","LTEJ"' -E TRUE \
		-c `pwd`/doc/compare.txt \
		-o `pwd`/fig1/1subspecies/beta_LN_ -h 3 -w 5
	grep -P -v '^\s*#' script/beta_pcoa.R > fig1/script/beta_pcoa_fieldI.R


### alpha多样性
	
	箱线图、稀释取线、样品稀释取线

	alpha_boxplot.sh -i `pwd`/result/alpha/index.txt -m '"chao1","richness","shannon_e"' \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' \
	-o `pwd`/fig1/1/alpha_ -h 3 -w 5
	# 稀释曲线
	alpha_rare.sh -i `pwd`/result/alpha/rare.txt \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' \
	-o `pwd`/fig/1/alpha_ -h 3 -w 5

### Taxonomy 门+变形菌纲
	cut -f 3-4 result/taxonomy_8.txt|sort|uniq|grep 'Proteobacteria' # 为什么会有这么多结果，只选5类继续分析
	cat <(grep -v 'Proteobacteria' result/tax/sum_p.txt) <(grep 'proteobacteria' result/tax/sum_c.txt) > result/tax/sum_pc.txt
tax_stackplot.sh -i `pwd`/result/tax/sum_ -m '"pc"' -n 10 \
	-d `pwd`/doc/design.txt -A groupID -B '"LTEJ","LIND","LSoil1","HTEJ","HIND","HSoil1"' -O FALSE \
	-o `pwd`/fig1/1/tax_pc_ -h 3 -w 5


### 门及纲水平差异


## 图2. 随机森林分类
	
	# 用family水平建模，用HN数据training，用LN验证。randomForest_family.R
	randomForest_class.R
	1. 纲水平建模，展示贡献度，和样品中热图
	使用高HN和HN下籼粳稻纲水平0.3%丰度的15个Feature机器学习；保存预测结果confusion.txt，整理16.4%错误率，TEJ 37.4%错误；
	2018/5/21 删除三个澳大亚利(纬度为负)粳稻, D4032, D4038, F4053; 标记A50/ZH11为TEJ，而IR24为IND，各分为Hn/Ln两种情况；
	错误率降低为15.3%，TEJ为36%;Top1 feature也变为了Nitrospira，Deltaproteobacteria
	2. Top feature：用各组柱状图/箱线图分类展示，再加梯度排序
	3. 在nrt和时间序列中验证

	
## 图3. 亚种差异与氮相关

	1. 差异OTUs曼哈顿图，维恩图
	由筛选组，改为筛选亚组(品种)中位数的OTUs: 原万5为343个OTUs，万一为942个；最终丰度为0.2%
	compare_sub.R # 修改丰度筛选group为groupID2，接下来 rm plot_volcano ; make plot_venn; make rmd

grep -P -v '^\s*#' script/compare_sub.R > fig1/script/compare.R

	差异OTUs在两块地曼哈顿、韦恩图;	曼哈顿图要写标颜色为门、纲,	plot_manhattan_pc.r
	plot_manhattan.sh -i result/compare/LTEJ-LIND_all.txt
	# 我们重点是突出IND，大多数是IND特异的，添加IND vs TEJ的组，重画曼哈顿图，让IND为向上实心三角

	# 曼哈顿图的代码和数据
	grep -P -v '^\s*#' script/plot_manhattan_pc.r > fig1/script/plot_manhattan_pc.r
	cp result/tax/sum_pc.txt fig1/data/
	cp result/compare/*IND-*TEJ_all.txt fig1/data/

	# 维恩图的代码和数据
	cp result/compare/diff.list* fig1/data/
	diff.list.vennHTEJ_HIND_DLTEJ_LIND_DCDE.r

	2. 差异OTUs在时间序列中变化
	alpha_boxplot.sh -i result/tax/sum_c.txt -d `pwd`/doc/design.txt -A groupID -B '"A50Cp0","A50Cp1","A50Cp2","A50Cp3","A50Cp7","A50Cp10","A50Cp14","A50Cp21","A50Cp28","A50Cp35","A50Cp42","A50Cp49","A50Cp63","A50Cp70","A50Cp77","A50Cp84","A50Cp91","A50Cp98","A50Cp112","A50Cp119"' \
	-m '"Deltaproteobacteria","Actinobacteria","Alphaproteobacteria","Clostridia","Betaproteobacteria","Nitrospira"' -t TRUE -o result/randomForest/time_ -n TRUE # 42以后没有数据呢？改用alpha_boxplot.sh

	3. 绘制维恩图的共有饼图
	# fig1/3compare.rmd 绘制时间序列的图，最后添加共有的成份
	# 筛选HL/LN下TEJ-IND共同下调的菌
	cat result/compare/?TEJ-?IND_sig.txt | grep 'Depleted' | cut -f 1 | sort | uniq -d > fig1/3compare/otu_IND_common_specific.txt
	cat result/compare/?TEJ-?IND_sig.txt | grep 'Enriched' | cut -f 1 | sort | uniq -d > fig1/3compare/otu_TEJ_common_specific.txt
	# 两块地保守上调、下调的OTUs
	cat fig1/3compare/otu_IND_common_specific.txt fig1/3compare/otu_TEJ_common_specific.txt > fig1/3compare/otu_common.txt
	# 并用faprotax注释
	# 绘制nrt和A50差异与籼粳稻共有
	tail -n 70 ~/rice/xianGeng/fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.txt | cut -f 1 > fig1/4nrt/venn_nrt_indiaHL.txt
	tail -n 6 ~/rice/xianGeng/fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D.txt | cut -f 1 > fig1/4nrt/venn_NRTsnp_indiaHL.txt

	5. 差异菌功能有无热图 plot_heatmap_timecourse.R
	filter_otus_by_sample.sh -f result/faprotax/element_tab.txt -o result/faprotax/xiaogeng -d doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"'
	# 结果可用STAMP进一步探索
	
	# 筛选时间序列中上/下调两大类的OTU进行功能分析
	mkdir -p fig/3
	awk '$2>0' fig/2/otu_IND_common_specific_time_cor6.txt | cut -f 1 > fig/3/timecournse_increase.id
	awk '$2<0' fig/2/otu_IND_common_specific_time_cor6.txt | cut -f 1 > fig/3/timecournse_decrease.id
	# 以Incease为例, decrease
	type=decrease
	filter_otus_from_otu_table.py -i result/otutab_norm_tax.biom -o timecourse/${type}.biom --otu_ids_to_exclude_fp fig/3/timecournse_${type}.id --negate_ids_to_exclude
	/usr/bin/python2.7 /mnt/bai/yongxin/software/FAPROTAX_1.1/collapse_table.py -i timecourse/${type}.biom -o timecourse/${type}.faprotax -g /mnt/bai/yongxin/software/FAPROTAX_1.1/FAPROTAX.txt --collapse_by_metadata 'taxonomy' -v --force # --out_report result/faprotax/report 
	filter_otus_by_sample.sh -f timecourse/${type}.faprotax -o fig/3/timecournse_faprotax_${type} -d doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"'
	# increase的差异功能类型均为IND>TEJ>soil，且以芳香、氮 相关
	compare.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
	-p 0.01 -q 0.01 -F 1.2 -t 0.0005 \
	-d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1"' \
	-o `pwd`/result/compare/
	# decrease的差异，stamp打开报错，但有时成功；下调无N循环相关，有
	
	# 可视化菌的功能有无
	## 筛选report为功能有无表
	grep 'OTU_' -B 1 result/faprotax/report | grep -v -P '^--$' > result/faprotax/report.clean
	faprotax_report_sum.pl -i result/faprotax/report.clean -o result/faprotax/report
	#OTU功能注释列表：result/faprotax/report.otu_func
	#功能包含OTU列表：result/faprotax/report.func_otu
	#OTU功能有无矩阵：result/faprotax/report.mat
	# plot_heatmap_timecourse.R 绘制时间序列的图，再添加相应菌的主要功能，
	# 同时对时间序列中不表达的也可视化功能:IND的功能绘制于 fig/2/otu_IND_common_specific_time_faprotax_noabundance.txt，TEJ单一条目录为 fig/2/otu_TEJ_common_specific_time_faprotax_noabundance.txt"
	# 再对时间序列中0点去掉重新计算，发现分为了4组，在原文件基础上添加-0标志


	# 菌种功能注释整体差异
	rm result/compare_far/diff.list
	compare.sh -i `pwd`/result/faprotax/element_tab.txt -c `pwd`/doc/compare.txt -m "wilcox" \
	-p 0.01 -q 0.05 -F 1.2 -t 0.0005 \
	-d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1","V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50LnCp6","A56LnCp6"' \
	-o `pwd`/result/compare_far/ -N FALSE
	batch_venn.pl -i doc/venn.txt -d result/compare_far/diff.list
	# 注释比较结果
	rm result/compare_far/diff.list.venn*.xls.*
	batch2.pl -i 'result/compare_far/diff.list.venn*.xls' -d result/compare_far/database.txt -o result/compare_far/ -p vennNumAnno.pl
	# 绘制箱线图
	# 确定要展示的Features，两组共有
	tail -n 39 result/compare_far/diff.list.vennHTEJ_HIND_DLTEJ_LIND_D.xls.xls|cut -f 1 > result/compare_far/IND.list
	tail -n 9 result/compare_far/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls|cut -f 1 > result/compare_far/TEJ.list
	make plot_fa_barplot # 绘制单个功能的箱线图
	# 修改alpha_boxplot.R为alpha_boxplot_far.R

    grep -P -v '^\s*#' script/alpha_boxplot_far.R > fig1/script/alpha_boxplot_far.R

## 菌群与nrt关系

### 差异菌/功能与氮相关基因显著相关
	
	# 来自胡斌整理的氮相关基因doc/N-related genes in rice.docx共9个基因，先在doc/rice_nitrogen_list.xlsx中惠惠相关ID，保存为doc/rice_nitrogen_list.txt, 其中第5列RAP_SNP的ID与SNP注释文件对应
	dos2unix doc/rice_nitrogen_list.txt # windows转换为linux
	# 整理出SNP数据中对应的基因型、提取相应的位点
	mkdir -p result/nitrogen_cor
	cut -f 5 doc/rice_nitrogen_list.txt | tail -n+2 | tr '\n' '|' # 提取ID并替换为|分隔
	grep -P 'OS10G0554200|OS08G0155400|OS02G0112100|OS02G0595900|OS01G0704100|OS01G0547600|OS03G0687000|OS04G0509600|OS06G0706400|OS06G0706500' /mnt/bai/yongxin/rice/miniCore/180319/gemma/snp.anno > result/nitrogen_cor/all_snp.list # 筛选到10个基因在miniCore中存在502个相关位点
	cut -f 3 result/nitrogen_cor/all_snp.list | sort | uniq -c # 8个MODERATE，467个MODIFIER和27个LOW
	grep -P 'HIGH|MODERATE' result/nitrogen_cor/all_snp.list > result/nitrogen_cor/good_snp.list # 其中重要SNP仅有8个，来自4个基因
	# 提SNP对应基因型
	cat /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_*.hmp.txt > /tmp/temp
	cut -f 1 result/nitrogen_cor/good_snp.list|tr '\n' '|' # 获取列表
	grep -P '2m655515\t|2m657013|6m29839102|6m29839240|8m3183208|10m21759092|10m21761740|10m21761997' /tmp/temp > /tmp/temp1
	cat <(head -n1 /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_1.hmp.txt) /tmp/temp1 > result/nitrogen_cor/good_snp.geno # 添加标题
	# 用excel转置 result/nitrogen_cor/good_snp.geno.t，添加注释
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$0} NR>FNR{print $0,a[$1]}' ../miniCore/doc/minicore_list.txt result/nitrogen_cor/good_snp.geno.t > result/nitrogen_cor/good_snp.geno.t.txt


	# 统计SNP基因型作为分组信息，来统计氮功能丰度组间P值
	nitrogen_cor.r # 保存实验设计+基因型，方便识别SNP不同基因型在籼粳稻中区别
	# 整理这4个重要SNP信息表，见SNP_list.xlsx
	# 统计基因型与亚种分布
	sed -i '1 s/^/SampleID\t/' result/nitrogen_cor/design.txt
	head -n1 result/nitrogen_cor/design.txt|tr '\t' '\n'|awk '{print NR,$1}'
	# NRT2.1 - 17; 1.1A - 21; 1.1B - 14，统计每个亚种内基因型的数量，可看到亚种内主要的SNP类型
	cut -f 2,8,20 result/nitrogen_cor/design.txt|sort|uniq|cut -f 2,3|sort|uniq -c
	grep -P '10m21759092|2m655515\t|8m3183208' /mnt/bai/yongxin/rice/miniCore/180319/gemma/T2.ann.vcf # 查询SNP变化位置、碱基和AA详细

	grep -P -v '^\s*#' script/nitrogen_cor.r > fig1/script/nitrogen_cor.R


### 关键氮高效基因NRT不同形态、突变体可部分解析亚种差异

	"A50LnCp6","A56LnCp6","A50LnCp7","A56LnCp7","A50LnSz7","A56LnSz7","A50HnCp6","A56HnCp6","A50HnCp7","A56HnCp7","A50HnSz7","A56HnSz7","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","nrtHnCp7","ZH11HnCp7","ZH11LnCp7","nrtLnCp7","nrtHnSz7","ZH11HnSz7","nrtLnSz7","ZH11LnSz7"
	# nrt vs ZH11(TEJ): 主图："V3703HnCp6","ZH11HnCp6", 附图："V3703LnCp6","ZH11LnCp6",
	# 近等基因系：主图："A50LnCp7","A56LnCp7", 附图："A50LnCp6","A56LnCp6",
	# 主图/附图各分析一次：alpha, beta, 差异OTUs, venn: HTEJ_HIND_D LTEJ_LIND_D V3703HnCp6_ZH11HnCp6; HTEJ_HIND_D LTEJ_LIND_D A50HnCp7_A56HnCp7

主图4组: http://bailab.genetics.ac.cn/report/16Sv2/xiangeng_wilcoxon_main/

附图4组: http://bailab.genetics.ac.cn/report/16Sv2xiangeng_wilcoxon_supp

主图3.10.2和3.10.3中共有菌数量但不，但丰度很高，计算其占全部菌的相对丰度。详见fig1/venn_overlap_abundance.R



### 宏基因组KO注释

	从金桃处获得KO表和丰度，获取KO的功能描述；在kegg中没找到，picurst中没找到(输出结果有注释但不完整)，google搜索KO description download，找到biostar解答；https://www.genome.jp/kegg-bin/get_htext?ko00001.keg 中的Download htext/jason下载KO和描述；htext方便检索，而jason方便在线分析，如jason2table
	
	cd ~/rice/xianGeng/meta
	grep 'D      K' ko00001.keg | cut -c 8- | sed 's/  /\t/' > KO_description.txt #|  cat -A | less -S
	# 比较两组差异
    compare.sh -i ko.txt -c compare.txt -m "wilcox"         -p 0.01 -q 0.05 -F 1.2 -t 0 -N FALSE         -d design.txt -A group -B '"HnNrt","HnZH11"'         -o compare/ -N FALSE -U 100 # Pvalue和FDR并不显著，可能秩和检验需要较多的样本数
	# 注释KO
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' KO_description.txt ko.txt > ko.anno

	# 绘制单个KO和箱线图
    cd ~/rice/xianGeng/meta/nitrogen/
	alpha_boxplot.sh -i ko.txt -m '"K02568","K02567","K00363","K10535","K00362"' \
	-d design.txt -A group -B '"HnZH11","HnNrt"' \
	-o ./ -h 2 -w 2.5 -t TRUE -n TRUE -U 1000000
    K02567 NapA;K02568 NapB;K00363 NirD

    # 2018/11/14 完整KO表
	cd ~/rice/xianGeng/meta/compare/
    # 制作kotab.txt, design.txt, compare.txt
    rm -r compare
    compare.sh -i kotab.txt -d design.txt -c compare.txt \
        -m "ttest" -p 0.05 -q 1 -F 1.2 -t 0 -N FALSE -A GroupID -B '"HnZH11","HnNrt","LnZH11","LnNrt"' \
        -o compare/ -N FALSE -U 100 
    cat compare/summary.txt
        # Pvalue和FDR并不显著，可能秩和检验需要较多的样本数，HN要比LN明显特别多

    ## https://www.kegg.jp 搜索Nitrogen metabolism ，打开PATHWAY map00910，点击Ortholog table打开KO与物种对应表，保存第一行KO和对应描述为nitrogene/ko.id
    # 设置分析本批关键字，同类分析替换即可 nitrogen/sulfur 有代谢通路, potassium 只有M00454中两个KO 	K07646 K07667, phosphorus 相关KO由张鹏帆贡献
    type=sulfur
    mkdir -p ${type}
    # 转换为KO ID对应基因名称，并添加差异比较对应表头
    sed ":a;N;s/\n/ /g;ta" ${type}/ko_raw.id |sed "s/\t/\n/g"| sed "s/ /\t/;s/(//;s/)//;" | cut -f 1 -d '['| sed '1 i HnZH11_HnNrt'|less > ${type}/ko.id
    # 筛选KO对应比较结果
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' compare/HnZH11-HnNrt_all.txt ${type}/ko.id | grep -v -P '^$' > ${type}/HnNrt_all.txt
    # 添加注释
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$7]=$8} NR>FNR{print $0,a[$1]}' ~/github/Metagenome/denovo1/kegg/ko00001.tsv  ${type}/HnNrt_all.txt > ${type}/HnNrt_all_anno.txt
    # 查看正对照K02568 K02567 K00363
    # 筛选显著差异并标颜色，用于在线展示
    # awk '$4<0.05' ${type}/HnNrt_all.txt | cut -f 1-4,7-8 | awk '{if ($2>0) {print $0"\tred"}else {print $0"\tgreen"}}' > ${type}/HnNrt_diff.txt
    # cat ${type}/HnNrt_diff.txt
    # KO添加描述
    # head -n1 ~/github/Metagenome/denovo1/kegg/ko00001.tsv|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$7]=$8} NR>FNR{print $0,a[$1]}' ~/github/Metagenome/denovo1/kegg/ko00001.tsv ${type}/HnNrt_diff.txt > ${type}/HnNrt_diff_anno.txt
    # N和S都有上升的基因，钾无差异，P有一个下调


## 图5. 可培养菌

    # 可培养比例与自然样品  ~/culture/jingmei171208.rice.med.soy/makefile.man # 2018/10/17 培养菌与自然样本比较
    # 统计科数量
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/rep_seqs_tax.txt culture_A50/otu_cultured.txt > culture_A50/otu_cultured.tax
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/rep_seqs_tax.txt culture_IR24/otu_cultured.txt > culture_IR24/otu_cultured.tax
    cat culture_A50/otu_cultured.tax culture_IR24/otu_cultured.tax | cut -f 6 | sort | uniq | wc -l
    # 分菌见 ~/culture/rice/makefile.man # 2018/10/19 水稻分菌整理
    # 总结：1098个单菌，去掉一个序列异常的1097，去冗余519，序列比对剩512条序列



## 图5. 微生物与表型关联——2018/11/19去掉

	1. 微生物、多样性、PCoA主轴与表型相关分析
	# 先使用胡斌整理表型数据+faprotax中氮通路相关
	# 方法1：script/phenotype_cor.R直接关联，spearman相关系数只有0-0.2，但能看到nitrogen_amonification正相关，而固
	
	# 方法2. 采用分组协变量关联，需要基因型的聚类/PCA信息
	# PCA信息 /mnt/bai/yongxin/rice/miniCore/180319/gemma/pca4.txt
	# 行名 head -n1 /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_1.hmp.txt
	head -n1 /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_1.hmp.txt|cut -f 12-|tr '\t' '\n' > temp.txt
	paste temp.txt /mnt/bai/yongxin/rice/miniCore/180319/gemma/pca4.txt | cut -f 1,3- | sed '1 i variety\tPC1\tPC2\tPC3\tPC4' | less > fig/4cor/genotype_pca4.txt
	
	# OTU/属水平(比OTU数量少至有描述)与张小宁整理表型关联(品种对应): HN/LN分开关联phenotype_cor.Rmd
	cp ~/rice/miniCore/mwas/phenotype/minicore?NPhenotype.txt doc/ # 准备miniCore中HN/LN表型
	
	# 表型-faprotax相关 phenotype_cor_faprotax.Rmd

	# 2018/6/5 相关热图+注释
	# LN条件下OTUs与表型关联，筛选>0.4相关的值进行注释物种，差异OTUs，和功能 phenotype_cor2.Rmd ；表型数据重新整理minicore低氮下为单株，excel整理
	sed -i '/\t0\t0\t0/d' doc/phenotype_sample_raw.txt # 删除缺失样品
	# 株的OTU对应株的表型，OTU_11与tiller相关仅为0.33(可能植株波动大，或不对应，会规律只有平均才能看出来)，而OTU_11与分均值还有0.44的相关。那OTUs的均值对应是否会更高呢？改为均值结果更好。

	# 2018/6/11 ## 筛选相关系数 > 0.4，# 添加faprotax有无 X OTU丰度，绘制泡泡图，line 310
	
	# 2018/6/22 用品种中全部OTUs计算相关系数，筛选相关系数>0.4，丰度大于均值0.1%的OTUs，物种名改为低级注释
	cp xiangeng4wilcox/result/compare/LTEJ-LIND_all.txt fig1/5pheno/ # 品种0.2%过滤OTU列表

    # 2018/11/16 株高、鲜重与OTU相关性 fig1\5pheno\LN_otu_mean_pheno_cor.r.txt，用echart可视化为LN_otu_mean_pheno_cor.r.pdf


## syncom重组体系表型分析
    # 数据整理为 wet/Xu-NbSycomFunctionData.txt, 分析代码见 fig1/5Syncom.Rmd



## 附图

1. 样本稀释曲线，平滑处理
alpha_rare_sample.R
2. 2块地所有样品IND/TEJ一起PCoA, CPCoA
3. 主图c/d的unifrac距离图
4. alpha diversity chao1 richness
5. 机器学习层级选择准确度柱状图
6. 共有OTUs饼图
7. 完整的OTUs时间序列和功能注释
8. 主图的功能差异分析完整版本
9. 图4a的，其它N相关功能与SNP关联
10. 图4b的beta其它距离， constraind
11. 附图专用4组多样性分析

## 附数据统计

	# 1. 品种数：68 IND
	cd ~/rice/xianGeng/fig1/ST
		27 TEJ
	sed -i 's/,/\t/g' 01.variety_geotable.csv # csv替换为tsv
	# 数据筛选95个，并添加95个
	awk '{FS=OFS="\t"} NR==FNR{a[$2]=$0} NR>FNR{print $0,a[$2]}' ~/rice/miniCore/doc/minicore_list.txt 01.variety_geotable.csv > 01.variety_geotable.txt 
	cut -f 14 01.variety_geotable.txt|tail -n+2|sort|uniq|wc -l # 44个国家

	# 2. 凌水海南样本数和测序量
	# 需要使用compare_sub.R中来原代码，在3compare中汇总样品和测序量

	# 3. 查看共有16个下调菌的丰度范围
	cd ~/rice/xianGeng
	awk '{FS=OFS="\t"} NR==FNR{a[$1]=$14} NR>FNR{print a[$1]}' xiangeng4wilcox/result/compare/HTEJ-HIND_all.txt fig1/3compare/otu_TEJ_common_specific.txt|sort -nr
	awk '{FS=OFS="\t"} NR==FNR{a[$1]=$14} NR>FNR{print a[$1]}' xiangeng4wilcox/result/compare/LTEJ-LIND_all.txt fig1/3compare/otu_TEJ_common_specific.txt|sort -nr
	awk '$2<0.05' fig1/ST/09.tiller_cor_p.txt|wc -l # 381个有214个显著P<0.05
	cd fig1/ST
	paste 09.tiller_cor.txt 09.tiller_cor_p.txt | cut -f 1,2,4> 09.tiller_cor_pr.txt

## 上传数据至NCBI PRJNA478068

	# 整理的最终上传实验设计 fig1/metadata.txt
	cut -f 16 fig1/metadata.txt|sort|uniq -c # 587个亚种数据来自miniCore，127+114=241来自nrt
	wc -l fig1/metadata.txt # 检查是否有样品重名
	cut -f 1 fig1/metadata.txt|sort|uniq|wc -l 
	mkdir -p seq/submit # 创建提交数据目录

	# 1. 籼粳稻样本，587个
	for RPM in `grep 'minicore' fig1/ST/02.design.txt|cut -f 1`; do
		cp ../miniCore/clean_data/sample/${RPM}.fq.gz seq/submit/
	done

	# 2. 时间序列己上传，见 https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA435900
	for RPM in `grep 'nrt' fig1/ST/02.design.txt|cut -f 1`; do
		cp /mnt/bai/yongxin/rice/zjj.nitrogen/180116/clean_data/sample/${RPM}.fq.gz seq/submit/
	done


## 上传数据至基因组所 PRJCA001214

# 宏基因组数据上传，共36个样品，其中只比较了HnZH11 vs HnNrt共6个样品，原始212GB，过滤后172GB [yongxin@meta:~/rice/nrt1.1b/submit]$

# 扩增子数据，828个样品上传

	# 整理的最终上传实验设计 fig1/metadata.txt
	cut -f 16 fig1/metadata.txt|sort|uniq -c # 587个亚种数据来自miniCore，127(随机森林)+114=241来自nrt
	wc -l fig1/metadata.txt # 检查是否有样品重名
	cut -f 1 fig1/metadata.txt|sort|uniq|wc -l 
	mkdir -p seq/submitGSA # 创建提交数据目录

	# 1. 籼粳稻样本，587个
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

+--------+-----------+----------------------------------------------------------------------------------+
| cra_id | accession | alias                                                                            |
+--------+-----------+----------------------------------------------------------------------------------+
|   1357 | CRA001362 | Rice nrt1.1b related root metagenome                                             |
|   1369 | CRA001372 | Rice subspecies indica and japonica, and NRT1.1b related 16S amplicon sequencing |
+--------+-----------+----------------------------------------------------------------------------------+

# 共享中间文件和分析流程代码 submit用于投稿，publish用于正式发表后共享

	## 在fig1中创建index.Rmd并生成网页 http://210.75.224.110/submit/rice_microbiome, username: rice, password: microbiome
	cd ~/rice/xianGeng/fig1 # 共享目录
	ln -sf `pwd` /var/www/html/submit/rice_microbiome # 链接至外网
	cp ~/github/Amplicon/16Sv2/rmd/.htaccess ./ # 加密
	htpasswd /mnt/bai/yongxin/bin/config/users rice # 添加新用户和密码

	## 获得分析流程
	pipeline=fig1/pipeline.sh
	make -n -B library_split_stat|grep -v '#' > $pipeline
	make -n -B fq_qc|grep -v '#' >> $pipeline
	make -n -B beta_calc|grep -v '#' >> $pipeline

	## 根据最新实验设计和OTU表整理结果
	cd fig1/
	cp /mnt/bai/yongxin/github/Amplicon/16Sv2/script/stat_plot_functions.R script/






## 补充结果2018/10/29 

	目前文章需要补充的图和表

	主图

	1.分菌详细流程图：我基金中的那个流程图吗？
	原图参考2018面上项目基金D:\life\grant\NSFC2018\figure\figure.pptx中图3，复制到fig1/20180610-组合图.pptx中：和白老师定方案，秦媛绘制

	2.菌保图：树吗？
	~/culture/rice/makefile.man 中 	# 绘制物种树，不带丰度和可培养 ，绘制进化树，无论科62，目22，都不能确保位于树的同一分枝上。
	改为graphlan绘制物种树，并添加标签，保存为fig1/graphlan.ai/pdf
	# graphlan绘制培养菌对应样本比例 /mnt/bai/yongxin/culture/jingmei171208.rice.med.soy/makefile.16s中sub分别为A50和IR24绘制；
	图5c中的圈图绘制代码 ~/culture/rice/verify	

	附图

	1. 无机氮功能实验结果
	2. A50实验结果（同IR24，数据）Figure6S.pdf

	Table

	1.分离菌占16S结果比例及taxonomy组成：自然样品，先不整理

	2.菌保信息对应表：verified & COTU，详见clutre/rice/makefile.man中# 制作isolate(鉴定结果)物种注释、来源和序列表 181024/TableS10.xls	# 制作COUT注释和序列信息 181024/TableS9.xls

	3.功能实验设计及表型结果（包括有机氮/无机氮，培养基成分，浓度，sycoms菌组成，表型数据等）



## 补充结果2018/11/17

1. IND/TEJ差异共有OTU在表型相关系数；

共有差异OTU列表，见： http://bailab.genetics.ac.cn/report/16Sv2/xiangeng_wilcoxon_main/result-otu.html#htej_hind_dltej_lind_dcde-5 http://bailab.genetics.ac.cn/report/16Sv2/xiangeng_wilcoxon_main/result-otu.html#htej_hind_eltej_lind_ecde-5

    wd=result/subspecies_otu_cor
    mkdir -p $wd
    # 筛选IND/TEJ富集的OTU
    tail -n 141 xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_D.xls.xls | awk '{print $0"\tIND_enriched"}'> $wd/IND_enriched_otu.id
    tail -n 16 xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_ELTEJ_LIND_E.xls.xls | awk '{print $0"\tTEJ_enriched"}'> $wd/TEJ_enriched_otu.id
    cat $wd/IND_enriched_otu.id $wd/TEJ_enriched_otu.id | sed '1 i OTUID\tAbundance\tTaxonomy\tType' > $wd/IND_TEJ_otu.id
    # 添加OTU相关
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' ./fig1/5pheno/LN_otu_mean_pheno_cor.r.txt $wd/IND_TEJ_otu.id | less -S > $wd/IND_TEJ_otu.cor # 添加亚种
    # 在ehcart绘制热图 http://www.ehbio.com/ImageGP/index.php/Home/Index/PHeatmap.html

2. Venn中共有比例

共有和特有OTU列表见: http://bailab.genetics.ac.cn/report/16Sv2/xiangeng_wilcoxon_main/result-otu.html#htej_hind_dltej_lind_dv3703hncp6_zh11hncp6_dde-5
    
    wd=result/nrt_venn_abundance
    mkdir -p $wd
    # 脚本见 fig1/venn_overlap_abundance.R
    cp xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D.xls.xls $wd/HTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D.xls
    cp xiangeng_wilcoxon_main/result/compare/diff.list.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls.xls $wd/HTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.xls
    # 替换名称
    # 整理为附表 TableSX.nrt_veen_list.xlsx
    

3. Figure S10统计

    # 检索原代码
    /rice/xianGeng/fig1]$grep 'Tiller number' *.Rmd
    # 5phenotype_cor.Rmd 中有代码和图，带wilcoxon检验的P值

4. Figure 3.I/J 统计注释菌的总丰度，注释到N的总丰度；
    # ../result/faprotax/report.mat在菌和功能对应表
    # ../result/41culture/otu_mean_tax.txt 丰度和注释
    # 查看注释菌的丰度
    awk '{a=a+$2} END {print a}' ../result/41culture/otu_mean_tax.txt # 总量99.752%
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1]}' ../result/41culture/otu_mean_tax.txt ../result/faprotax/report.mat | awk '{a=a+$1} END {print a}' # 注释的菌有62.972
    # 用R脚本统计每类的丰度，和几类总和的丰度 sum_faprotax_anno.R，补充为3comapre_en.rmd中3E段落中末属


5. 机器学习方法写清楚：使用了那个包，调用了哪个函数。
    /mnt/bai/yongxin/rice/xianGeng/fig1/随机森林分类方法描述.docx

6. Nitrogene related KO STAMP: Welch's t-test, Storey FDR


# 实验绘图

    otu=`cat wet/181220functional.bac.id|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'`

	alpha_boxplot.sh -i result/otutab.txt -d doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1","V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50LnCp6","A56LnCp6"' \
	-m $otu -t TRUE -o result/otu_boxplot/fuc_ -n TRUE 


# 2019/1/5 NBT审稿意义补充

文章结果： http://210.75.224.110/submit/rice_microbiome/

原主图：

## IND和TEJ富集菌在根中富集吗？分别比较HN和HN下水稻vsSoil
	
	# 修改实验设计，添加新比较组
rm -fr `pwd`/result/compare/
mkdir -p `pwd`/result/compare/
	# 丰度默认0.2，结果从381个OTU变为了224，代码与原来一样，可能compare_OTU.sh修改，改用compare.sh也不行？查也没找到原因；改为0.02也有很多无法对应；改有有参方法
compare_OTU_ref.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/compare_soil.txt -m "wilcox" \
        -p 0.01 -q 0.05 -F 1.2 -t 0.2 \
        -d `pwd`/doc/design.txt -A groupID -B '"HTEJ","HIND","HSoil1","LTEJ","LIND","LSoil1","V3703HnCp6","ZH11HnCp6","V3703LnCp6","ZH11LnCp6","A50LnCp7","A56LnCp7","A50LnCp6","A56LnCp6"' \
        -o `pwd`/result/compare/ -r xiangeng_wilcoxon_main/result/compare/LIND-LTEJ_all.txt
	# 添加LN下比较土壤LN下情况
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print a[$1],$0}' result/compare/LIND-LSoil1_all.txt xiangeng_wilcoxon_main/result/compare/LIND-LTEJ_sig.txt | sed '1 s/^/LIND-LSoil/' >  xiangeng_wilcoxon_main/result/compare/LIND-LTEJ_sig1.txt
	cut -f 1 xiangeng_wilcoxon_main/result/compare/LIND-LTEJ_sig1.txt|sort|uniq -c # LN下相比土壤133下调，111上调，28个不显著
	cut -f 7 xiangeng_wilcoxon_main/result/compare/LIND-LTEJ_sig1.txt|sort|uniq -c # LN下相比TEJ 34下调，244上调
	awk '$7=="Enriched"' xiangeng_wilcoxon_main/result/compare/LIND-LTEJ_sig1.txt|cut -f 1 |sort|uniq -c # LN下IND富集的244个菌，139下调，82上调，23个缺失

	# 追加特异的 fig1/ST/05.HIND-HTEJ_all.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print a[$1],$0}' result/compare/HIND-HSoil1_all.txt xiangeng_wilcoxon_main/result/compare/HIND-HTEJ_sig.txt | sed '1 s/^/HIND-HSoil1/' >  xiangeng_wilcoxon_main/result/compare/HIND-HIND-HTEJ_sig1.txt
	cut -f 1 xiangeng_wilcoxon_main/result/compare/HIND-HIND-HTEJ_sig1.txt|sort|uniq -c # HN下79下调，73上调，16个不显著
	# 筛选IND富集的
	cut -f 7 xiangeng_wilcoxon_main/result/compare/HIND-HIND-HTEJ_sig1.txt|sort|uniq -c # 147上调，71下调
	awk '$7=="Enriched"' xiangeng_wilcoxon_main/result/compare/HIND-HIND-HTEJ_sig1.txt|cut -f 1 |sort|uniq -c # HN下IND富集的147个菌，79下调，56上调，12个缺失

    # 查看共有部分土壤中富集情况
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print a[$1],$0}' result/compare/HIND-HSoil1_all.txt fig1/ST/05.vennHTEJ_HIND_ELTEJ_LIND_E.txt > fig1/ST/05.vennHTEJ_HIND_ELTEJ_LIND_E_HNsoil.txt
    tail -n 16 fig1/ST/05.vennHTEJ_HIND_ELTEJ_LIND_E_HNsoil.txt | cut -f 1 |sort|uniq -c # IND共有的16个14个富含，2个不显著
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print a[$1],$0}' result/compare/HIND-HSoil1_all.txt fig1/ST/05.vennHTEJ_HIND_DLTEJ_LIND_D.txt > fig1/ST/05.vennHTEJ_HIND_DLTEJ_LIND_D_HNsoil.txt
    tail -n 141 fig1/ST/05.vennHTEJ_HIND_DLTEJ_LIND_D_HNsoil.txt | cut -f 1 | sort | uniq -c # TEJ共有的51个IND富集，79个土壤富集，11个不显著

    # 查看IND特征与NRT调控共有的菌 fig4, ts8
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print a[$1],$0}' result/compare/HIND-HSoil1_all.txt fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D.txt > fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D_HNsoil.txt
    tail -n 6 fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DA50LnCp7_A56LnCp7_D_HNsoil.txt | cut -f 1 |sort|uniq -c # IND共有的6个5个富含，1个不显著
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$6} NR>FNR{print a[$1],$0}' result/compare/HIND-HSoil1_all.txt fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D.txt > fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D_HNsoil.txt
    tail -n 70 fig1/ST/07.vennHTEJ_HIND_DLTEJ_LIND_DV3703HnCp6_ZH11HnCp6_D_HNsoil.txt | cut -f 1 | sort | uniq -c # TEJ共有的37个IND富集，27个土壤富集，6个不显著

# 2019/1/30 水稻序列拼接

cd ~/rice/xianGeng/wet
cp -r /mnt/bai/haoran/20190111_16sBlast/seq ./

# 安装cap3 在线版http://doua.prabi.fr/software/cap3 本地版 http://seq.cs.iastate.edu/cap3.html 
conda install cap3

# 制作输入fa文件，将同一序列来源的多条测序结果合并，可以手动制作，也可以使用perl脚本自动批量
file=RiceP14C02
format_seq2fasta.pl -i "seq/${file}_*.seq" -o ${file}.fa
cap3 ${file}.fa # 拼接结果在屏幕上，且保存了一系列文件，*.fa.cap.contigs即可
sed -i "1 s/Contig1/${file}/" ${file}.fa.cap.contigs

# 批处理拼接
mkdir -p contigs
for file in `cut -f 2 16s_full_length_list.txt`; do
    echo $file
    format_seq2fasta.pl -i "seq/${file}_*.seq" -o ${file}.fa
    cap3 ${file}.fa > /temp/temp.txt
    # 改名
    sed -i "1 s/Contig1/${file}/" ${file}.fa.cap.contigs
    # 移至目录
    mv ${file}.fa.cap.contigs contigs/
    # 删除其它临时文件
    rm ${file}.*
done

# 序列合并
cat contigs/* > 16s_full_length_list.fa
format_fasta_1line.pl -i 16s_full_length_list.fa -o 16s_full_length_list1.fa # 还生成16s_full_length_list1.fa.tsv
# 追加菌ID
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$1} NR>FNR{print a[$1],$0}' 16s_full_length_list.txt 16s_full_length_list1.fa.tsv > 16s_full_length_list2.fa.tsv
# 追加到菌保表TableS11
sed -i 's/_//g' 16s_full_length_list2.fa.tsv
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3} NR>FNR{print $0,a[$3]}' 16s_full_length_list2.fa.tsv tableS11_1492.txt > tableS11_full.txt
# 有些ID没有对应上，统计数量
cut -f 11 tableS11_full.txt|sort|uniq -c|less # 99个空值

# 2019/2/18 根据编辑部分意见修改 统计、图点、具体P值和数值

# Talbe S10. 补充 13512 信息 CFU，原始数据：culture/rice/result/culture_bacteria.xls
# 添加每个CFU的物种注释：门、目、属
cd ~/rice/xianGeng/fig1/190218/table_raw
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3"\t"$5"\t"$7} NR>FNR{print $0,a[$4]}' tableS10_cotu.txt tableS10_cfu.txt > tableS10_cfu_anno.txt

# 与Indica enriched菌多少和N相关，3compare_en.rmd ### Stat 141 Indica commmon
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$5} NR>FNR{print $1,a[$1]}' result/compare/database.txt fig1/ST/S07.otu_TEJ_common_faprotax.txt > fig1/ST/S07.otu_TEJ_common_faprotax_mean.txt



# 2019/2/20 检验和整理文件和最终代码，提交github

cd ~/github
git clone git@github.com:YongxinLiu/Zhang2019NBT.git
cd ~/github/Zhang2019NBT
mkdir -p data fig1 fig2 fig3 fig4 fig5 fig6 script

# 通用数据
# cp ~/rice/xianGeng/ASV2/fig/design.txt data/
cp ~/rice/xianGeng/fig1/ST/02.design.txt data/design.txt
cp ~/rice/xianGeng/fig1/ST/02.otu.fa data/otu.fa
cp ~/rice/xianGeng/fig1/ST/02.otutab.txt data/otutab.txt
cp ~/rice/xianGeng/result/taxonomy_8.txt data/
cp ~/rice/xianGeng/result/faprotax/element_tab.txt data/
cp ~/rice/xianGeng/result/tax/sum_pc.txt data/
cp ~/rice/xianGeng/fig1/data/diff.list data/
cp ~/rice/xianGeng/fig1/data/HIND-HTEJ_all.txt data/
cp ~/rice/xianGeng/fig1/data/LIND-LTEJ_all.txt data/
cat data/?IND-?TEJ_all.txt | grep 'Enriched' | cut -f 1 | sort | uniq -d > data/otu_IND_common_specific.txt
cat data/?IND-?TEJ_all.txt | grep 'Depleted' | cut -f 1 | sort | uniq -d > data/otu_TEJ_common_specific.txt
cp ~/rice/xianGeng/result/faprotax/report.mat data/faprotax_report.mat
# 不考虑地块，亚种与土比较
cp ~/rice/xianGeng/fig1/data/*Soil* data/
# 基于OTU表群体丰度均值和物种注释
cp ~/rice/xianGeng/result//41culture/otu_mean_tax.txt data/

# 脚本
cp ~/rice/xianGeng/fig1/script/stat_plot_functions.R script/
cp ~/rice/xianGeng/fig1/script/compare.R script/
cp ~/rice/xianGeng/fig1/script/alpha_boxplot_far.R script/

# 图1.
cp data/TableS1varieties_geo.txt fig1/varieties_geo.txt
cp ~/rice/xianGeng/ASV2/fig/data/bray_curtis.txt fig1/
cp ~/rice/xianGeng/ASV2/fig/data/alpha.txt fig1/

# 图2.
cp /mnt/bai/yongxin/rice/xianGeng/result/tax/sum_f.txt fig2/
cp /mnt/bai/yongxin/rice/xianGeng/xiangeng0wilcox/result/compare_f/HTEJ-HIND_all.txt fig2/f_HTEJ-HIND_all.txt 

# 图3.
# 筛选用于分析和可视化的通路
cp ~/rice/xianGeng/doc/faprotax.id fig3/
# IND/TEJ在LN和HN中共同变化的通路
cp ~/rice/xianGeng/result/compare_far/???.list fig3/

# 图4
cp ~/rice/xianGeng/result/nitrogen_cor/good_snp.geno fig4/
cp ~/rice/xianGeng/fig1/4nrt/venn_nrt_indiaHL.txt fig4/
cp ~/rice/xianGeng/fig1/4nrt/venn_NRTsnp_indiaHL.txt fig4/
cp ~/rice/xianGeng/meta/nitrogen/ko.txt fig4/
cp ~/rice/xianGeng/meta/nitrogen/design.txt ~/github/Zhang2019NBT/fig4/meta_design.txt

# 图5.


# 图6.
cp ~/rice/xianGeng/wet/Xu-NbSycomFunctionData.txt fig6/Sycom.txt

# 附表，整理最终排版后附表更新到github
mkdir table && cd table 

# 二次检查文中来自附表的数值
## 表2国家数量，只有坐标，没有国家，可以查询。
cp ~/rice/xianGeng/fig1/ST/01.variety_geotable.txt ST02.VarietiesInfoAll.txt
cut -f 14 ST02.VarietiesInfoAll.txt|tail -n+2|sort|uniq|wc -l

# 表11目数量
cut -f 7 ST11.txt|tail -n+3|sort|uniq|wc -l # 22个目
# 表11科数量
cut -f 8 ST11.txt|tail -n+3|sort|uniq|wc -l # 57个科
# 筛选四大菌门
grep -P 'Actinobacteria|Bacteroidetes|Firmicutes|Proteobacteria' ST11.txt|cut -f 8|tail -n+3|sort|uniq|wc -l # 56个科

# 统计图中科
cd ~/github/Zhang2019NBT/fig5
# 原始为55个和30个科
wc -l family_fig5*
# sort uniq确保没有重复，复制出错
sort family_fig5a.txt|uniq|wc -l
sort family_fig5b.txt|uniq|wc -l
# 去除Unname
grep -v 'Unname' family_fig5a.txt|wc -l # 5a有34个
grep -v 'Unname' family_fig5b.txt|wc -l # 5a有29个
# 去除Unname并集，共45个
cat family_fig5?.txt|grep -v 'Unname'|sort|uniq # 人工查看是否有相近拼写错误
cat family_fig5?.txt|grep -v 'Unname'|sort|uniq|wc -l
# 去除Unname交集，共18个
cat family_fig5?.txt|grep -v 'Unname'|sort|uniq -d|wc -l
# 以上统计的是图中千分之一以上的的科，排除末命名的数量，而不是可培养的数量

# 依照图添加门和可培养Yes标记于5a/b文件
# 筛选四大均门，非Unname，且可培养
cat family_fig5?.txt|grep -P 'Actinobacteria|Bacteroidetes|Firmicutes|Proteobacteria'|grep -v 'Unname'|grep 'Yes'|cut -f 1|sort|uniq |  wc -l  # 27个科
cat family_fig5?.txt|grep -v 'Unname'|grep 'Yes'|cut -f 1|sort|uniq |  wc -l  # 27个科


# 统计四大菌门中可培养的科数量
cd ~/culture/jingmei171208.rice.med.soy
cd culture_A50 # 图5b
# 筛选可培养的，再看物种注释属于四大菌门和科
cat culture.sum
grep 'ring_alpha' 5_annotation.txt|wc -l # 筛选所有OTU 127
grep 'ring_shape' 5_annotation.txt|wc -l # 筛选可培养OTU 91
wc -l  otu_cultured.tax # 91个可培养的物种注释
cp ../0_ha_otu_culture.txt ./ # 原始物种注释在这里
cut -f 7 0_ha_otu_culture.txt|tail -n+2|sort|uniq -c|wc -l # 30科，与fig5/family_fig5b.txt一致
# 筛选四大菌门127个，可培养91个，科中非命名的90个，再非冗余21个
grep -P 'Actinobacteria|Bacteroidetes|Firmicutes|Proteobacteria'  0_ha_otu_culture.txt | grep 'YES' | cut -f 7| grep -v 'unnamed' | sort|uniq| wc -l




# push
git add .
git commit -m "final submit"
git push origin master



# 2019/4/5
# 错误1. OTU序列数量与OTU表不一样，otu.fa序列多于OTU表的行数
# 解读方法：依照OTU表筛选
grep -c '>' result/otu.fa
wc -l <(tail -n+2 result/otutab.txt)
# 备份原始数据
cp result/otu.fa result/otu.fa.bak
# 提取OTU表ID
cut -f 1 result/otutab.txt | tail -n+2 > temp/otuid.2
# 输出丢失的ID是什么
grep '>' result/otu.fa|sed 's/>//' > temp/otuid.1
cat temp/otuid.1 temp/otuid.2 | sort | uniq -u
# 可直接删除这两个ID，下面是提供新的
usearch10 -fastx_getseqs result/otu.fa.bak -labels temp/otuid.2 -fastaout result/otu.fa
grep -c '>' result/otu.fa
# 备份到github和fig1
cp result/otu.fa fig1/data/
cp result/otu.fa ~/github/Zhang2019NBT/data/

# 错误2. 附表10单菌CPU表，统计忘去表头多一行culture/rice/result/culture_bacteria.xls
# culture/rice/result/culture_bacteria.xls 手工加一行，并添加来源
sort -k1,1 ~/culture/rice/result/culture_bacteria.xls|less #确认缺失
# 手动加一行 L1P01A3 ，添加11667复制上一行，具体为添加 L1P01A3	38	49	COTU_68	Actinobacteria	Actinomycetales	Nocardioides	IR24


# 检查21个氮相关菌占总量的丰度
# 21个菌在列表在 D:\home\github\Zhang2019NBT\fig3\faprotax5.txt
fig3.Rmd permutation_test_faprotax 处生成 faprotax5_nitr.txt

# 筛选HN/LN中这些OTU
cd ~/github/Zhang2019NBT/data
mkdir -p temp
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' data/LIND-LTEJ_all.txt fig3/faprotax5_nitr.txt | cut -f 1-15 > temp/faprotax5_nitr_LIND-LTEJ.txt
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' data/HIND-HTEJ_all.txt fig3/faprotax5_nitr.txt | cut -f 1-15 > temp/faprotax5_nitr_HIND-HTEJ.txt

# 筛选141个OTUs
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' data/LIND-LTEJ_all.txt data/otu_IND_common_specific.txt | cut -f 1-15 > temp/IND141_LIND-LTEJ.txt
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' data/HIND-HTEJ_all.txt data/otu_IND_common_specific.txt | cut -f 1-15 > temp/IND141_HIND-HTEJ.txt

# 统计IND中21个氮相关菌在141个菌中的相对丰度
# (10.21+10.966)/(42.779+39.099) = 0.25862869146779354649600625320599

# 2019/4/5 统计IND_TEJ中全部平均值中21个氮相关菌在141个菌中的相对丰度
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$9} NR>FNR{print $1,a[$1]}' data/IND_TEJ_all.txt fig3/faprotax5_nitr.txt | awk '{a=a+$2} END {print a}' # 9.257
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$9} NR>FNR{print $1,a[$1]}' data/IND_TEJ_all.txt data/otu_IND_common_specific.txt | awk '{a=a+$2} END {print a}' # 27.678
# 9.257/ 27.678 = 25.67%
# 计算japonica的氮相关的OTU丰度 # OTU_86在16个中比例为6.25%，相对丰度为0.065   0.276，0.126   0.182
echo 'OTU_86' > fig3/faprotax5_nitr_TEJ.txt
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$9} NR>FNR{print $1,a[$1]}' data/IND_TEJ_all.txt fig3/faprotax5_nitr_TEJ.txt | awk '{a=a+$2} END {print a}' # 0.133494
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$9} NR>FNR{print $1,a[$1]}' data/IND_TEJ_all.txt data/otu_TEJ_common_specific.txt | awk '{a=a+$2} END {print a}' # 17.0228
# 0.133494 / 17.0228 = 0.78%

# 基于figureS8d挑选结果所有样品均值
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' table/ST08d.txt fig3/faprotax5_nitr.txt | awk '{a=a+$2} END {print a}' # 9.21975
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' table/ST08d.txt data/otu_IND_common_specific.txt | awk '{a=a+$2} END {print a}' # 35.9249
# 9.21975/35.9249=33.445%

# 添加table7b/d丰度
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$9} NR>FNR{print $1,a[$1]}' data/IND_TEJ_all.txt fig3/faprotax_IND.txt > fig3/faprotax_IND_abundance.txt
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$9} NR>FNR{print $1,a[$1]}' data/IND_TEJ_all.txt fig3/faprotax_IND.txt| awk '{a=a+$2} END {print a}'

awk '{FS=OFS="\t"} NR==FNR{a[$1]=$9} NR>FNR{print $1,a[$1]}' data/IND_TEJ_all.txt fig3/faprotax_TEJ.txt > fig3/faprotax_TEJ_abundance.txt
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$9} NR>FNR{print $1,a[$1]}' data/IND_TEJ_all.txt fig3/faprotax_TEJ.txt | awk '{a=a+$2} END {print a}'


# 2019/4/7 统计氮相关的ID的物种注释各级别数量
cd ~/github/Zhang2019NBT
awk '{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' table/ST03.c.taxonomy.txt fig3/faprotax_IND_enriched_nitr.txt > fig3/faprotax_IND_enriched_nitr_tax.txt
for i in `seq 1 6`;do cut -f $i fig3/faprotax_IND_enriched_nitr_tax.txt | sort | uniq | wc -l ; done

awk '{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' table/ST03.c.taxonomy.txt fig3/faprotax_TEJ_enriched_nitr.txt > fig3/faprotax_TEJ_enriched_nitr_tax.txt
for i in `seq 1 6`;do cut -f $i fig3/faprotax_TEJ_enriched_nitr_tax.txt | sort | uniq | wc -l ; done


# 2019/4/7 实验筛选菌与OTU相似度验证
cd ~/rice/xianGeng/wet
#菌保编号来自 表12b，调整格式为 选菌实验验证-20实验菌 /mnt/bai/yongxin/rice/xianGeng/wet/verify1_20bac.ID
#对应序列来自 表11，调整格式为 选菌实验验证-1079seq /mnt/bai/yongxin/rice/xianGeng/wet/verify2_1079seq.fa
grep 'Indica-enriched' verify1_20bac.ID | cut -f 1 > verify1_20bac.idonly
usearch10 -fastx_getseqs verify2_1079seq.fa -labels verify1_20bac.idonly -fastaout verify1_20bac.fa

#OTU ID来自附表8c中最后三部分，来自46+70+3的，即indica富集(任意一块地)且受nrt1.1b调控，/mnt/bai/yongxin/rice/xianGeng/wet/verify3_119.id
#OTU序列来自 附表3b 5141条OTU序列 /mnt/bai/yongxin/rice/xianGeng/wet/verify3_otu_5141.fa
usearch10 -fastx_getseqs verify3_otu_5141.fa -labels verify3_119.id -fastaout verify3_119.fa # 119 found

# 比较相似度：培养菌比对OTU
makeblastdb -in verify3_119.fa -dbtype nucl
blastn -query verify1_20bac.fa -db verify3_119.fa -out verify1_20bac_119.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 9 
awk '$3>97' verify1_20bac_119.blastn | wc -l # 只有R1889一个96.29%，不到97%，可能是因为经COTU中转导致
# 截取V5-V7也是R1889为96%
blastn -query verify1_20bac_v57.fa -db verify3_119.fa -out verify1_20bac_119.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 9 
awk '$3>97' verify1_20bac_119.blastn | wc -l # 只有R1889一个96.29%，不到97%，可能是因为经COTU中转导致
# R1889对应 “# 1vs3 实验菌vsOTU”中~/rice/xianGeng/select_bac/syncom_otu.blastn为I12吗？有相同的相似度，I12     OTU_166 96.29，最好为OTU_371，但其不存在于差异OTU列表
# 使用旧版序列，也是有一个小于100%
blastn -query ../select_bac/syncom_1492r.fa -db verify3_119.fa -out verify1_20bac_119.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 9 
awk '$3>97' verify1_20bac_119.blastn | wc -l # 只有R1889一个96.29%，不到97%，可能是因为经COTU中转导致

