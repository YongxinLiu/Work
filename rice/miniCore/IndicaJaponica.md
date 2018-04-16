<!-- TOC -->

- [1. 基于OTU分析表型](#1-基于otu分析表型)
    - [HN](#hn)
    - [LN](#ln)
    - [(HN-LN)/HN](#hn-lnhn)
- [2. 籼粳差异微生物组与氮吸收NRT](#2-籼粳差异微生物组与氮吸收nrt)
    - [2.1. GWAS关联分析](#21-gwas关联分析)
    - [2.2. OTUs](#22-otus)
    - [2.3. 元素循环Faprotax](#23-元素循环faprotax)

<!-- /TOC -->


# 1. 基于OTU分析表型
OTU表分为三类：高氮 HN，低氮LN，高氮 与低氮之差比例(HN-LN)/HN

## HN

## LN

基于HN的方法制作LN相关表型

    mkdir -p LN

#### 筛选低氮的样品

    # a. 获得所有LN的样品名和组名，筛选LN下的样品
    grep 'L$' doc/design.txt | cut -f 1,5> LN/design.txt
    Rscript scripts/otutab_sample_subset.r -i result/otu_table_filter.txt -d LN/design.txt -o LN/otutab0.txt
    #cut -f 1 LN/design.txt | grep -v -P 'M4111Le|S4161Hb|S4161Hc|L4110Hb' > LN/design.sampleID
    #usearch10 -otutab_sample_subset result/otu_table_filter.txt -labels LN/design.sampleID -output LN/otutab0.txt
    
    # b. 计算bray_curtis距离，筛选每组内的距离，最小Top3输出，并筛选
    usearch10 -beta_div LN/otutab0.txt -tree result/rep_seqs.fa.tree -filename_prefix LN/ -metrics bray_curtis
    Rscript scripts/beta_pcoa_group.r -i LN/bray_curtis.txt -d doc/design.txt -n GroupID -o LN/pcoa_bray
    #tail -n+2 LN/pcoa_bray_samples_top3.txt | cut -f 1 > LN/pcoa_bray_samples_top3.id
    #usearch10 -otutab_sample_subset LN/otutab0.txt -labels LN/pcoa_bray_samples_top3.id -output LN/otutab1.txt
    Rscript scripts/otutab_sample_subset.r -i LN/otutab0.txt -d LN/pcoa_bray_samples_top3.txt -o LN/otutab1.txt

# c. 按组合并
cut -c1-5 LN/pcoa_bray_samples_top3.id > LN/pcoa_bray_samples_top3.grp
paste LN/pcoa_bray_samples_top3.id LN/pcoa_bray_samples_top3.grp > LN/pcoa_bray_samples_top3.group
usearch10 -otutab_group LN/otutab1.txt -labels LN/pcoa_bray_samples_top3.group -output LN/otutab2.txt
# d. 筛选LN有表型的196个
head -n1 otu_filter/7otu_table_ratio_LN_top100.txt|cut -f 2-|sed 's/\t/\n/g' > result/variety196.id
usearch10 -otutab_sample_subset LN/otutab2.txt -labels result/variety196.id -output LN/otutab3.txt
usearch10 -otutab_stats LN/otutab3.txt -output LN/otutab3.stat
cat LN/otutab3.stat # 从77M至19.8M

### 筛选LN OTU表

    cp HN/otutab3.txt HN/otutab.txt
    # 标准化为比例
    usearch10 -otutab_counts2freqs HN/otutab.txt -output HN/otutab_freq.txt # 标准化为1

### a. alpha多样性分析
usearch10 -otutab_norm HN/otutab.txt -sample_size 10000 -output HN/otutab_norm.txt 
usearch10 -alpha_div HN/otutab_norm.txt -output HN/alpha.txt 
head -n1 HN/alpha.txt|tr '\t' '\n'|awk '{print NR,$0}' #共14种多样性指数，其中chao1, richness和shannon_e分别在4，10和13列
Rscript scripts/gemma_fam.r -i HN/alpha.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam -t FALSE

### b. beta多样性
mkdir HN/beta
usearch10 -cluster_agg result/rep_seqs.fa -treeout result/rep_seqs.fa.tree
usearch10 -beta_div HN/otutab.txt -tree result/rep_seqs.fa.tree -filename_prefix HN/beta/
# 制作样品对品种的表：组太多，没法观察颜色，只保留Admix
cut -f 1 -d '(' doc/design_group.txt > doc/design_group7.txt
# 结果HN/pcoa_bray15.txt即为1-5轴
Rscript scripts/beta_pcoa_dist.r -i HN/beta/bray_curtis.txt -d doc/design_group7.txt -n Region -o HN/pcoa_bray_curtis_ -w 8 -e 5
# 对三种距离和PCoA进行计算，在 unifrac_binary PCoA 3/4轴上，R4159异常，检查alpha中robbins为零，richness倒数第4(sort -k10,10n HN/alpha.txt|less)
for i in bray_curtis unifrac unifrac_binary; do
	Rscript scripts/beta_pcoa_dist.r -i HN/beta/${i}.txt -d doc/design_group7.txt -t ${i} -n Region -o HN/pcoa_${i}_ -w 8 -e 5
done
### Constrained PCoA by group，PC1-5 HN/pcoa_bray_constrained_15.txt
Rscript scripts/beta_cpcoa.r -i HN/otutab.txt -d doc/design_group7.txt -n Region -o HN/pcoa_bray_constrained -w 8 -e 5
for i in bray_curtis unifrac unifrac_binary bray_constrained; do
	Rscript scripts/gemma_fam.r -i HN/pcoa_${i}_15.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam -t FALSE
done


# c. 制作各分类级结果 /mnt/bai/yongxin/rice/miniCore/180319/HN/tax*
for i in p c o f g;do
	usearch10 -sintax_summary temp/otus_no_host.tax -otutabin HN/otutab.txt -rank ${i} -output HN/tax_${i}.txt
	Rscript scripts/gemma_fam.r -i HN/tax_${i}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam
done

# d. 制作为功能分类级，依赖picrust结果
for l in 1 2 3 4;do 
	# 筛选样品
	usearch10 -otutab_sample_subset function/metagenome_predictions.L${l}.txt -labels HN/pcoa_bray_samples_top3.id -output temp/temp_${l}_1
	# 按组合并
	usearch10 -otutab_group temp/temp_${l}_1 -labels HN/pcoa_bray_samples_top3.group -output temp/temp_${l}_2
	# 筛选HN有表型的196个
	usearch10 -otutab_sample_subset temp/temp_${l}_2 -labels result/variety196.id -output HN/KO${l}.txt
	# 统计结果
	usearch10 -otutab_stats HN/KO${l}.txt -output HN/KO${l}.stat
	cat HN/KO${l}.stat
	Rscript scripts/gemma_fam.r -i HN/KO${l}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam
done

# e. 元素循环farprotax
make faprotax_calc # 默认使用gg97注释
# Assigned 2735 records to groups, 9831 records were leftovers
# 修改makefile中gg13.5为rdp_v16s，可注释的物种提高了10%
# Assigned 3060 records to groups, 9506 records were leftovers
# 筛选HN下样品
# 筛选样品
sed -i '1 s/group/#OTU ID/' faprotax/otu_table_tax.faprotax
usearch10 -otutab_sample_subset faprotax/otu_table_tax.faprotax -labels HN/pcoa_bray_samples_top3.id -output temp/temp_far_1
# 按组合并
usearch10 -otutab_group temp/temp_far_1 -labels HN/pcoa_bray_samples_top3.group -output temp/temp_far_2
# 筛选HN有表型的196个
usearch10 -otutab_sample_subset temp/temp_far_2 -labels result/variety196.id -output HN/faprotax.txt
usearch10 -otutab_counts2freqs HN/faprotax.txt -output HN/faprotax_freq.txt # 标准化为1
# 统计结果
usearch10 -otutab_stats HN/faprotax.txt -output HN/faprotax.stat
cat HN/faprotax.stat

# f. 转换微生物组数据为Gemma表型
# OTU表转置
#awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++){a[FNR,i]=$i}} END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' HN/otutab_norm.txt | sed 's/\t$//' > HN/otutab_norm.txtt
#awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++){a[FNR,i]=$i}} END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' HN/tax_p.txt | sed 's/\t$//' > HN/tax_p.txtt
## 基因型数据位于/mnt/zhou/chulab/miniCore/snp1.5x/T2，其中bam文件追加表型，有不存在缺失应该补NA，目前需要用Excel中添加并托整行，再补齐所有缺失
#awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' HN/tax_p.txtt /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam | cut -f 1-6,9- > HN/tax_p.gemma
# 用R脚本实现合并，并补NA，输出为输入*.fam用于gemma分析，表头位于*.fam.header中
Rscript scripts/gemma_fam.r -i HN/tax_p.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam	
# 循环处理所有：alpha己转置需单独处理
for i in otutab_freq faprotax; do
	Rscript scripts/gemma_fam.r -i HN/${i}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam
done


## (HN-LN)/HN





# 2. 籼粳差异微生物组与氮吸收NRT

    cd rice/miniCore/180319/
    mkdir IndTej

## 2.1. GWAS关联分析

## 2.2. OTUs

## 2.3. 元素循环Faprotax

    # 设置变量
    mv IndTej/HN/faprotax_freq.txt IndTej/HN/fapro.txt
    i=fapro
    
    # 转换表型为fam文件
	Rscript scripts/gemma_fam.r -i IndTej/HN/${i}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam
	cp IndTej/HN/${i}.txt.fam T2.fam
    
	# 生成表型列名对应表
    # gemma表型文件fam不允许表头，单独保存关联结果后改名
    # 去掉前5行，取前99行，删除()"等符号，添加行号为gemma找列和输出文件名
	tail -n+6 IndTej/HN/${i}.txt.fam.header|head -n99|sed 's/[\(\)\"]//g'|awk '{print NR"\t"NR$0}'| sed "s/\t/\t${i}/"|sed '1 s/V6/Wax/' > IndTej/HN/${i}.list
    
    # 并行gemma关联基因型和表型
    # 使用paralle进行48线程并行
	parallel -j 48 "gemma -bfile T2 -k gemma/kinship.txt -c gemma/pca4.txt -lmm 4 -n {1} -o {1}" ::: `cut -f 1 IndTej/HN/${i}.list`
    
	# 批量按列编号改名
	awk 'BEGIN{OFS=FS="\t"}{system("mv output/"$1".assoc.txt output/"$2".assoc.txt");system("mv output/"$1".log.txt output/"$2".log.txt")}' <(cat IndTej/HN/${i}.list)
	
    # 筛选极显著的位点P<0.001
    # 注意并行时变量$符需要转义
	parallel -j 33 "cat <(head -n 1 output/{1}.assoc.txt) <(awk '\$14<0.001' output/{1}.assoc.txt) > output/{1}.assoc.sig" ::: `cut -f 2 IndTej/HN/${i}.list`
    
	# qqman批量绘制manhanttan plot
	parallel -j 33 "Rscript gemma/qqman.R -i output/{1}.assoc.sig -t 14" ::: `cut -f 2 IndTej/HN/${i}.list`
    
	# 添加SNP注释
    # parallel有太多$需要转义不方便，改用循环
	for RPM in `cut -f 2 IndTej/HN/${i}.list`; do
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6} NR>FNR{print $0,a[$2]}' <(sed '1 s/ID/rs/' gemma/snp.anno) output/${RPM}.assoc.sig |sed 's/\t$/\tNA/'|cut -f 1-3,14- > output/${RPM}.assoc.anno &
	done

    # 检索NRT1.1b基因的P值，无关键SNP
    grep 'OS10G0554200' output/*.assoc.anno|awk '$4<0.001'|sort -k4,4g|less -S
    # 检索SNP无结果
    grep '10m21761740' output/fapro12aerobic_ammonia_oxidation.assoc.txt # 无结果，SNP丢了吗？
    wc -l T2.bim # 原SNP 2288867
    wc -l output/fapro12aerobic_ammonia_oxidation.assoc.txt # 分析结果 2153792，少了135075

    
