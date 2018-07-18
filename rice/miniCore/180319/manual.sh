#!/bin/bash

# 16S扩增子操作手册 16S Amplicon pipeline manual, v1.4 2017/12/23

# 按如下命令执行，并修改makefile中的参数
# 进入工作目录
cd ~/rice/miniCore/180319



# 0. 配置实验设计参数

# 复制最新版流程配置文件+操作手册至工作目录
cp ~/rice/timecourse/ma* ./

# 修改工作目录wd
pwd

# 初始化工作目录
make init


# 从这里开始
ln ../temp/seqs_unique.fa temp/
ln ../temp/seqs_usearch.fa temp/
cp ../doc/design.txt doc/
cp ../doc/group_compare.txt doc/
touch merge_library
touch derep

make cluster_otu
make otu_table
make assign_tax # 基于gg13.8注释物种
make tree # 构建进化树，用于多样性分析
make alpha # 计算alpha多样性
make alpha_usearch_rare # (可选)，采用usearch结果绘制
make alpha_qiime_rare # (可选)，采用qiime进行稀释曲线分析，大量抽样比较费时；可继续运行其它步骤
make beta # 计算beta多样性
make graphlan # 绘制高丰度物种树
make ggtree # 绘制高丰度进化树
make diversity # 多样性绘制箱线和散点

make filter # 过滤高丰度菌分析，默认0.05%。如lefse只有1%，ruben文章常用0.5-0.1%
make rediv # 新表重新多样性分析
make draw_div # 多样性绘图

# 功能预测：采和Picrust将OTU转换为KO，可以在level1-4及PCoA轴与基因型与表型关联
make picrust_calc # 修改usearch10为vsearch提高比较速度



## 基因型和SNP相关分析 2018/7/3

	# 来自 ~/rice/xianGeng/manual.sh # 图4. 菌群与nrt关系

	## 2. miniCore菌与相关基因相关分析
	cd ~/rice/miniCore/180319
	# 小宁整理了16个候选基因
	less /mnt/bai/xiaoning/past/software/GEMMA/gwas_candidate_gene.txt
	cp /mnt/bai/xiaoning/past/software/GEMMA/gwas_candidate_gene.txt doc/gwas_candidate_gene.txt
	# 在 http://rapdb.dna.affrc.go.jp 上查找每个基因的symbol, name, description and links
	# 来自胡斌整理的氮相关基因doc/N-related genes in rice.docx共9个基因，先在doc/rice_nitrogen_list.xlsx中惠惠相关ID，保存为doc/rice_nitrogen_list.txt, 其中第5列RAP_SNP的ID与SNP注释文件对应
	dos2unix doc/gwas_candidate_gene.txt # windows转换为linux
	# 整理出SNP数据中对应的基因型、提取相应的位点
	mkdir -p result/gene_cor
	cut -f 3 doc/gwas_candidate_gene.txt | tail -n+2 | tr '\n' '|' # 提取ID并替换为|分隔
	grep -P 'OS10G0554200|OS08G0155400|OS03G0203200|OS06G0110000|OS11G0104300|OS01G0746400|OS04G0550600|OS11G0587000|OS09G0441900|OS05G0407500|OS05G0333200|OS01G0927000|OS05G0196500|OS04G0509600|OS01G0909200|OS02G0699000|OS01G0883800' /mnt/bai/yongxin/rice/miniCore/180319/gemma/snp.anno > result/gene_cor/all_snp.list # 筛选到10个基因在miniCore中存在689个相关位点
	cut -f 3 result/gene_cor/all_snp.list | sort | uniq -c # 1 HIGH, 38个MODERATE，615个MODIFIER和35个LOW
	grep -P 'HIGH|MODERATE' result/gene_cor/all_snp.list > result/gene_cor/good_snp.list # 筛选高和中等
	cut -f 4 result/gene_cor/good_snp.list|sort|uniq -c # 其中重要SNP仅有39个，来自9个基因
	# 提SNP对应基因型
	cat /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_*.hmp.txt > /tmp/temp
	cut -f 1 result/gene_cor/good_snp.list|tr '\n' '|' # 获取列表
	grep -P '1m39596400|1m39596699|1m39596713|1m39596808|1m39600096|1m39602780|1m39605678|1m40659706|1m40659949|1m40663740|2m28751006|2m28751738|2m28752575|4m27567935|4m27568586|4m27568855|5m5941257|5m5941396|5m5943888|5m5944557|5m5944851|5m5944944|8m3183208|9m16411594|9m16414341|9m16414429|9m16415032|9m16415203|9m16415254|10m21759092|10m21761740|10m21761997|11m22230369|11m22230384|11m22230417|11m22230594|11m22230660|11m22230777|11m22230789' /tmp/temp > /tmp/temp1
	cat <(head -n1 /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_1.hmp.txt) /tmp/temp1 > result/gene_cor/good_snp.geno # 添加标题
	# 转置
	awk 'BEGIN{FS=OFS="\t"}{for(i=1;i<=NF;i++){a[FNR,i]=$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' <(cut -f 1,12- result/gene_cor/good_snp.geno) | sed 's/\t$//' > result/gene_cor/good_snp.geno.t
	# 添加品种信息到SNP
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$0} NR>FNR{print $0,a[$1]}' ../doc/minicore_list.txt <(sed '1 s/rs/EFD/' result/gene_cor/good_snp.geno.t) > result/gene_cor/good_snp.geno.t.txt

	# 统计SNP基因型作为分组信息，来统计LN条件下alpha多样性相关
	cp ~/rice/xianGeng/script/nitrogen_cor.r script/snp_microbiome_wilcox.R

	# 添加P值和基因注释和描述
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4} NR>FNR{print $0,a[$1]}' result/gene_cor/all_snp.list <(sed 's/X//' result/gene_cor/alpha/pvalue.txt|cut -f 1,4,10,13) |less -S > result/gene_cor/alpha/pvalue.txt2
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$3]=$0} NR>FNR{print $0,a[$5]}' doc/gwas_candidate_gene.txt result/gene_cor/alpha/pvalue.txt2 |less -S > result/gene_cor/alpha/pvalue.txt3

	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4} NR>FNR{print $0,a[$1]}' result/gene_cor/all_snp.list <(sed 's/X//' result/gene_cor/beta/pvalue.txt|cut -f 1-4) |less -S > result/gene_cor/beta/pvalue.txt2
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$3]=$0} NR>FNR{print $0,a[$5]}' doc/gwas_candidate_gene.txt result/gene_cor/beta/pvalue.txt2 |less -S > result/gene_cor/beta/pvalue.txt3

	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4} NR>FNR{print $0,a[$1]}' result/gene_cor/all_snp.list <(sed 's/X//' result/gene_cor/genus/pvalue.txt) |less -S > result/gene_cor/genus/pvalue.txt2
	head -n 1 result/gene_cor/genus/pvalue.txt2 | sed 's/\t/\n/g'|awk '{print NR,$1}'
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$3]=$1} NR>FNR{print a[$513],$0}' doc/gwas_candidate_gene.txt result/gene_cor/genus/pvalue.txt2 |less -S > result/gene_cor/genus/pvalue.txt3

	 awk '{for(i=1;i<=NF;i++){a[FNR,i]=$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' result/gene_cor/genus/pvalue.txt3 > result/gene_cor/genus/pvalue.txt4 # 转置
	# 按丰度排序




	# 整理结果为D:\work\rice\miniCore\180319\result\gene_cor\

	gene_cor.r # 保存实验设计+基因型，方便识别SNP不同基因型在籼粳稻中区别
	# 整理这4个重要SNP信息表，见SNP_list.xlsx
	# 统计基因型与亚种分布
	sed -i '1 s/^/SampleID\t/' result/gene_cor/design.txt
	head -n1 result/gene_cor/design.txt|tr '\t' '\n'|awk '{print NR,$1}'
	# NRT2.1 - 17; 1.1A - 21; 1.1B - 14，统计每个亚种内基因型的数量，可看到亚种内主要的SNP类型
	cut -f 2,8,20 result/gene_cor/design.txt|sort|uniq|cut -f 2,3|sort|uniq -c
	grep -P '10m21759092|2m655515\t|8m3183208' /mnt/bai/yongxin/rice/miniCore/180319/gemma/T2.ann.vcf # 查询SNP变化位置、碱基和AA详细







## 1 数据筛选：

# OTU表先去除叶绿体和线粒体，HN和LN分开；标准化，求均值；再按门-属合并

make filter_otu # 排除宿主
make filter_otu1 # 按组指定1/2列筛选，测试用:基因型名+LN

# 按beta多样性筛选样品
otutab=otu_filter/2otu_table_filter.txt
# 评估原始OTU表
usearch10 -otutab_stats ${otutab} -output ${otutab}.sum 
cat ${otutab}.sum 
# 建树
usearch10 -cluster_agg result/rep_seqs.fa -treeout result/rep_seqs.fa.tree # 与result/rep_seqs.tree类不同，每个括号换行
#usearch10 -otutab_counts2freqs ${otutab} -output ${otutab}.norm # 标准化为1
#usearch10 -otutab_stats ${otutab}.norm -output ${otutab}.norm.sum
#cat ${otutab}.norm.sum
## 转换为总和数值10000的数值
#usearch10 -otutab_freqs2counts ${otutab}.norm -output ${otutab}.rare
#usearch10 -otutab_stats ${otutab}.rare -output ${otutab}.rare.sum
#cat ${otutab}.rare.sum
# 抽样至2W
usearch10 -otutab_norm ${otutab} -sample_size 20000 -output ${otutab}.norm
usearch10 -otutab_stats ${otutab}.norm -output ${otutab}.norm.sum 
cat ${otutab}.norm.sum 


# 筛选样品中最稳定的重复：Beta多样性分析
# GWAS分析每品种基因型数据只能对应单样品数据，
# 1. 使用otutab_counts2freqs标准化的结果，发现距离无结果
# 2. 采用原始数据计算，发现结果可能偏大
# 3. 采用rarefraction标准化的结果，与原始序列相同，可能它也进行了标准化
sub="LN/beta_norm"
mkdir -p ${sub}
# 计算距离矩阵，用count2freqs无结果，尝试用原始文件
usearch10 -beta_div ${otutab} -tree result/rep_seqs.fa.tree -filename_prefix ${sub}/ -metrics bray_curtis,unifrac
# 绘制图：基于otutab_freqs2counts的标准化结果，大量点坐标重合?原因为bray_curtis矩阵矩阵大部分全为0，而unifrac全为-1
beta_pcoa.r -i ${sub}/bray_curtis.txt -t bray_curtis -d doc/design.txt -n GroupID -w 16 -e 10 -o ${sub}/pcoa_bray # 整体看不清，需要每个小组查看
# 按小组展示并筛选数据，B4017L为例
beta_pcoa.r -i ${sub}/bray_curtis.txt -t bray_curtis -d doc/design_B4017L.txt -n GroupID -w 16 -e 10 -o ${sub}/pcoa_bray # 发现组内波动极大
# 按功能分类着色
beta_pcoa.r -i ${sub}/bray_curtis.txt -t bray_curtis -d ../doc/design1.txt -n Region -w 16 -e 10 -o ${sub}/pcoa_bray_region 

# 写脚本scripts/beta_pcoa_group.r，筛选每组的距离，选前3输出
Rscript scripts/beta_pcoa_group.r -i ${sub}/bray_curtis.txt -t bray_curtis -d doc/design.txt -n GroupID -o ${sub}/pcoa_bray 
wc -l ${sub}/pcoa_bray_samples_* # 1118个样，筛选为579个样

# 筛选这579个样品
# 备份子集OTU表
cp otu_filter/2otu_table_filter.txt otu_filter/2otu_table_filter_LN.txt
tail -n+2 LN/beta_norm/pcoa_bray_samples_top3.txt | cut -f 1 > LN/beta_norm/sample_ids.txt
usearch10 -otutab_sample_subset otu_filter/2otu_table_filter_LN.txt -labels LN/beta_norm/sample_ids.txt -output otu_filter/3otu_table_best_LN.txt
usearch10 -otutab_stats otu_filter/3otu_table_best_LN.txt -output otu_filter/3otu_table_best_LN.txt.sum 
cat otu_filter/3otu_table_best_LN.txt.sum

sub="LN/beta_filter"
mkdir -p ${sub}
usearch10 -beta_div otu_filter/3otu_table_best_LN.txt -tree result/rep_seqs.fa.tree -filename_prefix ${sub}/ -metrics bray_curtis
beta_pcoa.r -i ${sub}/bray_curtis.txt -t bray_curtis -d doc/design.txt -n GroupID -w 16 -e 10 -o ${sub}/pcoa_bray # 579个样与原来看不出差异
# 绘制按地区着色
beta_pcoa.r -i ${sub}/bray_curtis.txt -t bray_curtis -d ../doc/design1.txt -n Region -w 16 -e 10 -o ${sub}/pcoa_bray_region # 579个样与原来看不出差异



# OTU表按组合并，是数值直接相加，而不是标准化求平均，这样可能受样本权重影响 http://www.drive5.com/usearch/manual/cmd_otutab_group.html
awk 'BEGIN{OFS=FS="\t"} NR===FNR {a[$$1]=$$0} NR>FNR {print a[$$1]}
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$6} NR>FNR {print $0,a[$1]}' doc/design.txt LN/beta_norm/sample_ids.txt > LN/beta_norm/sample_ids_group.txt
usearch10 -otutab_group otu_filter/3otu_table_best_LN.txt  -labels LN/beta_norm/sample_ids_group.txt -output otu_filter/4otu_table_group_LN.txt
usearch10 -otutab_stats otu_filter/4otu_table_group_LN.txt -output otu_filter/4otu_table_group_LN.txt.sum 
cat otu_filter/4otu_table_group_LN.txt.sum 

sub="LN/beta_merge"
mkdir -p ${sub}
usearch10 -beta_div otu_filter/4otu_table_group_LN.txt -tree result/rep_seqs.fa.tree -filename_prefix ${sub}/ -metrics bray_curtis
# 制作按组的实验设计
cut -f 5- ../doc/design1.txt | tail -n+2 | sort| uniq |sed '1 s/^/SampleID\tvariety\tDescription\tRegion\n/' > doc/design_group.txt
cut -f 6,8 ../doc/design1.txt | tail -n+2 | sort| uniq |sed '1 s/^/SampleID\tRegion\n/' > doc/design_group.txt
# 绘制按地区着色
beta_pcoa.r -i ${sub}/bray_curtis.txt -t bray_curtis -d doc/design_group.txt -n Region -w 16 -e 10 -o ${sub}/pcoa_bray_region # 579个样与原来看不出差异



# 筛选高丰度:按丰度排序
usearch10 -otutab_sortotus otu_filter/4otu_table_group_LN.txt -output otu_filter/5otu_table_sort_LN.txt

# 标准化为百分比
usearch10 -otutab_counts2freqs otu_filter/5otu_table_sort_LN.txt -output otu_filter/6otu_table_ratio_LN.txt # 标准化为1

# 筛选前Top 1000个OTUs
head -n 1000 otu_filter/6otu_table_ratio_LN.txt > otu_filter/7otu_table_ratio_LN_top100.txt



# 2018/4/10 筛选OTU表，按OTU行和样品列分别筛选。再按注释合并

# 1.1 OTU表按OTU注释筛选
# 1.1.1去嵌合体
cp result/rep_seqs.fa temp/rep_seqs.fa
grep -c '>' temp/rep_seqs.fa # 12566
usearch10 -uchime2_ref temp/rep_seqs.fa \
	-db /mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.fa -strand plus -mode balanced \
	-chimeras temp/otus_chimeras.fa -threads 24
# 获得非嵌合体序列ID
cat temp/rep_seqs.fa temp/otus_chimeras.fa|grep '>' | sort | uniq -u | sed 's/>//' > temp/no_chimeras.id
# 筛选非嵌合体
usearch10 -fastx_getseqs temp/rep_seqs.fa -labels temp/no_chimeras.id -fastaout temp/otus_no_chimeras.fa
grep -c '>' temp/otus_no_chimeras.fa # 9950
# 1.1.2 去非细菌序列(线粒体、叶绿体、细菌)
usearch_silva=/mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.udb
usearch10 -sintax temp/otus_no_chimeras.fa \
	-db ${usearch_silva} -sintax_cutoff 0.8 -strand both \
	-tabbedout temp/otus_no_chimeras.tax -threads 24
grep -P -v 'Mitochondria|Chloroplast|Eukaryota|\t$' temp/otus_no_chimeras.tax | cut -f 1 > temp/otus_no_host.id
wc -l temp/otus_no_host.id # 9625
# 获得最终的廒和OTU表
usearch10 -fastx_getseqs temp/otus_no_chimeras.fa -labels temp/otus_no_host.id -fastaout temp/otus_no_host.fa
# 再用rdp分类注释，查看是否有末分类结果
usearch_rdp=/mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb
usearch10 -sintax temp/otus_no_host.fa \
	-db ${usearch_rdp} -sintax_cutoff 0.8 -strand both \
	-tabbedout temp/otus_no_host.tax -threads 24
grep -P 'Mitochondria|Chloroplast|Eukaryota|\t$' temp/otus_no_host.tax # 无Mitochondria，有大量Chloroplast来自p:Cyanobacteria/Chloroplast为正常，无末注释
# Usearch筛选OTU失败，改为手动筛选
# usearch10 -otutab_otu_subset result/otu_table.txt -labels temp/otus_no_host.id -output result/otu_table_filter.txt
# ---Fatal error--- ../otutab.cpp(30) assert failed: SampleIndex < m_SampleCount
head -n1 result/otu_table.txt > result/otu_table_filter.txt
awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/otu_table.txt temp/otus_no_host.id >> result/otu_table_filter.txt
usearch10 -otutab_stats result/otu_table.txt -output result/otu_table.stat
cat result/otu_table.stat
usearch10 -otutab_stats result/otu_table_filter.txt -output result/otu_table_filter.stat
cat result/otu_table_filter.stat # 从117M至77M，最大值从256K至130K，主要是OTU_1线粒体导致



# 1.2 按样品筛选

# 先按高低氮分开，以HN为例

# 1.2.1 筛选高氮的样品
# a. 获得所有HN的样品名和组名，筛选HN下的样品
grep 'H$' doc/design.txt | cut -f 1,5> HN/design.txt
cut -f 1 HN/design.txt | grep -v -P 'F4056Hf|S4161Hb|S4161Hc|L4110Hb' > HN/design.sampleID
usearch10 -otutab_sample_subset result/otu_table_filter.txt -labels HN/design.sampleID -output HN/otutab0.txt
# b. 计算bray_curtis距离，筛选每组内的距离，最小Top3输出，并筛选
usearch10 -beta_div HN/otutab0.txt -tree result/rep_seqs.fa.tree -filename_prefix HN/ -metrics bray_curtis
Rscript scripts/beta_pcoa_group.r -i HN/bray_curtis.txt -d doc/design.txt -n GroupID -o HN/pcoa_bray
tail -n+2 HN/pcoa_bray_samples_top3.txt | cut -f 1 > HN/pcoa_bray_samples_top3.id
usearch10 -otutab_sample_subset HN/otutab0.txt -labels HN/pcoa_bray_samples_top3.id -output HN/otutab1.txt
# c. 按组合并
cut -c1-5 HN/pcoa_bray_samples_top3.id > HN/pcoa_bray_samples_top3.grp
paste HN/pcoa_bray_samples_top3.id HN/pcoa_bray_samples_top3.grp > HN/pcoa_bray_samples_top3.group
usearch10 -otutab_group HN/otutab1.txt -labels HN/pcoa_bray_samples_top3.group -output HN/otutab2.txt
# d. 筛选HN有表型的196个
head -n1 otu_filter/7otu_table_ratio_LN_top100.txt|cut -f 2-|sed 's/\t/\n/g' > result/variety196.id
usearch10 -otutab_sample_subset HN/otutab2.txt -labels result/variety196.id -output HN/otutab3.txt
usearch10 -otutab_stats HN/otutab3.txt -output HN/otutab3.stat
cat HN/otutab3.stat # 从77M至19.8M

## OTU表及相关数据

### a. 筛选HN最终表
#最终筛选的结果归为otutab.txt
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



# 1.3 Gemma分析
cd /mnt/bai/yongxin/rice/miniCore/180319/
mkdir -p gemma
cd gemma
# bed和bim文件
cp /mnt/zhou/chulab/miniCore/snp1.5x/T2.b* ./
# kinship文件
cp /mnt/bai/xiaoning/past/software/GEMMA/output/T2.cXX.txt kinship.txt
# pca文件作为co variation
cp /mnt/bai/xiaoning/past/software/gcta_1.91.3beta/sativa.pca.eigenvec pca4.txt
# 表型数据，以属为例 HN/tax_g.txt.fam
cp ../HN/tax_g.txt.fam T2.fam
# 以第6列正对照wax为例
gemma -bfile T2 \ # 读取bed bim fam文件
	-k kinship.txt \ # 品种间关系
	-lmm 4 \ # 4种计算模型
	-n 1 \ # 分析表型所有fam文件中的列，第6列为1
	-o 1 \ # 输出文件名
	-c pca4.txt # covariates文件
# 显示各列编号和名称
head -n1 output/6.assoc.txt|tr '\t' '\n'|awk '{print NR,$0}'
# 按14列筛选p<0.001
cat <(head -n 1 output/6.assoc.txt) <(awk '$14<0.001' output/6.assoc.txt) > output/6.assoc.txt.sig
# 绘制曼哈顿图
cp /mnt/bai/xiaoning/project/jingying/Unoise_NRT/gwas_plot.R ../scripts/ # 备份源码
cp /mnt/bai/xiaoning/project/jingying/Unoise_NRT/gwas_plot.R ./ # 当前文件夹内修改
Rscript qqman.R -i output/6.assoc.txt.sig -t 14 -o output/6.assoc.txt.sig.pdf



# 并行对每个表型处理，在alpha为例
#生成列与名称对应的表，去掉前5行，取前99行，删除()"等符号，添加行号为gemma找列和输出文件名
tail -n+6 ../HN/alpha.txt.fam.header|head -n99|sed 's/[\(\)\"]//g'|awk '{print NR"\talpha"NR$0}' > alpha.txt
# parallel -j 33 "gemma -bfile T2 -k kinship.txt -c pca4.txt -lmm 4 -n {1} -o {2}" ::: `cut -f 1 alpha.txt` ::: `cut -f 2 alpha.txt` # 此种方法为组全式，即14x14，不可取
	parallel -j 33 "gemma -bfile T2 -k kinship.txt -c pca4.txt -lmm 4 -n {1} -o {1}" ::: `cut -f 1 alpha.txt`
#批量改名和筛选
	awk 'BEGIN{OFS=FS="\t"}{system("mv output/"$1".assoc.txt output/"$2".assoc.txt");system("mv output/"$1".log.txt output/"$2".log.txt")}' <(cat alpha.txt)
	# 注意变量$符需要转义
	parallel -j 33 "cat <(head -n 1 output/{1}.assoc.txt) <(awk '\$14<0.001' output/{1}.assoc.txt) > output/{1}.assoc.txt.sig" ::: `cut -f 2 alpha.txt`
#批量绘制manhanttan plot
	parallel -j 33 "Rscript qqman.R -i output/{1}.assoc.txt.sig -t 14" ::: `cut -f 2 alpha.txt`



# 并行对每个表型处理 alpha, pcoa_unifrac_15, order, genus, faprotax
i=pcoa_unifrac_15
#生成列与名称对应的表，去掉前5行，取前99行，删除()"等符号，添加行号为gemma找列和输出文件名
	cp ../HN/${i}.txt.fam T2.fam
	tail -n+6 ../HN/${i}.txt.fam.header|head -n99|sed 's/[\(\)\"]//g'|awk '{print NR"\t"$0}'| sed "s/\t/\t${i}/" > ${i}.txt
	parallel --xapply -j 33 "gemma -bfile T2 -k kinship.txt -c pca4.txt -lmm 4 -n {1} -o {2}" ::: `cut -f 1 ${i}.txt` ::: `cut -f 2 ${i}.txt`
#批量改名
#	awk 'BEGIN{OFS=FS="\t"}{system("mv output/"$1".assoc.txt output/"$2".assoc.txt");system("mv output/"$1".log.txt output/"$2".log.txt")}' <(cat ${i}.txt)
	# 筛选1e-3显著的绘图，注意变量$符需要转义
	parallel -j 33 "cat <(head -n 1 output/{1}.assoc.txt) <(awk '\$14<0.001' output/{1}.assoc.txt) > output/{1}.assoc.txt.sig" ::: `cut -f 2 ${i}.txt`
#批量绘制manhanttan plot
	parallel -j 33 "Rscript qqman.R -i output/{1}.assoc.txt.sig -t 14" ::: `cut -f 2 ${i}.txt`
# 添加注释
#	parallel -j 33 "awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[\$1]=\$2"\t"\$3"\t"\$4"\t"\$5} NR>FNR{print \$0,a[\$2]}' <(sed '1 s/ID/rs/' snp.anno) output/{1}.assoc.txt.sig |sed 's/\t$/\tNA/'|cut -f 1-3,14- > output/{1}.assoc.txt.sig.anno" ::: `cut -f 2 ${i}.txt`
for RPM in `cut -f 2 ${i}.txt`; do
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3"\t"$4"\t"$5} NR>FNR{print $0,a[$2]}' <(sed '1 s/ID/rs/' snp.anno) output/${RPM}.assoc.txt.sig |sed 's/\t$/\tNA/'|cut -f 1-3,14- > output/${RPM}.assoc.txt.sig.anno
done

# 筛选修改基因列表
# 显示所有候选基因
mkdir -p list
grep -P 'OS08G0155400|OS03G0203200|OS06G0110000|OS11G0104300|OS01G0746400|OS04G0550600|OS11G0587000|OS09G0441900|OS05G0407500|OS05G0333200|OS01G0927000|OS05G0196500|OS04G0509600|OS01G0909200|OS02G0699000|OS01G0883800' output/*.anno > list/candidate.txt
# 显示所有小于1e-5的重要SNP
grep -P 'HIGH|MODERATE' output/*.anno | awk '$4<0.00001' | sort -k4,4g |less > list/VIP_e-5.txt



# SNP注释，3列为SNPID，8列为注释
cp /mnt/bai/xiaoning/test_1/variant_calling/data/test1/snpEff/T2.ann.vcf gemma/
snp_anno=/mnt/bai/yongxin/rice/miniCore/180319/gemma/T2.ann.vcf
# 制作可读的注释
#提取ID，注释，效应和基因ID
grep -v "#" $snp_anno|cut -f 3,8|cut -f 1 -d ','|tr "|" "\t"|cut -f 1,3,4,6|sed '1 s/.*/ID\tAnnotation\tAnnotation_impact\tLocus_ID/'>temp
#rapdb注释,2列Locus_ID，3列Description
gene_anno=IRGSP-1.0_representative_annotation_2018-03-29.tsv
wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/${gene_anno}.gz 
gunzip ${gene_anno}.gz
# 追加Description至temp上，转换大小写snpeff一致
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$4]}' <(cut -f 2-3 ${gene_anno}|sed 's/Os/OS/;s/g/G/') temp > snp.anno1
# 添加RAPG ID
wget http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU_2018-03-29.txt.gz
gunzip RAP-MSU_2018-03-29.txt.gz
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$4]}' <(cut -f 1 -d '.' RAP-MSU_2018-03-29.txt|sed 's/Os/OS/;s/g/G/'|sed '1 i Locus_ID\tRGAP') snp.anno1 > snp.anno
# 追加到sig上的annotation
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3"\t"$4"\t"$5} NR>FNR{print $0,a[$2]}' <(sed '1 s/ID/rs/' snp.anno) output/otutab_freq10OTU_4.assoc.txt.sig |sed 's/\t$/\tNA/'|cut -f 1-3,14- > output/otutab_freq10OTU_4.assoc.txt.sig.anno
sort -k4,4g output/otutab_freq10OTU_4.assoc.txt.sig.anno|less # 按P值排序，找到10号上的峰
# 查询 是否有相关基因
grep -P 'OS08G0155400|OS03G0203200|OS06G0110000|OS11G0104300|OS01G0746400|OS04G0550600|OS11G0587000|OS09G0441900|OS05G0407500|OS05G0333200|OS01G0927000|OS05G0196500|OS04G0509600|OS01G0909200|OS02G0699000|OS01G0883800' output/otutab_freq10OTU_4.assoc.txt.sig.anno



# 计算OTU表与PCoA轴的相关性，筛选与主轴相关的差异菌

# 籼粳稻分化基因 IndTej

	# 研究籼粳稻分化基因与菌群的关系，重点关注NRT1.1b
	dir=IndTej
	mkdir -p ${dir}

## 筛选籼粳稻的样本

	#SNP共有样品添加分类注释，并筛选IND/TEJ
	awk 'NR==FNR {a[$1]=$2} NR>FNR {print $0"\t"a[$1]}' doc/design_group7.txt result/variety196.id |grep -P 'IND|TEJ' > ${dir}/variety.id
	cut -f 2 ${dir}/variety.id|sort|uniq -c # 97, 67:30

## Top100 OTUs in HN
	
	mkdir -p ${dir}/HN
	cut -f 1 ${dir}/variety.id>${dir}/variety.id1
	# 标准化的结果筛选会变为0，用R实现
#	usearch10 -otutab_sample_subset HN/otutab_freq.txt -labels ${dir}/variety.id1 -output ${dir}/HN/otutab_freq.txt
#	usearch10 -otutab_sample_subset HN/faprotax_freq.txt -labels ${dir}/variety.id1 -output ${dir}/HN/faprotax_frep.txt
	# 筛选样品，再标准化
	usearch10 -otutab_sample_subset HN/otutab.txt -labels ${dir}/variety.id1 -output ${dir}/HN/otutab.txt
	usearch10 -otutab_counts2freqs ${dir}/HN/otutab.txt -output ${dir}/HN/otutab_freq.txt
	usearch10 -otutab_sample_subset HN/faprotax.txt -labels ${dir}/variety.id1 -output ${dir}/HN/faprotax.txt
	usearch10 -otutab_counts2freqs ${dir}/HN/faprotax.txt -output ${dir}/HN/faprotax_freq.txt

	# 循环处理两个变量
	for i in otutab_freq faprotax_freq; do Rscript scripts/gemma_fam.r -i ${dir}/HN/${i}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam; done
	# 准备gemma文件
	#ln gemma/T2.bim ./
	#ln gemma/T2.bed ./
	cp ${dir}/HN/otutab_freq.txt.fam T2.fam
	i=otutab_freq
	#生成列与名称对应的表，去掉前5行，取前99行，删除()"等符号，添加行号为gemma找列和输出文件名
	tail -n+6 ${dir}/HN/${i}.txt.fam.header|head -n99|sed 's/[\(\)\"]//g'|awk '{print NR"\t"NR$0}'| sed "s/\t/\t${i}/" > ${dir}/HN/${i}.list
	parallel -j 33 "gemma -bfile T2 -k gemma/kinship.txt -c gemma/pca4.txt -lmm 4 -n {1} -o {1}" ::: `cut -f 1 ${dir}/HN/${i}.list`
	#批量改名和筛选
	awk 'BEGIN{OFS=FS="\t"}{system("mv output/"$1".assoc.txt output/"$2".assoc.txt");system("mv output/"$1".log.txt output/"$2".log.txt")}' <(cat ${dir}/HN/${i}.list)
	# 注意变量$符需要转义
	parallel -j 33 "cat <(head -n 1 output/{1}.assoc.txt) <(awk '\$14<0.001' output/{1}.assoc.txt) > output/{1}.assoc.txt.sig" ::: `cut -f 2 ${dir}/HN/${i}.list`
	#批量绘制manhanttan plot
	parallel -j 33 "Rscript gemma/qqman.R -i output/{1}.assoc.txt.sig -t 14" ::: `cut -f 2 ${dir}/HN/${i}.list`
	# 添加SNP注释,parallel有太多""需要转义不方便，改用循环
#	parallel -j 33 "awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3"\t"$4"\t"$5} NR>FNR{print $0,a[$2]}' <(sed '1 s/ID/rs/' gemma/snp.anno) output/{1}.assoc.txt.sig |sed 's/\t$/\tNA/'|cut -f 1-3,14- > output/{1}.assoc.txt.sig.anno" ::: `cut -f 2 ${dir}/HN/${i}.list|head -n3`
	for RPM in `cut -f 2 ${dir}/HN/${i}.list`; do
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3"\t"$4"\t"$5} NR>FNR{print $0,a[$2]}' <(sed '1 s/ID/rs/' gemma/snp.anno) output/${RPM}.assoc.txt.sig |sed 's/\t$/\tNA/'|cut -f 1-3,14- > output/${RPM}.assoc.txt.sig.anno &
	done
#	# 移至分类目录
#	mv output/* ${dir}/HN/
#	# 按顺序查看显著的基因
#	for RPM in `cut -f 2 ${dir}/HN/${i}.list`; do
#	sort -k4,4g ${dir}/HN/${RPM}.assoc.txt.sig.anno|sed "1 s/$/${RPM}/"|less -S # sort by pvalue, focus on chr10 21.7M end
#	done


	# 查看OTU的描述：重点关注OTU25,7,13,23
	less result/rep_seqs_tax.txt
	# 检索NRT1.1b，最高p<e-5
	grep 'OS10G0554200' IndTej/HN/*.assoc.txt.sig.anno|awk '$4<0.001'|sort -k4,4g|less -S 
	grep 'OS08G0155400' IndTej/HN/*.assoc.txt.sig.anno|awk '$4<0.001'|less # NRT1.1A无显著，也无非同义SNP
	grep 'OS03G0707600' IndTej/HN/*.assoc.txt.sig.anno|awk '$4<0.001'|sort -k4,4g| less # della


2018/4/25 查找高低氮下Ind氨化效率高品种

	~/rice/miniCore/180319/IndTej/script/Indica_ammoniation_variety.R
	筛选了每个品种氨化的品种排序
	manual/nitrate_ammonification_IND_[HL]N.txt # 合并排序，发现HN/LN中高氨化能力的品种前30个有24个共有

