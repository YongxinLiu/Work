#!/bin/bash

# 16S扩增子操作手册 16S Amplicon pipeline manual, v1.4 2017/12/23

# 按如下命令执行，并修改makefile中的参数
# 进入工作目录
cd /mnt/bai/yongxin/wheat/profile



# 0. 配置实验设计参数

# 复制最新版流程配置文件+操作手册至工作目录
cp ~/rice/timecourse/ma* ./

# 修改工作目录wd
pwd

# 初始化工作目录
make init

# 人工编写实验设计，按design.xlsx中要求每个表另存为txt格式上传于doc/文件夹中
# 编辑并上传实验设计至doc/，必选(mappingfiles/summary/material/library/group_compare), 可选(group_venn, group_tern)

# 如文库较多，可将所有库保存于design_raw.txt中，使用split_design.pl命令生成按库分割的mappingfile
dos2unix doc/*
split_design.pl -i doc/design_raw.txt
# 基于mappingfile生成实验设计，需要指定一个存在文件提取header，默认为L1
cat <(head -n1 doc/L1.txt|sed 's/#//g') doc/L* |grep -v '#' > doc/design.txt
# 检查实验设计首行是否名称唯一
wc -l doc/design.txt && cut -f 1 doc/design.txt | sort | uniq | wc -l

# 上传原始数据*.fq.gz至clean_data/
for i in AGTCAA AGTTCC ATGTCA CCGTCC GTAGAG GTCCGC GTGAAA GTGGCC GTTTCG CGTACG;do
	ln ~/seq/171121.lane2.wheat.rice/L1-Index-${i}/FCHY555BCXY_L1_CWHPEPI00001607_Index-${i}_1.fq.gz clean_data/${i}_1.fq.gz
	ln ~/seq/171121.lane2.wheat.rice/L1-Index-${i}/FCHY555BCXY_L1_CWHPEPI00001607_Index-${i}_2.fq.gz clean_data/${i}_2.fq.gz
done


# 简化文件名(可选) 
rename 's/Med-//g;s/_HLY73BCXY_L1//g;s/clean\.//g' *.gz 
# 原始数据文件名重命名为L1/2/3_1/2.fq.gz，并格式化文本信息
make rename
# 修改文库列表list
# ll clean_data/*.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|sort -u|tr "\n" " "
tail -n+2 doc/library.txt |cut -f 1|tr "\n" " " # 另一种方法，但需要有library.txt
# 显示实验设计,修改实验的主要和次要分组类型g1/g2，如groupID, genotype, compartment, soiltype，没有g2可使用默认batch
head -n 3 doc/design.txt
# 获取主、次、第三分组信息：具体根据实验设计进一步筛选
tail -n+2 doc/design.txt|cut -f 5|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
# 最简单的实验只需一组信息，第二组填batch和1；复杂实验需要次、三分组信息
# 仔细全文阅读配置文件，修改所有与本项目相关参数



# 1. 数据质控
make qc # 原始数据质控
make merge # 双端数据合并再质控
make extract_barcodes # 提取样本的barcodes
make split # 拆分样品
make cutadapt # 切除双端接头
make stat # 统计文库各步操作stat each library

make multiqc # 汇总数据合并前后的评估报告

grep -v -c '#' doc/L* | sed 's/doc\/L//;s/.txt:/\t/' | sort -k1,1n # 统计各库中样品数据
cut -f 1,20 clean_data/multiqc_data/multiqc_fastqc.txt|grep '_1'|sed 's/L//;s/_1//;s/\.0//'|sort -k1,1n # 展示库数据量，数据列可变，有可能为13，18，19，20


# 2. 标准流程
make merge_library # 合并所有文库，并简单统计绘图
make derep # 数据去冗余
make unoise # 鉴定生物特征序列
make rm_chimeras # (可选)去嵌合体，执行完Unoise后，不执行此步可直接跳过
make otu_table # 生成OTU表
make assign_tax # 基于gg13.8注释物种
make tree # 构建进化树，用于多样性分析
make alpha # 计算alpha多样性
make alpha_usearch_rare # (可选)，采用usearch结果绘制
make alpha_qiime_rare # (可选)，采用qiime进行稀释曲线分析，大量抽样比较费时；可继续运行其它步骤
make beta # 计算beta多样性
make graphlan # 绘制高丰度物种树
make ggtree # 绘制高丰度进化树
make diversity # 多样性绘制箱线和散点

cat result/qc.sum # 显示数据量统计


# 3. 个性分析
make filter # 过滤高丰度菌分析，默认0.05%。如lefse只有1%，ruben文章常用0.5-0.1%
make rediv # 新表重新多样性分析
make draw_div # 多样性绘图
make draw_tax # 物种门-属级别比较
make draw_otu # OTUs水平差异比较
make draw_ter # OTUs三元图比较

make alpha_usearch_rare # usearch rarefraction
make alpha_qiime_rare # alpha rarefraction


make rmd # 生成可读性报告 write report 



# 4. 可选高级分析
make culture_graphlan # 绘制高丰度菌可培养比例物种树

# 5. 新功能测试 
# rarefraction
alpha_rare_usearch.sh -d /mnt/bai/yongxin/wheat/NP/doc/design.txt -m FALSE -A groupID -B '"BSNP","RsBDHNP","BDHNP"' -C compartment -D '"rhizosphere","root","soil"' -o result -g TRUE -h 5 -w 8 -s 7



# 手动分析部分
# 1. 参考时间序列分析
git clone https://github.com/bulksoil/LifeCylceAnalysis # 只有几行代码和部分数据

# 2. 分析不同基因型、时间、深度PCoA result_k1-c/01diversity.r

# 3. 分析taxonomy门水平，两个基因型，soil, L1, L2; result_k1-c/02taxonomy_egr.r
awk '$6=="soil" && $8=="1"' doc/design.txt|sort -k7,7n|cut -f 5|uniq|awk '{print "\""$1"\""}'|tr "\n" "," # 筛选第一土层，时间排序
awk '$6=="XY54" && $8=="1"' doc/design.txt|sort -k7,7n|cut -f 5|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
awk '$6=="XY54" && $8=="2"' doc/design.txt|sort -k7,7n|cut -f 5|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
awk '$6=="J411" && $8=="1"' doc/design.txt|sort -k7,7n|cut -f 5|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
awk '$6=="J411" && $8=="2"' doc/design.txt|sort -k7,7n|cut -f 5|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
make test_part3 # 下载并修改文件名 tax_stack_phylumpro_abc.pdf -- tax_stack_phylumpro_J411_1_group.pdf, tax_stack_phylumpro_sample.pdf -- tax_stack_phylumpro_J411_1_sample.pdf

# 4. 按root、XY54、J411分别出报告，XY54中样品在需要剔除
# 筛选某一基因型，按时间和深度分别排序
awk '$6=="soil"' doc/design.txt|sort -k7,7n -k8,8n|cut -f 5|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
# a. 基于上行结果修改分组，b. doc改为doc/soil，version=wheat_profile_soil0
mkdir -p doc/soil
sed 's/XY54/soil/g' doc/group_compare.txt > doc/soil/group_compare.txt
sed '1 s/$/soil/g' doc/summary.txt > doc/soil/summary.txt
rm diversity draw_div
make diversity
make rmd

awk '$6=="XY54"' doc/design.txt|sort -k7,7n -k8,8n|cut -f 5|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
# a. 基于上行结果修改分组，b. doc改为doc/XY54，version=wheat_profile_XY540
mkdir -p doc/XY54
sed 's/XY54/XY54/g' doc/group_compare.txt > doc/XY54/group_compare.txt
sed '1 s/$/XY54/g' doc/summary.txt > doc/XY54/summary.txt
rm diversity draw_div
make diversity
make rmd

awk '$6=="J411"' doc/design.txt|sort -k7,7n -k8,8n|cut -f 5|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
# a. 基于上行结果修改分组，b. doc改为doc/J411，version=wheat_profile_J4110
mkdir -p doc/J411
sed 's/XY54/J411/g' doc/group_compare.txt > doc/J411/group_compare.txt
sed '1 s/$/J411/g' doc/summary.txt > doc/J411/summary.txt
rm diversity draw_div
make diversity
make rmd

# 5. 绘制组内bray-curtis distance箱线图分析 beta_group_inner.r
开始按各genotype分为三组绘制组内箱线图，没规律。再按土层分开，有明显季节周期规律。

# 6. 绘制几个主要OTUs的时间曲线 DAOTU_OTUs.r
"OTU_2","OTU_6","OTU_7","OTU_12","OTU_13"
绘制OTU_2, OTU以上全部在同一基因型L1下变化 规律

# 7. 绘制组间相关性 DAOTU_BC_corplot.r
采用按组合并，vegdist计算距离，再使用corplot可视化，并与pearson相关系数 比较


梯度上色，按组等间隔着色；讲什么要先说明它的意义和应用。

7天根系有土套；

XY54 228L天时有Tenericutes感染





# 常见问题
1. split_libraries_fastq.py: 
a. incorrect value for phred_offset, : 修改phred_score 33/64
b. 右端引物中检中大量N
复制makefile_16s至doc目录，程序中新增--sequence_max_n 3，
cp temp/L1_split/split_library_log.txt temp/L1_split/split_library_log.txt.0
rm L1.split
make split # N序列从80%降至10%，可接受

2. 批处理文库
for i in L1 L2 L3 L4 L5 L6 L7 L8 L9 L10; do
	unzip ${i}_1_fastqc.zip
	unzip ${i}_2_fastqc.zip
done

2. miniCore LN的文库中clean数据极少，将一个数据级
检查qc.sum文件发现是merge中绝大多数无法合并
检查右端质量L21与L1比较，发现L21从135以后，质量进差，甚至中位数全接近最小值。

3. awk NR FNR ref
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$2"\t"$3"\t"$2$3"\t"$2$3} NR>FNR {print $0,a[$6]$7}' design0_materialID.txt design1_raw.txt | cut -f 11,1-3,10,8-9 | sed 's/#SampleID/GroupID/' | sed 's/#Sample\tID/variety\tDescription/' |sed 's/#SampleIDreplicate/#SampleID/' > design2_clean.txt

4. 删除部分更新文库的的完成标准 qc merge
for i in `seq 1 20`;do
#	rm L${i}.merge
	fastqc L${i}_1.fq.gz --quiet --extract &
done

5. 判断文件是否存在
venn="/mnt/bai/yongxin/rice/strigolactone.LiJY/doc/group_venn.txt"
if [ ! -f "${venn}" ]; then
	echo "hello"
fi
if test ! -f "${venn}" ; then
	echo "hello"
fi
# 但makefile中无法运行
http://blog.csdn.net/zhangxinrun/article/details/19561567 使用makefile中的语法也无法运行
if test $(shell if [ -f $(venn) ]; then echo "TRUE"; else echo "FALSE"; fi;); then ;echo "hello";fi
