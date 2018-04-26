<!-- TOC -->

- [1. 基于OTU分析亚种差异](#1-基于otu分析亚种差异)
    - [1.1. HN](#11-hn)
            - [**硝酸盐氨化在Indica中富集**](#硝酸盐氨化在indica中富集)
    - [1.2. LN](#12-ln)
        - [1.2.1. 筛选低氮的样品](#121-筛选低氮的样品)
        - [1.2.2. 筛选LN OTU及相关表](#122-筛选ln-otu及相关表)
    - [1.3. (HN-LN)/HN](#13-hn-lnhn)
- [2. 籼粳差异微生物组与氮吸收NRT](#2-籼粳差异微生物组与氮吸收nrt)
    - [2.1. 籼粳多样性分析](#21-籼粳多样性分析)
        - [2.1.1. Alpha多样性](#211-alpha多样性)
        - [2.1.2. Beta多样性](#212-beta多样性)
        - [gemma新pcoa](#gemma新pcoa)
    - [2.2. GWAS关联分析](#22-gwas关联分析)
    - [2.3. OTUs](#23-otus)
    - [2.4. PCoA](#24-pcoa)
    - [2.5. 元素循环Faprotax](#25-元素循环faprotax)
- [3. 附录](#3-附录)
    - [3.1. 检查Gemma丢弃SNP原因](#31-检查gemma丢弃snp原因)
    - [3.1. 检查 10m21759092 SNP基因](#31-检查-10m21759092-snp基因)
    - [3.2 使用GAPI分析LN WAX/PC2和模拟数据](#32-使用gapi分析ln-waxpc2和模拟数据)
    - [3.2 使用GAPI分析部分SNP](#32-使用gapi分析部分snp)
- [文章图表](#文章图表)
    - [Fig1. 实验设计+多样性](#fig1-实验设计多样性)
    - [Fig2. 差异菌](#fig2-差异菌)
        - [Fig2.1 差异菌 Lefse展示](#fig21-差异菌-lefse展示)
        - [Fig2.2 机器学习找marker菌](#fig22-机器学习找marker菌)
        - [Fig2.2 氮利用上存在差异](#fig22-氮利用上存在差异)
    - [Fig3. 网络分析](#fig3-网络分析)
    - [Fig4. NRT品种与基因](#fig4-nrt品种与基因)

<!-- /TOC -->

**本课题在rice miniCore项目下，只以籼稻TEJ和粳稻IND两个亚种进行分析**
    
    # 工作目录为IndTej子目录
    cd ~/rice/miniCore/180319/
    wd=/mnt/bai/yongxin/rice/miniCore/180319
    
# 1. 基于OTU分析亚种差异
OTU表分为三类：高氮HN，低氮LN，高氮 与低氮之差比例(HN-LN)/HN

## 1.1. HN
    
    cd IndTej/HN/  

    # 准备makefile文件
    cp ~/rice/rice.epi/180409/ma* ./

    # 修改工作目录，然后
    make init

    # 准备预分析结果: 实验设计，OTU表和代表序列
    cat $wd/doc/design_group7.txt | sed '1 s/Region/groupID/' > doc/design.txt
    cp $wd/HN/otutab.txt temp/otu_table_raw.txt
    cp $wd/result/rep_seqs.fa result/rep_seqs.fa
    # 编写实验设计：doc/group_compare.txt # IND TEJ
    
    # 运行流程
    touch merge_librar
    touch derep
    touch cluster_otu
    touch otu_table

    # 从这里开始
    make otu_table_filter
    make assign_tax
    make tree
    make alpha
    make beta
    make graphlan
    make ggtree # 报错是因为实验设计只有两列，去除行名仅剩一列无法成为数据框，保证实验设计至少3列。
    make diversity
    
    make rmd # draw_tax报错因为行名重复，可跳过继续

    # 计算HN下KEEG
    mkdir -p function
    touch picrust_calc
    cp ../../HN/KO*.txt function/
    # 修改makefile中输入他们的为KO1.txt
    make picrust_draw # 分别做2/3级差异，无氮相关
    
    # 计算元素循环faprotax
    mkdir -p faprotax
    touch faprotax_calc
    cp ../../HN/faprotax.txt faprotax/
    # 修改makefile中输入他们的为KO1.txt
    make faprotax_draw # 分别做2/3级差异，无氮相关   


#### **硝酸盐氨化在Indica中富集**
    
按组分别进行T-test/Wilcoxon检查所有氮相关功能，DAKO_egr_boxplot.r中。
发现11个氮相关功能中，Indica的nitrate_ammonification和nitrite_ammonification显著高，nitrogen_respiration, nitrite_respiration, nitrate_respiration，且在HN和LN中一致。而nitrogen_fixation更低。
    
查看相关菌  
    
    grep -P -A 14 '^# nitrate_ammonification' ../../faprotax/report

最要是OTU_9, OTU_102等，但在差异OTU中没有显著。但查看原始数据比较结果的均值差异明显。改用wilcoxon检验，其为上调的第 1/2 个。
    
秩合检验更靠谱，负二项在菌群 研究中经常出现异常结果
    
## 1.2. LN

基于HN的方法制作LN相关表型

    cd IndTej/LN/  

    # 准备makefile文件
    cp ../HN/ma* ./
    # 修改工作目录，然后
    make init
    # 准备预分析结果: 实验设计，OTU表和代表序列
    cp $wd/LN/otutab.txt temp/otu_table_raw.txt
    cp $wd/result/rep_seqs.fa result/rep_seqs.fa
    cp ../HN/doc/* doc/
    
    # 运行流程, 从这里开始
    touch merge_library derep cluster_otu otu_table
    make diversity
    make rmd # draw_tax报错因为行名重复，可跳过继续

    # 计算LN下KEEG
    mkdir -p function
    touch picrust_calc
    cp ../../LN/KO*.txt function/
    # 修改makefile中输入他们的为KO1.txt
    make picrust_draw # 分别做2/3级差异，无氮相关
    
    # 计算元素循环faprotax
    mkdir -p faprotax
    touch faprotax_calc
    cp ../../LN/faprotax.txt faprotax/
    # 修改makefile中输入他们的为KO1.txt
    make faprotax_draw # 分别做2/3级差异，无氮相关


### 1.2.1. 筛选低氮的样品

    # a. 获得所有LN的样品名和组名，筛选LN下的样品
    grep 'L$' doc/design.txt | cut -f 1,5> LN/design.txt
    Rscript scripts/otutab_sample_subset.r -i result/otu_table_filter.txt -d LN/design.txt -o LN/otutab0.txt
    
    # b. 计算bray_curtis距离，筛选每组内的距离，最小Top3输出，并筛选
    usearch10 -beta_div LN/otutab0.txt -tree result/rep_seqs.fa.tree -filename_prefix LN/ -metrics bray_curtis
    Rscript scripts/beta_pcoa_group.r -i LN/bray_curtis.txt -d doc/design.txt -n GroupID -o LN/pcoa_bray
    Rscript scripts/otutab_sample_subset.r -i LN/otutab0.txt -d LN/pcoa_bray_samples_top3.txt -o LN/otutab1.txt

    # c. 按组合并
    paste <(tail -n+2 LN/pcoa_bray_samples_top3.txt|cut -f 1) <(tail -n+2 LN/pcoa_bray_samples_top3.txt|cut -c1-5) > LN/pcoa_bray_samples_top3.group
    usearch10 -otutab_group LN/otutab1.txt -labels LN/pcoa_bray_samples_top3.group -output LN/otutab2.txt

    # d. 筛选LN有表型的196个
    usearch10 -otutab_sample_subset LN/otutab2.txt -labels result/variety196.id -output LN/otutab3.txt
    usearch10 -otutab_stats LN/otutab3.txt -output LN/otutab3.stat
    cat LN/otutab3.stat

### 1.2.2. 筛选LN OTU及相关表

    cp LN/otutab3.txt LN/otutab.txt
    # 标准化为比例
    usearch10 -otutab_counts2freqs LN/otutab.txt -output LN/otu.txt # 标准化为1

    # a. alpha多样性分析
    usearch10 -otutab_norm LN/otutab.txt -sample_size 10000 -output LN/otutab_norm.txt 
    usearch10 -alpha_div LN/otutab_norm.txt -output LN/alpha.txt 
    head -n1 LN/alpha.txt|tr '\t' '\n'|awk '{print NR,$0}' #共14种多样性指数，其中chao1, ricLNess和shannon_e分别在4，10和13列
    Rscript scripts/gemma_fam.r -i LN/alpha.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam -t FALSE

    ### b. beta多样性
    mkdir -p LN/beta
    usearch10 -beta_div LN/otutab.txt -tree result/rep_seqs.fa.tree -filename_prefix LN/beta/
    # 对三种距离和PCoA进行计算，在 unifrac_binary PCoA 3/4轴上 # R4159异常，检查alpha中robbins为零，ricLNess倒数第4(sort -k10,10n LN/alpha.txt|less)
    for i in bray_curtis unifrac unifrac_binary; do
        Rscript scripts/beta_pcoa_dist.r -i LN/beta/${i}.txt -d doc/design_group7.txt -t ${i} -n Region -o LN/pcoa_${i}_ -w 8 -e 5
    done
    ### Constrained PCoA by group，PC1-5 LN/pcoa_bray_constrained_15.txt
    Rscript scripts/beta_cpcoa.r -i LN/otutab.txt -d doc/design_group7.txt -n Region -o LN/pcoa_bray_constrained -w 8 -e 5
    
    # 转换四4种距离主坐标轴1-5为表型 
    for i in bray_curtis unifrac unifrac_binary bray_constrained; do
        Rscript scripts/gemma_fam.r -i LN/pcoa_${i}_15.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam -t FALSE
    done


    # c. 制作各分类级结果 /mnt/bai/yongxin/rice/miniCore/180319/LN/tax*
    for i in p c o f g;do
        usearch10 -sintax_summary temp/otus_no_host.tax -otutabin LN/otutab.txt -rank ${i} -output LN/tax_${i}.txt
        Rscript scripts/gemma_fam.r -i LN/tax_${i}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam
    done

    # d. 制作为功能分类级，依赖picrust结果
    for l in 1 2 3 4;do 
        # 筛选样品
        usearch10 -otutab_sample_subset function/metagenome_predictions.L${l}.txt -labels LN/pcoa_bray_samples_top3.id -output temp/temp_${l}_1
        # 按组合并
        usearch10 -otutab_group temp/temp_${l}_1 -labels LN/pcoa_bray_samples_top3.group -output temp/temp_${l}_2
        # 筛选LN有表型的196个
        usearch10 -otutab_sample_subset temp/temp_${l}_2 -labels result/variety196.id -output LN/KO${l}.txt
        # 统计结果
        usearch10 -otutab_stats LN/KO${l}.txt -output LN/KO${l}.stat
        cat LN/KO${l}.stat
        Rscript scripts/gemma_fam.r -i LN/KO${l}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam
    done

    # e. 元素循环farprotax
 
    Rscript scripts/otutab_sample_subset.r -i faprotax/otu_table_tax.faprotax -d LN/pcoa_bray_samples_top3.txt -o temp/temp_far_1 # 筛选样品
    usearch10 -otutab_group temp/temp_far_1 -labels LN/pcoa_bray_samples_top3.group -output temp/temp_far_2 # 按组合并
    # 筛选LN有表型的196个
    usearch10 -otutab_sample_subset temp/temp_far_2 -labels result/variety196.id -output LN/faprotax.txt
    usearch10 -otutab_counts2freqs LN/faprotax.txt -output LN/fapro.txt # 标准化为1

    # f. 转换otu和farpro微生物组数据为Gemma表型
    for i in otu faprox; do
        Rscript scripts/gemma_fam.r -i LN/${i}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam
    done


## 1.3. (HN-LN)/HN





# 2. 籼粳差异微生物组与氮吸收NRT

    cd rice/miniCore/180319/

    # 工作目录IndTej
    mkdir -p IndTej/
    
    # 筛选籼粳品种ID
    awk 'NR==FNR {a[$1]=$2} NR>FNR {print $0"\t"a[$1]}' doc/design_group7.txt result/variety196.id |grep -P 'IND|TEJ' > IndTej/variety.id
    # 统计数量67：30
	cut -f 2 IndTej/variety.id|sort|uniq -c
    # 建立标准实验设计，方便绘图
    awk '{print $1"\t"$1"\t"$2}' IndTej/variety.id|sed '1 i SampleID\tSampleID\tgroupID'|less>IndTej/design.txt

    # 先以低氮LN为例，初步分析低氮的结果差异大
    mkdir -p IndTej/LN
     # 筛选OTU表
    Rscript scripts/otutab_sample_subset.r -i LN/otutab.txt -d IndTej/variety.id -o IndTej/LN/otutab.txt


## 2.1. 籼粳多样性分析

    # 统计OTU表
    usearch10 -otutab_stats IndTej/LN/otutab.txt -output IndTej/LN/otutab.stat
    cat IndTej/LN/otutab.stat

    # 抽样标准化至3万
    usearch10 -otutab_norm IndTej/LN/otutab.txt -sample_size 30000 -output IndTej/LN/otutab_norm.txt 
    usearch10 -otutab_stats IndTej/LN/otutab_norm.txt -output IndTej/LN/otutab_norm.stat
    cat IndTej/LN/otutab_norm.stat

### 2.1.1. Alpha多样性

    # 计算14种多样性指数
    usearch10 -alpha_div IndTej/LN/otutab_norm.txt -output IndTej/LN/alpha.txt 
    # 绘制箱线图，脚本位于script目录中，结果在文件以方法开头；居然全有显著差异3.8e-8
	alpha_boxplot.sh -i IndTej/LN/alpha.txt -m '"chao1","richness","shannon_e"' -d IndTej/design.txt -A groupID -o IndTej/LN/
   
### 2.1.2. Beta多样性

    # 计算标准化OTU表的距离
    mkdir -p IndTej/LN/beta/
    usearch10 -beta_div IndTej/LN/otutab.txt -tree result/rep_seqs.fa.tree -filename_prefix IndTej/LN/beta/
    # 绘图及输出主坐标轴PC1-5
	beta_pcoa.sh -i IndTej/LN/beta/ -m '"jaccard_binary","bray_curtis","unifrac_binary","unifrac"' -d IndTej/design.txt -A groupID -o IndTej/LN/beta/ -c doc/compare.txt # compare里IND TEJ

### gemma新pcoa

    IndTej/LN/beta/15.txt

    # 转换表型为fam文件
    Rscript scripts/gemma_fam.r -i IndTej/LN/beta/15.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam -t FALSE
    cp IndTej/LN/beta/15.txt.fam T2.fam
    # gemma计算单个SNP位点：不显著
	gemma -bfile T2 -k gemma/kinship.txt -c gemma/pca4.txt -lmm 4 -n 3 -o 3 -snps temp/OS10G0554200.snp # p值0.1，挑单个位点分析p值与全体一致
    # anova统计表型存在显著差异
	alpha_boxplot.sh -i IndTej/LN/beta/15.txt -m '"PC1","PC2","PC3"' -d IndTej/design.txt -A groupID -o IndTej/LN/
    
    # 改用chr10
    cp /mnt/zhou/chulab/miniCore/snp1.5x/Chr10.bed ./
    cp /mnt/zhou/chulab/miniCore/snp1.5x/Chr10.bim ./
    cp T2.fam Chr10.fam
    gemma -bfile Chr10 -k gemma/kinship.txt -c gemma/pca4.txt -lmm 4 -n 3 -o 3 -snps temp/OS10G0554200.snp # p值0.069，减少SNP确实增加了p值
    cat output/3.assoc.txt
    # 挑先单个位点
    plink --bfile /mnt/zhou/chulab/miniCore/snp1.5x/Chr10 --snp 10m21759092  -make-bed # 默认为 plink.*
    cp T2.fam plink.fam   
    gemma -bfile plink -k gemma/kinship.txt -c gemma/pca4.txt -lmm 4 -n 3 -o 3 -snps temp/OS10G0554200.snp # p值仍为0.069，减少SNP不再增加p值显著了
    cat output/3.assoc.txt

    # 模拟数据
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3} NR>FNR{print $0,a[$1]}' IndTej/design.txt plink.fam|sed 's/IND/1/;s/TEJ/100/;s/\t$/\tNA/'|less > Chr10.fam
    gemma -bfile Chr10 -k gemma/kinship.txt -c gemma/pca4.txt -lmm 4 -n 1 -o 1 -snps temp/OS10G0554200.snp # p值仍为0.069，减少SNP不再增加p值显著了
    # 使用阳性对照仍然不显著
    
    

## 2.2. GWAS关联分析

## 2.3. OTUs

## 2.4. PCoA
    
    # 发现LN下PCoA bray_cutis第二轴和unifrac第一轴IND和TEJ明显分析
    
    # 筛选97个籼粳道PCoA bray and unifrac
    Rscript scripts/otutab_sample_subset.r -i LN/pcoa_bray_curtis_15.txt -d IndTej/variety.id -o IndTej/LN/pcoa_bray.txt -t FALSE

    i=pcoa_bray
    
    # 转换表型为fam文件
	Rscript scripts/gemma_fam.r -i IndTej/LN/${i}.txt -d /mnt/zhou/chulab/miniCore/snp1.5x/T2.fam -t FALSE
	cp IndTej/LN/${i}.txt.fam T2.fam
    
	# 生成表型列名对应表
    # gemma表型文件fam不允许表头，单独保存关联结果后改名
    # 去掉前5行，取前99行，删除()"等符号，添加行号为gemma找列和输出文件名
	tail -n+6 IndTej/LN/${i}.txt.fam.header|head -n99|sed 's/[\(\)\"]//g'|awk '{print NR"\t"NR$0}'| sed "s/\t/\t${i}/"|sed '1 s/V6/Wax/' > IndTej/LN/${i}.list
    
    # 并行gemma关联基因型和表型
    # 使用paralle进行48线程并行
	parallel -j 48 "gemma -bfile T2 -k gemma/kinship.txt -c gemma/pca4.txt -lmm 4 -n {1} -o {1} -snps temp/OS10G0554200.snp" ::: `cut -f 1 IndTej/LN/${i}.list`
    
	# 批量按列编号改名
	awk 'BEGIN{OFS=FS="\t"}{system("mv output/"$1".assoc.txt output/"$2".assoc.txt");system("mv output/"$1".log.txt output/"$2".log.txt")}' <(cat IndTej/LN/${i}.list)
	
    # 筛选极显著的位点P<0.001
    # 注意并行时变量$符需要转义
	parallel -j 33 "cat <(head -n 1 output/{1}.assoc.txt) <(awk '\$14<0.001' output/{1}.assoc.txt) > output/{1}.assoc.sig" ::: `cut -f 2 IndTej/LN/${i}.list`
    
	# qqman批量绘制manhanttan plot
	parallel -j 33 "Rscript gemma/qqman.R -i output/{1}.assoc.sig -t 14" ::: `cut -f 2 IndTej/LN/${i}.list`
    
	# 添加SNP注释
    # parallel有太多$需要转义不方便，改用循环
	for RPM in `cut -f 2 IndTej/LN/${i}.list`; do
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6} NR>FNR{print $0,a[$2]}' <(sed '1 s/ID/rs/' gemma/snp.anno) output/${RPM}.assoc.sig |sed 's/\t$/\tNA/'|cut -f 1-3,14- > output/${RPM}.assoc.anno &
	done
    
    # 查看10m21759092位点结果
    grep '10m21759092' output/pcoa_bray* | less >temp/10m21759092_pcoa_bray_all.txt
    # 第14列，pscore: pc2为0.1，其它为0.3-0.98；13列lrt为0.09，12列p_wald
    grep '10m21759092' output/pcoa_bray* | less >temp/10m21759092_pcoa_bray_1.txt

    

## 2.5. 元素循环Faprotax

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
    wc -l gemma/output/otutab_freq10OTU_4.assoc.txt # 2287362
    wc -l output/fapro*.assoc.txt # 只有wax为2287400，其它为2153792，可能是我们有太少样品分析

    gemma -bfile T2 -k gemma/kinship.txt -c gemma/pca4.txt -notsnp -lmm 4 -n 2 -o 2 # minor allele frequency cutoff is not used, 2153792
    nohup gemma -bfile T2 -k gemma/kinship.txt -c gemma/pca4.txt -miss 0.001 -maf 0.001 -lmm 4 -n 2 -o 3 &bg # miss and maf to 0.001, 2153792
    gemma -bfile T2 -k gemma/kinship.txt -c gemma/pca4.txt -miss 0.000001 -maf 0.000001 -lmm 4 -n 2 -o 4 # miss and maf to 0.00001, 2153792
    gemma -bfile T2 -k gemma/kinship.txt -lmm 4 -n 2 -o 5 # 去掉covariates也无效，仍然选择丢了一些SNP
    # 指定十号染色体的SNP，看是否还扔
    grep 'OS10G0554200' gemma/snp.anno|cut -f 1 |grep '10m21759092' > temp/OS10G0554200.snp # 90，应该为 10m21759092 ，不是10m21761740

    time gemma -bfile T2 -k gemma/kinship.txt -c gemma/pca4.txt -miss 0.000001 -maf 0.000001 -lmm 4 -n 2 -o 6 -snps temp/OS10G0554200.snp -pmin 0 # 90个SNP只分析了47个，从61起跳过了8行。添加pmin-100或0也无影响


# 3. 附录

## 3.1. 检查Gemma丢弃SNP原因

    # 难道我这97个表型该位点SNP无差异吗？检查 10m21761740
    # 筛选10m21761740行和表头
    cat <(head -n 7 /mnt/bai/yongxin/rice/miniCore/mwas/genotype/Chr10.vcf|tail -n 1) <(grep '10m21761740' /mnt/bai/yongxin/rice/miniCore/mwas/genotype/Chr10.vcf) > temp/10m21761740.vcf

    # 替换0/0为G($4)，1/1为A($5)
    sed -i 's/#//;s/0\/0/G/g;s/1\/1/A/g' temp/10m21761740.vcf 

    # 在excel中转置，并选择所有基因型
    awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++){a[FNR,i]=$i}} END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' temp/10m21761740.vcf  | sed 's/\t$//' |  cut -f 2 -d '_' > temp/10m21761740.vcft
    # 添加亚种分类
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$2]=$13} NR>FNR {print $0,a[$1]}' /mnt/bai/yongxin/rice/miniCore/doc/minicore_list.txt temp/10m21761740.vcft | sort -k3,3 -k2,2 | sed 's/\t$/\tNA/' > temp/10m21761740.subspecie # 查看按基因型分组还不错
    grep -P 'TEJ$|IND$' temp/10m21761740.subspecie|less # 在粳籼稻间无差异

## 3.1. 检查 10m21759092 SNP基因
    
   # 筛选10m21759092行和表头
    cat <(head -n 7 /mnt/bai/yongxin/rice/miniCore/mwas/genotype/Chr10.vcf|tail -n 1) <(grep '10m21759092' /mnt/bai/yongxin/rice/miniCore/mwas/genotype/Chr10.vcf) > temp/10m21759092.vcf

    # 替换0/0为G($4)，1/1为A($5)
    sed -i 's/#//;s/0\/0/G/g;s/1\/1/A/g' temp/10m21759092.vcf 

    # 在excel中转置，并选择所有基因型
    awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++){a[FNR,i]=$i}} END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' temp/10m21759092.vcf  | sed 's/\t$//' |  cut -f 2 -d '_' > temp/10m21759092.vcft
    # 添加亚种分类
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$2]=$13} NR>FNR {print $0,a[$1]}' /mnt/bai/yongxin/rice/miniCore/doc/minicore_list.txt temp/10m21759092.vcft | sort -k3,3 -k2,2 | sed 's/\t$/\tNA/' > temp/10m21759092.subspecies # 查看按基因型分组还不错
    grep -P 'TEJ$|IND$' temp/10m21759092.subspecies|less # 在粳籼稻间无差异
    # 发现 L4105 N4126 两个IND中 G 应为A才完美，手动修改Chr10.fam后再测试

## 3.2 使用GAPI分析LN WAX/PC2和模拟数据
    
    cut -f 1,6,8,12 Chr10.fam|grep -v -P 'NA\tNA'|sed '1 i Taxa\tWax\tPC2\tStimulate'|less>~/rice/miniCore/mwas/GAPIT/simulate_1_2_PC2_WAX_IND_TEJ/10m21759092.subspecies
    # 运行GAPIT，需要输入文件有列名，无行名，列名第一列为样品。
    cd ~/rice/miniCore/mwas/GAPIT/simulate_1_2_PC2_WAX_IND_TEJ
    grep 10m21759092 GAPIT..*GWAS.Results.csv
    # 结果PC2仅为0.17，模拟数据1e-15次方，
    # 查看WAX的数据分布？

## 3.2 使用GAPI分析部分SNP

    cd ~/rice/miniCore/mwas/GAPIT
    mkdir -p genotype
    # 筛选/mnt/bai/xiaoning/past/software/tassel-5-standalone目录中sum_geno_11.hmp.txt文件
    # 筛选SNP中引起移码和非同义突变
    cut -f 3 ~/rice/miniCore/180319/gemma//snp.anno|sort|uniq -c # 共228万，其中6万+重点
    grep -P 'HIGH|MODERATE' ~/rice/miniCore/180319/gemma/snp.anno |sed '1 i rs' > genotype/snp.vip
    # 筛选hmp格式
    cat /mnt/bai/xiaoning/past/software/tassel-5-standalone/sum_geno_*.hmp.txt > temp
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' temp genotype/snp.vip > genotype/sum_vip_hmp.txt

    # simulate_1_2_PC2_WAX_IND_TEJ.r 使用6万个SNP与Wax, LN bray curtis PC2和1：100关联，几分钟可以出结果
    # 发现Wax有唯一结果过e-8, PC2仅有0.097， 模拟数据1:100有最高值e-15；改为1：10，sd=100%，p=0.0005, fdr=0.21; 1:2 sd=100%，p=0.3, fdr=1; 1:2 sd=50%, p=0.3, fdr=1;
    # 1:2, sd=0, p=1e-15, fdr=1e-12; 1:2, sd=10%, p=1e-7, fdr=1e-4; 1:2, sd=20%, p=1e-7, fdr=1e-4;
    grep 10m21759092 GAPIT..*GWAS.Results.csv

A_mean | B_mean | sd | p_value |fdr
-|-|-|-|-
1 | 100 | 0 | 1e-15 | 1e-12
1 | 2 | 0 | 1e-15 | 1e-12
1 | 2 | 10% | 1e-7 | 1e-4
1 | 2 | 20% | 3e-3 | 0.6
1 | 2 | 30% | 0.6 | 1
1 | 1.1 | 0 | 1e-15 | 1e-12
1 | 1.1 | 10% | 0.4 | 1
-0.1 | 0.1 | 100% | 4.8e-4 | 0.13

结论：GWAS分析中，对于discrete类型数据，如果与基因型一致，可以获得极高的p值，如100:0和1：0结果没区别；对于连续型数据，扰动到结果影响很大，如1：2时2倍差异，波动sd=10%即从e-15降为-7，20%只剩e-3，30%即趋近1，无任何显著。

**关联IND/TEJ与氨化作用**

simulate_1_2_PC2_WAX_IND_TEJ.r

# 文章图表

## Fig1. 实验设计+多样性
地图，实验设计展示- 秦媛

    cd ~/rice/IndTej
    cp /mnt/bai/qinyuan/rice/minicore/1.html figure

## Fig2. 差异菌

### Fig2.1 差异菌 Lefse展示

    cd /mnt/bai/yongxin/rice/miniCore/180319/IndTej/HN
    # 基于OTU表，实验设计和物种注释，生成门-属水平汇总
    Rscript ~/github/Amplicon/16S/R/lefse1.1.r # 过滤1/万 - 1/100，千1最适合出cladogram
    
    # lefse用docker(结果权限为bennyyu，需要当前和输出目录有o+w权限)
    
    # 1. 格式转换
    chmod o+w . temp
	docker run --rm -v `pwd`:/tmp/ --name=lefse biobakery/lefse format_input.py /tmp/result/lefse.txt /tmp/temp/lefse.in -c 1 -s 2 -u 3 -o 1000000
    # 2. 统计分析
	docker run --rm -v `pwd`:/tmp/ --name=lefse biobakery/lefse run_lefse.py /tmp/temp/lefse.in /tmp/temp/lefse.res -a 0.2 -w 0.2 -l 1
    # 无显著结果？默认anova, wilcoxon和LDR为0.05, 0.05, 和2 。也没结果？可以尝试只保留高丰度的菌(根本找不到差异)
    # 改用wilcoxon检验的结果，手动对应差异菌添加：把-正常为组类型、LDA值和P值；改为分组，AB/RF贡献和P值；但贡献太小改为相对丰度。添加两个继续绘图。
    
    # 3. 绘制柱状图
    docker run --rm -v `pwd`:/tmp/ --name=lefse biobakery/lefse plot_res.py /tmp/temp/lefse.res /tmp/temp/lefse.res.png 
    docker run --rm -v `pwd`:/tmp/ --name=lefse biobakery/lefse plot_res.py /tmp/temp/lefse.res /tmp/temp/lefse.res.pdf --format pdf 
    # 4. 绘制树状图
    docker run --rm -v `pwd`:/tmp/ --name=lefse biobakery/lefse plot_cladogram.py /tmp/temp/lefse.res /tmp/temp/lefse.cladogram.pdf --format pdf 
    # 单个OTU不好看，还需要仔细修修饰。

热图：http://bailab.genetics.ac.cn/report/16s/miniCore_IndTej_HN_wilcox_p01/result-otu.html#ind-vs-tej-6 基于wilcoxon检验的结果，没有热图聚类不够好。
修改脚本DAOTU_egr.r改为DAOTU_egr_sortHeatmap.r，绘制不聚类的结果heatmap_noclust.pdf，仍然规律不明显。还是柱状图合理。

### Fig2.2 机器学习找marker菌
    
    


### Fig2.2 氮利用上存在差异

    Faprotax结果表明硝酸盐/亚硝酸盐在高低氮下均有显著差别。详见HN/LN - faprotax中


    
## Fig3. 网络分析

菌与PCoA轴相关

网络 —— 

## Fig4. NRT品种与基因
    
    cd ~/rice/zjj.nitrogen/180116
    # 查看OTU_9在CP 2016年中箱线图
    alpha_boxplot.sh -i result/otu_table_norm.txt -d doc/design.txt -A groupID -B '"A50HnCp6","A56HnCp6","A58HnCp6","IR24HnCp6","V3703HnCp6","ZH11HnCp6","A50LnCp6","A56LnCp6","A58LnCp6","IR24LnCp6","V3703LnCp6","ZH11LnCp6"' -m '"OTU_9"' -t TRUE -o temp/CP6_ -n TRUE
    alpha_boxplot.sh -i faprotax/otu_table_tax.faprotax -d doc/design.txt -A groupID -B '"A50HnCp6","A56HnCp6","A58HnCp6","IR24HnCp6","V3703HnCp6","ZH11HnCp6","A50LnCp6","A56LnCp6","A58LnCp6","IR24LnCp6","V3703LnCp6","ZH11LnCp6"' -m '"nitrate_ammonification"' -t TRUE -o temp/CP6_ -n TRUE
    # 查看OTU_9在CP 2017年中箱线图
    alpha_boxplot.sh -i result/otu_table_norm.txt -d doc/design.txt -A groupID -B '"A50HnCp7","A56HnCp7","A58HnCp7","IR24HnCp7","nrtHnCp7","ZH11HnCp7","A50LnCp7","A56LnCp7","A58LnCp7","IR24LnCp7","nrtLnCp7","ZH11LnCp7"' -m '"OTU_9"' -t TRUE -o temp/CP7_
    alpha_boxplot.sh -i faprotax/otu_table_tax.faprotax -d doc/design.txt -A groupID -B '"A50HnCp7","A56HnCp7","A58HnCp7","IR24HnCp7","nrtHnCp7","ZH11HnCp7","A50LnCp7","A56LnCp7","A58LnCp7","IR24LnCp7","nrtLnCp7","ZH11LnCp7"' -m '"nitrate_ammonification"' -t TRUE -o temp/CP7_ -n TRUE
    # 查看OTU_9在Sz 2017年中箱线图
    alpha_boxplot.sh -i result/otu_table_norm.txt -d doc/design.txt -A groupID -B '"A50HnSz7","A56HnSz7","A58HnSz7","IR24HnSz7","nrtHnSz7","ZH11HnSz7","A50LnSz7","A56LnSz7","A58LnSz7","IR24LnSz7","nrtLnSz7","ZH11LnSz7"' -m '"OTU_9"' -t TRUE -o temp/SZ7_ 
    alpha_boxplot.sh -i faprotax/otu_table_tax.faprotax -d doc/design.txt -A groupID -B '"A50HnSz7","A56HnSz7","A58HnSz7","IR24HnSz7","nrtHnSz7","ZH11HnSz7","A50LnSz7","A56LnSz7","A58LnSz7","IR24LnSz7","nrtLnSz7","ZH11LnSz7"' -m '"nitrate_ammonification"' -t TRUE -o temp/SZ7_ -n TRUE


    cd /mnt/bai/yongxin/rice/zjj.nitrogen/CP2016
    alpha_boxplot.sh -i result/otu_table.txt -d doc/design.txt -A groupID -B '"A50HnCp5","A50LnCp5","A56HnCp5","A56LnCp5","A58HnCp5","A58LnCp5","D77HnCp5","D77LnCp5","D86HnCp5","D86LnCp5","IR24HnCp5","IR24LnCp5","soilHnCp5","soilLnCp5","V3703HnCp5","V3703LnCp5","Zh11HnCp5","Zh11LnCp5"' -m '"OTU_32"' -t TRUE -o temp/CP5_ -n TRUE




    
    




    
