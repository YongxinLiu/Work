准备工作见云笔记：水稻miniCore与分蘖 ### 按亚种合并

	# 参考 ~/rice/miniCore/180319/manual.sh L165继续
    cd ~/rice/integrate16s/v2
	wd=LN2/result
	mkdir -p ${wd}
	# 按之前筛选的样品编号提取样品
	awk '$5=="L"' ../v2OTU/fig210928/data/miniCore_metadata.txt|cut -f1|less -S > ${wd}/sample_ids.txt
	usearch10 -otutab_sample_subset result/otutab.txt -labels ${wd}/sample_ids.txt -output ${wd}/otutab_raw.txt
	# 按组合并，usearch group直接求合
	awk '$5=="L"' ../v2OTU/fig210928/data/miniCore_metadata.txt|cut -f1,4|less -S > ${wd}/sample_ids_group.txt
	usearch10 -otutab_group ${wd}/otutab_raw.txt -labels ${wd}/sample_ids_group.txt -output ${wd}/otutab.txt
	# 统计，最小17656万，最大226363万
	usearch10 -otutab_stats ${wd}/otutab.txt -output ${wd}/otutab.txt.sum
    cat ${wd}/otutab.txt.sum
	# 采用标准流程生成alpha, beta多样性，按最小样本量抽平
    # 详见 /mnt/bai/yongxin/rice/integrate16s/v2/LN
    
# 从OTU表开始生成多样性结果
    
    cd ~/rice/integrate16s/v2/LN2
    # 使用initial_project.sh 创立新项目
    make init
	touch otutab_create
	touch otutab_filter
    biom convert -i result/otutab.txt -o result/otutab.biom --table-type="OTU table" --to-json
    # 统计OTU表
    biom summarize-table -i result/otutab.biom > result/otutab.biom.sum
    head -n 30 result/otutab.biom.sum # 18712
    # 检查表格
    csvtk -t stat result/otutab.txt
    # 统计每样总量
    # cat result/otutab.txt|datamash --header-in --header-out sum 2-203|datamash transpose|sort -k2,2n|less
	# 根据OTU表统计 cat result/otutab.txt.sum ，修改抽样量19318，只选择两个极小的5千和1.3万
	make otutab_norm
	# 设置sintax的cutoff
    ln -sf `pwd`/../result/otu.fa result/
	make tax_assign
	make tax_sum
	# 生成树时间长，可touch tree_make，再ln旧文件
	touch tree_make
    ln -sf `pwd`/../result/otu.tree result/otu.tree
	# 多样性计算
	make alpha_calc
	# Beta多样性计算
	make beta_calc

    # 准备数据
    cp result/tax/sum_g.txt ../fig/data/cultivar_LN_genus.txt



    # 添加实验设计继续
    cut -f 2- ~/rice/miniCore/doc/minicore_list.txt > doc/design.txt
    # 原始实验设计，只有分组
    mv doc/design.txt doc/design_raw.txt
    cp ~/rice/integrate16s/v2OTU/fig/191031/LNmean/metadata.txt doc/design.txt

	rm -rf alpha_boxplot 
	rm -rf DA_compare 
	make DA_compare # 绘制alpha、beta、taxonomy和差异OTU比较

    # 指定OTU进行比对，采用Measurable OTU
    compare_OTU_ref.sh -i `pwd`/result/otutab.txt -c `pwd`/doc/"tiller"/compare.txt -m "wilcox" \
        -p 0.05 -q 0.2 -F 1.5 -t 0.3 \
        -d `pwd`/doc/design.txt  -A tiller_cat -B `cat doc/"tiller"/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'` \
        -o `pwd`/result/compare/  -r ../miniCore/otutabK1.txt

	rm -f plot_volcano # 删除OTU差异比较可化标记
	make plot_manhattan # 绘制差异比较的火山图、热图、曼哈顿图
	make plot_venn # 绘制OTU差异共有/特有维恩图
	make DA_compare_tax # 高分类级差异比较，维恩图绘制，2为reads count负二项分布统计
	make rmd # 生成网页报告，必须依赖的只有alpha, beta, taxonomy

    # 绘制差异菌的箱线图
    mkdir -p result/compare/boxplot
    cat <(tail -n+2 ../fig/191031/LNmean/DA_CPM_FC.txt) <(tail -n+2 ../fig/191031/HNmean/DA_CPM_FC.txt) | cut -f 1 | sort | uniq > temp/OTUID.txt
    alpha_boxplot.sh -i result/otutab.txt -d doc/design.txt -A tiller_cat -B '"T6","T1"' -m `cat temp/OTUID.txt|awk '{print "\""$1"\""}'|tr "\n" ","|sed 's/,$//'` -t TRUE -o result/compare/boxplot/ -n TRUE -h 2 -w 3

# GWAS分析

## 1. 制作emmax表型文件

    mkdir -p pheno
    # alpha
    awk 'NR==FNR{a[$1]=$10} NR>FNR {print $1,$2,a[$2]}' result/alpha/index.txt ../emmax/snp.tfam | sed 's/ $/ NA/;s/ /\t/g' > pheno/alpha_richness
        # 对数2转换，log(x)/log(2)
    awk 'NR==FNR{a[$1]=$10} NR>FNR {print $1,$2,log(a[$2])/log(2)}' result/alpha/index.txt ../emmax/snp.tfam | sed 's/-inf$/NA/' > pheno/log2.alpha_richness
    # beta2-4(PC1-3), beta多样性有负数无法log2转换
    for i in `seq 2 4`;do
        awk 'NR==FNR{a[$1]=$2} NR>FNR {print $1,$2,a[$2]}' <(cut -f 1,${i} result/beta/bray_curtis14.txt) ../emmax/snp.tfam | sed 's/ $/ NA/;s/ /\t/g' > pheno/beta_bc${i}; done
    cut -f 1 ../fig/191031/otutab_mean.txt | tail -n+2 > pheno/otuid
    # 实现特征表筛选(可选)，并根据fam顺序制作表型文件和对数转换结果
    Rscript ~/github/Amplicon/16Sv2/script/gwas_otutab2emmax.R --input result/otutab_norm.txt \
        --list pheno/otuid --fam /mnt/bai/yongxin/rice/integrate16s/v2OTU/genotype/emmax/filtered.tfam \
        --output pheno/


## 2. GWAS关联

    # 必须有输出目录，否则 Segmentation fault (core dumped)
    # -c ../emmax/snp.cov \ 群体结构中有NA，无法分析，使用小宁gcta计算结果，但结果均为e-100次方明显错误
    awk '{print $1"\t"$0}' /mnt/bai/xiaoning/past/software/gcta_1.91.3beta/T2.sativaPca > ../emmax/snp.pca
    # 改用志文计算结果Plink计算
    time plink --vcf /mnt/zhou/zhiwen/rice_GWAS/rice_snp/T2.vcf --pca -out snp.pca
    rm snp.pca
    ln snp.pca.eigenvec snp.pca

    cd ~/rice/integrate16s/v2OTU/LN
    mkdir -p emmax
    i=alpha_richness
    snp=/mnt/bai/yongxin/rice/integrate16s/v2OTU/genotype/emmax/filtered
    time emmax -v -d 5 -t ${snp} \
        -p pheno/${i} \
        -k ${snp}.aBN.kinf \
        -c ${snp}.pca.eigenvec \
        -o emmax/${i}

    mkdir -p emmax_log2
    for i in `cat pheno/otuid|head -n1`; do
    echo $i
    time emmax -v -d 5 -t ${snp} \
        -p pheno/log2.${i} \
        -k ${snp}.aBN.kinf \
        -c ${snp}.pca.eigenvec \
        -o emmax_log2/${i}
    done

## 3. 结果筛选并可视化

    dir=emmax
    # 按P值过滤<1e-3，输出ID两次和P值，格式化ID为染色体
    i=alpha_richness
    awk '$4<1e-100' emmax/${i}.ps | awk '{print $1"\t"$1"\t"$4}' | sed 's/m/\t/' | awk '{print $1"\t"$3"\t"$2"\t"$4}' | sed '1 i CHR\tSNP\tBP\tP' > emmax/${i}.qqman
    # 结果全为上百次方，全显著

# alpha_richness alpha_chao1 alpha_shannon_e beta_bc2 beta_bc3 beta_bc4
for i in `cut -f 2 pheno/otu.id|tail -n+30`; do
# i=alpha_shannon_e
awk '$4<0.001' ${dir}/${i}.ps | awk '{print $1"\t"$1"\t"$4}' | sed 's/chr//;s/.s_/\t/' | awk '{print $1"\t"$3"\t"$2"\t"$4}' | sed '1 i CHR\tSNP\tBP\tP' > ${dir}/${i}.qqman
qqman.R -i ${dir}/${i}.qqman
done

## 4. QTL定位

    







# 1. 处理序列 Processing sequences

	# 0. 准备工作 Preparation

	## 0.1 准备流程配置文件

	# 设置工作目录
	wd=rice/miniCore
	# 创建环境代码见~/github/Work/initial_project.sh

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt
	sed -i 's/\t/\tL171121_/' doc/library.txt # time check SeqLibraryList.xlsx
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz

	# 标准多文库实验设计拆分，保存模板中design页为doc/design_raw.txt
	split_design.pl -i doc/design_raw.txt
	# 从其它处复制实验设计
	cp ~/ath/jt.HuangAC/batch3/doc/L*.txt doc/
	# 删除多余空格，windows换行符等(苹果用户勿用)
	sed -i 's/ //g;s/\r//' doc/*.txt 
	head -n3 doc/L1.txt
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
	# 输入为temp/filtered.fa，输出为temp/uniques.fa
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

	# OTU table filter samples and OTU
	# 推荐过滤低测序量<5000的样本，筛选大于1RPM的OTU
	make otutab_filter 


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

## 1.17. 有参考构建OTU表-非16S不可用

	# Reference based OTU table
	# otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
	make otutab_gg



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

