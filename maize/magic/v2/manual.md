
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
	wd=maize/magic/v2/
	# 创建环境代码见~/github/Work/initial_project.sh

	## 准备实验设计

	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt
	sed -i 's/\t/\tL171121_/' doc/library.txt # time check SeqLibraryList.xlsx
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/L190220/"$2"_1.fq seq/"$1"_1.fq");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/L190220/"$2"_2.fq seq/"$1"_2.fq");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz

	# 标准多文库实验设计拆分，保存模板中design页为doc/design_raw.txt
	split_design.pl -i doc/design_raw.txt
	# 从其它处复制实验设计
	cp ~/ath/jt.HuangAC/batch3/doc/L*.txt doc/
	# 删除多余空格，windows换行符等
	sed -i 's/ //g;s/\r/\n/' doc/*.txt 
	head -n3 doc/L1.txt
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt | sed 's/#//g') <(cat doc/L* |grep -v '#'|grep -v -P '^SampleID\t') > doc/design.txt
	# 检查是否相等
	wc -l doc/design.txt
	cut -f 1 doc/design.txt|sort|uniq|wc -l

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


## 1.2. 按实验设计拆分文库为样品


	# 拆分样品
	head -n3 doc/L1.txt
	# 按L1/2/3...txt拆分library为samples
	# 输入为seq/L*.fq，输出为seq/sample/*.fq
	make library_split
	make library_split_stat
	# 统计结果见result/split有txt/pdf/png，推荐看png方便快速查看每张位图，主要缺失来自1，2号库
	# 查看样本量排序
	sort -k2,2n result/sample_split.log|less # 有40多个样本少于10000

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

# 获得了新lane的数据，重命名，并合并叶丝期数据
mv temp/filtered.fa temp/filtered_JT.fa # 拔节期jointing state, JT
ln /mnt/bai/yongxin/maize/magic/temp/filtered.fa temp/filtered_Sk.fa # 叶丝期 silking, SK
cat temp/filtered_JT.fa temp/filtered_Sk.fa > temp/filtered.fa
# 合并实验设计
mv doc/design.txt doc/design_JT.txt
ln /mnt/bai/yongxin/maize/magic/doc/design.txt doc/design_Sk.txt # 叶丝期 silking, SK
cat doc/design_Sk.txt <(tail -n+2 doc/design_JT.txt) > doc/design.txt 

    # (第一阶段结束，获得纯净扩增子序列temp/filtered.fa，可提供此文件从下面开始)


## 1.6. 序列去冗余

	# Remove redundancy, get unique reads
	# 输入为temp/filtered.fa，输出为temp/uniques.fa
	make fa_unqiue
    # mini=300, Unique 27761

## 1.7. 挑选OTU

	# Pick OTUs
	# unoise3速度比cluster_otus慢上百倍，更精细但结果也更多
	# 输入为temp/uniques.fa，输出为temp/Zotus.fa
	make otu_pick
    # 再次修改mini=500, Amplicon为9655

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

## 1.17. 有参考构建OTU表

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
	
	# Group compareing by edgeR or wilcox
	# 可选负二项分布，或wilcoxon秩和检验
	make DA_compare
	make DA_compare_tax
	make plot_volcano
	make plot_heatmap
	make plot_manhattan

# 3. 高级分析

## 3.9 培养菌注释

	# 默认为水稻，包括相似度、覆盖度、丰度和物种注释，请修改参数处菌库位置和注释文件
	make culture


# 4. 个性分析

## 4.1 基因型数据 

    mkdir -p snp/
    # 链接所有表型数据
    ln ../snp/cubic_* snp/
    # snp/cubic_*共5个文件，hmp2plink代表hmp格式转换为plink格式，maf0.02为min allel frequency 0.02，bed/bim/fam为一组，分别为二进制的样本与SNP表，snp列表，样本对应的表型列表；_1404_Kinship.txt为kinship相似度矩阵，最相近即本身为2；_PopStructure.txt为群体结构数据，即前10个主成分
    # snp和样本数量，样式分别为chr1.s_717和MG_49
    wc -l snp/cubic_1404_hmp2plink_maf0.02.bim snp/cubic_1404_hmp2plink_maf0.02.fam # 1404个样品，11,825,030个SNP
    # 1404个样本的列表
    cut -f 2 snp/cubic_1404_hmp2plink_maf0.02.fam -d ' ' > snp/cubic_1404.id
    # 注释文件
    cp -r ../snp/variant_effect/ snp/
    gunzip snp/variant_effect/cubic_1404_chr*.gz

	# 注：gapit在R中太慢跑不动；gemma一直报数据格式不对；tassel跑的慢，拆分成100份还要2天，对PCoA报错；只有emmax跑通且快，详见笔记“emmax GWAS分析流程”
	mkdir -p emmax
    time plink --bfile snp/cubic_1404_hmp2plink_maf0.02 \
        --recode 12 --output-missing-genotype 0 --transpose \
        --out emmax/snp 
    # 创建矩阵，v输出过程，
    time emmax-kin -v snp/snp 
    # 制作协变量，要求与表型一致，第三列为截距1，后面PCA前10轴
    awk 'NR==FNR{a[$1]=$0} NR>FNR {print $1,a[$2]}' snp/cubic_PopStructure.txt snp/snp.tfam |sed 's/\t/ 1 /' | sed 's/\t/ /g' > snp/snp.cov

## 4.2 微生物组数据

    # 备份原始OTU表和实验设计
    mv result/otutab.txt result/otutab.txt.190311
    mv doc/design.txt doc/design.txt.190311 
    # 制作新的实验设计
    head -n1 result/otutab.txt|sed 's/\t/\n/g'|awk '{print $1"\t"$1"\t"$1}'|sed 's/#OTUID\t#OTUID\t#OTUID/SampleID\tgroupID\tgroup/'|less > doc/design.txt

    mkdir -p 16S/
    # otutab与fam列表对应样本：与基因型对应的ID为ID_MG(Z列，26列)
    # 筛选新八叶期(JT开头)样品，样本名、组名和基因型对应ID
    cut -f 1,3,26 doc/design.txt | grep -v -P '^Sk' > doc/design_JT.id
    # 对应OTU样本名至基因型，生成以基因型为样本名的新OTU表
    script/otutab2genotype.Rmd
    # 重新计算物种和多样性
    usearch10 -otutab_stats result/otutab.txt -output result/otutab.stat
    # 删除了2个样品小于10000
    usearch11 -otutab_rare result/otutab.txt -sample_size 10000 -output result/otutab_norm.txt 
    make tax_sum # 重新计算物种组成
    make beta_calc # alpha, beta

    # Alpha
    mkdir -p pheno
    awk 'NR==FNR{a[$1]=$10} NR>FNR {print $1,$2,a[$2]}' result/alpha/index.txt emmax/snp.tfam | sed 's/ $/ NA/' > pheno/alpha_richness.txt
    awk 'NR==FNR{a[$1]=$4} NR>FNR {print $1,$2,a[$2]}' result/alpha/index.txt emmax/snp.tfam | sed 's/ $/ NA/' > pheno/alpha_chao1.txt
    awk 'NR==FNR{a[$1]=$13} NR>FNR {print $1,$2,a[$2]}' result/alpha/index.txt emmax/snp.tfam | sed 's/ $/ NA/' > pheno/alpha_shannon_e.txt
    
    # Beta 手动生成PC1-4 "bray_curtis","unweighted_unifrac","weighted_unifrac"
    beta_pcoa.sh -i `pwd`/result/beta/ -m '"weighted_unifrac"' \
        -d `pwd`/doc/design.txt  -A groupID -E FALSE -o `pwd`/result/beta/ -h 5 -w 8 
    # 制作表型文件
    for j in bray_curtis unweighted_unifrac weighted_unifrac; do
    for i in `seq 2 4`;do
    awk 'NR==FNR{a[$1]=$2} NR>FNR {print $1,$2,a[$2]}' <(cut -f 1,${i} result/beta/${j}14.txt) emmax/snp.tfam | sed 's/ $/ NA/' > pheno/beta_${j}${i}.txt
    done
    done
    # 编写alpha, beta表型列表pheno/list_alpha.txt

    # 制作 Top 30 OTU
    # 按丰度排序
    usearch11 -otutab_sortotus result/otutab_norm.txt -output result/otutab.sort.txt
    # 取前99
    head -n 100 result/otutab.sort.txt > result/otutab.sort.top.txt
    # 转置
    awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++){a[FNR,i]=$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]" "}print ""}}' result/otutab.sort.top.txt > result/tout.txt
    # 制作OTU与列号索引文件
    sed 's/^#OTU ID/ID/;s/ $//' <(head -n1 result/tout.txt) |tr ' ' '\n'|awk '{print NR"\t"$0}' > pheno/otu.id
    # 批量生成OTU表型文件， parellel中有$需转义
    parallel --xapply -j 3 \
    "awk 'NR==FNR{a[\$1]=\$2} NR>FNR {print \$1,\$2,a[\$2]}' <(cut -f 1,{1} -d ' ' result/tout.txt) emmax/snp.tfam | sed 's/ \$/ NA/' > pheno/{2}.txt " \
    ::: `tail -n+2 pheno/otu.id | cut -f 1` \
    ::: `tail -n+2 pheno/otu.id | cut -f 2`

## 4.3 基因型与微生物组关联

    # emmax 表型+基因型+协变量
    mkdir -p emmax
    for i in `cat pheno/list_alpha.txt`; do
    # i=alpha_shannon_e
    time emmax -t snp/snp \
    -p pheno/${i}.txt \
    -k snp/snp.aBN.kinf \
    -c snp/snp.cov \
    -o emmax/${i} 
    done
    # Single phenotyp 69m, 在biocloud上运行12个表型过夜，OTU在meta上运行

    # meta上批量运行
    mkdir -p emmax
    scp -r 192.168.0.110:~/maize/magic/v2/pheno .
    ln ../emmax/snp.* emmax/
    for i in `cut -f 2 pheno/otu.id|tail -n+2`; do
    # i=alpha_shannon_e
    time emmax -t emmax/snp \
        -p pheno/${i}.txt \
        -k emmax/snp.aBN.kinf \
        -c emmax/snp.cov \
        -o emmax/${i} 
    done
    dir=emmax
    # alpha_richness alpha_chao1 alpha_shannon_e beta_bc2 beta_bc3 beta_bc4
    for i in `cut -f 2 pheno/otu.id|tail -n+2`; do
    # i=alpha_shannon_e
    awk '$4<0.001' ${dir}/${i}.ps | awk '{print $1"\t"$1"\t"$4}' | sed 's/chr//;s/.s_/\t/' | awk '{print $1"\t"$3"\t"$2"\t"$4}' | sed '1 i CHR\tSNP\tBP\tP' > ${dir}/${i}.qqman
    qqman.R -i ${dir}/${i}.qqman
    done


    # 绘制曼哈顿图
    dir=emmax
    # alpha_richness alpha_chao1 alpha_shannon_e beta_bc2 beta_bc3 beta_bc4
    for i in `cat pheno/list_alpha.txt`; do
    # i=alpha_shannon_e
    awk '$4<0.001' ${dir}/${i}.ps | awk '{print $1"\t"$1"\t"$4}' | sed 's/chr//;s/.s_/\t/' | awk '{print $1"\t"$3"\t"$2"\t"$4}' | sed '1 i CHR\tSNP\tBP\tP' > ${dir}/${i}.qqman
    ~/github/Metagenome/denovo1/script/qqman.R -i ${dir}/${i}.qqman
    done



## 4.4 关联结果注释和分析

### 注释显著SNP对应SNP类型注释

	# 顺序注释1-10号染色体的单个表型大约2m，而使用合并10条染色体注释文件则需30m
#	for i in `seq 1 10`; do \
#	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{if ($4<0.00001){print $4,a[$2]}}' snp/variant_effect/cubic_1404_chr${i}_vepOut.txt \
#	emmax/alpha_shannon_e.qqman | grep 'chr' >> emmax/alpha_shannon_e.1e5.txt; done

	# 批量注释表型 cat pheno/list_alpha.txt     cut -f 2 pheno/otu.id|tail -n+2
	#for j in `cat pheno/list_alpha.txt`; do \
	for j in `cut -f 2 pheno/otu.id|tail -n+2`; do \
	for i in `seq 1 10`; do \
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{if ($4<0.00001){print $4,a[$2]}}' ../snp/variant_effect/cubic_1404_chr${i}_vepOut.txt \
	emmax/${j}.qqman | grep 'chr' >> emmax/1e5.${j}.txt; done; done




### 添加基因功能注释



## 附录

### 两批OTU对应编号

    mkdir -p idlink
    tail -n+2 pheno/otu.id | cut -f 2 > idlink/v2.id
	usearch10 -fastx_getseqs result/otu.fa -labels idlink/v2.id -fastaout idlink/v2.fa
    sed -i 's/>/>v2/' idlink/v2.fa

    tail -n+2 ../pheno/otu.id | cut -f 2 > idlink/v1.id
	usearch10 -fastx_getseqs ../result/otu.fa -labels idlink/v1.id -fastaout idlink/v1.fa
    makeblastdb -dbtype nucl -in idlink/v1.fa

	blastn -query idlink/v2.fa -db idlink/v1.fa -out idlink/v2-v1.txt -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads 9
    grep '100.00' idlink/v2-v1.txt | sed '1 i qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' | sed 's/ /\t/g' > idlink/v2-v1_perfectmatch.txt

    # gemma
    mkdir -p gemma
    cp result/alpha/index.txt.fam snp.fam
    # 样本顺序为字母顺序，不是fam中顺序
    i=alpha
	tail -n+6 result/alpha/index.txt.fam.header|head -n99|sed 's/[\(\)\"]//g'|awk '{print NR"\t"NR$0}'| sed "s/\t/\t${i}/" > gemma/${i}.list
    # 计算alpha多样性中第9列的richness关联
    gemma -bfile snp -k snp.kin -c snp.pca -lmm 4 -n 9 -o gemma/alpha9richness # le

    parallel -j 33 "gemma -bfile T2 -k gemma/kinship.txt -c snp.pca -lmm 4 -n {1} -o gemma/{1}" ::: `cut -f 1 gemma/${i}.list`

	parallel -j 33 "gemma -bfile snp/cubic_1404_hmp2plink_maf0.02 -k snp/cubic_1404_Kinship.txt -c snp/cubic_PopStructure.txt -lmm 4 -n {1} -o {1}" ::: `cut -f 1 ${dir}/HN/${i}.list`

    # gapit
    script/gapit.r # 使用内存超800G，手动终止

    # tassel
    # 转换ped格式为hmp
    time ~/bin/TASSEL5/run_pipeline.pl -plink -ped  cubic_1404_hmp2plink_maf0.02.ped -map cubic_1404_hmp2plink_maf0.02.map  -export cubic_1404_hmp_maf0.02 -exportType Hapmap -Xmx50g
    # 
    # 生成tassel格式
    mkdir -p tassel/16s
    cut -f 1,4,10,13 result/alpha/index.txt|sed '1 s/Sample/<Trait>/' > tassel/16s/alpha.txt
    cut -f 1,3  tassel/16s/alpha.txt >  tassel/16s/alpha_richness.txt
    cut -f 1,2  tassel/16s/alpha.txt >  tassel/16s/alpha_chao1.txt
    cut -f 1,4  tassel/16s/alpha.txt >  tassel/16s/alpha_shannon.txt
    cut -f 1,2,3 result/beta/bray_curtis14.txt|sed '1 s/Samples/<Trait>/' > tassel/16s/beta_bc.txt
    cut -f 1,2  tassel/16s/beta_bc.txt >  tassel/16s/beta_bc1.txt
    cut -f 1,3  tassel/16s/beta_bc.txt >  tassel/16s/beta_bc2.txt


# meta-gwas 在meta上分析GWAS
    
    [meta@meta:]scp -r yongxin@210.75.224.110:~/software/tassel-5-standalone/ ~/soft/tassel5
    screen -R gwas
    ssh yongxin@210.75.224.32
    cd ~/maize/magic/
    scp -r yongxin@210.75.224.110:~/maize/magic/snp ./
    scp -r yongxin@210.75.224.110:~/maize/magic/tassel ./
    scp -r yongxin@210.75.224.110:~/maize/magic/result ./


### 拆分snp批量运行tassel

    cd ~/maize/magic/snp
    # 11M行，输出1位数字，每个1.2M SNP拆为10份，一夜tassel只运行20%，改为100份
    split -a 2 -d -l 120000 cubic_1404_hmp_maf0.02.hmp.txt split_hmp
    # 提取ID
    ls split_hmp??|cut -c10- > split_list.txt
    # 添加表头
    header=`head -n1 split_hmp00`
    for i in `tail -n+2 split_list.txt`;do sed -i "1 i $header" split_hmp$i; done
    cd ..

/home/meta/soft/tassel5/run_pipeline.pl \
        -Xms10g -Xmx30g \
        -fork1 -h snp/split_hmp00 \
        -fork2 -r tassel/16s/${file}.txt \
        -fork3 -q snp/cubic_PopStructure.txt \
        -fork4 -k snp/cubic_1404_Kinship.txt \
        -combine5 -input1 -input2 -input3 -intersect \
        -combine6 -input5 -input4 \
        -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
        -mlmOutputFile $out/mlmOut00 \
        -export $out/mlmExport00 \
        -runfork1 -runfork2 -runfork3 -runfork4 \
        > $out/out00
# Data sets will not be joined because both phenotypes have values for PC1

    out=tassel/beta
    mkdir -p $out
    file=beta_bc1
    mkdir -p $out
    parallel --xapply -j 50 \
    "/home/meta/soft/tassel5/run_pipeline.pl \
        -Xms10g -Xmx30g \
        -fork1 -h snp/split_hmp{1} \
        -fork2 -r tassel/16s/${file}.txt \
        -fork3 -q snp/cubic_PopStructure.txt \
        -fork4 -k snp/cubic_1404_Kinship.txt \
        -combine5 -input1 -input2 -input3 -intersect \
        -combine6 -input5 -input4 \
        -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
        -mlmOutputFile $out/mlmOut{1} \
        -export $out/mlmExport{1} \
        -runfork1 -runfork2 -runfork3 -runfork4 \
        > $out/out{1}" \
        ::: `cat snp/split_list.txt`
    # 结果过滤与合并
    echo -e "CHR\tSNP\tBP\tP" > tassel/$file.txt
    for i in `cat snp/split_list.txt`;do
    awk 'BEGIN{FS=OFS="\t"} {if($7>=0 && $7<=1){print $3,$2,$4,$7}}' ${out}/mlmOut${i}_split_hmp${i}_+_${file}_+_cubic_PopStructure_stats.txt >> tassel/$file.txt
    echo -e "CHR\tSNP\tBP\tP" > tassel/$file.sig
    awk '$4<0.001' tassel/$file.txt >> tassel/$file.sig
    done
    # 绘图
    # cp /home/meta/soft/Metagenome/denovo1/script/qqman.R ./
    /usr/bin/Rscript qqman.R -i tassel/${file}.sig

    # 类似的结果文件 tassel/alpha/mlmOut00_split_hmp00_+_alpha_richness_+_cubic_PopStructure_stats.txt

    # 结果整理和可视化manhattan plot


## Reference

### Tassel参考流程

    TRAIT=$1
    ID=$2
    perl ~/software/tassel3-standalone/run_pipeline.pl \
        -Xms1g -Xmx3g \
        -fork1 -h ~/hjliu/cubic/geno/cubic_1404_maf0.02_ID$ID.hmp.gz \
        -fork2 -r ~/hjliu/cubic/AllTrait/${TRAIT}.txt \
        -fork3 -q ~/hjliu/cubic/GWAS/cubic_PopStructure.txt \
        -fork4 -k ~/hjliu/cubic/GWAS/cubic_1404_Kinship.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
        -mlmOutputFile mlmOut_cubic_${TRAIT}_ID${ID} \
        -export mlmExport_${TRAIT}_ID${ID} \
        -runfork1 -runfork2 -runfork3 -runfork4 \
        > out_${TRAIT}_ID${ID}


### emma参考流程
    
    i=$1
    ~/software/emmax-beta-07Mar2010/emmax -v -d 10 -t genotype/Plink/cubic_391_maf0.05.plk -p phenotype/${i}.rld.txt -k genotype/cubic_391_Kinship_emmax.txt -c genotype/cubic_PopStructure_391_peer_emmax.txt.bak -o GWAS/${i}.rld.peer
    gzip GWAS/${i}.rld.peer.ps

### 基因组v3.25    ensemble plants 上可以找之前版本下载 ，是用的  ensemble 的vep 注释的
阈值可以适当调整    1e-8 是非常严格的了        零星过阈值的不考虑
我们是20Kb范围内至少要有两个显著位点  才算
