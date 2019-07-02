
	# 快速分析 Quick Start(所需文件准备好)
	find . -name "*" -type f -size 0c | xargs -n 1 rm -f # 清理零字节文件，用于从头重新分析项目清空makefile点位文件
	make init
	make validate_mapping
	make sample_merge
	make extract_barcodes
	make split_libraries # 样本拆分
	make split_libraries_stat
	make fq_trim
	make fa_unqiue
	make otu_pick
	make chimera_ref
	make host_rm
	make otutab_create
	make otutab_filter
	make otutab_norm
	make tax_assign
	make tax_sum
	make tree_make
	make identify_isolate

# 与实验菌建立联系
## 筛选纯菌并建索引
tail -n+2 result/culture_select.xls | cut -f 1,2|sed 's/^/OTU_/g;s/;/\t/g'|less>result/culture_select.tax
filter_fasta.py -f result/otu.fa -o result/culture_select.fa -s result/culture_select.tax
sed -i 's/OTU/COTU/' result/culture_select.fa
makeblastdb -dbtype nucl -in result/culture_select.fa





# 每个库分别找菌，分为A50/IR24 H/L氮
cat doc/L1.txt <(tail -n+2 doc/L2.txt) > doc/A50L.txt # 按分类合并
cat doc/L4.txt <(tail -n+2 doc/L5.txt) > doc/IR24H.txt
ln doc/L3.txt doc/A50H.txt
ln doc/L6.txt doc/IR24L.txt
temp=temp
result=result
for lib in A50L A50H IR24L IR24H; do
filter_samples_from_otu_table.py -i ${result}/otu_table.biom -o ${result}/${lib}_otu_table.biom --sample_id_fp doc/${lib}.txt
biom convert -i ${result}/${lib}_otu_table.biom -o ${result}/${lib}_otu_table.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU //g' ${result}/${lib}_otu_table.txt
identify_isolate.sh -f ${lib}_otu_table.txt -o ${lib}
tail -n+2 result/${lib}culture_select.xls | cut -f 1,2|sed 's/^/OTU_/g;s/;/\t/g'|less>result/${lib}culture_select.tax
filter_fasta.py -f result/rep_seqs.fa -o result/${lib}culture_select.fa -s result/${lib}culture_select.tax
makeblastdb -dbtype nucl -in result/${lib}culture_select.fa
done


# 1. 处理序列 Processing sequences

	# 0. 准备工作 Preparation

	## 0.1 准备流程配置文件

	# 创建环境代码见~/github/Work/initial_project.sh
	# 设置工作目录
	wd=culture/medicago/190626
	## 准备实验设计
	cd ~/$wd
	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt，参考SeqLibraryList.xlsx，与/mnt/bai/yongxin/seq/amplicon目录中文件对应
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz


	## 写mappingfile, s为物种，p为板数；多个library需要多个mappingfile
	# 可指定品种、部位和培养基类型
	# 单个文库
	write_mappingfile_culture2.pl -o doc/L1.txt -s medicago -L L1 -v A17 -c Rhizosphere -m R2A -B 1 -p 1
	# 批量相同属性文库
	for i in `seq 1 10`; do write_mappingfile_culture2.pl -o doc/L${i}.txt -s medicago -L L${i} -v A17 -c Rhizosphere -m R2A -B 1 -p 48; done
	# 按Library信息批量生成
	# awk '{if(NR>2){system("echo "$1" "$6" "$7" "$8)}}' doc/library.txt
	awk '{if(NR>2){system("write_mappingfile_culture2.pl -o doc/"$1".txt -s medicago -L "$1" -v "$6" -c "$7" -m "$8" -B "$9" -p "$10)}}' doc/library.txt
	# L9, L10手动修改个性化数据，在Excel中手动修改 


	# 删除多余空格，windows换行符等(MAC用户禁用)
	sed -i 's/ //g;s/\r//' doc/*.txt
	# 查看数据格式、列名，了解和检查实验设计
	head -n3 doc/L1.txt
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt | sed 's/#//g') <(cat doc/L* |grep -v '#'|grep -v -P '^SampleID\t') > doc/design.txt
	# 检查是否相等
	wc -l doc/design.txt
	cut -f 1 doc/design.txt|sort|uniq|wc -l
	# 查看冗余的列(仅上方不等时使用)
	cut -f 1 doc/design.txt|sort|uniq -c| less

## 1.2. 按实验设计拆分文库为样品


	# 拆分样品
	head -n3 doc/L1.txt
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
