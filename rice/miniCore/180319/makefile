SHELL:=/bin/bash


# 16S扩增子流程配置文件，请根据项目具体情况修改 v1.4 2017/12/23
# Config file of 16S Amplicon pipeline, please modify according to experiment design. v1.4 2017/12/23



# 1. 测序文库参数 Standard pipeline parameter

# 工作目录 working directory
# pwd # 修改wd为新工作目录
wd=/mnt/bai/yongxin/rice/miniCore/180319
# make init # 建立流程所需目录

# 文库列表 library list
# ll clean_data/*.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|sort -u|tr "\n" " "
# tail -n+2 doc/library.txt |cut -f 1|tr "\n" " "
list=L1 L2

# 文库建库方法
# 类型，单左侧barcode选择barcode_single_end，单右端和双端均选择barcode_paired_stitched。barcode_paired_stitched for barcode in end and both
lib_type=barcode_single_end
# 正向barcode长度 forword barcode length
bc1=10
# 反向barcode长度 reverse barcode length
bc2=0
# Barcode类型，值为前两个barcode长度加和 barcode type, usually length equal barcode 1 plus barcode 2
bt=10
# 质量值类型，分33或64两种；determine_phred-score.pl clean_data/L1_1.fq.gz
phred=64
# OTU鉴定方法：默认 unoise3 ，适合小数据，样本没有多批次和背景差异不大； cluster_otus 大数据，背景差异较大
cluster=cluster_otus


# 2. 分组信息 Group

# 主要分组列名 Primary design group column
g1=variety
# 主要分组筛选，不筛选可为空
# 获取组信息 cat mwas/microbiome/genotype.list|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
# 主要分组总表：
g1_list='"A4001","A4002","A4003","A4004","A4005","A4006","A4007","A4008","A4009","A4010","Aro33","Aro34","Aro35","Aro36","Aro37","Aro38","B4011","B4012","B4013","B4014","B4015","B4016","B4017","B4018","B4019","B4020","C4021","C4022","C4023","C4024","C4025","C4026","C4027","C4028","C4029","C4030","D4031","D4032","D4033","D4034","D4035","D4036","D4037","D4038","D4039","D4040","E4041","E4042","E4043","E4044","E4045","E4046","E4047","E4048","E4049","E4050","F4051","F4052","F4053","F4054","F4055","F4056","F4057","F4058","F4059","F4060","G4061","G4062","G4063","G4064","G4065","G4066","G4067","G4068","G4069","G4070","H4071","H4072","H4073","H4074","H4075","H4076","H4077","H4078","H4079","H4080","J4081","J4082","J4083","J4084","J4085","J4086","J4087","J4088","J4089","J4090","K4091","K4092","K4093","K4094","K4096","K4097","K4098","K4099","K4100","L4101","L4103","L4104","L4105","L4106","L4107","L4108","L4109","L4110","M4111","M4112","M4113","M4114","M4115","M4116","M4117","M4118","M4119","N4121","N4122","N4123","N4124","N4125","N4126","N4127","N4128","N4129","N4130","P4131","P4132","P4133","P4134","P4135","P4136","P4137","P4138","P4139","P4140","Q4141","Q4142","Q4143","Q4144","Q4145","Q4146","Q4147","Q4148","Q4149","Q4150","R4151","R4152","R4153","R4154","R4155","R4156","R4157","R4158","R4159","R4160","S4161","S4162","S4163","S4165","S4166","S4167","S4168","S4169","S4170","T4171","T4172","T4173","T4174","T4175","T4176","T4177","T4178","T4179","W4193","W4194","W4195","W4196","W4197","W4198","W4199","W4200","X4201","X4202","X4203","X4204","X4205","X4206","X4207","X4208","X4209","X4210","Y4211","Y4212","Y4213","Y4214","Y4215","Y4216","Y4217"'
# 次要分组列名 Secondary design group column，没有先真batch
g2=Description
# 次要分组筛选，不筛选可为空
# 次要分组总表：tail -n+2 doc/design.txt|cut -f 7 |sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
# "H","L"
g2_list='"L"'

# 第三分组，可按此分组分别画Constrained PCoA，本示例是在不同土壤类型下画品种间差异
# 第三分组总表：tail -n+2 doc/design.txt|cut -f 8 |sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
g3=Description
g3_list='"H","L"'

# 合并主要和次要分组 default FALSE, if merge_group=TRUE, must group_order=FALSE
merge_group=FALSE
# 图例是否按主要组筛选样式排序，default TRUE ，但必须与merge_group保持相反
group_order=TRUE
# 成对比较，TRUE为默认比对group_compare.txt，而FALSE则自动两两比对
pair_compare=TRUE
# 仅批次按形状显示，默认按分组形状 Only Shape batch, default FALSE
batch=FALSE



# 3. 报告参数
# 报告输入信息目录 A50Cp A50Sz IR24Cp IR24Sz
sub=""
doc=doc/${sub}
# 报告输出目录
version=rice_TC_${sub}_1
# 报告输出是否精简版 report elite report, if FALSE report all figure and table, TRUE report frequently used figure
elite_report=TRUE
# 图片长宽和字体大小，7组以下用默认，7组以上改为8x5或更大； figure size, recommend 4x2.5, 5x3(default), 8x5, 16x10, text_size 6, 7(default), 8
width=8
height=5
text_size=7
# 图中显示taxonomy的数量，5，8(default)，10
tax_number=10
# 按丰度和分类单元过滤OTU OTU taxonomy and abundance filter parameter
# 丰度按万分之五过滤 # threshold of filter low abundance OTU，OTU太多计算过慢可改为万一
thre=0.0005
# 物种目前只去除叶绿体(k__Bacteria     p__Cyanobacteria        c__Chloroplast )和线粒体(k__Bacteria     p__Proteobacteria       c__Alphaproteobacteria  o__Rickettsiales        f__mitochondria)。不过滤p_xxx # filter some phylum	p__Actinobacteria,p__Bacteroidetes,p__Firmicutes,p__Proteobacteria p__Cyanobacteria,p__Chloroflexi c__Chloroplast,f__mitochondria
taxonomy=c__Chloroplast,f__mitochondria
# 显著性P值过滤 # threshold of filter differentially abundance OTU
pvalue=0.05
# 统计检验方式fdr, in edgeR have fdr or nonw
fdr="fdr"
logFC=1.3
ellipse=TRUE


# 4. 不常用参数
## 输入输出目录文件 Input and output directory and files
## 可变配置文件目录，包括6个文本文件，主要个性group_*.txt来设置比较组、维恩图和三元图；可在doc在建子目录，复制并编辑，修改此处目录
seq=clean_data
summary=${wd}/${doc}/summary.txt
library=${wd}/doc/library.txt
design=${wd}/doc/design.txt
compare=${wd}/${doc}/group_compare.txt
venn=${wd}/${doc}/group_venn.txt
tern=${wd}/${doc}/group_tern.txt
temp=temp
result=result
## 过滤OTU表结果目录 result based on filter OTU table
result_f=result_k1-c

## 日志文件，记录数据量整体过滤和OTU过滤 log file for basic statistics
log_reads=result/log_reads.txt
log_otus=result/log_otus.txt
log_usearch=result/log_usearch.txt

## 过滤序列质量>19为99%准确度 base quality, accurate > 99%; 29 means 99.9%
quality=19
## 过滤N的数字，默认0，在引物中有N时最多3
N=3
## 16S primers F799 and R1192 
# 5` primer used for 16S
primer5=AACMGGATTAGATACCCKG
# 3` primer used for 16S, must reverse compliment
primer3=GGAAGGTGGGGATGACGT 
# 引物匹配错误率，建议0.15，可调到最高0.25
er=0.25

## 保留扩增子的最小长度，细菌799-1192用300，真菌ITS1用220 # min length, recommend 300 for bacterial 16S and 220 for ITS
min_len=300
## 最小样本量 # sample min count, filter samples less than thre_count
thre_count=5000
## 用于聚类的序列最低8，可选样品重复数量(replications numbers)、3xrep，甚至 1/1M # min count of unique reads, reconmend 1/1000000?
minuniquesize=15
## 比对序列的相似度阈值，默认0.97 # similarity of cluster OTU
sim=0.97
## 最大使用计算机线程数，主要给clustero多序列比对使用 # threads number used: 32
p=48
## 用于筛选绘制圈图ggtree和igraphlan的OTU # filter OTU percentage > 0.5% for draw taxonomy and phylogenetic tree, 0.1% about 150 OTU is too much to show
tax_per=0.005
## OTU物种注释的方法 # rdp, blast, rtax, mothur, uclust, sortmerna , default=uclust, recommend rdp is better
method=rdp
## Alpha多样性分析的抽样数据 # alpha rarefaction count, recoomend 10000, at least 5000
rarefaction=10000
## OTU表达丰度样式，默认为百分比percentage，可选css, rpm # add css or percentage mean add normlized value, default no sample data
otu_stat_style=percentage

# 数据库 database；目前RDP虽然全，但是出现分类明显不准确且缺少Firmicute问题？暂时用gg13.8
## silva 128 99%, 492M uchime2_ref 中建议不去，如果去用最大的数据库，替换原29M rdp为新492M silva 128 99%
rdp=/mnt/bai/public/ref/silva/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/99/99_otus_16S.fasta
## 绿色基因细菌16S数据库多序列比对文件，用于建立多序列比对和进化树 97% 763M, 99% 1.5G, SILVA 128 99% 13G, RDP 11.5仅细菌比对文件有78G，过G计算会很困难
gg_align=/mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta
## RDP 11.5 16S细菌和古菌序列: 注释比例高，但出现分类时间过长、进化树聚类不一致、缺少Firmicute等问题？改为Roboert处下载的RDP trainset 16
gg_seq=/mnt/bai/public/ref/rdp/rdp_16s_v16_sp_qiime.fa
## RDP 11.5 16S细菌和古菌物种注释信息
gg_tax=/mnt/bai/public/ref/rdp/rdp_16s_v16_sp_qiime.tax.full

### RDP数据库用于去除嵌合体 rdp gold database, for remove chimera
#rdp=/mnt/bai/public/ref/rdp_gold.fa
### 绿色基因细菌16S数据库多序列比对文件，用于建立多序列比对和进化树  greengene bacterial 16S database
#gg_align=/mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta 
## 绿色基因细菌16S数据库 greengene bacterial 16S database，虽然旧、不完整，但快、准。
#gg_seq=/mnt/bai/public/ref/gg_13_8_otus/rep_set/97_otus.fasta
### 绿色基因细菌16S数据库物种注释信息 greengene bacterial 16S database
#gg_tax=/mnt/bai/public/ref/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt


# culture_graphlan
# 筛选指定组样品并鉴定培养比例,且样品也要对应
type=""
filter=filter_${type}_k1
thre2=0.0005
otu_table=${wd}/${result_f}/otu_table.txt
cluture_db=/mnt/bai/yongxin/culture/rice/result/${type}culture_select.fa


# ## 不同列计算PCoA, compartment(第一轴分开), genotype(1/2/3/4分不开), site(第三轴分开), day
group_color=day
time_course=TRUE



# OTU过滤参数
# 过滤结果输出目录
of=otu_filter

# 载入16S扩增子主流程，修改请保证向前兼容，最后只添加新分枝流程
# Loading amplicon 16S main pipeline. Don't change the following file unless necessary.
include /mnt/bai/yongxin/ref/amplicon/16s/makefile_16s

