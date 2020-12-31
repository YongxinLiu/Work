# 水稻菌保库

## 2020/1/3版本，新增近400个新菌保

    cut -f 1,2 rice_stock_200604_full.txt|grep -v -P '\t0$'|tail -n+2|sed 's/^/>/;s/\t/\n/'|less > sequence.fa
    makeblastdb -in sequence.fa -dbtype nucl
    # 物种注释
    usearch10 --sintax sequence.fa --db /mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb  -strand both \
      --tabbedout temp/otus.sintax --sintax_cutoff 0
    cut -f 1,2 temp/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
      > sequence.tax
    
    # 提取V5-V7区
	cutadapt -g AACMGGATTAGATACCCKG -e 0.2 sequence.fa -o temp/cut5.fa
	# 反向时先切反向引物，1193R ACGTCATCCCCACCTTCC，正常为1492R GGTTACCTTGTTACGACTT找不到，可选 1391R GACGGGCGGTGTGTRCA F TGYACACACCGCCCGTC
	cutadapt -a GGAAGGTGGGGATGACGT -e 0.2 temp/cut5.fa -o sequenceV5-7.fa
	# 长度分布，370-385，多数为379
	grep -v '>' sequenceV5-7.fa|awk '{print length($0)}'|sort -n|uniq -c
	# 冗余度
	grep -c '>' sequenceV5-7.fa # 1626
	grep -v '>' sequenceV5-7.fa | sort|uniq -c|wc -l # 806



## 2020/1/3版本，新增近400个新菌保

## 2019/7/31版本，手动修正一些菌
    cut -f 1,2 rice_stock_190731_full.txt|grep -v -P '\t0$'|tail -n+2|sed 's/^/>/;s/\t/\n/'|less>sequence.fa
    makeblastdb -in sequence.fa -dbtype nucl

## 2019/7/12版本

    # 菌保列表 rice_stock_190712.txt，补充序列，和物种注释
    # 先组装3端序列，没有再被充之前单端的，最后再物种注释

    # 原始序列 /mnt/bai/haoran/20190111_16sBlast
    # 分析位置 ~/rice/xianGeng/wet
    # 发表文章整理 /var/www/html/culture_collection/data ，完整的不全，只有和反向

    # 追加列表全长序列
    cd ~/culture/rice/stock/
    id=190712
    mkdir -p temp
    dos2unix /mnt/bai/yongxin/rice/xianGeng/wet/tableS11_full.txt rice_stock_${id}.txt
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' /mnt/bai/yongxin/rice/xianGeng/wet/tableS11_full.txt rice_stock_${id}.txt > temp/add_full
    # 追加单端
    sed 's/\r//' /var/www/html/culture_collection/data/16S_rice_culture_collection.fasta > temp/nbt.fa
    format_fasta_1line.pl -i temp/nbt.fa -o temp/nbt1.fa
    sed -i '1 i CNGB\tSeq1492R' temp/nbt1.fa.tsv
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$10]}' temp/nbt1.fa.tsv temp/add_full > rice_stock_${id}_full.txt

    # 生成菌保序列库
    # Excel修改编辑rice_stock_${id}_full.txt，删除重复的ID列，条件最后列Final =IF(Q2="",R2,Q2)
    dos2unix rice_stock_190712_full.txt
    # 只有1023条序列，需要手动检查缺失的原因
    cut -f 1,19 rice_stock_190712_full.txt|tail -n+2|grep -v -P '\t0$'|sed 's/^/>/;s/\t/\n/' | less > sequence.fa # wc -l

    # 以后手动添加序列，用cap3拼接


# 附录

## 1. 2019/1/30 水稻序列拼接 —— 参考代码

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
    # 编号转换
    sed 's/_//g' 16s_full_length_list.txt > 16s_full_length_list1.txt
    # 追加菌ID
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$1} NR>FNR{print a[$1],$0}' 16s_full_length_list1.txt 16s_full_length_list1.fa.tsv > 16s_full_length_list2.fa.tsv
    # 追加到菌保表TableS11，已经替换过啦！
    sed -i 's/_//g' 16s_full_length_list2.fa.tsv
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3} NR>FNR{print $0,a[$3]}' 16s_full_length_list2.fa.tsv tableS11_1492.txt > tableS11_full.txt
    # 有些ID没有对应上，统计数量
    cut -f 11 tableS11_full.txt|sort|uniq -c|less # 99个空值

## 2. 2020/1/3 水稻新增序列序列拼接

    # 设置环境变量
    wd=/mnt/bai/yongxin/culture/rice/wet
    batch=S200103
    cd ${wd}

    # 准备文件
    rm -rf ${batch}
    rm -rf temp &&  mkdir -p temp
    cp -r /mnt/bai/haoran/20191223_R2ABacteriaCulture16sSeq/Seq ./
    mv Seq ${batch}
    cp /mnt/bai/haoran/20191223_R2ABacteriaCulture16sSeq/R2A分菌纯化信息.xlsx ${batch}.xlsx 
    # 保存样本元数据为${batch}_metadata.tsv

    # 文件重命名，批量修改规则，无规则需手动修改
    rename 's/P\d_//' ${batch}/*.seq # 删除P+数字_
    rename 's/_\w\d\d_/_/' ${batch}/*.seq # 删除中间段孔号，会出重名问题？
    rename 's/_\w\d\d_/_/' ${batch}/*.seq 2>&1 | wc -l # 统计标准错误，有16个文件重复

    # 统计样本名与序列名是否一致
    # 样本数量及是否有重名
    tail -n+2 ${batch}_metadata.tsv | wc -l # 统计样本数，431
    tail -n+2 ${batch}_metadata.tsv | sort | uniq | wc -l # 统计样本是否存在重名，数值不等则存在重名，使用 uniq -u 输出重名
    # 制作测序结果列表
    ls $batch/ | wc -l # 统计序列文件数，1185
    ls $batch/ | grep 'seq' | cut -f 1 -d '_' | sort | uniq > ${batch}.ID # 保存序列ID
    wc -l ${batch}.ID # 431条
    # 比较两者是否一致，完全一致无输出，输出结果为其中单个文件中结果
    cat <(tail -n+2 ${batch}_metadata.tsv|cut -f 1) ${batch}.ID | sort | uniq -u
    dup=`cat <(tail -n+2 ${batch}_metadata.tsv|cut -f 1) ${batch}.ID | sort | uniq -u|tr '\n' '|'|sed 's/|$//'`
    grep -P ${dup} ${batch}_metadata.tsv ${batch}.ID
    # 有错误返回修改

    # 制作输入fa文件，将同一序列来源的多条测序结果合并，可以手动制作，也可以使用perl脚本自动批量
    file=`cat ${batch}.ID | head -n1`
    format_seq2fasta.pl -i "${batch}/${file}_*.seq" -o temp/${file}.fa
    cap3 temp/${file}.fa # 拼接结果在屏幕上，且保存了一系列文件，*.fa.cap.contigs即可
    sed -i "1 s/Contig1/${file}/" temp/${file}.fa.cap.contigs

    # 批处理拼接
    mkdir -p contigs
    for file in `cat ${batch}.ID`; do
        echo $file
        format_seq2fasta.pl -i "${batch}/${file}_*.seq" -o temp/${file}.fa
        cap3 temp/${file}.fa > temp/${file}.print
        # 改名
        sed -i "1 s/Contig1/${file}/" temp/${file}.fa.cap.contigs
        # 移至目录
        mv temp/${file}.fa.cap.contigs contigs/
        # 删除其它临时文件
        rm temp/${file}.*
    done

    # 添加各引物序列至ID
    cut -f 1 ${batch}_metadata.tsv|tail -n+2 > temp/ID.txt
    # 序列合并
    cat contigs/* > temp/${batch}.fa
    format_fasta_1line.pl -i temp/${batch}.fa -o temp/${batch}_1.fa # *.tsv表格用于补家
    # 追加菌ID，样本信息表保存为${batch}.txt
    # awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/${batch}_1.fa.tsv ${batch}_metadata.tsv > ${batch}_seq.txt
    wc -l ${batch}*
    # 发现样本列表中的529个，测序结果中有431个ID，拼接的只有377个

    dos2unix ${batch}/*
    for primer in 27F 515F 1492R; do
    # primer=27F
    format_seq2fasta.pl -i "${batch}/*_${primer}.seq" -o temp/${primer}.fa
    format_fasta_1line.pl -i temp/${primer}.fa -o temp/${primer}_1.fa # *.tsv表格用于补家
    sed -i "s/_${primer}//" temp/${primer}_1.fa.tsv
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/${primer}_1.fa.tsv temp/ID.txt > temp/${primer}.tsv
    done
    # 合并序列
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/${batch}_1.fa.tsv temp/ID.txt > temp/merged.tsv # 小于merge，存在ID错误的列表编号
    # 合并所有
    paste temp/merged.tsv temp/27F.tsv temp/515F.tsv temp/1492R.tsv | cut -f 1,2,4,6,8 | sed '1 i ID\tMerged\t27F\t515F\t1492R' > temp/${batch}_seq.tsv
    # 筛选代表序列Final，优先级为1. Merge $2，2. 515F $4，3. 1492R $5，4. 27F $3，全空补位
    awk 'BEGIN{OFS=FS="\t"}{if($2!=""){print $0,$2;} \
        else if($4!=""){print $0,$4;} \
        else if($5!=""){print $0,$5;} \
        else if($3!=""){print $0,$3;} \
        else{print $0"\t";}}' temp/${batch}_seq.tsv \
        | sed 's/Merged$/Final_Merge_515F_1492R_27F/' \
        > ${batch}_seq5.tsv
    # 提取最终序列作为库序列，cut1,5列；tail去标题；grep去空行；sort排序；sed转换为fa格式
    cut -f 1,6 ${batch}_seq5.tsv|tail -n+2|grep -v -P '\t$'|sort -h|sed 's/^/>/;s/\t/\n/'|less -S>${batch}.fa
    # 筛选有序列但表中没有的ID，即ID不对应的序列
    # cat ${batch}.ID <(grep '>' ${batch}.fa|sed 's/>//')|sort|uniq -u > ID.error

    # 与之前菌库比较
    blastn -query ${batch}.fa -db /mnt/bai/yongxin/culture/rice/stock/sequence.fa -out ${batch}.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 

    # 追加至菌库
    cat /mnt/bai/yongxin/culture/rice/stock/rice_stock_190731.fa ${batch}.fa > /mnt/bai/yongxin/culture/rice/stock/sequence.fa
    makeblastdb -in /mnt/bai/yongxin/culture/rice/stock/sequence.fa -dbtype nucl
