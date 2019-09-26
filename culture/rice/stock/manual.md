# 水稻菌保库

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

## 2019/7/31版本，手动修正一些菌
    cut -f 1,2 rice_stock_190731_full.txt|grep -v -P '\t0$'|tail -n+2|sed 's/^/>/;s/\t/\n/'|less>sequence.fa
    makeblastdb -in sequence.fa -dbtype nucl

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
