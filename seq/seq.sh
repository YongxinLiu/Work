# amplicon 诺禾致源测序

## novagene汇总

    # novagene 下载客户端安装和使用 2019/9/11
    # 在bailab远程桌面上，从报告系统下载客户端至download目录并解压，添加可执行权限和移动至环境变量
    chmod +x Downloads/linuxnd/linuxnd
    sudo mv Downloads/linuxnd/linuxnd /usr/local/bin/
    # 登陆和查看
    linuxnd login -u X101SC19070595-Z01-F003 -p hn5hypc7
    linuxnd list # 查看用户根目录，确定是否登陆成功
    linuxnd list oss://gxxkeyan@126.com/H101SC19070595/KY_kehu_JK/X101SC19070595-Z01/X101SC19070595-Z01-F003/ # 查看邮件路径文件
    linuxnd cp -d oss://gxxkeyan@126.com/H101SC19070595/KY_kehu_JK/X101SC19070595-Z01/X101SC19070595-Z01-F003/2.cleandata/ ~/seq/L${id}/ # 查看指定路径下的2.cleandata

    # 根据反向标准获得完整信息
    # 准备Index列表拆分，需要IndexRC列，可由ID或Index列检索到
	# 情况1：如保存ID列为index.ID，获取index和IndexRC，并补充到实验设计
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' \
        /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt index.ID > index.txt
    # 情况2：如保存Index序列，保存为index，获取indexID和IndexRC，并补充到实验设计
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$2]=$0} NR>FNR {print a[$1]}' \
        /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt index > index.txt
    cat index.txt
    # 保存至库文件和登记 SeqLibraryList.xlsx
    
    # 批量整理和改名
    find 2.cleandata/ .gz|grep '_1.clean.fq.gz'|sort|sed s/1.clean.fq.gz// # 填入第一列，制作成表格rename.txt
    awk 'BEGIN{OFS=FS="\t"}{system("ln "$1"1.clean.fq.gz "$2"_1.fq.gz")}' rename.txt
    awk 'BEGIN{OFS=FS="\t"}{system("ln "$1"2.clean.fq.gz "$2"_2.fq.gz")}' rename.txt

    # 质控和报告汇总
    fastqc -t 48 L*.fq.gz
    multiqc ./ # 详见multiqc_report.html
    # 提取各样品数据量，不同批数据量的列会有变化，可能要更改列的数值5、6等
    l=`head -n1 multiqc_data/multiqc_fastqc.txt|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep 'Total Sequences'|cut -f 1`
    cut -f 1,${l} multiqc_data/multiqc_fastqc.txt | sed 's/_.\t/\t/' | uniq > datasize.txt
    cat datasize.txt
    # 分双端统计md5值
    md5sum L*_1.fq.gz > /tmp/md5sum1.txt
    md5sum L*_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt
    # 汇总至amplicon目录
    ls L*.gz
    ls L*.gz|wc
    ln L*.gz ../amplicon/

## 190921 郭晓璇 8个库
    id=190921
    mkdir ~/seq/L${id}/
    cd ~/seq/L${id}/
    linuxnd login -u X101SC19070595-Z01-F005 -p 4e4pjpug 
    linuxnd cp -d oss://gxxkeyan@126.com/H101SC19070595/KY_kehu_JK/X101SC19070595-Z01/X101SC19070595-Z01-F005/2.cleandata/ ~/seq/L${id}/ # 查看指定路径下的2.cleandata

## 190909 宝原7
    id=190909
    mkdir ~/seq/L${id}/
    cd ~/seq/L${id}/
    linuxnd login -u X101SC19070595-Z01-F004 -p 3gy397n0
    linuxnd cp -d oss://gxxkeyan@126.com/H101SC19070595/KY_kehu_JK/X101SC19070595-Z01/X101SC19070595-Z01-F004/2.cleandata/ ~/seq/L${id}/ # 查看指定路径下的2.cleandata


## 190904 王超玉米时间序列5个库
    id=190904
    mkdir ~/seq/L${id}/
    cd ~/seq/L${id}/
    linuxnd login -u X101SC19070595-Z01-F003 -p hn5hypc7
    linuxnd cp -d oss://gxxkeyan@126.com/H101SC19070595/KY_kehu_JK/X101SC19070595-Z01/X101SC19070595-Z01-F003/2.cleandata/ ~/seq/L${id}/ # 查看指定路径下的2.cleandata

## 160914 16S+ITS
    cd ~/seq/160914xx.rice.TF.16S.ITS/clean_data
    zcat GB_HNVN7BCXX_L1_1.clean.fq.gz GB_H2NLWBCXY_L1_1.clean.fq.gz | pigz -p 30 > ~/seq/novagene/L160914_GB_1.fq.gz &
    zcat GB_HNVN7BCXX_L1_2.clean.fq.gz GB_H2NLWBCXY_L1_2.clean.fq.gz | pigz -p 30 > ~/seq/novagene/L160914_GB_2.fq.gz &
    zcat GF_H2NLWBCXY_L1_1.clean.fq.gz GF_HNVN7BCXX_L1_1.clean.fq.gz | pigz -p 30 > ~/seq/novagene/L160914_GF_1.fq.gz &
    zcat GF_H2NLWBCXY_L1_2.clean.fq.gz GF_HNVN7BCXX_L1_2.clean.fq.gz | pigz -p 30 > ~/seq/novagene/L160914_GF_2.fq.gz &

## 170122 xx.absolute.quantification
    cd ~/seq/170122xx.absolute.quantification/clean_data
    rename 's/clean.//' *gz
    rename 's/AQ44_HCHTTBCXY_L1/L170122_AQ44/' *.gz
    rename 's/AQ44-2_H552MBCXY_L2/L170122_AQ44-2/' *.gz
    ln AQ44_HCHTTBCXY_L1_1.clean.fq.gz ~/seq/novagene/L170122_AQ44_1.fq.gz

## 170302
    # 删除共有和冗余名字
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    # 提取文件列表
    ls *_1.fq.gz|cut -f 1 -d '_' > index.txt
    # 修改添加新名另存为library.txt
    dos2unix index.txt
    # 批量链接
    awk 'BEGIN{OFS=FS="\t"}{system("ln "$1"_1.fq.gz ~/seq/novagene/"$2"_1.fq.gz")}' index.txt
    awk 'BEGIN{OFS=FS="\t"}{system("ln "$1"_2.fq.gz ~/seq/novagene/"$2"_2.fq.gz")}' index.txt

## 170327
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^/L170327_/' *.gz
    ln *.gz ~/seq/novagene/

## 170424
    cd ~/seq/170424.ath/clean_data
    rename 's/clean.//' *gz
    rename 's/JT201703283101_HCJ2YBCXY_L2/L170424_T1/' *.gz
    rename 's/JT201703283202_HCJ2YBCXY/L170424/' *.gz
    ln *.gz ~/seq/novagene/

## 170425
    cd ~/seq/170425.clover/clean_data
    rename 's/clean.//' *gz
    awk 'BEGIN{OFS=FS="\t"}{system("ln "$1"_1.fq.gz ~/seq/novagene/"$2"_1.fq.gz")}' index.txt
    awk 'BEGIN{OFS=FS="\t"}{system("ln "$1"_2.fq.gz ~/seq/novagene/"$2"_2.fq.gz")}' index.txt

## 170623
    cd ~/seq/170623rice.nitrogen/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/yoyo-/L170623_/' *.gz
    ln *.gz ~/seq/novagene/

## 170703.GuoXX
    cd ~/seq/170703.GuoXX.CTK/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^/L170703_/' *.gz
    ln *.gz ~/seq/novagene/

## 170703.jt.ath.terpene
    cd ~/seq/170703.jt.ath.terpene/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    ls *_1.fq.gz|cut -f 1 -d '_' > index.txt
    awk 'BEGIN{OFS=FS="\t"}{system("ln "$1"_1.fq.gz ~/seq/novagene/"$2"_1.fq.gz")}' index.txt
    awk 'BEGIN{OFS=FS="\t"}{system("ln "$1"_2.fq.gz ~/seq/novagene/"$2"_2.fq.gz")}' index.txt

## 170707.ZhouJM.GaoCL
    cd ~/seq/170707.ZhouJM.GaoCL/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/JT201706032/L170707_/' *gz
    ln *.gz ~/seq/novagene/

## 170718.JT.terpene3
    cd ~/seq/170718.JT.terpene3/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/JT201706133234-/L170718_/' *gz
    ln *.gz ~/seq/novagene/

## 170813JY.rice.epi
    cd ~/seq/170813JY.rice.epi/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/YoYo-2-/L170813_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|cut -c1-14|uniq

## 170814.ZN.rice.sd1
    cd ~/seq/170814.ZN.rice.sd1/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/ZN-/L170814_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|cut -c1-14|uniq

## 170816.JT.wheat.NP
    cd ~/seq/170816.JT.wheat.NP/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/JT201706244047-/L170816_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|cut -c1-14|uniq

## 170823JY.clover
    cd ~/seq/170823JY.clover/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/Med-/L170823_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|cut -c1-14|uniq

## 170831.XX.AQ
    cd ~/seq/170831.XX.AQ/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^/L170831_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|sed 's/_[12].fq.gz//'|uniq

## 170919.XX_AQ
    cd ~/seq/170919.XX_AQ/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^/L170919_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|sed 's/_[12].fq.gz//'|uniq

## 170920.culture_med
    cd ~/seq/170920.culture_med/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^/L170920_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|sed 's/_[12].fq.gz//'|uniq

## 170926.culture_rice
    cd ~/seq/170926.culture_rice/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^/L170926_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|sed 's/_[12].fq.gz//'|uniq

## 170928.culture_ath
    cd ~/seq/170928.culture_ath/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^JT201708100/L170928_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|sed 's/_[12].fq.gz//'|uniq

## 170929.culture_wheat
    cd ~/seq/170929.culture_wheat/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^/L170929_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|sed 's/_[12].fq.gz//'|uniq

## 171013-20.GXX.AQ
    cd ~/seq/171013-20.GXX.AQ/clean_data
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz
    rename 's/^/L171013_/' *gz
    ln *.gz ~/seq/novagene/
    ls *.gz|sed 's/_[12].fq.gz//'|uniq

## 190724
    cd ~/seq/L190724/
    mv 2.cleandata/DYM1902_FKDL190748303-1a/* ./
    md5sum -c MD5_DYM1902_FKDL190748303-1a.txt
    rename 's/DYM1902_FKDL190748303-1a/L190724/;s/.clean//' *.gz
    fastqc -t 9 *.gz
    ln *.gz ~/seq/amplicon/



# amplicon 华大基因扩增子包lane测序数据

## 171018 lane1 钱景美苜蓿真菌
    wd=L171018 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 171121 lane2
    wd=L171121 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 171122 lane3
    wd=L171122 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 171212 lane4
    wd=L171212 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 171213 lane5
    wd=L171213 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 171221 华大lane6 重测
    wd=L171221 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 171220 华大lane7
    wd=L171220 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 180108 华大lane8
    
    wd=L180108 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 180210.lane9.ath3T

    wd=L180210 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 180301 lane10

    wd=L180301 && cd ~/seq/$wd
    # 转换64为33、准备index、拆lane、统计

## 180528 lane11 定量
    
    wd=L180528
    cd ~/seq/$wd
    cd Clean/CWHPEPI00001823/
    # 确定为64位格式，先转换，详见附见2
    rename 's/FCH7GG7BCX2_L2_CWHPEPI00001823/lane_64/' *.gz
    # 有预拆好的Library移动并改名即可
    mv Clean/CWHPEPI00001823/seq/*.fq ./
    rename "s/^/${wd}_/" *.fq
    # 质控和压缩上传

## 180719 lane12 nrt1.1a

    wd=L180719
    cd ~/seq/$wd
    # 附录：查看质量报告和评估文件完整性
    # 此lane数据被华大拆分为5个文件夹？分别处理
    # 有三个文件夹中只包括一个样品，直接移动并改名即可
    mv Clean/CW541/FCH7F2JBCX2_L1_CWHPE18070019_* ./
    rename 's/FCH7F2JBCX2_L1_CWHPE18070019/L180719_CACTCA/' FCH7F2JBCX2_L1_CWHPE18070019_*.gz
    mv Clean/TARDF16/FCH7F2JBCX2_L1_CWHPE18070017-A_* ./
    rename 's/FCH7F2JBCX2_L1_CWHPE18070017-A/L180719_CCGTCC/' FCH7F2JBCX2_L1_CWHPE18070017-A_*.gz
    mv Clean/TARDB28/FCH7F2JBCX2_L1_CWHPE18070018_* ./
    rename 's/FCH7F2JBCX2_L1_CWHPE18070018/L180719_CAAAAG/' FCH7F2JBCX2_L1_CWHPE18070018_*.gz
	# TOM文件夹包括两个library吗？检查只有 CAGGCG
    mv Clean/TOM/FCH7F2JBCX2_L1_CWHPE18070020_* ./
    rename 's/FCH7F2JBCX2_L1_CWHPE18070020/L180719_CAGGCG/' FCH7F2JBCX2_L1_CWHPE18070020_*.gz
    # AC目录正好多一个TCTGAG，逐个移动并改Index为反向互补，再移至根目录
    awk 'BEGIN{OFS=FS="\t"}{system("mv Clean/AC/"$2"_1.fq.gz Clean/L180719_"$3"_1.fq.gz ")}' index.txt
    awk 'BEGIN{OFS=FS="\t"}{system("mv Clean/AC/"$2"_2.fq.gz Clean/L180719_"$3"_2.fq.gz ")}' index.txt
    mv Clean/*.gz ./
    # 附录：质控与压缩

## 180816 lane13 sd1
    
    # 设置批次编号
    wd=L180816
    cd ~/seq/$wd
    # 附录：查看质量报告和评估文件完整性
    # 统一lane文件为当前目录的Lane_1/2.fq.gz
    mv Clean/D20170718/FCHFLH2BCX2_L2_CWHPE18070298_* ./
    rename 's/FCHFLH2BCX2_L2_CWHPE18070298/lane/' *.fq.gz
	# 附录：检查数据质量，如果64跳到附录2
    # 附录：准备Index列表拆分，需要IndexRC列，可由ID或Index列检索到
    tail -n+2 ~/rice/zn.sd1/b2/doc/library.txt |awk '{print "Index"$5"\t"$3"\t"$2}' > index.txt
    # 附录：并行拆库，注意必须有wd变量，且indexRC位于第三列
    # 附录：质控与压缩

## 181009.lane14

    # 移动全部文件至根目录
    cd ~/seq/181009.lane14
    for f in `find . -maxdepth 3 | grep 'fq.gz'`; do
        mv $f ./
    done
    # 改名
    rename 's/FCHFLMTBCX2_L1_Index-/L181009_/' *.gz


## 181023WangGD

    # 简化目录和文件名
    mv Clean/CCPM/* ./
    rename 's/FCHL5W2BCX2_L2_CWHPE18090036/lane/' *.fq.gz
    # 并行拆库
    prefix=L181023
    

## lane16 181024WangChao

    # 简化目录和文件名
    mv Clean/DYMSLR18/FCHFLMNBCX2_L2_CWHPE18100007_* ./
    rename 's/FCHFLMNBCX2_L2_CWHPE18100007/lane/' *.fq.gz
    # 并行拆库
    prefix=L181024


## 181025.lane15

    # 简化目录和文件名
    mv Clean/SDRS20180912/FCHL72LBCX2_L2_CWHPEPI00001955_* ./
    rename 's/FCHL72LBCX2_L2_CWHPEPI00001955/lane/' *.fq.gz
    # 按附录代码操作

    ## 郭晓璇的Index13没有结果，沟通后尝试Index11的RC GGCTAC在lane中检索发现较多，修改index列表最后一列，重新分析
    rm L181025_AGTCAA* # 删除错误Index相关文件
    idx=GGCTAC
    prefix=L181025


## 181122.lane17
    
    wd=L181122
    cd ~/seq
    mv F18FTSNCKF3152_LIBagiR $wd
    cd $wd
    # 简化目录和文件名
    mv Clean/G5WH_WWAT/FCHHH57BCX2_L2_CWHPEPI00001961_* ./
    rename 's/FCHHH57BCX2_L2_CWHPEPI00001961/lane/' *.fq.gz

## 190123.lane18

    wd=L190123
    cd ~/seq/$wd
    mv Clean/W18NP_N8R_WTCBAC_ems541/FCHNCTYBCX2_L2_CWHPEPI00001985_* .
    rename 's/FCHNCTYBCX2_L2_CWHPEPI00001985/lane/' *.fq.gz
    sed -i 's/Clean\/W18NP_N8R_WTCBAC_ems541\/FCHNCTYBCX2_L2_CWHPEPI00001985/lane/' md5.txt

## 190220.lane19

    wd=L190220
    cd ~/seq/$wd
    mv Clean/DYMJTR18-WTCBAC2/FCHT7YTBCX2_L1_CWHPEPI00001995_* ./
    rename 's/FCHT7YTBCX2_L1_CWHPEPI00001995/lane/' *.fq.gz
    sed -i 's/Clean\/DYMJTR18-WTCBAC2\/FCHT7YTBCX2_L1_CWHPEPI00001995/lane/' md5.txt

## 190516.lane20

    wd=L190516
    cd ~/seq/$wd
    md5sum -c Clean.md5.txt 
    mv Clean/CWHPEPI00002033/FCHVH3FBCX2_L1_CWHPEPI00002033_*.gz ./
    rename 's/FCHVH3FBCX2_L1_CWHPEPI00002033/lane/' *.fq.gz

## 190620.lane21

    wd=L190620
    cd ~/seq/$wd
    md5sum -c Clean.md5.txt 
    mv Clean/CWHPEPI00002058/FCHVH3FBCX2_L1_CWHPEPI00002033_*.gz ./
    rename 's/FCHVH3FBCX2_L1_CWHPEPI00002033/lane/' *.fq.gz

## 190802.lane22

    wd=L190802
    cd ~/seq/$wd
    md5sum -c Clean.md5.txt 
    mv Clean/CWHPEPI00002079/FCH2NF7BCX3_L2_CWHPEPI00002079_*.gz ./
    rename 's/FCH2NF7BCX3_L2_CWHPEPI00002079/lane/' *.fq.gz
    # 按下方`扩增子lane处理通用代码代码`操作


## 扩增子lane处理通用代码

    # 下载方法：远程桌面 yongxin至bailab
    cd ~/software/GUI-linux/GUI-linux-3.0
    sh bgionline.sh # 复制用户名，密码需手动输入

### 1. 华大lane返回通用处理代码

    # 查看质量报告和评估文件完整性
    tar xvzf upload.tar.gz # 解压报告
    firefox upload/index.html # firefox、下载或远程桌面查看
    md5sum -c md5.txt # 检查下载数据准确性，一致输出OK

    # 准备Index列表拆分，需要IndexRC列，可由ID或Index列检索到
	# 情况1：如保存ID列为index.ID，获取index和IndexRC，并补充到实验设计
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' \
        /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt index.ID > index.txt
    # 情况2：如保存Index序列，保存为index，获取indexID和IndexRC，并补充到实验设计
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$2]=$0} NR>FNR {print a[$1]}' \
        /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt index > index.txt
    cat index.txt
    # 保存至库文件和登记 SeqLibraryList.xlsx

	# 检查数据质量格式，如果64跳转换，否则跳至下节
	determine_phred-score.pl lane_1.fq.gz
	# 如果为64，改原始数据为33
	rename 's/lane/lane_64/' lane_*
	# 关闭质量控制，主要目的是格式转换64至33，而不是过滤序列，否则扩增子必须成对不满足则usearch无法合并
	time fastp -i lane_64_1.fq.gz -I lane_64_2.fq.gz -o lane_1.fq.gz -O lane_2.fq.gz -6 -A -G -Q -L -w 8
	# 1lane 80GB, 2 threads, 102min
    # 查看fastp报告，数据问题为双端之和，且单位为百万(M)
    # firefox upload/index.html # fastp.html

    # 并行拆库，按目录名见附件编号A18xxxx，j可按文库数量调整，推荐默认24 # |head -n2 只处理前两个新增的
    parallel -j 24 "zcat lane_1.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > ${wd}_{1}_1.fq" ::: `cut -f 3 index.txt | grep  -v '^$'`
    parallel -j 24 "zcat lane_2.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > ${wd}_{1}_2.fq" ::: `cut -f 3 index.txt | grep  -v '^$'`

    # 质控和报告汇总
    fastqc -t 48 L*.fq
    multiqc ./ # 详见multiqc_report.html
    # 提取各样品数据量，不同批数据量的列会有变化，可能要更改列的数值5、6等
    l=`head -n1 multiqc_data/multiqc_fastqc.txt|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep 'Total Sequences'|cut -f 1`
    cut -f 1,${l} multiqc_data/multiqc_fastqc.txt | sed 's/_.\t/\t/' | uniq > datasize.txt
    cat datasize.txt

    # 压缩空间及上传
    pigz *.fq
    # 分双端统计md5值
    md5sum L*_1.fq.gz > /tmp/md5sum1.txt
    md5sum L*_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt

    # 汇总至amplicon目录
    ls L*.gz
    ls L*.gz|wc
    ln L*.gz ../amplicon/


## 汇总为amplicon目录
    cd ~/seq
    mkdir -p amplicon && cd amplicon
    ln /mnt/bai/yongxin/seq/novagene/*.gz ./ # 共242个文件，121个库 L160914-L171013
    # 第二批目前20个库，171018-190220
    ls|grep -P '^L' > lane_list.txt # 可选 ls -d L*
    wc -l lane_list.txt
    for i in `ls -d L*`; do ln ${i}/${i}*.gz amplicon/;done
    ls|wc # 库数据为1006个文件，503个库





# meta 

    meta@meta:/mnt/m2/data/meta/

## 2017/12/11 ath 2.5T
    
    cd /mnt/m2/data/meta/ath/2.5T
    # 只有去宿主后的clean_data 137G，而且重复了两份clean_data(tar包含有fq, fastqc, log和out有总体宿主或剩余数量)和remove_host_fq，没有去宿主前和统计信息
    # 样本已经为单样一个文件，直接改名
    mkdir -p seq
    awk 'BEGIN{OFS=FS="\t"}{system("mv remove_host_fq/"$3"/"$3".rmhost.clean.1.fq.gz seq/"$1"_1.fq.gz ");system("mv remove_host_fq/"$3"/"$3".rmhost.clean.2.fq.gz seq/"$1"_2.fq.gz ");}' <(tail -n+2 metadata.txt)
    # 统计样本md5值和数据量见附录代码


## 2017/12/27 rice nrt1.1b
    
    cd /mnt/m2/data/meta/rice/nrt1.1b
    # 只有去宿主后的clean_data 137G，而且重复了两份clean_data(tar包含有fq, fastqc, log和out有总体宿主或剩余数量)和remove_host_fq，没有去宿主前和统计信息
    # 样本已经为单样一个文件，直接改名
    mkdir -p seq
	# 按实验设计按文件夹批量合并再改名，需要输入文件每个样本一个目录
	p=30
    for i in `tail -n+2 metadata.txt|cut -f3`; do
        zcat `find 1.clean_reads/${i}/ -name *.gz | grep '_1.fq'` | pigz -p ${p} > seq/${i}_1.fq.gz &
        zcat `find 1.clean_reads/${i}/ -name *.gz | grep '_2.fq'` | pigz -p ${p} > seq/${i}_2.fq.gz &
    done
    awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$3"_1.fq.gz seq/"$1"_1.fq.gz ");system("mv seq/"$3"_2.fq.gz seq/"$1"_2.fq.gz ");}' <(tail -n+2 metadata.txt)
    # 统计样本md5值和数据量见附录代码


## 2018/11/19 ath 3T
    
    cd /mnt/m2/data/meta/ath/3T
    # 整理代码见 /mnt/m1/yongxin/ath/3T/manual.md
    for i in `tail -n+2 design.txt|cut -f3`; do
        zcat `find 01.filter/${i}/ -name *.gz | grep '_1.fq'` | pigz -p ${p} > seq/${i}_1.fq.gz &
        zcat `find 01.filter/${i}/ -name *.gz | grep '_2.fq'` | pigz -p ${p} > seq/${i}_2.fq.gz &
    done
    awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$3"_1.fq.gz seq/"$1"_1.fq.gz ");system("mv seq/"$3"_2.fq.gz seq/"$1"_2.fq.gz ");}' <(tail -n+2 design.txt)


## 2018/11/21 rice miniCore

    cd /mnt/m2/data/meta/rice/miniCore
    # 了解文件格式
    less -S 00.rawdata1/data.info # 样本和180个原始数据对应表
    find 01.cleandata/*/* -name *.gz | wc -l # 整理好的180个对应文件，己去除宿主
    # 整理原始序列按样本合并、上传、改名
    find 00.rawdata*/*/* -name *.gz|wc # 725文件，单数吗？
    find 00.rawdata*/*/* -name *.gz|cut -f 3 -d '/'|sort|uniq|wc -l # 720，有5个重复，应该是重复多次
    # 非冗余整理到同一文件夹中
    mkdir -p temp
    for i in `find 00.rawdata*/*/* -name *.gz`; do
        ln -sf `pwd`/$i temp/
    done
    # 制作实验设计批量两文件合并并重命名
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$2"_1.fq.gz temp/"$3"_1.fq.gz | pigz -p 30 > seq/"$1"_1.fq.gz")}' <(tail -n+2 metadata.txt|cut -f 1,9,10)
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$2"_2.fq.gz temp/"$3"_2.fq.gz | pigz -p 30 > seq/"$1"_2.fq.gz")}' <(tail -n+2 metadata.txt|cut -f 1,9,10)


## 2018/11/28 rice wetdry
    
    # 制作实验设计改名并移动至seq目录中
    cd /mnt/m2/data/meta/rice/wetdry
    mkdir -p seq
    awk 'BEGIN{OFS=FS="\t"}{system("mv 2.cleandata/"$2"/"$2"_1.clean.fq.gz seq/"$1"_1.fq.gz");system("mv 2.cleandata/"$2"/"$2"_2.clean.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 design.txt)

## 2018/12/21 wheat rootrot
    cd /mnt/m2/data/meta/wheat/rootrot
    ossutil config -e oss.aliyuncs.com -i LTAIOvpGVPaADctY -k xGt7LX4KrReiKpEDgFKgSSthPPJfZg
    ossutil cp oss://novo-data-nj/customer-DIDtwmIC/ ./ -r -f --jobs 3 --parallel 2
    # 根据下载文件编写metadata.txt
    mkdir -p seq
    awk 'BEGIN{OFS=FS="\t"}{system("mv 2.cleandata/"$2"/"$2"_1.clean.fq.gz seq/"$1"_1.fq.gz");system("mv 2.cleandata/"$2"/"$2"_2.clean.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 metadata.txt)

## 2019/4/18 rice miniCore2 非粳籼的66个样品

    cd /mnt/m2/data/meta/rice/miniCore2
    # 了解文件格式
    less -S metadata.txt # 66样本芯片、Lane编号
    # 整理原始序列按样本合并、上传、改名
    find raw/*/* -name *.gz|wc # 264文件
    find raw/*/* -name *.gz|cut -f 3 -d '/'|sort|uniq|wc -l # 264非冗余文件
    # 非冗余整理到同一文件夹中
    mkdir -p temp
    for i in `find raw/*/* -name *.gz`; do
        ln -sf `pwd`/$i temp/
    done
    # 制作实验设计批量两文件合并并重命名
    mkdir -p seq
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$2"_1.fq.gz temp/"$3"_1.fq.gz | pigz -p 30 > seq/"$1"_1.fq.gz")}' <(tail -n+2 metadata.txt|cut -f 1,6,7)
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$2"_2.fq.gz temp/"$3"_2.fq.gz | pigz -p 30 > seq/"$1"_2.fq.gz")}' <(tail -n+2 metadata.txt|cut -f 1,6,7)


## 2019/8/18 medicago wt+lyr4 10个苜蓿宏基因组样品

    # 只复制报告和cleandata
    cd /mnt/m2/data/meta/med/lyr4
    # md5校验，15G的文件要1分半
    ls|grep 'QJM' > cat dir.list
    for i in `cat dir.list`; do
        cd $i
        time md5sum -c MD5*.txt
        cd ..
    done
    # 提取文件名列表，编辑为metadata.txt
    find ./ *.gz|grep '_1'|sort|sed 's/1.clean.fq.gz//'
    # 移动改名
    awk 'BEGIN{OFS=FS="\t"}{system("mv "$2"1.clean.fq.gz "$1"_1.fq.gz")}' <(tail -n+2 metadata.txt)
    awk 'BEGIN{OFS=FS="\t"}{system("mv "$2"2.clean.fq.gz "$1"_2.fq.gz")}' <(tail -n+2 metadata.txt)

## 附录: 宏基因组通用操作代码
    
    # 质控和报告汇总
    cd seq
    fastqc -t 48 *.fq.gz
    multiqc ./ # 详见multiqc_report.html
    # 提取各样品数据量，不同批数据量的列会有变化，可能要更改列的数值5、6等
    l=`head -n1 multiqc_data/multiqc_fastqc.txt|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep 'Total Sequences'|cut -f 1`
    cut -f 1,${l} multiqc_data/multiqc_fastqc.txt | sed 's/_.\t/\t/' | uniq > datasize.txt
    cat datasize.txt

    # 分双端统计md5值
    md5sum *_1.fq.gz > /tmp/md5sum1.txt
    md5sum *_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt



# rna

    meta@meta:/mnt/m2/data/rna/

## 2017/11/27 王鑫拟南芥无菌苗 /mnt/m2/data/rna/171127sterility/ 诺禾致源测序


## 2018/11/20 王鑫拟南芥时间序列 /mnt/m2/data/rna/181120TimeCourse 美格基因测序
    
    # 数据位于meta服务器meta用户的/mnt/m2/data/rna/181120TimeCourse目录[meta@meta:/mnt/m2/data/rna/181120TimeCourse]$
    sed -i 's/\r//' md5sum.txt # 文件尾有windows换行符，替换后才能运行md5sum
    # 整理文件名与样本名对照表
    rename 's/_\w+_R/_/' *.gz
    rename 's/_001.fastq/.fq/' *.gz # 两次删除和替换与库名一致
    # FTP上传GSA/181120RNASeq


## 附录. RNA-Seq通用操作代码

    # 测试文件完整性
    md5sum -c md5sum.txt
    # 根据实验设计批量改名
    mkdir -p seq·
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$9"_1.fq.gz temp/"$10"_1.fq.gz | pigz -p 16 > seq/"$1"_1.fq.gz")}' <(tail -n+2 metadata.txt) &
    awk 'BEGIN{OFS=FS="\t"}{system("zcat temp/"$9"_2.fq.gz temp/"$10"_2.fq.gz | pigz -p 16 > seq/"$1"_2.fq.gz")}' <(tail -n+2 metadata.txt) &

    # 质控和报告汇总
    fastqc -t 48 *.fq.gz
    multiqc ./ # 详见multiqc_report.html
    # 提取各样品数据量，不同批数据量的列会有变化，可能要更改列的数值5、6等
    l=`head -n1 multiqc_data/multiqc_fastqc.txt|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep 'Total Sequences'|cut -f 1`
    cut -f 1,${l} multiqc_data/multiqc_fastqc.txt | sed 's/_.\t/\t/' | uniq > datasize.txt
    cat datasize.txt

    # 分双端统计md5值
    md5sum L*_1.fq.gz > /tmp/md5sum1.txt
    md5sum L*_2.fq.gz > /tmp/md5sum2.txt
    paste /tmp/md5sum1.txt /tmp/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt



# 其它

## 180115.BSA.wx BSA测序筛选突变体位点
    mv clean_data/*.gz ./
    # 通用诺禾改名
    # 质控和报告汇总
    # 分双端统计md5值

## 180601.ath.WGS 基因组测序
    mv customer-sMIXvUoG/02.cleanData/*.gz ./
    # 质控和报告汇总
    # 分双端统计md5值


## 诺合数据下载工具
	# 下载和安装 https://helpcdn.aliyun.com/document_detail/50452.html
	cd ~/bin
	wget http://docs-aliyun.cn-hangzhou.oss.aliyun-inc.com/assets/attach/50452/cn_zh/1524643963683/ossutil64?spm=a2c4g.11186623.2.6.TfLBVN -O ossutil
	chmod +x ossutil
    ## 数据下载指定目录
	# 进入screen环境防断网
	screen -d -r seq
	# 设置下载目录，每次不同数据请修改此处
	wd=~/seq/180601.ath.WGS/
	mkdir -p $wd
	# 创建配置文件，不同用户名请修改用户名和密码
	ossutil config -e oss.aliyuncs.com -i LTAIhcpcybDqnbxI -k c3oJcIIDzc9z3H7f8CENdC4taSqUoG
	# 下载文件
	ossutil cp oss://novo-data-nj/customer-sMIXvUoG/ $wd -r -f --jobs 3 --parallel 2
    
    # 通用诺禾改名
    rename 's/clean.//' *gz
    rename 's/_\w+_L\d_/_/' *.gz

    # 质控和报告汇总
    fastqc -t 48 *.fq.gz
    multiqc ./ # 详见multiqc_report.html
    # 提取各样品数据量，不同批数据量的列会有变化，可能要更改列的数值5、6等
    l=`head -n1 multiqc_data/multiqc_fastqc.txt|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep 'Total Sequences'|cut -f 1`
    cut -f 1,${l} multiqc_data/multiqc_fastqc.txt | sed 's/_.\t/\t/' | uniq > datasize.txt
    cat datasize.txt

    # 分双端统计md5值
    md5sum *_1.fq.gz > md5sum1.txt
    md5sum *_2.fq.gz > md5sum2.txt
    paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > md5sum.txt
    cat md5sum.txt

# 备份

## 1amplicon / 5Xiaoning 扩增子备份
    
    # 远程桌面 yongxin@210.75.224.110(biocloud)，挂载需使用bennyyu密码zpp...，可用sudo passwd重置

    # novagene
    cp -r ~/seq/novagene /media/yongxin/1amplicon/ # 121 samples, 232.4 GB, 也可选程桌面托拖拽复制

    # BGI
    cd ~/seq
    # 链接所有数据至BGI目录
    for i in `ls -d L*`; do ln ${i}/${i}* BGI/;done
    ls BGI/*_1.fq.gz|wc # 331 文件
    # 统计每组中样品
    ls BGI/*_1.fq.gz | cut -f 1 -d '_' | uniq -c
    # 质控汇总BGI
    cp -r ~/seq/novagene /media/yongxin/1amplicon/ # 331 samples, 861.1 GB, 也可选程桌面托拖拽复制


## 2meta5T / 6QinYuan 宏基因组去宿主后数据备份
        
    # 复制 yongxin@meta 中 ath 目录中的 2.5T/3T ， rice 目录中的 miniCore/nrt1.1b/wetdry ，这些目录中的submit目录，包括质控和去宿主后的clean文件，统计文件(质控和宿主比例)，实验设计metadata.txt和md5值。
    # 将复制完成的2meta5T内容全盘复制到6QinYuan中

## 3RNABSA / 7Zhiwen 其它：转录组、BSA测序结果备份
    # yongxin@biocloud上
    target=/media/yongxin/3RNABSA/
    scp -r meta@192.168.0.32:~/data/rna/* ${target}
    cd ~/seq
    cp -r 180115.BSA.wx ${target}
    cp -r 180601.ath.WGS ${target}


## 附录代码：

	# 按indexRC提取完整index信息，库中为index.ID，获取index和IndexRC，并补充到实验设计
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$3]=$0} NR>FNR {print a[$2]}' \
        /mnt/bai/yongxin/ref/culture/IlluminaIndex48.txt doc/library.txt