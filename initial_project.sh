# 16S/ITS v2 version project initial

## 设置工作目录
    # 2018/7/18 rice/integrate16s 重新分析miniCore类群和SL分蘖等9套水稻数据
    # 2018/7/24 rice/nrt1.1a 氮吸收利用基因nrt1.1a与菌群
    # 2018/7/24 ath/integrate16s 整理拟南芥萜类二半萜和三萜数据
    # 2018/8/7 ath/wx.16s 王鑫16S功能基因验证
    # 2018/8/29 rice/zn.sd1/b2 水稻第二批SD1
    # 2018/10/11 ath/jt.HuangAC/syncom 三萜重组
    # 2018/10/23 culture/ath/gaochulei/control 培养菌对照
    # 2018/10/26 ath/CCPM 拟南芥代谢物突变体
    # 2018/10/31 maize/magic 玉米GWAS
    # 2018/12/1 medicago/AMF/its 苜宿AMF ITS数据
    # 2018/12/1 ehbio/qianxubo 肠道测试数据
    # 2018/12/7 ath/myb28 基因对叶际微生物的影响
    # 2018/12/8 wheat/FHB/perithecium 赤霉病子囊壳微生太
    # 2018/12/13 culture/ath/starting ath分菌起始样品比较 starting
    # 2018/12/13 culture/ath/leaf ath分菌起始叶样品比较分离菌 leaf
    # 2018/12/31 rice/timecourse/v2 水稻时间序列基于新数据格式和新版流程分析
    # 2019/1/21 ath/jt.HuangAC/coevolve 回答审稿人是否
    # 2019/2/1 ath/SA 水杨酸定量课题
    # 2019/2/26 data/temp/3tb3 三萜第三批，拆分原始数据提交GSA
    # 2019/2/27 wheat/profile 小麦时间序列重复分析
    # 2019/2/28 maize/magic/v2 玉米拔节期MAGIC群体数据
    # 2019/3/7 rice/zn.sd1/v3 水稻SD1第三批
    # 2019/3/11 rice/Gprotein/v2 水稻G蛋白第二批
    # 2019/3/25 medicago/culture_start 苜蓿分菌起始样品
    # 2019/4/9 maize/magic/v2salt 玉米亲本耐盐
    # 2019/5/14 medicago/AMF2 苜蓿按新思路从头分析
    # 2019/6/3 rice/hinge1/its 水稻hinge1真菌
    # 2019/7/1 rice/hinge1/16s 水稻hinge1细菌
    # 2019/7/6 rice/integrate16s/v2/LN 水稻GWAS LN
    # 2019/7/6 rice/integrate16s/v2/HN 水稻GWAS HN
    # 2019/9/10 medicago/AMF3 苜蓿按新思路并修改实验设计
    # 2019/10/8 rice/hinge1/16s2 水稻hinge1细菌重测数据替换部分文库
    # 2019/10/24 rice/alkali 水稻耐盐碱能力
    # 2019/10/25 rice/timecourse 水稻时间序列
    # 2019/11/20 ath/pi.absort.GaoCL 拟南芥磷吸引
    # 2019/11/20 rice/integrate16s/v2OTU/pubGB2019 水稻公共GB数据95个
    # 2019/12/13 rice/blast 水稻稻瘟病 - 卢训丽
    # 2020/2/17 wheat/FHB/sourcetrack 小麦赤霉病 - 陈云
    # 2020/5/29 rice/rnaSL - 张婧赢
    # 2020/11/8 rice/miniCore2 - miniCore样本重拆库为双端
    # 2020/11/9 rice/strigolactone.LiJY2 - strigolactone.LiJY样本重拆库为双端
    # 2020/12/3 medicago/AMF4 苜蓿按新思路并修改实验设计
    # 2020/12/10 medicago/AMF4/SynCom2Sand 苜蓿SynCom2沙子体系实验有参分析
    # 2020/12/16 medicago/AMF4/SynCom1Sand 苜蓿SynCom1沙子体系实验有参分析
    # 2020/12/24 medicago/AMF4/SynCom1SandDenovo 苜蓿SynCom1沙子体系实验
    # 2020/12/27 medicago/AMF4/SynCom2SandDenovo 苜蓿SynCom2沙子体系实验
    # 2020/12/31 medicago/AMF4/SynCom1FlowPotDenovo 苜蓿SynCom1在FlowPot体系中实验



## 项目建立代码
    wd=medicago/AMF4/SynCom1FlowPotDenovo
    # 创建github目录，用于备份流程、文档
    cd ~
    mkdir -p ~/github/Work/$wd
    cd ~/github/Work/$wd
    cp ~/github/Amplicon/16Sv2/parameter.md ./
    cp ~/github/Amplicon/16Sv2/manual.md ./
    # 链接代码至工作区
    mkdir -p ~/$wd
    ln -s `pwd`/parameter.md ~/$wd/makefile
    ln -s `pwd`/manual.md ~/$wd/manual.sh
    cd ~/$wd



# 16S cutlure v2 version project initial

## 设置工作目录
    # 2019/6/26 culture/medicago/190626 苜蓿分菌测序
    # 2019/6/26 culture/rice/190626 水稻分菌测序R2A和厌氧
    # 2019/7/4 culture/wheat1907 小麦分菌测序
    # 2019/8/10 culture/maize/190810 玉米分菌3个index 1Nova+2BGI

## 项目建立代码
    wd=culture/maize/190810
    # 创建github目录，用于备份流程、文档
    cd ~
    mkdir -p ~/github/Work/$wd
    cd ~/github/Work/$wd
#    ln ~/github/Amplicon/16Sculture2/parameter.md ./
#    ln ~/github/Amplicon/16Sculture2/manual.md ./
    cp ~/github/Amplicon/16Sculture2/parameter.md ./
    cp ~/github/Amplicon/16Sculture2/manual.md ./
    # 链接代码至工作区
    mkdir -p ~/$wd
    ln -s `pwd`/parameter.md ~/$wd/makefile
    ln -s `pwd`/manual.md ~/$wd/manual.sh
    cd ~/$wd

# 其它项目：菌库序列和物种注释整理

## 设置工作目录
    # 2019/7/12 culture/rice/stock 水稻分菌测序
    # 2019/9/28 other/baiyang/190626AgrobacteriumMatchStock 鉴定农杆菌是否存在于菌库中

## 项目建立代码
    wd=other/baiyang/190626AgrobacteriumMatchStock
    # 创建github目录，用于备份流程、文档
    cd ~
    mkdir -p ~/github/Work/$wd
    cd ~/github/Work/$wd
    touch manual.md
    # 链接代码至工作区
    mkdir -p ~/$wd
    ln -s `pwd`/manual.md ~/$wd/manual.sh
    cd ~/$wd



# 备份项目实验设计
    
    cd ~/github/Work
    # 项目列表 project_list.txt
    for wd in `cat project_list.txt | cut -f 2 `; do
    #wd=medicago/AMF3
    cp -r ~/$wd/doc/ ~/github/Work/$wd/doc/
    done



# 备份代码

    # 文件加入缓冲区
    cd ~/github/Work
    find . -size +50M -print # 检查大文件
    git add . 
    # 提交修改
    git commit -m "v2019" 
    # 推送到github
    git push origin master
