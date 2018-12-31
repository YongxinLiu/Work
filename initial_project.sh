# 16S v2 version project initial

# 设置工作目录
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

wd=rice/timecourse/v2
# 创建github目录，用于备份流程、文档
cd 
mkdir -p ~/github/Work/$wd
cd ~/github/Work/$wd
cp ~/github/Amplicon/16Sv2/parameter.md ./
cp ~/github/Amplicon/16Sv2/manual.md ./
# 链接代码至工作区
mkdir -p ~/$wd
ln -s `pwd`/parameter.md ~/$wd/makefile
ln -s `pwd`/manual.md ~/$wd/manual.sh
cd ~/$wd

# 备份代码
# 文件加入缓冲区
cd ~/github/Work
git add . 
# 提交修改
git commit -m "2018/12/7" 
# 推送到github
git push origin master
