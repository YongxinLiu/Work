# 16S v2 version project initial

# 设置工作目录
# 2018/7/18 rice/integrate16s 重新分析miniCore类群和SL分蘖等9套水稻数据
# 2018/7/24 rice/nrt1.1a 氮吸收利用基因nrt1.1a与菌群
# 2018/7/24 ath/integrate16s 整理拟南芥萜类二半萜和三萜数据
# 2018/8/7 ath/wx.16s 王鑫16S功能基因验证
# 2018/8/29 rice/zn.sd1/b2 水稻第二批SD1
# 2018/10/11 ath/jt.HuangAC/syncom 三萜重组



wd=ath/jt.HuangAC/syncom
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
git commit -m "2018 new project Oct." 
# 推送到github
git push origin master
