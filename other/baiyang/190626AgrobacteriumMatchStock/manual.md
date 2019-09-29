# 查找某农杆菌是否存在于菌库中 —— 黄安诚

    cd ~/other/baiyang/190626AgrobacteriumMatchStock/

## 菌名查找Agrobacterium rhizogenes strains-可信度不高

    grep 'Agrobacterium' ~/culture/all/culture_select.tax
    # 在拟南芥和番茄中有一些属，但没有种

    # 检查数据库是否物种注释中来源存在
    # 查询种
    grep 'rhizogenes' /db/usearch/*.fa > species.db
    wc -l species.db # 71条记录
    cut -f 1 -d ':' species.db | uniq -c # # 只在silva和rdp中有
    less species.db # RDP和Silva中有Rhizobium_rhizogenes
    # 查询属
    grep 'Agrobacterium' /db/usearch/*.fa > genus.db
    wc -l genus.db # 2801条记录
    cut -f 1 -d ':' genus.db | uniq -c # 只在gg和silva中有
    less genus.db # RDP和Silva中有Rhizobium_rhizogenes


## 序列查找

    # 参考基因组来自 Valdes Franco JA, Collier R, Wang Y, Huo N, Gu Y, Thilmony R, Thomson JG. 2016. Draft genome sequence of Agrobacterium rhizogenes strain NCPPB2659. Genome
Announc 4(4):e00746-16. doi:10.1128/genomeA.00746-16.
    # 下载地址：https://www.ncbi.nlm.nih.gov/nuccore/LYBK00000000 ，点WGS LYBK01000001-LYBK01000018，Download，fasta 
    wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/LY/BK/LYBK01/LYBK01.1.fsa_nt.gz
    gunzip LYBK01.1.fsa_nt.gz
    # 与所有菌库比
    blastn -query LYBK01.1.fsa_nt -db ~/culture/all/culture_select.fa -out LYBK01_all.blastn -outfmt 6  -num_alignments 100 -evalue 1 -num_threads 9 
    # 存在99.97相近的 tomato_9, med_21, ath_94
    # ~/culture/rice/stock
    #blastn -query LYBK01.1.fsa_nt -db ~/culture/rice/stock/sequence.fa -out LYBK01_rice.blastn -outfmt 6  -num_alignments 100 -evalue 1 -num_threads 9 
    # ath_94，王鑫找A907.fa
    makeblastdb -dbtype nucl -in A907.fa
    blastn -query LYBK01.1.fsa_nt -db A907.fa -out LYBK01_ath.blastn -outfmt 6  -num_alignments 1 -evalue 1 -num_threads 9 # 全长有5个mismatch


