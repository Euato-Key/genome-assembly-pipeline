# 数据下载教程

如果您的 `sequence_files` 目录中数据不完整，或者需要下载数据库，请参考以下指南。

## 1. 基因组数据库下载

### Kraken2 数据库 (用于污染去除)
Kraken2 需要预构建的数据库。标准库体积很大 (>50GB)，建议下载 Standard-8 或 PlusPF-8 等精简版。

```bash
mkdir -p database/kraken_db
cd database/kraken_db

# 下载 Standard-8 数据库 (约 8GB)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz
tar -xzvf k2_standard_08gb_20231009.tar.gz
```
*下载地址参考: [BenLangmead/aws-indexes](https://benlangmead.github.io/aws-indexes/k2)*

### BUSCO 数据库 (用于质量评估)
Busco 会在运行时自动下载所需的 lineage 数据集，但建议手动下载并指定路径以避免网络问题。

```bash
mkdir -p database/busco_db
cd database/busco_db

# 查看可用数据集
busco --list-datasets

# 下载特定数据集 (例如哺乳动物 mammalia_odb10)
wget https://busco-data.ezlab.org/v5/data/lineages/mammalia_odb10.2024-01-08.tar.gz
tar -xzvf mammalia_odb10.2024-01-08.tar.gz
```

## 2. 公共测试数据集下载

如果您没有自己的测序数据，可以使用公共数据库中的数据进行练习。推荐使用 **E. coli (大肠杆菌)** 数据，因为体积小，跑得快。

### HiFi 数据
```bash
# E. coli K12 MG1655 PacBio HiFi reads
wget -O test_data/HiFi.fasta.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/Ecoli/Ecoli_K12_MG1655_PacBio_HiFi.fastq.gz
gunzip test_data/HiFi.fasta.gz
# 转换为 fasta (如果需要)
sed -n '1~4s/^@/>/p;2~4p' test_data/HiFi.fastq > test_data/HiFi.fasta
```

### Hi-C 数据
Hi-C 数据通常很大，练习时可以只下载前几百万条 reads。
```bash
# 示例链接 (需替换为真实 SRA 链接)
# 使用 SRA Toolkit 下载
fastq-dump --split-files --gzip -X 1000000 SRR12345678
```

## 3. 常用数据下载工具

- **SRA Toolkit**: 下载 NCBI SRA 数据。
- **Aspera (ascp)**: 高速下载 EBI/NCBI 数据。
- **Wget / Curl**: 下载直链文件。

## 4. 验证数据完整性
下载完成后，务必检查 md5sum：
```bash
md5sum -c md5sum.txt
```
