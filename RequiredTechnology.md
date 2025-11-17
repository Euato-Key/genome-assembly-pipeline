### **一、 计算环境 (Computing Environment)**

这是项目的基础，不可或Rack。

*   **操作系统 (OS):** **Linux** (推荐 CentOS 7, Rocky Linux 8/9, 或 Ubuntu 20.04/22.04 LTS)。
*   **核心要求:** 熟练使用**命令行界面 (CLI)** 和 **SSH** 远程登录。
*   **硬件配置 (Hardware):**
    *   **CPU:** **64核** 或以上。
    *   **内存 (RAM):** **512GB** 起步，推荐 **1TB** 或更高（Hifiasm非常消耗内存）。
    *   **存储 (Storage):** **10TB** 以上的高速存储空间 (SSD/NVMe)。
*   **作业调度系统 (Job Scheduler):** **Slurm** 或 SGE (如果你在HPC集群上工作)。

### **二、 核心软件栈 (Core Software Stack)**

按照流程图的顺序，为你指定每个步骤的核心工具。

1.  **软件包管理 (Package Management):**
    *   **Conda / Mamba:** **必须掌握**。使用它来创建独立环境并安装以下所有软件，以解决复杂的依赖关系。

2.  **数据质控与预处理 (QC & Pre-processing):**
    *   HiFi数据 `ccs提取`: **PacBio SMRT Link** 套件中的 `ccs` 命令。
    *   HiFi/ONT数据 `去潜在污染`: **Kraken2** (用于物种分类) + **`seqtk subseq`** (用于根据Kraken2结果提取或剔除reads)。
    *   ONT数据 `过滤`: **pychopper**。
    *   ONT数据 `质控`: **fastplong**。
    *   Hi-C/WGS数据 `质控`: **fastp**。

3.  **基因组评估 (Genome Survey):**
    *   K-mer计数: **jellyfish**。
    *   基因组大小估计: **GenomeScope2.0** (Genoscope的现代继任者)。

4.  **基因组组装 (Genome Assembly):**
    *   Contig级别组装: **Hifiasm** (v0.16.1或更高版本)。这是整个流程的核心，利用HiFi和ONT数据生成高质量的contigs。

5.  **染色体挂载 (Scaffolding):**
    *   Hi-C数据比对: **Chromap**。
    *   染色体级别组装: **yaHs** (通常与Hifiasm一同分发)。

6.  **手动校正 (Manual Curation):**
    *   可视化与编辑: **Juicebox Assembly Tools (JBAT)**。这是进行精细手动调整的黄金标准工具。

7.  **质量评估 (Quality Assessment):**
    *   连续性评估: **QUAST**。
    *   完整性评估: **BUSCO** (v5或更高版本)。

### **三、 编程与流程管理 (Programming & Workflow Management)**

这是将所有工具粘合起来，实现自动化的关键。

*   **必备脚本语言:** **Bash Shell Scripting**。用于编写执行命令、管理文件和构建简单流程的脚本。
*   **高级脚本语言:** **Python 3**。
    *   **核心库:** `Biopython` (序列处理), `pandas` (数据分析), `matplotlib`/`seaborn` (绘图)。用于解析复杂的输出文件、生成统计报告和图表。
*   **流程管理系统 (推荐):** **Nextflow**。
    *   这是构建可重复、可扩展、可移植的生物信息学流程的工业标准。当你需要反复运行或与他人协作时，它的价值无与伦比。

---

### **技术栈清单总结 (Checklist)**

**如果你想从零开始搭建这个能力，请按以下顺序学习和准备：**

1.  **[√] 基础:**
    *   [ ] 熟练掌握Linux命令行。
    *   [ ] 学会使用Conda/Mamba管理软件环境。
    *   [ ] 掌握Bash脚本编程基础。

2.  **[√] 核心工具:**
    *   [ ] **质控:** `fastp`, `fastplong`
    *   [ ] **组装:** `Hifiasm`
    *   [ ] **挂载:** `Chromap`, `yaHs`
    *   [ ] **评估:** `QUAST`, `BUSCO`
    *   [ ] **校正:** `Juicebox`

3.  **[√] 进阶能力:**
    *   [ ] 掌握Python和相关生信库。
    *   [ ] 学习使用Nextflow构建自动化流程。

4.  **[√] 知识背景:**
    *   [ ] 理解不同测序技术（HiFi, ONT, Hi-C）的原理。
    *   [ ] 了解基因组组装的关键概念（N50, BUSCO得分等）。