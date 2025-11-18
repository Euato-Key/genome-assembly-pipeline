# 使用示例

本文档提供了基因组组装流程的各种使用场景示例。

## 🎯 场景 1: 完整流程运行

适用于：首次运行完整的基因组组装流程

### 步骤 1: 准备环境

```bash
# 安装conda环境
bash setup_environment.sh

# 激活环境
conda activate genome_assembly

# 验证环境
bash scripts/validate.sh
```

### 步骤 2: 配置参数

```bash
# 复制配置模板
cp config.example.yaml config.yaml

# 编辑配置文件
nano config.yaml
```

修改关键参数为你的实际路径。

### 步骤 3: 运行流程

```bash
# 运行完整流程
bash assembly_pipeline.sh run

# 在另一个终端监控进度
tail -f logs/assembly_*.log
```

---

## 🎯 场景 2: 跳过某些步骤

适用于：部分数据已处理，想跳过某些步骤

修改 config.yaml：

```yaml
steps:
  run_qc: false              # 跳过质控（数据已质控）
  run_genome_survey: false   # 跳过基因组调查
  run_kraken: false          # 跳过污染检测
  run_hifiasm: true          # 直接开始组装
```

运行：
```bash
bash assembly_pipeline.sh run
```

---

## 🎯 场景 3: 从断点续跑

适用于：流程中断后继续运行

```bash
# 查看上次完成的步骤
cat logs/status.log

# 从断点继续
bash assembly_pipeline.sh resume
```

---

## 🎯 场景 4: 仅运行质控

```yaml
steps:
  run_qc: true
  run_fastp: true
  # 其他步骤设为 false
```

质控结果在：`01_qc/` 目录

---

## 🎯 场景 5: 测试运行

```bash
# 运行测试脚本
bash test.sh

# 使用测试配置运行
CONFIG_FILE=config.test.yaml bash assembly_pipeline.sh run
```

---

## 📊 查看结果

### BUSCO结果
```bash
cat 05_busco/MyProject/short_summary.txt
```

### QUAST结果
```bash
cat 06_quast/report.txt
```

### 组装统计
```bash
grep -c ">" 04_hifiasm/MyProject.contigs.fa
```

---

## 🧹 清理和重新运行

### 清理中间文件
```bash
bash assembly_pipeline.sh clean
```

### 完全重新运行
```bash
rm -rf 0[1-9]_* logs/
bash assembly_pipeline.sh run
```

---

更多信息请参考：
- `README.md` - 完整文档
- `QUICKSTART.md` - 快速开始
- `PROJECT_STRUCTURE.md` - 项目结构
