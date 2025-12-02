# 参数调节指南 (Parameter Tuning Guide)

本指南帮助您根据不同的基因组特性、测序平台和硬件资源，对 **config/config.sh** 中的关键参数进行合理调节，以获得最佳的组装质量和运行效率。

---

## 1. 资源配置

### `THREADS`
- **默认**: `12`
- **建议**: 将其设置为服务器可用 CPU 核心数的 70%~80%。
  - 64 核服务器 → `48`
  - 32 核服务器 → `24`
  - 个人电脑 (8 核) → `4`（保留系统资源）
- **注意**: 过高的线程数在 I/O 密集型任务（如 fastp）可能导致磁盘争用，适当调低可提升整体稳定性。

---

## 2. Fastp（数据质控）

| 参数 | 含义 | 推荐值 | 说明 |
|------|------|--------|------|
| `FASTP_MIN_LEN` | 最小读长 | `145`（Illumina 150bp） | 若数据质量较差，可降低至 `100` 或 `75`，但会增加噪声。 |
| `FASTP_THREAD` | 线程数 | `${THREADS}` | 与全局线程保持一致。 |
| `FASTP_QUALIFIED_QUALITY` | 质量阈值 | `20` | 过滤低质量碱基。 |
| `FASTP_CUT_MEAN_QUALITY` | 滑动窗口平均质量 | `20` | 适用于大多数数据。 |

---

## 3. Jellyfish & GenomeScope（基因组调查）

### `KMER_SIZE`
- **默认**: `19`
- **调节原则**:
  - 小基因组 (<100 Mb) → `17` 或 `19`
  - 中等基因组 (100 Mb‑1 Gb) → `19`‑`21`
  - 大基因组 (>1 Gb) → `21`‑`23`
- **原因**: 较大的基因组需要更长的 k‑mer 以避免碰撞；但 k‑mer 过大时覆盖度下降。

### `GENOME_SIZE_EST`
- **默认**: `1g`
- **调节**: 设为预估基因组大小（如 `3g`、`500m`），影响 Jellyfish 哈希表大小。若设置过小，Jellyfish 会自动扩容但消耗更多内存。

---

## 4. Kraken2（污染过滤）

| 参数 | 含义 | 推荐值 | 说明 |
|------|------|--------|------|
| `KRAKEN_CONFIDENCE` | 置信度阈值 (0‑1) | `0.5` | 对 HiFi 长读建议较高阈值，降低误删。 |
| `KRAKEN_MIN_HIT_GROUPS` | 最小命中组数 | `300` | HiFi 数据建议保持较高，以确保分类可靠。 |

---

## 5. Hifiasm（基因组组装）

> 这些参数在 `scripts/run_assembly.sh` 中的 `run_assembly` 函数里可自行修改。

| 参数 | 含义 | 推荐值 | 说明 |
|------|------|--------|------|
| `-l` (purge level) | 单倍型纯化级别 (0‑3) | `1` | `0` 适用于近交系，`2` 适用于高杂合度物种。 |
| `--primary` | 生成主要 contig | `--primary` | 保持默认即可。 |
| `--h1` / `--h2` | Hi‑C 辅助分型 | 视需求而定 | 若想得到两套 haplotype，加入此参数。 |

---

## 6. BUSCO（完整性评估）

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `BUSCO_MODE` | 评估模式 | `genome` |
| `BUSCO_DB` | 数据库路径 | 参考 `config.sh` 中的 `BUSCO_DB`，建议使用对应物种的 lineage（如 `mammalia_odb10`）。 |

---

## 7. Chromap（Hi‑C 比对）

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `CHROMAP_PRESET` | 预设模式 | `hic` |
| `CHROMAP_REMOVE_DUP` | 去除 PCR 重复 | `--remove-pcr-duplicates` |

---

## 8. YAHS（Scaffolding）

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `-l` (min contig length) | 最小 contig 长度 | `0`（不过滤） |
| `-e` (酶切位点) | 限制性酶 | 若已知建库酶可填入，如 `GATC`，否则留空自动检测 |

---

## 9. Juicer（可视化）

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `JUICER_JVM_MEM` | Java 最大堆内存 | 小基因组 `10g`‑`20g`，大型基因组 `80g`‑`120g`（需对应机器内存）。 |

---

## 10. 常见场景预设

| 场景 | 关键参数调整 |
|------|-------------------|
| **细菌（5 Mb）** | `THREADS=4`，`KMER_SIZE=15`，`GENOME_SIZE_EST=10m`，`JUICER_JVM_MEM=4g` |
| **哺乳动物（3 Gb）** | `THREADS=48`，`KMER_SIZE=21`，`GENOME_SIZE_EST=3g`，`KRAKEN_CONFIDENCE=0.5`，`JUICER_JVM_MEM=120g` |
| **高杂合植物（>5 Gb）** | `THREADS=64`，`KMER_SIZE=21`，`-l 2`（hifiasm）或使用 `--h1 --h2`，`JUICER_JVM_MEM=150g` |

---

### 使用方式
1. 打开 `config/config.sh`。
2. 根据上述建议修改对应变量的值。
3. 保存后运行 `bash scripts/run_assembly.sh`。

如需进一步细化参数（例如不同测序平台的特定设置），欢迎在 `docs/parameter_tuning.md` 中补充或直接联系我。祝您组装顺利！
