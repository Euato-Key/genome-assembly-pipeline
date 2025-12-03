# 解决 hifiasm CPU 指令集兼容性问题

## 问题背景

`hifiasm` 是一个高性能的基因组组装工具，官方预编译版本通常使用 AVX2 等高级 CPU 指令集进行优化。如果您的 CPU 不支持这些指令集，运行时会出现以下错误：

```bash
Illegal instruction (core dumped)
```

## 本项目的解决方案

本项目已经集成了自动检测和适配机制，使用 **Bioconda** 作为主要的工具管理方式，并在必要时自动从源码编译通用版本。

### 🔍 自动检测流程

当您运行 `bash scripts/setup_tools.sh` 时，脚本会：

1. **检测 CPU 指令集**
   - 检查 `/proc/cpuinfo` 中的 `flags` 字段
   - 判断是否支持 `avx2`、`sse4_2` 等指令集

2. **智能选择安装方式**
   - ✅ **支持 AVX2**：直接从 Bioconda 安装优化版（性能最佳）
   - ⚠️ **不支持 AVX2**：自动从源码编译通用版（兼容性最佳）

3. **验证安装**
   - 运行 `bash scripts/validate_environment.sh` 检查所有工具
   - 显示 CPU 指令集支持情况

## 使用方法

### 方法 1：自动安装（推荐）

```bash
# 1. 运行安装脚本（会自动检测 CPU 并选择合适的安装方式）
bash scripts/setup_tools.sh

# 2. 验证环境
bash scripts/validate_environment.sh
```

**预期输出示例：**

```
[INFO] 正在检测 CPU 指令集...
[INFO] 检测到 SSE4.2 指令集支持（不支持 AVX2）
[INFO] CPU 不支持 AVX2，将从源码编译通用版 hifiasm...
[INFO] 克隆 hifiasm 源码...
[INFO] 编译 hifiasm（使用通用 x86-64 指令集）...
[INFO] 通用版 hifiasm 编译并安装完成
```

### 方法 2：手动检查 CPU 指令集

如果您想先了解自己的 CPU 支持情况：

```bash
# 查看 CPU 型号
lscpu | grep -i "model name"

# 查看支持的指令集
grep -E "^flags" /proc/cpuinfo | head -1

# 检查是否支持 AVX2
grep -q avx2 /proc/cpuinfo && echo "支持 AVX2" || echo "不支持 AVX2"
```

### 方法 3：强制重新编译通用版

如果您已经安装了 Bioconda 版本但遇到问题，可以强制重新编译：

```bash
# 1. 激活 conda 环境
conda activate genome_assembly

# 2. 卸载现有的 hifiasm
conda remove -y hifiasm

# 3. 安装编译依赖
conda install -y gcc_linux-64 gxx_linux-64 make zlib git

# 4. 克隆并编译
cd ~/code/genome-assembly-pipeline_docker/tools
git clone https://github.com/chhylp123/hifiasm.git
cd hifiasm
make clean
make CFLAGS="-g -O3 -march=x86-64 -Wall" -j$(nproc)

# 5. 复制到 conda 环境
CONDA_PREFIX=$(conda info --base)/envs/genome_assembly
cp hifiasm $CONDA_PREFIX/bin/
chmod +x $CONDA_PREFIX/bin/hifiasm

# 6. 验证
hifiasm --version
```

## 性能对比

| CPU 指令集 | 安装方式 | 相对性能 | 兼容性 |
|-----------|---------|---------|--------|
| **AVX2** | Bioconda 预编译版 | 100% (最快) | 仅限支持 AVX2 的 CPU |
| **SSE4.2** | 源码编译 (x86-64) | ~85% | 兼容大多数现代 CPU |
| **通用版** | 源码编译 (x86-64) | ~70-85% | 兼容所有 x86-64 CPU |

> **注意**：即使使用通用版，hifiasm 仍然是一个非常高效的组装工具，性能损失主要体现在大规模数据集上。

## 常见问题

### Q1: 如何确认当前使用的是哪个版本？

```bash
# 激活环境
conda activate genome_assembly

# 查看 hifiasm 位置
which hifiasm

# 查看版本
hifiasm --version

# 检查是否为编译版本（如果在 tools/hifiasm_build 目录下则为编译版）
ls -la $(which hifiasm)
```

### Q2: 编译失败怎么办？

如果自动编译失败，可能的原因：

1. **缺少 git**：`conda install -y git`
2. **网络问题**：手动下载源码
   ```bash
   wget https://github.com/chhylp123/hifiasm/archive/refs/heads/master.zip
   unzip master.zip
   cd hifiasm-master
   make CFLAGS="-g -O3 -march=x86-64 -Wall"
   ```
3. **编译器问题**：确保安装了 `gcc_linux-64` 和 `gxx_linux-64`

### Q3: 可以在 Docker 中使用吗？

可以！本项目支持 Docker 部署。在 Dockerfile 中，脚本会自动检测容器内的 CPU 指令集并选择合适的安装方式。

### Q4: 为什么不直接使用通用版？

- **性能考虑**：如果您的 CPU 支持 AVX2，使用优化版可以显著提升速度（特别是处理大型基因组时）
- **灵活性**：自动检测机制确保在不同硬件上都能正常运行

## 技术细节

### 编译选项说明

```makefile
# 原始 Makefile（使用本机指令集）
CFLAGS = -g -O3 -march=native -Wall

# 通用版 Makefile（兼容所有 x86-64）
CFLAGS = -g -O3 -march=x86-64 -Wall
```

- `-march=native`：使用本机所有可用的指令集（可能不兼容其他机器）
- `-march=x86-64`：仅使用基础 x86-64 指令集（兼容性最好）
- `-march=core-avx2`：使用 AVX2 指令集（性能好但需要 CPU 支持）

### Bioconda 包说明

Bioconda 提供的 `hifiasm` 包通常使用 `-march=x86-64-v2` 或更高级别编译，这在大多数现代 CPU 上都能运行，但在一些老旧或虚拟化环境中可能会遇到问题。

## 相关资源

- [hifiasm GitHub](https://github.com/chhylp123/hifiasm)
- [Bioconda hifiasm 包](https://anaconda.org/bioconda/hifiasm)
- [CPU 指令集参考](https://en.wikipedia.org/wiki/X86-64#Microarchitecture_levels)

## 更新日志

- **2025-12-03**：添加自动 CPU 指令集检测和智能安装机制
- **2025-12-03**：集成 Bioconda 作为主要包管理工具
- **2025-12-03**：添加源码编译备选方案

---

如有其他问题，请查看项目主 README 或提交 Issue。
