# 腾讯云 IDE 快速部署指南

## 问题说明

在腾讯云 IDE 等云环境中，使用 `conda` 安装大量包时可能会遇到：
- 依赖解析超时（卡在 `Solving environment` 步骤）
- 网络连接不稳定
- 镜像源同步延迟

## 解决方案

本项目提供了**三种**安装方式，按推荐顺序：

### 方案 1：快速安装脚本（推荐）⭐

使用优化的快速安装脚本，自动使用 `mamba` 加速：

```bash
# 在腾讯云 IDE 终端中运行
cd /workspace  # 或您的项目目录

# 运行快速安装脚本
bash scripts/setup_tools_fast.sh
```

**特点：**
- ✅ 使用 `mamba` 替代 `conda`（速度提升 5-10 倍）
- ✅ 分批安装，避免依赖解析超时
- ✅ 仅安装核心工具，可选工具失败不影响主流程
- ✅ 自动检测 CPU 指令集

**预计时间：** 5-10 分钟

---

### 方案 2：标准安装（已优化）

如果快速脚本遇到问题，使用改进后的标准脚本：

```bash
bash scripts/setup_tools.sh
```

**改进点：**
- 分 4 个步骤安装，每步单独处理
- 添加错误重试机制
- 支持多个镜像源备选

**预计时间：** 10-20 分钟

---

### 方案 3：手动安装（最灵活）

如果自动脚本都失败，可以手动逐步安装：

```bash
# 1. 创建环境
conda create -n genome_assembly_env -y python=3.9

# 2. 激活环境
conda activate genome_assembly_env

# 3. 安装 mamba（加速后续安装）
conda install -y -n base -c conda-forge mamba

# 4. 使用 mamba 安装核心工具
mamba install -y -c bioconda -c conda-forge \
    fastp jellyfish samtools chromap yahs wget pigz

# 5. 安装 hifiasm（根据 CPU 选择）
# 如果 CPU 支持 AVX2：
mamba install -y -c bioconda hifiasm

# 如果 CPU 不支持 AVX2（会报 Illegal instruction 错误）：
# 需要从源码编译，运行：
bash scripts/compile_hifiasm.sh

# 6. 可选工具（失败可跳过）
mamba install -y -c bioconda genomescope2 kraken2 busco
```

---

## 常见问题

### Q1: 卡在 "Solving environment" 怎么办？

**方法 1：** 强制中断（Ctrl+C）并使用快速脚本
```bash
# 按 Ctrl+C 中断
bash scripts/setup_tools_fast.sh
```

**方法 2：** 安装 mamba 后重试
```bash
conda install -y -n base -c conda-forge mamba
# 然后编辑脚本，将所有 conda 替换为 mamba
```

**方法 3：** 减少一次性安装的包数量（手动安装）

---

### Q2: 网络超时怎么办？

```bash
# 增加 conda 超时时间
conda config --set remote_read_timeout_secs 600

# 或使用国内镜像（清华源）
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
```

---

### Q3: hifiasm 报 "Illegal instruction" 错误？

这是 CPU 指令集不兼容问题，需要编译通用版本：

```bash
# 1. 激活环境
conda activate genome_assembly_env

# 2. 安装编译依赖
conda install -y -c conda-forge gcc_linux-64 gxx_linux-64 make git

# 3. 运行编译脚本
bash scripts/compile_hifiasm.sh

# 4. 验证
hifiasm --version
```

---

### Q4: 如何验证安装是否成功？

```bash
# 激活环境
conda activate genome_assembly_env

# 运行验证脚本
bash scripts/validate_environment.sh
```

**预期输出：**
```
[INFO] Docker found: Docker version 20.10.x
[INFO] Git found: git version 2.x.x
[INFO] Python3 found: Python 3.9.x
[INFO] hifiasm found: hifiasm 0.19.x
[INFO] CPU does not support AVX2 - using generic build
[INFO] All required tools are available. Environment validation passed.
```

---

## 腾讯云 IDE 特殊注意事项

### 1. 持久化存储

腾讯云 IDE 的 `/workspace` 目录是持久化的，建议将项目放在这里：

```bash
cd /workspace
git clone <your-repo>
cd genome-assembly-pipeline_docker
```

### 2. 内存限制

云 IDE 可能有内存限制，建议：
- 使用快速安装脚本（占用内存更少）
- 避免同时运行多个大型任务
- 编译时使用 `-j2` 而不是 `-j$(nproc)`

### 3. 网络代理

如果云环境配置了代理，可能需要：

```bash
# 查看代理设置
env | grep -i proxy

# 如果需要，为 conda 配置代理
conda config --set proxy_servers.http http://your-proxy:port
conda config --set proxy_servers.https https://your-proxy:port
```

---

## 推荐工作流程

### 首次使用

```bash
# 1. 进入项目目录
cd /workspace/genome-assembly-pipeline_docker

# 2. 快速安装（推荐）
bash scripts/setup_tools_fast.sh

# 3. 验证环境
conda activate genome_assembly_env
bash scripts/validate_environment.sh

# 4. 配置参数
vim config/config.sh

# 5. 运行测试
bash scripts/run_assembly.sh
```

### 后续使用

```bash
# 每次打开 IDE 后
cd /workspace/genome-assembly-pipeline_docker
conda activate genome_assembly_env

# 运行分析
bash scripts/run_assembly.sh
```

---

## 性能对比

| 安装方式 | 预计时间 | 成功率 | 适用场景 |
|---------|---------|--------|---------|
| **快速脚本** | 5-10 分钟 | ⭐⭐⭐⭐⭐ | 云 IDE、网络不稳定 |
| **标准脚本** | 10-20 分钟 | ⭐⭐⭐⭐ | 本地环境、网络良好 |
| **手动安装** | 15-30 分钟 | ⭐⭐⭐ | 自定义需求、调试 |

---

## 故障排除

### 完全重置环境

如果安装出现严重问题，可以完全重置：

```bash
# 1. 删除环境
conda deactivate
conda env remove -n genome_assembly_env -y

# 2. 清理缓存
conda clean --all -y

# 3. 重新安装
bash scripts/setup_tools_fast.sh
```

### 查看详细日志

```bash
# 启用 conda 调试模式
conda install --debug -c bioconda hifiasm

# 查看 conda 配置
conda config --show

# 查看已安装的包
conda list
```

---

## 获取帮助

如果遇到问题：

1. **查看文档**：`docs/hifiasm_cpu_compatibility.md`
2. **运行验证**：`bash scripts/validate_environment.sh`
3. **查看日志**：`logs/` 目录下的日志文件
4. **提交 Issue**：附上错误信息和系统信息

---

**最后更新：** 2025-12-03  
**适用版本：** 腾讯云 IDE、阿里云 IDE、VS Code Remote 等云环境
