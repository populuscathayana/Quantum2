# 环境与复现说明

本仓库复现的是**已有计算结果的核验、数值处理和绘图**。所需计算输出已经提供，不需要安装量子化学计算软件。

## 1. 环境

验证使用 Python 3.12.14、NumPy 2.3.5、Pillow 12.3.0，以及 R 4.3.1 / ggplot2 3.4.3。Python 依赖版本固定在 `requirements.txt`。建议使用 Python 3.12；R 绘图需要支持 `linewidth` 的 ggplot2（3.4 或更新版本）。

Python 和 R 可分别从 [Python 官方下载页](https://www.python.org/downloads/) 与 [CRAN](https://cran.r-project.org/) 安装。只处理数值数据时不需要 R。

在仓库根目录执行：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Windows PowerShell 中创建和启用环境：

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

如需重绘谱图，在 R 中安装 ggplot2：

```r
install.packages("ggplot2", repos = "https://cloud.r-project.org")
```

验证 R 可用：

```bash
Rscript -e 'library(ggplot2); print(packageVersion("ggplot2"))'
```

PNG 输出自动选择可用的 Quartz 或 Cairo 后端；无图形桌面的环境需使用支持 Cairo 的 R。PDF 使用核心 Helvetica 字体；结构和轨道 PNG 依次尝试 Arial、Liberation Sans、DejaVu Sans。若系统没有这些字体，可将 `QUANTUM2_FONT` 和 `QUANTUM2_FONT_BOLD` 分别设为普通与粗体 TrueType 字体文件的完整路径。字体差异可能改变字形和布局，不改变数值结果。

## 2. 三种复现范围

### 核验计算文件并重建数值表格

```bash
python reproduce.py --data-only --verify-raw
```

`--verify-raw` 对仓库中的 16 项计算逐一检查：

- 输入方法、4 核 / 8 GB 设置、正常结束与文件完整性。
- 组成、电荷、自旋、原子顺序、输出与 fchk 几何一致性、E/E 构型和成键关系。
- 优化频率与 fchk 一致性，并由笛卡尔 Hessian 独立重算 42 / 90 个内部模式。
- 20 / 40 个单重激发态与既有数据记录一致；由 fchk 跃迁偶极核验振子强度。
- NMR 屏蔽张量的迹、完整 J 矩阵及 FC / SD / PSO / DSO 分量之和。

这是科学数据核验，不依赖原来的服务器账号、调度器、归档包或传输日志。原始工作流程的核验记录保留在相应 JSON 中；其 `source_directory`、原始文件名通过 `data/manifest.json` 映射到公开副本。

数据重建以已提供的科学 JSON 为输入，重新生成基础跃迁/屏蔽/耦合 CSV，以及理论和实验比较表。`--verify-raw` 先检查这些 JSON 与计算文件是否相符，但不覆盖它们。检查到不一致会退出，不自动改写数据。

### 重绘结构与谱图

```bash
python reproduce.py
```

重新生成单体结构图、基准 UV / NMR 图、五个体系的理论比较图和实验比较图，以及 R 绘图所用的展宽曲线。所有图输出到原有 `analysis/` 子目录，同时保留 PNG 与矢量 PDF（结构图为 PNG）。

轨道图默认使用已提交版本。完整重建轨道图：

```bash
python reproduce.py --verify-raw --orbitals
```

轨道由 `data/calculations/conf-02/uv.fchk` 中的基函数和 MO 系数重建，包含解析重叠矩阵核验；这一步比普通绘图慢，不执行 SCF 计算。

如果 `Rscript` 不在命令搜索路径中：

```bash
python reproduce.py --rscript /完整路径/Rscript
```

### 单独检查

```bash
python scripts/verify_public_data.py
python -m unittest discover -s scripts -p 'test_experiment.py' -v
python scripts/check_report.py
```

复现入口也会自动运行实验数据转换测试和报告数值检查。运行时不要加 Python 的 `-O` 参数，以保留科学断言检查。

## 3. 输入与输出的界限

**稳定输入：** `data/` 中的计算副本、`references/experimental/` 中的来源摘录，以及下列科学数据记录：

- `analysis/stage1/audit.json`
- `analysis/dimers/*/audit.json`
- `analysis/{stage2,conf01}/verified-results.json`
- `analysis/dimer-properties/*/{uv,nmr}/verified-results.json`

**重建输出：** 基础 CSV、`analysis/comparison/` 与 `analysis/experimental/` 的比较数据、展宽曲线、PNG / PDF 图。启用 `--orbitals` 时另更新轨道 PNG 和 `orbital-validation.json`。

命令会覆盖这些派生文件。修改过图或表时，请先提交或另存改动。复现不会修改输入文件、科学 JSON、报告 TeX 或报告 PDF。

## 4. 预期结果

- 16 项计算通过公开副本核验；单体和二聚体分别有 42、90 个正振动模式。
- conf-02 最强 UV 跃迁：257.44 nm，振子强度 0.9375。
- 五个体系展宽峰位依次为：256.20、257.44、264.32、263.53、265.31 nm。
- NIST UV 数值摘录共 126 点，峰顶平台为 255.56–256.41 nm。
- conf-02 的两个实验耦合比较：15.6217 对 15.2 Hz；11.5586 对 10.0 Hz。

浮点运算末位、图像字体和 PDF 元数据可能随运行库变化。复现检查以数值一致性、图表内容和报告引用为准，不要求生成图像逐字节相同。

## 5. 报告排版（可选）

已提供报告 PDF；重绘数据不要求 TeX 环境。需要修改排版时，使用带中文支持和 Fandol 字体的 XeLaTeX 环境。从 `experiment-report/` 目录编译两次：

```bash
cd experiment-report
xelatex -interaction=nonstopmode -halt-on-error sorbic_acid_uv_nmr.tex
xelatex -interaction=nonstopmode -halt-on-error sorbic_acid_uv_nmr.tex
```

需保留上一级 `analysis/`，因为 TeX 使用相对路径引用图表。源码会回退到 Fandol 与 TeX Gyre 字体；校徽源图未分发，缺少时自动省略。自行编译的封面和字形因此可能与现有 PDF 略有不同。

这些命令会更新报告 PDF。复现入口本身不调用 TeX。源码数值检查不等同于 PDF 编译或版面检查。

## 6. 实验来源

离线复现使用本仓库的数值摘录，不需要下载文献全文或访问数据库。核对原文时，请使用 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) 中的原始链接与 JSON 中记录的图表位置。

UV 为文献数值化参考曲线；NMR 为明确列出的峰位及耦合常数，不是作者自行测量的数据。引用、使用和再分发这些资料时应保留来源并遵守原始条件。
