# 数据来源与第三方内容

## 紫外实验参考

[NIST Chemistry WebBook：山梨酸，CAS 110-44-1](https://webbook.nist.gov/cgi/cbook.cgi?ID=C110441&Mask=400)。

- 原始文献：Timmons C. J., *UV atlas of organic compounds*, 1966, 1, B7/1。
- 数据标识：NIST / RAS UV No. 19947。
- [JCAMP 来源](https://webbook.nist.gov/cgi/cbook.cgi?Index=0&JCAMP=C110441&Type=UVVis)；访问日期 2026-10-09。
- 原文件标注数据来源为 INSTITUTE OF ENERGY PROBLEMS OF CHEMICAL PHYSICS, RAS；所有者字段为 INEP CP RAS, NIST OSRD，并附有 2007 年美国商务部长代表美国的汇编版权及保留权利声明。

仓库提供此次比较所用的 126 个数值点（波长与 log10 ε）和相应派生曲线，未重新分发完整 JCAMP 文件或原图集。这是历史图谱的数值化数据；条目未给出溶剂、浓度和温度。NIST 收录不等于本仓库有权将第三方谱图改授 MIT 许可。再次使用或分发应遵循来源的使用与引用条件。

## 氢核实验参考

大槻崇，*NMRを用いた食品中の食品添加物分析法の開発に関する研究*，国立医薬品食品衛生研究所，2011 年 4 月。

[厚生劳动科学研究成果数据库报告 PDF](https://mhlw-grants.niph.go.jp/system/files/2010/104031/201033054A/201033054A0001.pdf)，访问日期 2026-10-09。山梨酸峰位位于 Fig. 1（印刷页 11 / PDF 第 13 页），实验条件见印刷页 3 和 Table 1。溶剂 DMSO-d6，600 MHz，25 °C，DSS-d6 内标。

`nmr-h-reference.json` 保留六组明确标注的峰位及其原子对应关系。仓库不分发报告全文或原谱图，不从图像估算线宽、强度或多重峰形。

## 碳核与耦合常数实验参考

Nishiyama Y., Tamura Y., Ibe A., et al. *Study on Isomers of Sorbic Acid Produced in Food*. Eisei Kagaku **37**(2), 89–96 (1991)。

[论文页面与 DOI](https://doi.org/10.1248/jhs1956.37.89)，访问日期 2026-10-09。数据取自 Table II / III 的 SOA 列（印刷页 93 / PDF 第 5 页）；实验使用 CDCl3，TMS 内标。

仓库只保留六个碳位移及两个已分辨的耦合常数。文献从羧基开始给碳编号，本项目从甲基开始编号，两者对应关系为“项目碳号 = 7 − 文献碳号”。保留原文归属，不根据理论峰序重写实验标签。论文全文和原谱图不随仓库分发。

各 JSON 中的 `local_source` 是原始文献档案名，供来源追溯，并非运行依赖。实验比较图和相关 CSV 使用上述第三方数据，不能仅凭本仓库许可证忽略原始权利。

## 本次计算、软件与标识

Gaussian 输入和计算结果用于记录方法与数值结果；输出中的软件版本、引用和原有声明保留。公开 `.out` 副本仅替换一行程序及临时目录路径；`.gjf` 和 `.fchk` 保持原样。本仓库不包含 Gaussian 程序、许可证或二进制检查点。已有结果的后处理不需要该软件。

Python、NumPy、Pillow、R、ggplot2 和 TeX 及其字体分别遵循各自许可；仓库不重新分发这些软件或字体。

报告中的浙江大学名称及封面标识用于学术报告署名，不表示机构对本仓库背书；校徽源图未分发。该标识不属于 MIT 许可的原创资源。

## 许可范围

根目录 [MIT License](LICENSE) 与 Quantum1 一致，适用于作者有权许可的原创代码和文档。第三方谱图数据、引文、软件声明和标识仍由各自权利人保有权利，本仓库不对它们作额外授权。
