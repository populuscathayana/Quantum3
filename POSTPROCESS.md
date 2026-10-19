# 后处理范围

`reproduce.py --check-only` 核对文件清单、Gaussian 能量与频率、ORCA 数值 Hessian、热化学参考态和报告数据。默认运行还重绘完整电子能阶梯、NEB/CI-NEB/IRC 对照、Gibbs 图、IRC 曲线与关键结构，并生成完整 IRC 动画。

## 能量与路径

- Gaussian 与 ORCA 各用自身参考态计算能差。
- 电子能路径与离散 RRHO Gibbs 自由能分别作图；热化学为 298.15 K、1 atm。
- 双分子活化自由能以两个分离反应物的自由能之和为基准。
- 原环加成 IRC 含 54 个接受点和 TS；有限 IRC 末点与另外优化的极小值分列。
- 构象分支采用已完成的 NEB、CI-NEB、精修鞍点和频率结果；构象 IRC 已停止。

## 结构与轨道

IRC 结构和动画使用真实接受点，保持原子编号，只作平移和正旋转配准。前线轨道等值面和虚频图随报告提供；默认复现保留这些成图。原始轨道 cube 网格不纳入发布包。

## 输出

计算文件在 `data/calculations/`，核验和绘图数据在 `analysis/`，文件清单为 `data/manifest.json`。派生图写到 `build/figures/`；可用 `--output` 指定其他目录。后处理不修改计算文件或报告，也不调用量子化学程序。
