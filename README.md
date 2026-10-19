# 1,3-丁二烯与乙烯的环加成反应

B3LYP-D3(BJ)/def2-SVP 气相计算：环加成过渡态、双向 IRC、NEB/CI-NEB 对照，以及环己烯 A/B 构象转化。

- **[实验报告 PDF（15 页）](experiment-report/cycloaddition_mechanism.pdf)**
- [报告 TeX 源文件](experiment-report/cycloaddition_mechanism.tex)
- [环境与离线复现](ENVIRONMENT.md)
- [数据索引](data/README.md) · [来源及第三方说明](THIRD_PARTY_NOTICES.md)

PDF 与 TeX 均注明本仓库地址。报告版本：2026-10-19。

## 主要结果

| 结果 | 数值 |
|---|---:|
| 环加成普通 NEB 电子势垒，以弱复合物为参考 | 15.7360 kcal/mol |
| 环加成 CI-NEB 电子势垒，以弱复合物为参考 | 15.7381 kcal/mol |
| Gaussian 环加成 TS 电子势垒，以弱复合物为参考 | 15.7407 kcal/mol |
| Gaussian 双分子活化自由能，298.15 K、1 atm | 26.4693 kcal/mol |
| ORCA 环扭转精修鞍点相对 A / B 的电子能高度 | 5.5907 / 0.0272 kcal/mol |
| ORCA 环扭转精修鞍点唯一虚频 | 67.56i cm⁻¹ |

环加成 IRC 连接弱结合反应物复合物与环己烯 B。A 是更稳定的独立构象；构象 CI-NEB 已收敛，双向构象 IRC 未完成并已停止。各程序分别使用自身能量参考，自由能以两个分离反应物为基准。

## 能量与路径

![完整电子能阶梯](experiment-report/figures/electronic-energy-ladder-complete.png)

![NEB、CI-NEB 与 IRC 对照](experiment-report/figures/neb-cineb-irc-comparison.png)

![Gibbs 自由能](experiment-report/figures/gibbs-free-energy-profile.png)

[完整 IRC 动画](experiment-report/figures/irc-full-path.gif) · [环加成 TS 虚频](experiment-report/figures/ts-imaginary-mode.gif) · [CI-NEB 峰顶虚频](experiment-report/figures/cineb-ci400-imaginary-mode.gif)

## 内容与复现

| 目录 | 内容 |
|---|---|
| `experiment-report/` | 报告 PDF、TeX、11 幅报告图及 3 个动画 |
| `sources/gaussian/` | 驻点、环加成 IRC 的输入/输出；驻点格式化检查点 |
| `sources/orca/` | 环加成峰顶及构象 A/B/TS 的数值频率、Hessian、梯度和 TS 优化结果 |
| `data/` | 热化学表、完整 IRC 点、NEB/CI-NEB 路径、结构及对照数据 |
| `scripts/` | 数据核对与移植后的绘图程序 |
| `MANIFEST.json` | 科学文件的大小及 SHA256 清单 |

```sh
python -m pip install -r requirements.txt
python reproduce.py --check-only
python reproduce.py
```

默认将重绘结果写到 `build/figures/`：完整能量阶梯、NEB/CI-NEB/IRC 对照、Gibbs 图、完整 IRC 曲线/结构/动画。报告现有图保持原样。复现使用仓库内的数据和输出，不运行 Gaussian 或 ORCA。

前线轨道等值面以成图形式提供；其原始波函数网格不纳入此轻量发布包。软件安装包、机器配置、调度记录和中间力调用文件保留在本地。
