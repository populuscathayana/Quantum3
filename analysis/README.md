# 数据索引

| 文件 | 物理量与来源 |
|---|---|
| `thermochemistry.json` | Gaussian 各驻点 E、ZPE、H、原始及对称数修正 G；298.15 K、1 atm |
| `energy-ladder.json` | 完整六能级及构象 NEB/CI/TS 高度，单位 kcal/mol |
| `path-comparison.json` | 原环加成 NEB52、CI400、Gaussian TS/IRC 的能量、成键距离与力 |
| `irc.json` | 原环加成 54 个 IRC 点及 TS；含坐标、能量、来源编号与反应坐标 |
| `conformer-audit.json` | ORCA A/B/TS 的能量、振动、残余力及与 Gaussian 鞍点的配准 |
| `conformer-bands.json` | 构象 NEB35、CI50 的完整末态图像、能量、真实力及有效力 |
| `cycloaddition-ci400/` | 环加成 CI 峰顶结构、模式及路径数据 |
| `structures/` | 共同原子编号下的 XYZ 结构 |

Gaussian 原始驻点文件在 `data/calculations/gaussian/`；ORCA 数值频率、Hessian 和梯度在 `data/calculations/orca/`。`.fchk` 是文本格式化检查点，`.hess` 是文本 Hessian，均可在不安装量子化学软件的情况下读取。输出副本隐去了主机名、个人工作目录及临时目录，计算数值保持原样。`data/manifest.json` 核对发布文件的字节。

能量阶梯以两个分离反应物的 Gaussian 电子能之和为零点。构象对照各程序采用自身 A/B 极小值；ORCA 各行统一采用独立端点频率计算的电子能。路径记录中的端点重新求能会有极小的网格旋转差异，保留各自原始值。

原环加成 IRC 的来源编号为 3099（前 28 点）、3286（后 6 点）、3100（反向 20 点）、3025（TS）。3344、3303 分别是另外优化的反应物复合物和产物 B；3011 是独立构象 A，3346 是 Gaussian 环扭转鞍点。

构象 CI-NEB 已收敛，A/B 和精修 TS 的数值频率已完成；构象双向 IRC 未完成并已停止。这里没有加入未计算的路径点或构象混合权重。
