# 环境与离线复现

## Python

建议 Python 3.12。依赖为 NumPy、Pillow、ASE、Matplotlib 和 pypdf；Matplotlib 提供可移植的 DejaVu 字体。

```sh
git clone https://github.com/populuscathayana/Quantum3.git
cd Quantum3
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py --check-only
python reproduce.py
```

Windows 可用 `.venv\Scripts\python.exe` 替代上述激活后的 `python`。

`--check-only` 核对文件清单、原始能量/频率/Hessian、参考态、路径数据和报告链接。默认运行另外重绘能量阶梯、路径对照、Gibbs 图和完整 IRC 曲线、结构及动画，输出至 `build/figures/`。`--output` 可指定另一个派生图目录。原始计算文件和已提交的报告图不被修改。

仅首次安装依赖需要网络。数据核对与重绘均离线执行，不需要安装量子化学软件。

## 报告 PDF 与 TeX

[现成 PDF](experiment-report/cycloaddition_mechanism.pdf) 包含全部 15 页和可点击的仓库链接；[TeX](experiment-report/cycloaddition_mechanism.tex) 包含同一链接及相对图片路径。

重新排版使用 XeLaTeX，需要 `fontspec`、`xeCJK`、`graphicx`、`amsmath`、`amssymb`、`array`、`tabularx`、`booktabs`、`float`、`caption`、`titlesec`、`geometry`、`xcolor`、`listings`、`xurl`、`hyperref`、`bookmark`、`placeins`。中文字体优先采用 Songti/Heiti/Kaiti SC，否则使用 Fandol；英文回退字体为 TeX Gyre Termes。

在另行安装的 TeX 环境中，可将输出放到独立目录：

```sh
mkdir -p build/tex
cd experiment-report
xelatex -interaction=nonstopmode -halt-on-error -output-directory=../build/tex cycloaddition_mechanism.tex
xelatex -interaction=nonstopmode -halt-on-error -output-directory=../build/tex cycloaddition_mechanism.tex
```

仓库地址列于报告正文末尾的中文“数据与代码”节；封面和科学内容保持原样。
