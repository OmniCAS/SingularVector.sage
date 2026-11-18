# SingularVector.sage

在 Lie 超代数的 Verma 模里精确计算奇奇异向量（Liu–Luo–Wang, *Odd singular vector formula for general linear Lie superalgebras*, arXiv:1903.04683）。
2026-11 重写：基于 SageMath 库（passagemath 10.8），全部是有理数与 ℚ 上多项式的精确运算；2022 年的工具保留在 [`legacy/`](legacy/README.md)。

这里的结果都是**算出来的**（精确计算），不是形式化证明；标“猜测”的是猜测。

## 安装

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python tests/test_small.py
```

`requirements.txt` 钉住 `passagemath-singular==10.8.14`、`passagemath-repl==10.8.14`（Singular、FLINT、PARI 的二进制 wheel，macOS arm64 / Linux x86_64，Python 3.11–3.14）。

**注意**：不要把 `passagemath-symbolics` 装进跑实验的同一个环境。2026-11-18 在 macOS arm64 上观察到：装了它之后，9 个变量的多项式环上的运算
（例如 osp(6\|12)、osp(8\|10)、osp(10\|8)）会在 libsingular 里段错误；卸掉之后同样的计算正常。`legacy/compare_with_2026.py` 需要符号环，
请单独建一个环境：`uv pip install passagemath-singular==10.8.14 passagemath-repl==10.8.14 passagemath-symbolics==10.8.14`。

## 结构

| 路径 | 内容 |
| --- | --- |
| `singular_vector/algebra.py` | 由矩阵给出的李超代数：括号（矩阵超交换子）、根、单根、高度、PBW 序、U(n⁻) 的左乘 |
| `singular_vector/osp.py` | osp(M\|2N)（M 奇偶都行），生成元 `T/P/Q/Rp/Rm` 与 Lean 形式化的 `Defs/OSP.lean` 相同；论文 §2.4 的根向量矩阵 `paper_matrix` |
| `singular_vector/gl.py` | gl(m\|n)（`Defs/GL.lean` 的约定）；论文公式 (3)–(7) 的 `S_formula` |
| `singular_vector/verma.py` | Verma 模；奇异向量方程组 D(λ)v = 0（对 λ 仿射）|
| `singular_vector/singular.py` | 符号判据、按 PBW 次数自上而下的精确解、“对每个 λ”的完整结论、随机点核对 |
| `singular_vector/export.py` | JSON 输出（多项式带 LaTeX） |
| `experiments/` | 具体实验，每个写一份 `data/experiments/<id>.json` |
| `data/` | 实验结果与逐例数据（见下） |
| `tests/test_small.py` | 回归测试（约 1 秒） |

约定：权坐标 λ = (a_N, …, a_1 \| b_1, …, b_M)，a_i 对应 δ_i、b_k 对应 ε_k；(δ, δ) = −1、(ε, ε) = 1；正根系是字典序（distinguished Borel）；
f_γ 是权为 −γ 的根向量（就是对应的生成元）；PBW 序把偶根向量排在前、奇根向量在后。

## 实验

在仓库根目录运行（`python` 指装好 passagemath 的那个）。每个实验输出 `data/experiments/<id>.json`，字段与姊妹站点“探索记录”相同：
`question`、`method`、`scale`、`result`、`outcome`（`confirmed` = 在所述范围内由精确计算证明；`evidence` = 支持性证据；`refuted`）、`reproduce`、`rows`（逐例数据）。

| id | 问题 | 命令 | 目前的结果 |
| --- | --- | --- | --- |
| `conjecture-8` | Conjecture 8（β = δ_s + ε_t）：权 λ−β 的奇异向量是否恰好一维？ | `python experiments/conjecture8.py --all` | 50 个代数、401 个（代数, 根），最大到 osp(14\|2)（2625 维）：对**每个** λ，维数在 (λ+ρ,β)=0 时为 1，否则为 0（confirmed） |
| `llw-formula` | 公式 (7) 对 gl(m\|n)（Theorem 1）以及用 `T a b` 写出时对 osp 的 δ_s − ε_t（Proposition 9）是否给出奇异向量？ | `python experiments/llw_formula.py` | m+n ≤ 6 的 gl 70 例、osp 41 例，随机点上全部成立，对照全部被拒 |
| `example-10` | 论文 Example 10 的 15 项对不对？ | `python experiments/example10.py` | 第 3、6 项单项式印错（权不对），按唯一的单因子修正后 15 个系数全对 |
| `singular-vectors` | 奇异向量逐项长什么样？ | `python experiments/singular_cases.py` | 逐例数据在 `data/cases/`（见下）；每个 PBW 单项式的系数都非零 |
| `gl-so-factorisation` | 对 N、s 的依赖 | `python experiments/factorization.py` | δ_s+ε_t 的向量 = gl(N) 部分（一次因子之积）⊗ osp(M\|2) 的 δ1+ε_t 向量（8 例） |
| `fbeta-coefficient` | **猜测 C1**：f_β 那一项的系数 | `python experiments/fbeta.py` | 是 c·∏_{γ∈Δ₀⁺,(β,γ^∨)>0} ∏_{j<(β,γ^∨)} ((λ+ρ,γ^∨) − j)，c = (−1)^{m+1}（M = 2m）或 (−1)^{m+1}/2（M = 2m+1）；35 例成立，其中 20 例是猜出之后才算的（evidence，不是证明） |

逐例数据 `data/cases/<id>.json`（索引 `data/cases/index.json`）：代数、β、ρ、单根、超平面（文本与 LaTeX）、PBW 序（每个根向量的 LaTeX 名）、
权空间维数、三条结论、归一化后奇异向量的每一项（单项式的 LaTeX、系数的展开式、因式分解式与 LaTeX、是否一次因子之积）。

## English summary

Exact computations (QQ and polynomial rings over QQ, SageMath library via passagemath 10.8) of odd singular vectors in Verma
modules of gl(m|n) and osp(M|2N), in the conventions of a Lean formalisation of Liu–Luo–Wang (arXiv:1903.04683).
Experiments write `data/experiments/<id>.json` (question, method, scale, result, outcome, reproduce, rows) and case files
`data/cases/<id>.json` with LaTeX, ready for a website.  Results are computed, not formally proved; C1 is a conjecture.
The 2022 tools are in `legacy/` (unchanged); `legacy/compare_with_2026.py` shows they give the same singular vectors.

License: MIT.
