# 数据格式（给网站用）

全部由 `experiments/` 里的脚本生成，可以重跑覆盖；数字是精确计算的结果。字符串里的公式是 LaTeX（KaTeX 可直接渲染），`text` 字段是 Sage 的纯文本。

## `experiments/<id>.json`：一条实验记录

| 字段 | 含义 |
| --- | --- |
| `id`、`title` | 实验名 |
| `question`、`method`、`scale`、`result` | 问题、方法、规模、结论（英文短句，网站可译） |
| `outcome` | `confirmed`（在所述范围内由精确计算证明）、`evidence`（支持性证据，如猜测 C1）、`refuted` |
| `reproduce` | 复现命令（在仓库根目录） |
| `generated_at`、`software` | 生成时间、passagemath 与 Python 版本 |
| `rows` | 逐例结果（每个实验各自的字段，见脚本文件头） |

现有：`conjecture-8`、`llw-formula`、`example-10`、`singular-vectors`、`gl-so-factorisation`、`fbeta-coefficient`；
`legacy-comparison.jsonl` 是 2022 工具与 2026 版的逐例比对。

## `cases/<id>.json`：一个（代数, 根）的奇异向量

`id` 形如 `osp6-2-d1pe1`（osp(6|2)，δ1+ε1；`p` = 加号，`m` = 减号）。索引 `cases/index.json` 列出全部案例的摘要。

| 字段 | 含义 |
| --- | --- |
| `algebra`、`algebra_latex`、`beta`、`beta_latex` | 代数与根 |
| `rho`、`rho_latex`、`simple_roots_latex` | ρ（坐标 (a_N…a_1 \| b_1…b_M)）与单根 |
| `hyperplane` | 超平面 (λ+ρ, β) = 0 解出 a_s 的式子（`text`、`latex`）；`free_variables` 是其余坐标 |
| `pbw_order` | 负根向量的 PBW 次序（`name`、`latex`、`odd`） |
| `weight_space_dim`、`equations`、`top_degree`、`symbol_blocks` | 权空间维数、方程数、顶次、各次数的符号块（列数、行数、秩） |
| `statements` | `symbol_criterion`（≤ 1 对每个 λ）、`exactly_one_on_H`、`every_lambda`（对每个 λ：一维 ⇔ 在超平面上） |
| `normalisation` | 单根向量单项式的系数为 1 |
| `terms` | 每一项：`word`、`latex`（单项式）、`pbw_degree`、`odd_factor`、`coefficient`（`text` 展开式、`factored`、`latex` 因式分解、`degree`、`linear_factors_only`）；项数超过 60 时为 `null`（`terms_omitted` 为真） |
| `f_beta_coefficient` | 根向量 f_β 本身那一项的系数（猜测 C1 说的就是它） |
