# 2022 年的工具（保留作参考）

`legacy/2022/src/` 是 2022-05 的原始代码，**一字未改**，只是从 `src/` 移到这里。2026-11 的重写在仓库根目录（`singular_vector/`、`experiments/`）。

## 当年的做法

- `spo.sage`：`SpO(m, n)` = spo(2m|2n)（m 是 sp 的秩，n 是 o 的秩），矩阵放在 gl(2m+2n) 里，sp 块在前；奇根矩阵手写，与论文 §2.4 的表逐项相同，δ_i − δ_j 用的是正确的减号。
- `AlgebraBSC.sage`、`MatSp.sage`：由矩阵基解出结构常数。
- `pbw.sage`、`pbws.sage`：在 **U(g)** 里按 PBW 序正规排序；`to_verma` 丢掉含 e 的项、把 h 换成权的值。
- `Lie.sage`：根的换算、`positive_dec`（枚举 Kostant 分拆）。
- `rex.py`：把方程导到 MATLAB 求解、生成 LaTeX 表格。

## 与 2026 版的对应（`compare_with_2026.py`，需要 `passagemath-symbolics`，请装在单独的环境里，见根目录 README）

- `SpO(m, n)` 是 2026 版的 `OSP(2n, 2m)`。下标：旧 k < m → u(k)；m ≤ k < 2m → v(k−m)；2m ≤ k < 2m+n → u(k−m)；k ≥ 2m+n → v(k−m−n)。旧的 δ_i 是新的 δ_{m+1−i}。
- 旧的根向量是新生成元的 ±1 倍（±2δ 为 1/2 倍）。
- 用旧代码自己的 `positive_dec`、`SPBWElement`、`to_verma` 建方程、在 Sage 里精确求解（不用 MATLAB），奇异向量与 2026 版**相同**（整体差 −1）：
  spo(2|4)、spo(2|6)、spo(4|4)、spo(2|8)、spo(4|6) 的 11 个（代数, 根），结果在 `data/experiments/legacy-comparison.jsonl`。

## 为什么重写

1. 两个相同的奇元相邻判为 0，只对 isotropic 奇根成立；osp(2M+1|2N) 的 ±δ_i 不对。`Oth` 只有 D 型（n ≥ 2），所以奇数维正交部分和 osp(2|2N) 都做不了。
2. `positive_dec` 会给出奇根重数 ≥ 2 的分拆（在 U(n⁻) 里是 0），当成未知数会多出自由变量。
3. `AlgebraBSC.check_lie_identity` 两个分支都返回 `True`（Jacobi 检查等于没做）；`Element.__eq__` 少参数、`Element.__mul__` 用了未定义的全局 `N`、`checke_lie_closure` 少 `@staticmethod`。
4. 权在 Sage 的符号环里，在 U(g) 里排序再投影，建方程组比 2026 版慢三到四个数量级（osp(8|2) 的 δ1+ε1：14.6 秒对 0.003 秒）；求解要经 MATLAB 与剪贴板，不可复现。

仍然有用的：论文记号下的矩阵表（2026 版的 `osp.paper_matrix` 是同一张表）；在 U(g) 里做 PBW 运算（以后若要处理 Šapovalov 元本身，可以从这里起步）。
