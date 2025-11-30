# solvers.py
# 数值分析·第三章：解线性方程组的直接法
# 实现方法：
#   1. 列主元高斯消元法
#   2. 全主元高斯消元法
#   3. 克劳特 LU 分解
#   4. 平方根法（Cholesky，适用于对称正定矩阵）
#   5. 追赶法（Thomas 算法，适用于三对角线性方程组）

import math


class SolverError(Exception):
    """求解过程中的数值或结构错误。"""
    pass


# ==================== 输入解析与线性性检查 ====================

def parse_matrix(entries):
    """
    entries: GUI 收集到的字符串矩阵（含右端项）
             形状：n 行，每行 n+1 个字符串
    返回: (A, b) 其中 A 为 n×n，b 为长度 n 的向量（float）
    若有非数值项，抛出 SolverError，视为“非线性或非法输入”。
    """
    A = []
    b = []
    for row in entries:
        if len(row) < 2:
            raise SolverError("每一行至少需要一个系数和一个常数项。")
        try:
            nums = [float(x) for x in row]
        except ValueError:
            # 不能转成 float，就把它当成“非线性 / 非法输入”
            raise SolverError("检测到非数值输入，无法构成线性方程组。")
        A.append(nums[:-1])
        b.append(nums[-1])
    n = len(A)
    # 简单检查维数一致
    for row in A:
        if len(row) != n:
            raise SolverError("系数矩阵不是方阵，无法作为标准线性方程组。")
    return A, b


# ==================== 格式化工具 ====================

def format_augmented_matrix(A, b, precision=4):
    """把增广矩阵 A|b 格式化成多行字符串，方便打印。"""
    n = len(A)
    if n == 0:
        return ""
    m = len(A[0])
    # 先把所有数字转成字符串，确定列宽
    str_cols = []
    for i in range(n):
        row_str = []
        for j in range(m):
            row_str.append(f"{A[i][j]:.{precision}g}")
        row_str.append(f"{b[i]:.{precision}g}")
        str_cols.append(row_str)

    col_widths = [0] * (m + 1)
    for row in str_cols:
        for j, s in enumerate(row):
            col_widths[j] = max(col_widths[j], len(s))

    lines = []
    for i in range(n):
        left = "  ".join(
            s.rjust(col_widths[j]) for j, s in enumerate(str_cols[i][:-1])
        )
        right = str_cols[i][-1].rjust(col_widths[-1])
        lines.append(f"[ {left} | {right} ]")
    return "\n".join(lines)


def format_matrix(M, name="M", precision=4):
    """把矩阵 M 格式化成字符串。"""
    n = len(M)
    if n == 0:
        return f"{name} = []"

    m = len(M[0])
    str_cols = []
    for i in range(n):
        row_str = []
        for j in range(m):
            row_str.append(f"{M[i][j]:.{precision}g}")
        str_cols.append(row_str)

    col_widths = [0] * m
    for row in str_cols:
        for j, s in enumerate(row):
            col_widths[j] = max(col_widths[j], len(s))

    lines = [f"{name} ="]
    for i in range(n):
        line = "  ".join(
            s.rjust(col_widths[j]) for j, s in enumerate(str_cols[i])
        )
        lines.append("  [ " + line + " ]")
    return "\n".join(lines)


def format_vector(v, name="x", precision=6):
    """向量打印。"""
    comp = ",  ".join(f"{val:.{precision}g}" for val in v)
    return f"{name} = [ {comp} ]^T"


def build_equation_strings(A, b, precision=4):
    """
    把 Ax = b 格式化成“1x1 + 2x2 - 3x3 = 1”这种形式。
    返回字符串列表（每方程一行）。
    """
    lines = []
    n = len(A)
    for i in range(n):
        terms = []
        for j, a in enumerate(A[i]):
            var = f"x{j+1}"
            if abs(a) < 1e-12:
                continue
            # 处理符号和系数显示
            if not terms:
                # 第一项
                coef_str = f"{a:.{precision}g}"
                terms.append(f"{coef_str}{var}")
            else:
                sign = "+" if a >= 0 else "-"
                coef_str = f"{abs(a):.{precision}g}"
                terms.append(f" {sign} {coef_str}{var}")
        if not terms:
            left = "0"
        else:
            left = "".join(terms)
        right = f"{b[i]:.{precision}g}"
        lines.append(f"{left} = {right}")
    return lines


# ==================== 高斯消元（列主元） ====================

def gaussian_partial_pivot(A, b):
    """
    列主元高斯消元，返回 (x, detail_text)。
    A, b 会被拷贝，不修改原数据。
    """
    n = len(A)
    # 深拷贝
    a = [row[:] for row in A]
    rhs = b[:]
    lines = []
    eps = 1e-14

    lines.append("【方法 1：列主元高斯消元法（Doolittle 视角）】")
    lines.append("")
    lines.append("初始增广矩阵 [A|b]：")
    lines.append(format_augmented_matrix(a, rhs))
    lines.append("")

    # 前向消元
    for k in range(n - 1):
        # 选列主元：在第 k 列，k..n-1 行中找绝对值最大的元素
        pivot_row = max(range(k, n), key=lambda i: abs(a[i][k]))
        if abs(a[pivot_row][k]) < eps:
            raise SolverError("矩阵接近奇异，无法继续高斯消元。")

        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            rhs[k], rhs[pivot_row] = rhs[pivot_row], rhs[k]
            lines.append(
                f"第 {k+1} 步：交换第 {k+1} 行和第 {pivot_row+1} 行（列主元选择）。"
            )
            lines.append(format_augmented_matrix(a, rhs))
            lines.append("")

        # 对 k+1..n-1 行进行消元
        for i in range(k + 1, n):
            factor = a[i][k] / a[k][k]
            a[i][k] = 0.0
            for j in range(k + 1, n):
                a[i][j] -= factor * a[k][j]
            rhs[i] -= factor * rhs[k]
        lines.append(f"第 {k+1} 步：用第 {k+1} 行消去第 {k+2}~{n} 行的 x{k+1}。")
        lines.append(format_augmented_matrix(a, rhs))
        lines.append("")

    # 回代
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        if abs(a[i][i]) < eps:
            raise SolverError("回代时遇到零主元，线性方程组可能无唯一解。")
        s = sum(a[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (rhs[i] - s) / a[i][i]

    lines.append("回代过程得到解向量：")
    lines.append(format_vector(x, name="x"))
    detail = "\n".join(lines)
    return x, detail


# ==================== 高斯消元（全主元：主元素法） ====================

def gaussian_full_pivot(A, b):
    """
    全主元高斯消元（主元素法），返回 (x, detail_text)。
    需要跟踪列交换，最终还原变量顺序。
    """
    n = len(A)
    a = [row[:] for row in A]
    rhs = b[:]
    col_perm = list(range(n))  # 记录列交换顺序
    lines = []
    eps = 1e-14

    lines.append("【方法 2：全主元高斯消元法（主元素法）】")
    lines.append("")
    lines.append("初始增广矩阵 [A|b]：")
    lines.append(format_augmented_matrix(a, rhs))
    lines.append("")
    lines.append("说明：每一步在未消元子矩阵中寻找绝对值最大的元素作为“主元素”，")
    lines.append("      允许做行交换和列交换，以提高数值稳定性。")
    lines.append("")

    for k in range(n - 1):
        # 在 a[i][j], i,j>=k 中寻找最大元素
        max_val = 0.0
        pivot_row = k
        pivot_col = k
        for i in range(k, n):
            for j in range(k, n):
                if abs(a[i][j]) > max_val:
                    max_val = abs(a[i][j])
                    pivot_row = i
                    pivot_col = j
        if max_val < eps:
            raise SolverError("矩阵接近奇异，无法继续全主元高斯消元。")

        # 行交换
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            rhs[k], rhs[pivot_row] = rhs[pivot_row], rhs[k]
            lines.append(f"第 {k+1} 步：交换第 {k+1} 行和第 {pivot_row+1} 行。")

        # 列交换
        if pivot_col != k:
            for i in range(n):
                a[i][k], a[i][pivot_col] = a[i][pivot_col], a[i][k]
            col_perm[k], col_perm[pivot_col] = col_perm[pivot_col], col_perm[k]
            lines.append(f"第 {k+1} 步：交换第 {k+1} 列和第 {pivot_col+1} 列（变量交换）。")

        lines.append("当前增广矩阵：")
        lines.append(format_augmented_matrix(a, rhs))
        lines.append("")

        # 消元
        for i in range(k + 1, n):
            factor = a[i][k] / a[k][k]
            a[i][k] = 0.0
            for j in range(k + 1, n):
                a[i][j] -= factor * a[k][j]
            rhs[i] -= factor * rhs[k]
        lines.append(f"用第 {k+1} 行消去第 {k+2}~{n} 行的第 {k+1} 列。")
        lines.append(format_augmented_matrix(a, rhs))
        lines.append("")

    # 回代，求的是“已重排变量”的解
    x_perm = [0.0] * n
    for i in range(n - 1, -1, -1):
        if abs(a[i][i]) < eps:
            raise SolverError("回代时遇到零主元，线性方程组可能无唯一解。")
        s = sum(a[i][j] * x_perm[j] for j in range(i + 1, n))
        x_perm[i] = (rhs[i] - s) / a[i][i]

    # 还原到原变量顺序
    x = [0.0] * n
    for j in range(n):
        x[col_perm[j]] = x_perm[j]

    lines.append("回代过程得到“重排变量”的解向量：")
    lines.append(format_vector(x_perm, name="x(重排)"))
    lines.append("")
    lines.append("根据列交换记录，还原到原变量顺序：")
    lines.append(f"变量重排顺序（索引从 1 开始）：{[i+1 for i in col_perm]}")
    lines.append(format_vector(x, name="x(原变量顺序)"))

    detail = "\n".join(lines)
    return x, detail


# ==================== 克劳特分解 ====================

def crout_lu(A):
    """
    克劳特分解：A = L U，其中 U 对角线为 1。
    返回 (L, U)。
    """
    n = len(A)
    eps = 1e-14
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]
    # U 对角线为 1
    for i in range(n):
        U[i][i] = 1.0

    for j in range(n):
        # 先算 L 的第 j 列
        for i in range(j, n):
            s = sum(L[i][k] * U[k][j] for k in range(j))
            L[i][j] = A[i][j] - s
        if abs(L[j][j]) < eps:
            raise SolverError("克劳特分解过程中出现零主元，分解失败。")
        # 再算 U 的第 j 行（j+1..n-1）
        for k in range(j + 1, n):
            s = sum(L[j][m] * U[m][k] for m in range(j))
            U[j][k] = (A[j][k] - s) / L[j][j]

    return L, U


def forward_substitution(L, b):
    """前代：解 L y = b（L 为下三角，一般对角）。"""
    n = len(L)
    y = [0.0] * n
    eps = 1e-14
    for i in range(n):
        if abs(L[i][i]) < eps:
            raise SolverError("前代时遇到零对角元。")
        s = sum(L[i][j] * y[j] for j in range(i))
        y[i] = (b[i] - s) / L[i][i]
    return y


def back_substitution(U, y):
    """回代：解 U x = y（U 为上三角）。"""
    n = len(U)
    x = [0.0] * n
    eps = 1e-14
    for i in range(n - 1, -1, -1):
        if abs(U[i][i]) < eps:
            raise SolverError("回代时遇到零对角元。")
        s = sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (y[i] - s) / U[i][i]
    return x


def crout_solver(A, b):
    """克劳特 LU 分解 + 前代 / 回代。返回 (x, detail_text)。"""
    lines = []
    lines.append("【方法 3：克劳特分解法（Crout LU）】")
    lines.append("")
    lines.append("系数矩阵 A：")
    lines.append(format_matrix(A, name="A"))
    lines.append("")
    lines.append("右端向量 b：")
    lines.append(format_vector(b, name="b"))
    lines.append("")

    L, U = crout_lu(A)

    lines.append("分解得到 A = L · U：")
    lines.append(format_matrix(L, name="L"))
    lines.append("")
    lines.append(format_matrix(U, name="U"))
    lines.append("")

    y = forward_substitution(L, b)
    x = back_substitution(U, y)

    lines.append("先解 L·y = b，得到中间变量：")
    lines.append(format_vector(y, name="y"))
    lines.append("")
    lines.append("再解 U·x = y，得到解向量：")
    lines.append(format_vector(x, name="x"))

    detail = "\n".join(lines)
    return x, detail


# ==================== 平方根法（Cholesky） ====================

def is_symmetric(A, tol=1e-12):
    """简易判断对称性。"""
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j] - A[j][i]) > tol:
                return False
    return True


def cholesky(A):
    """
    Cholesky 分解：A 必须为对称正定。
    返回下三角矩阵 L，使 A = L L^T。
    """
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    eps = 1e-14

    for i in range(n):
        for j in range(i + 1):
            s = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                val = A[i][i] - s
                if val <= eps:
                    raise SolverError("矩阵不是正定矩阵，平方根法失败。")
                L[i][j] = math.sqrt(val)
            else:
                if abs(L[j][j]) < eps:
                    raise SolverError("平方根法过程中遇到零对角元。")
                L[i][j] = (A[i][j] - s) / L[j][j]
    return L


def cholesky_solver(A, b):
    """平方根法（Cholesky）：适用于对称正定矩阵。"""
    if not is_symmetric(A):
        raise SolverError("系数矩阵不是对称矩阵，不能使用平方根法。")

    lines = []
    lines.append("【方法 4：平方根法（Cholesky 分解，适用于对称正定矩阵）】")
    lines.append("")
    lines.append("系数矩阵 A：")
    lines.append(format_matrix(A, name="A"))
    lines.append("")
    lines.append("右端向量 b：")
    lines.append(format_vector(b, name="b"))
    lines.append("")

    L = cholesky(A)
    LT = [[L[j][i] for j in range(len(L))] for i in range(len(L))]

    lines.append("分解得到 A = L · L^T：")
    lines.append(format_matrix(L, name="L"))
    lines.append("")
    lines.append(format_matrix(LT, name="L^T"))
    lines.append("")

    # 先解 L·y = b，再解 L^T·x = y
    y = forward_substitution(L, b)
    x = back_substitution(LT, y)

    lines.append("先解 L·y = b，得到中间变量：")
    lines.append(format_vector(y, name="y"))
    lines.append("")
    lines.append("再解 L^T·x = y，得到解向量：")
    lines.append(format_vector(x, name="x"))

    detail = "\n".join(lines)
    return x, detail


# ==================== 追赶法（Thomas 算法，三对角） ====================

def is_tridiagonal(A, tol=1e-12):
    """判断是否为三对角矩阵。"""
    n = len(A)
    for i in range(n):
        for j in range(n):
            if j == i or j == i - 1 or j == i + 1:
                continue
            if abs(A[i][j]) > tol:
                return False
    return True


def thomas_solver(A, b):
    """
    追赶法（Thomas 算法），适用于三对角线性方程组。
    返回 (x, detail_text)。
    """
    if not is_tridiagonal(A):
        raise SolverError("该系数矩阵不是三对角矩阵，不能使用追赶法。")

    n = len(A)
    # 提取三条对角线
    a = [0.0] * n  # 下对角线，a[0] 无用
    d = [0.0] * n  # 主对角线
    c = [0.0] * n  # 上对角线，c[n-1] 无用
    for i in range(n):
        d[i] = A[i][i]
        if i > 0:
            a[i] = A[i][i - 1]
        if i < n - 1:
            c[i] = A[i][i + 1]

    rhs = b[:]
    lines = []
    eps = 1e-14

    lines.append("【方法 5：追赶法（Thomas 算法，三对角线性方程组）】")
    lines.append("")
    lines.append("原始三对角矩阵的三条对角线：")
    lines.append(format_vector(a, name="下对角 a"))
    lines.append(format_vector(d, name="主对角 d"))
    lines.append(format_vector(c, name="上对角 c"))
    lines.append(format_vector(rhs, name="右端项 f"))
    lines.append("")

    # 前向“追赶”（修正系数）
    for i in range(1, n):
        if abs(d[i - 1]) < eps:
            raise SolverError("追赶法前向消去中出现零主对角元。")
        w = a[i] / d[i - 1]
        d[i] = d[i] - w * c[i - 1]
        rhs[i] = rhs[i] - w * rhs[i - 1]
        lines.append(f"第 {i} 步前向追赶：")
        lines.append(f"  w = a[{i}] / d[{i-1}] = {w:g}")
        lines.append(format_vector(d, name="当前主对角 d"))
        lines.append(format_vector(rhs, name="当前右端项 f"))
        lines.append("")

    # 回代
    x = [0.0] * n
    if abs(d[n - 1]) < eps:
        raise SolverError("追赶法回代中出现零主对角元。")
    x[n - 1] = rhs[n - 1] / d[n - 1]
    lines.append("回代：")
    lines.append(f"x_{n} = f_{n} / d_{n} = {x[n - 1]:g}")

    for i in range(n - 2, -1, -1):
        if abs(d[i]) < eps:
            raise SolverError("追赶法回代中出现零主对角元。")
        x[i] = (rhs[i] - c[i] * x[i + 1]) / d[i]
        lines.append(
            f"x_{i+1} = (f_{i+1} - c_{i+1} * x_{i+2}) / d_{i+1} = {x[i]:g}"
        )

    lines.append("")
    lines.append("最终解向量：")
    lines.append(format_vector(x, name="x"))

    detail = "\n".join(lines)
    return x, detail


# ==================== 结构分析：自动推荐方法 ====================

def analyze_matrix(A):
    """
    根据系数矩阵的结构特征，给出一个求解方法建议：
    - 三对角：追赶法
    - 近似对称：平方根法
    - 其他：列主元高斯消元
    返回 (method_key, description)
    """
    try:
        if is_tridiagonal(A):
            return (
                "thomas",
                "检测到系数矩阵为三对角矩阵，追赶法（Thomas 算法）可在 O(n) 时间内高效求解。"
            )
        if is_symmetric(A):
            return (
                "cholesky",
                "检测到系数矩阵近似对称，若进一步正定则适合平方根法（Cholesky），运算量约为普通 LU 的一半。"
            )
    except Exception:
        # 保底处理，不影响主流程
        pass

    return (
        "gauss",
        "未检测到明显结构特征，默认推荐列主元高斯消元法，一般情况下具有较好的数值稳定性。"
    )


# ==================== 统一调度函数 ====================

def solve_linear_system(A, b, method):
    """
    统一对外接口。
    method:
      - "gauss"        : 列主元高斯消元法
      - "gauss_full"   : 全主元高斯消元法（主元素法）
      - "crout"        : 克劳特分解法
      - "cholesky"     : 平方根法（对称正定矩阵）
      - "thomas"       : 追赶法（三对角）
    返回 (x, detail_text)
    """
    if method == "gauss":
        return gaussian_partial_pivot(A, b)
    elif method == "gauss_full":
        return gaussian_full_pivot(A, b)
    elif method == "crout":
        return crout_solver(A, b)
    elif method == "cholesky":
        return cholesky_solver(A, b)
    elif method == "thomas":
        return thomas_solver(A, b)
    else:
        raise SolverError(f"未知求解方法：{method}")
