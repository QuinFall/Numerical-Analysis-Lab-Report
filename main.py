# main.py
# 解线性方程组的直接法的编程实现 —— 王万淇（Quin）

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import csv

from solvers import (
    parse_matrix,
    build_equation_strings,  # 仍然保留，用于导出等场景（内部计算）
    solve_linear_system,
    analyze_matrix,
    SolverError,
)

# 颜色和样式
BG_MAIN = "#f5f7ff"
BG_PANEL = "#ffffff"
FG_TITLE = "#222244"
FG_TEXT = "#333333"
FG_HINT = "#666666"


# ---------- 辅助：把数字变成下标，生成更美观的方程字符串 ----------

SUBSCRIPT_DIGITS = {
    "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄",
    "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉",
}


def int_to_subscript(n: int) -> str:
    """把整数 n 转成下标数字，比如 12 -> '₁₂'。"""
    s = str(n)
    return "".join(SUBSCRIPT_DIGITS.get(ch, ch) for ch in s)


def format_equations_pretty(A, b, precision=4):
    """
    把 Ax = b 格式化成更美观的形式，例如：
    1.25·x₁ + 2.0·x₂ - 3·x₃ = 5
    区分系数和变量，用 '·' 连接，变量下标用 Unicode 下标。
    """
    lines = []
    n = len(A)
    for i in range(n):
        terms = []
        for j, a in enumerate(A[i]):
            if abs(a) < 1e-12:
                continue
            var = "x" + int_to_subscript(j + 1)
            coef = f"{abs(a):.{precision}g}"
            # 第一项单独处理符号
            if not terms:
                if a < 0:
                    terms.append(f"- {coef}·{var}")
                else:
                    terms.append(f"{coef}·{var}")
            else:
                sign = "+" if a >= 0 else "-"
                terms.append(f" {sign} {coef}·{var}")
        left = "0" if not terms else "".join(terms)
        right = f"{b[i]:.{precision}g}"
        lines.append(f"{left} = {right}")
    return "\n".join(lines)


class NAApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("解线性方程组的直接法 —— 上机实验")
        self.geometry("900x600")
        self.resizable(False, False)
        self.configure(bg=BG_MAIN)

        # 当前保存的输入/矩阵
        self.current_entries = None  # 原始字符串 entries
        self.current_A = None
        self.current_b = None

        # 最近一次求解的结果，用于导出
        self.last_detail_text = ""
        self.last_solution_x = None

        # 容器 frame，避免反复销毁窗口本身
        self.main_frame = tk.Frame(self, bg=BG_MAIN)
        self.main_frame.pack(fill="both", expand=True)

        # 初始界面
        self.show_welcome()

    # -------------- 公共方法：清空并切换页面 ----------------

    def clear_frame(self):
        for widget in self.main_frame.winfo_children():
            widget.destroy()

    # -------------- 页面 1：欢迎界面 ----------------

    def show_welcome(self):
        self.clear_frame()
        frame = self.main_frame

        title = tk.Label(
            frame,
            text="解线性方程组的直接法的编程实现\n—— 王万淇 2025年11月18日",
            font=("Microsoft YaHei", 22, "bold"),
            justify="center",
            bg=BG_MAIN,
            fg=FG_TITLE,
        )
        title.pack(pady=80)

        tip = tk.Label(
            frame,
            text="提示：点击键盘任意键或下方按钮开始实验",
            font=("Microsoft YaHei", 12),
            bg=BG_MAIN,
            fg=FG_HINT,
        )
        tip.pack(pady=10)

        btn_start = tk.Button(
            frame,
            text="开始实验",
            font=("Microsoft YaHei", 12, "bold"),
            width=15,
            command=self.show_dimension_input,
            bg="#4f7cff",
            fg="white",
            activebackground="#3c63d6",
            activeforeground="white",
            relief="flat",
        )
        btn_start.pack(pady=40)

        # 绑定任意键启动
        self.bind("<Key>", lambda event: self.show_dimension_input())

    # -------------- 页面 2：选择未知数个数 ----------------

    def show_dimension_input(self):
        # 取消欢迎界面的按键绑定
        self.unbind("<Key>")
        self.clear_frame()
        frame = self.main_frame

        # 左上角返回主页面按钮
        btn_home = tk.Button(
            frame,
            text="← 返回主页面",
            font=("Microsoft YaHei", 9),
            command=self.show_welcome,
            bg=BG_MAIN,
            fg=FG_TEXT,
            relief="flat",
        )
        btn_home.pack(anchor="nw", padx=10, pady=8)

        title = tk.Label(
            frame,
            text="步骤一：选择未知数（方程）个数 n",
            font=("Microsoft YaHei", 16, "bold"),
            bg=BG_MAIN,
            fg=FG_TITLE,
        )
        title.pack(pady=(10, 5))

        inner = tk.Frame(frame, bg=BG_MAIN)
        inner.pack(pady=10)

        lab = tk.Label(
            inner,
            text="请选择未知数个数 n（2 ~ 10）：",
            font=("Microsoft YaHei", 12),
            bg=BG_MAIN,
            fg=FG_TEXT,
        )
        lab.grid(row=0, column=0, padx=5, pady=5)

        if not hasattr(self, "var_n"):
            self.var_n = tk.StringVar(value="4")
        spin = ttk.Spinbox(
            inner,
            from_=2,
            to=10,
            textvariable=self.var_n,
            width=5,
            font=("Microsoft YaHei", 12),
        )
        spin.grid(row=0, column=1, padx=5, pady=5)

        btn_next = tk.Button(
            inner,
            text="生成方程输入模板",
            font=("Microsoft YaHei", 11, "bold"),
            command=self.show_equation_input,
            bg="#4f7cff",
            fg="white",
            activebackground="#3c63d6",
            activeforeground="white",
            relief="flat",
        )
        btn_next.grid(row=0, column=2, padx=10, pady=5)

    # -------------- 页面 3：输入方程系数 ----------------

    def show_equation_input(self):
        self.clear_frame()
        frame = self.main_frame
        frame.configure(bg=BG_MAIN)

        try:
            n = int(self.var_n.get())
        except Exception:
            messagebox.showerror("输入错误", "请先选择合法的未知数个数 n。")
            self.show_dimension_input()
            return

        # 左上角返回主页面
        btn_home = tk.Button(
            frame,
            text="← 返回主页面",
            font=("Microsoft YaHei", 9),
            command=self.show_welcome,
            bg=BG_MAIN,
            fg=FG_TEXT,
            relief="flat",
        )
        btn_home.pack(anchor="nw", padx=10, pady=8)

        title = tk.Label(
            frame,
            text="步骤二：输入线性方程组的系数与常数项",
            font=("Microsoft YaHei", 16, "bold"),
            bg=BG_MAIN,
            fg=FG_TITLE,
        )
        title.pack(pady=(5, 5))

        tip = tk.Label(
            frame,
            text="说明：在括号中输入系数，支持小数；未填视为 0。",
            font=("Microsoft YaHei", 11),
            bg=BG_MAIN,
            fg=FG_HINT,
        )
        tip.pack(pady=(0, 5))

        # 中间：滚动区域，容纳所有输入框（支持上下 + 左右滚动）
        work_area = tk.Frame(frame, bg=BG_MAIN)
        work_area.pack(fill="both", expand=True, padx=20, pady=5)

        canvas = tk.Canvas(work_area, bg=BG_PANEL, highlightthickness=0)
        scroll_y = tk.Scrollbar(work_area, orient="vertical", command=canvas.yview)
        scroll_x = tk.Scrollbar(work_area, orient="horizontal", command=canvas.xview)

        canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)

        canvas.grid(row=0, column=0, sticky="nsew")
        scroll_y.grid(row=0, column=1, sticky="ns")
        scroll_x.grid(row=1, column=0, sticky="ew")

        work_area.rowconfigure(0, weight=1)
        work_area.columnconfigure(0, weight=1)

        container = tk.Frame(canvas, bg=BG_PANEL)
        canvas.create_window((0, 0), window=container, anchor="nw")

        def _on_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        container.bind("<Configure>", _on_configure)

        # 生成 n 行 输入：()x1 + ()x2 + ... = ()
        self.entry_grid = []
        for i in range(n):
            row_frame = tk.Frame(container, bg=BG_PANEL)
            row_frame.pack(anchor="w", pady=4, padx=10)

            row_entries = []
            for j in range(n):
                e = tk.Entry(
                    row_frame,
                    width=6,
                    font=("Consolas", 11),
                    justify="right",
                    relief="solid",
                    bd=1,
                )
                e.insert(0, "")  # 空 = 0
                e.pack(side="left", padx=1)
                row_entries.append(e)

                lbl_var = tk.Label(
                    row_frame,
                    text=f"x{j+1}",
                    font=("Microsoft YaHei", 11),
                    bg=BG_PANEL,
                    fg=FG_TEXT,
                )
                lbl_var.pack(side="left", padx=(0, 2))

                if j < n - 1:
                    lbl_plus = tk.Label(
                        row_frame,
                        text=" + ",
                        font=("Microsoft YaHei", 11),
                        bg=BG_PANEL,
                        fg=FG_TEXT,
                    )
                    lbl_plus.pack(side="left")

            lbl_eq = tk.Label(
                row_frame,
                text=" = ",
                font=("Microsoft YaHei", 11),
                bg=BG_PANEL,
                fg=FG_TEXT,
            )
            lbl_eq.pack(side="left", padx=2)

            e_rhs = tk.Entry(
                row_frame,
                width=6,
                font=("Consolas", 11),
                justify="right",
                relief="solid",
                bd=1,
            )
            e_rhs.insert(0, "")
            e_rhs.pack(side="left", padx=1)
            row_entries.append(e_rhs)

            self.entry_grid.append(row_entries)

        # 底部按钮栏：居中放两个按钮
        btn_bar = tk.Frame(frame, bg=BG_MAIN)
        btn_bar.pack(side="bottom", pady=12)

        btn_done = tk.Button(
            btn_bar,
            text="数据输入完毕 → 检查并展示方程组",
            font=("Microsoft YaHei", 11, "bold"),
            command=self.check_and_show_methods,
            bg="#4f7cff",
            fg="white",
            activebackground="#3c63d6",
            activeforeground="white",
            relief="flat",
            padx=10,
            pady=3,
        )
        btn_done.grid(row=0, column=0, padx=8, pady=5)

        btn_back = tk.Button(
            btn_bar,
            text="返回上一步（重新选择 n）",
            font=("Microsoft YaHei", 10),
            command=self.show_dimension_input,
            bg="#ddddff",
            fg=FG_TEXT,
            activebackground="#c0c8ff",
            relief="flat",
            padx=8,
            pady=3,
        )
        btn_back.grid(row=0, column=1, padx=8, pady=5)

    # -------------- 检查线性性并展示方法选择 ----------------

    def check_and_show_methods(self):
        # 从 entry_grid 读取所有字符串
        raw_entries = []
        for row_entries in self.entry_grid:
            row_strs = []
            for e in row_entries:
                text = e.get().strip()
                if text == "":
                    text = "0"
                row_strs.append(text)
            raw_entries.append(row_strs)

        self.current_entries = raw_entries

        # 尝试解析为 A, b
        try:
            A, b = parse_matrix(raw_entries)
        except SolverError as e:
            messagebox.showerror(
                "输入错误",
                f"检测为非线性或非法输入：\n{e}\n\n将返回到系数输入界面。",
            )
            self.show_equation_input()
            return

        self.current_A = A
        self.current_b = b
        self.show_method_selection()

    # -------------- 页面 3：展示方程组 & 方法选择 ----------------

    def show_method_selection(self):
        self.clear_frame()
        frame = self.main_frame
        frame.configure(bg=BG_MAIN)

        # 左上角返回主页面
        btn_home = tk.Button(
            frame,
            text="← 返回主页面",
            font=("Microsoft YaHei", 9),
            command=self.show_welcome,
            bg=BG_MAIN,
            fg=FG_TEXT,
            relief="flat",
        )
        btn_home.pack(anchor="nw", padx=10, pady=8)

        title = tk.Label(
            frame,
            text="步骤三：确认线性方程组，并选择求解方法",
            font=("Microsoft YaHei", 16, "bold"),
            bg=BG_MAIN,
            fg=FG_TITLE,
        )
        title.pack(pady=(5, 10))

        # 显示方程组
        panel = tk.Frame(frame, bg=BG_PANEL, bd=1, relief="solid")
        panel.pack(pady=5, fill="x", padx=20)

        eq_frame = tk.Frame(panel, bg=BG_PANEL)
        eq_frame.pack(padx=10, pady=8, fill="x")

        eq_title = tk.Label(
            eq_frame,
            text="当前输入的线性方程组：",
            font=("Microsoft YaHei", 12, "bold"),
            anchor="w",
            bg=BG_PANEL,
            fg=FG_TEXT,
        )
        eq_title.pack(anchor="w")

        # 使用更美观的下标形式展示方程
        eq_text = format_equations_pretty(self.current_A, self.current_b)

        lbl_eq = tk.Label(
            eq_frame,
            text=eq_text,
            font=("Consolas", 11),
            justify="left",
            anchor="w",
            bg=BG_PANEL,
            fg=FG_TEXT,
        )
        lbl_eq.pack(anchor="w", pady=5)

        # 自动分析矩阵结构，推荐求解方法（用红色标亮）
        method_display = {
            "gauss": ("1", "列主元高斯消元法"),
            "gauss_full": ("2", "全主元高斯消元法"),
            "crout": ("3", "克劳特分解法"),
            "cholesky": ("4", "平方根法（Cholesky）"),
            "thomas": ("5", "追赶法（Thomas）"),
        }
        rec_key, rec_desc = analyze_matrix(self.current_A)
        self.recommended_method = rec_key
        rec_no, rec_name = method_display.get(rec_key, ("1", "列主元高斯消元法"))

        rec_label = tk.Label(
            eq_frame,
            text=f"系统分析结果：推荐优先使用 {rec_no}、{rec_name}\n理由：{rec_desc}",
            font=("Microsoft YaHei", 10, "bold"),
            justify="left",
            fg="red",
            bg=BG_PANEL,
            anchor="w",
        )
        rec_label.pack(anchor="w", pady=(0, 4))

        # 方法选择说明
        method_frame = tk.Frame(panel, bg=BG_PANEL)
        method_frame.pack(padx=10, pady=(0, 8), fill="x")

        lab_tip = tk.Label(
            method_frame,
            text="输入对应方法前的数字，来选择解以上线性方程组的方法：",
            font=("Microsoft YaHei", 12),
            anchor="w",
            justify="left",
            bg=BG_PANEL,
            fg=FG_TEXT,
        )
        lab_tip.pack(anchor="w", pady=5)

        methods_text = (
            "1、列主元高斯消元法（Gaussian Elimination with Partial Pivoting）\n"
            "2、全主元高斯消元法（主元素法，Full Pivoting）\n"
            "3、克劳特分解法（Crout LU 分解）\n"
            "4、平方根法（Cholesky 分解，适用于对称正定矩阵）\n"
            "5、追赶法（Thomas 算法，仅适用于三对角线性方程组）"
        )

        lab_methods = tk.Label(
            method_frame,
            text=methods_text,
            font=("Microsoft YaHei", 11),
            justify="left",
            anchor="w",
            bg=BG_PANEL,
            fg=FG_TEXT,
        )
        lab_methods.pack(anchor="w")

        choice_frame = tk.Frame(method_frame, bg=BG_PANEL)
        choice_frame.pack(anchor="w", pady=8)

        lab_choice = tk.Label(
            choice_frame,
            text="请输入方法编号（1~5）：",
            font=("Microsoft YaHei", 11),
            bg=BG_PANEL,
            fg=FG_TEXT,
        )
        lab_choice.pack(side="left", padx=5)

        self.var_method_choice = tk.StringVar(value=rec_no)
        entry_choice = tk.Entry(
            choice_frame,
            textvariable=self.var_method_choice,
            width=4,
            font=("Consolas", 11),
            justify="center",
            relief="solid",
            bd=1,
        )
        entry_choice.pack(side="left", padx=5)

        btn_solve = tk.Button(
            choice_frame,
            text="开始求解",
            font=("Microsoft YaHei", 11, "bold"),
            command=self.solve_with_selected_method,
            bg="#4f7cff",
            fg="white",
            activebackground="#3c63d6",
            activeforeground="white",
            relief="flat",
            padx=10,
            pady=3,
        )
        btn_solve.pack(side="left", padx=10)

        # 底部返回按钮
        btn_back = tk.Button(
            frame,
            text="返回重新输入方程",
            font=("Microsoft YaHei", 10),
            command=self.show_equation_input,
            bg="#ddddff",
            fg=FG_TEXT,
            activebackground="#c0c8ff",
            relief="flat",
        )
        btn_back.pack(pady=5)

    # -------------- 调用对应数值方法求解 ----------------

    def solve_with_selected_method(self):
        choice = self.var_method_choice.get().strip()
        method_map = {
            "1": "gauss",
            "2": "gauss_full",
            "3": "crout",
            "4": "cholesky",
            "5": "thomas",
        }
        if choice not in method_map:
            messagebox.showerror("输入错误", "请选择 1~5 之间的合法编号。")
            return
        method = method_map[choice]

        if self.current_A is None or self.current_b is None:
            messagebox.showerror("错误", "尚未获得有效的线性方程组。")
            return

        try:
            x, detail = solve_linear_system(self.current_A, self.current_b, method)
        except SolverError as e:
            messagebox.showerror("求解失败", str(e))
            return
        except Exception as e:
            messagebox.showerror("异常错误", f"求解过程中出现未知错误：\n{e}")
            return

        # 记录最近一次结果，便于导出
        self.last_detail_text = detail
        self.last_solution_x = x

        self.show_result()

    # -------------- 页面 4：展示求解过程与结果 ----------------

    def show_result(self):
        self.clear_frame()
        frame = self.main_frame
        frame.configure(bg=BG_MAIN)

        # 顶部栏：左上返回，右上导出按钮
        top_bar = tk.Frame(frame, bg=BG_MAIN)
        top_bar.pack(fill="x", pady=(5, 0))

        btn_home = tk.Button(
            top_bar,
            text="← 返回主页面",
            font=("Microsoft YaHei", 9),
            command=self.show_welcome,
            bg=BG_MAIN,
            fg=FG_TEXT,
            relief="flat",
        )
        btn_home.pack(side="left", padx=10)

        export_frame = tk.Frame(top_bar, bg=BG_MAIN)
        export_frame.pack(side="right", padx=10)

        btn_export_txt = tk.Button(
            export_frame,
            text="导出求解过程 TXT",
            font=("Microsoft YaHei", 9),
            command=self.export_txt,
            bg="#4f7cff",
            fg="white",
            activebackground="#3c63d6",
            activeforeground="white",
            relief="flat",
            padx=6,
            pady=2,
        )
        btn_export_txt.pack(side="left", padx=4)

        btn_export_csv = tk.Button(
            export_frame,
            text="导出解向量 CSV",
            font=("Microsoft YaHei", 9),
            command=self.export_csv,
            bg="#4f7cff",
            fg="white",
            activebackground="#3c63d6",
            activeforeground="white",
            relief="flat",
            padx=6,
            pady=2,
        )
        btn_export_csv.pack(side="left", padx=4)

        title = tk.Label(
            frame,
            text="步骤四：求解过程与最终结果展示",
            font=("Microsoft YaHei", 16, "bold"),
            bg=BG_MAIN,
            fg=FG_TITLE,
        )
        title.pack(pady=(0, 5))

        # 提示：标红标亮
        hint = tk.Label(
            frame,
            text="计算完成！请查看下方详细步骤，若需要重新实验，请点击左上角“返回主页面”按钮。",
            font=("Microsoft YaHei", 11, "bold"),
            bg=BG_MAIN,
            fg="red",
            wraplength=820,
            justify="center",
        )
        hint.pack(pady=(0, 5))

        # 再次展示方程组（用下标形式）
        eq_panel = tk.Frame(frame, bg=BG_PANEL, bd=1, relief="solid")
        eq_panel.pack(pady=5, fill="x", padx=20)

        eq_frame = tk.Frame(eq_panel, bg=BG_PANEL)
        eq_frame.pack(padx=10, pady=8, fill="x")

        eq_title = tk.Label(
            eq_frame,
            text="原线性方程组：",
            font=("Microsoft YaHei", 12, "bold"),
            anchor="w",
            bg=BG_PANEL,
            fg=FG_TEXT,
        )
        eq_title.pack(anchor="w")

        eq_text = format_equations_pretty(self.current_A, self.current_b)
        lbl_eq = tk.Label(
            eq_frame,
            text=eq_text,
            font=("Consolas", 11),
            justify="left",
            anchor="w",
            bg=BG_PANEL,
            fg=FG_TEXT,
        )
        lbl_eq.pack(anchor="w", pady=5)

        # 中央：大文本框，显示算法步骤 / 矩阵 / 解
        text_panel = tk.Frame(frame, bg=BG_MAIN)
        text_panel.pack(fill="both", expand=True, padx=20, pady=8)

        txt = tk.Text(
            text_panel,
            wrap="none",
            font=("Consolas", 11),
        )
        scroll_y = tk.Scrollbar(
            text_panel, orient="vertical", command=txt.yview
        )
        scroll_x = tk.Scrollbar(
            text_panel, orient="horizontal", command=txt.xview
        )
        txt.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)

        txt.grid(row=0, column=0, sticky="nsew")
        scroll_y.grid(row=0, column=1, sticky="ns")
        scroll_x.grid(row=1, column=0, sticky="ew")

        text_panel.rowconfigure(0, weight=1)
        text_panel.columnconfigure(0, weight=1)

        detail_text = self.last_detail_text or ""
        txt.insert("1.0", detail_text)
        txt.configure(state="disabled")  # 只读显示

        # 底部：继续使用当前 n
        btn_frame = tk.Frame(frame, bg=BG_MAIN)
        btn_frame.pack(pady=8)

        btn_again = tk.Button(
            btn_frame,
            text="继续使用当前 n 重新输入方程",
            font=("Microsoft YaHei", 10),
            command=self.show_equation_input,
            bg="#ddddff",
            fg=FG_TEXT,
            activebackground="#c0c8ff",
            relief="flat",
            padx=8,
            pady=3,
        )
        btn_again.pack(padx=10, pady=5)

    # -------------- 导出功能 ----------------

    def export_txt(self):
        """导出求解过程为 TXT 文件。"""
        if not self.last_detail_text:
            messagebox.showwarning("无内容", "当前没有可导出的求解过程。")
            return

        path = filedialog.asksaveasfilename(
            title="导出求解过程为 TXT",
            defaultextension=".txt",
            filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(self.last_detail_text)
            messagebox.showinfo("导出成功", f"求解过程已导出到：\n{path}")
        except Exception as e:
            messagebox.showerror("导出失败", f"导出 TXT 时出现错误：\n{e}")

    def export_csv(self):
        """导出解向量为 CSV 文件。"""
        if self.last_solution_x is None:
            messagebox.showwarning("无解向量", "当前没有可导出的解向量。")
            return

        path = filedialog.asksaveasfilename(
            title="导出解向量为 CSV",
            defaultextension=".csv",
            filetypes=[("CSV 文件", "*.csv"), ("所有文件", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["变量", "数值"])
                for i, val in enumerate(self.last_solution_x, start=1):
                    writer.writerow([f"x{i}", val])
            messagebox.showinfo("导出成功", f"解向量已导出到：\n{path}")
        except Exception as e:
            messagebox.showerror("导出失败", f"导出 CSV 时出现错误：\n{e}")


if __name__ == "__main__":
    app = NAApp()
    app.mainloop()
