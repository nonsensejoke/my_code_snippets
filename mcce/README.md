# MCCE 后处理脚本

这两个脚本用于 MCCE 计算完成后的后处理，都需要读取 MCCE 的 `step2_out.pdb`（构象坐标）。大多数情况下还会用到 `fort.38`（各构象在不同 pH 下的占比表）。

| 脚本 | 用途 |
| --- | --- |
| `patch_mcce_pdb.py` | 在指定 pH 下选出每个残基占比最高的状态，写回原始 PDB，得到质子化好的结构 |
| `export_residue_microstates.py` | 把单个残基的所有微观态排列在一个 PyMOL 会话里，方便对比 |

运行环境：Python 3.7+，只用到标准库。`export_residue_microstates.py` 还需要命令行可调用的 PyMOL。

---

## patch_mcce_pdb.py

在给定 pH 下，从 population 表中为每个残基选出占比最高的构象，用 `step2_out.pdb` 中该构象的坐标（骨架 `_000` 加上选中的侧链构象）替换原始 PDB 中的对应残基。

处理的残基包括：

- 可滴定残基：ASP、GLU、LYS、HIS、CYS
- 二硫键半胱氨酸：CYD
- 极性转子：ASN、GLN（默认处理，加 `--no-polar-rotamers` 可关闭）

如果 population 表里没有恰好等于目标 pH 的列，脚本会在相邻两个 pH 列之间线性插值。

### 用法

```bash
python patch_mcce_pdb.py \
  --input-pdb prot.pdb \
  --population fort.38 \
  --mcce-pdb step2_out.pdb \
  --ph 7.4
```

| 参数 | 说明 | 默认值 |
| --- | --- | --- |
| `--input-pdb` | 原始输入 PDB（必填） | |
| `--population` | MCCE population 表，如 `fort.38`，首行须以 `ph` 开头（必填） | |
| `--mcce-pdb` | MCCE `step2_out.pdb`（必填） | |
| `--ph` | 目标 pH，须在表格的 pH 范围内 | `7.4` |
| `--standard-out` | standard 版输出 | `patched_standard.pdb` |
| `--amber-out` | Amber 版输出 | `patched_amber.pdb` |
| `--report-out` | 占比报告 | `population_report.txt` |
| `--no-polar-rotamers` | 不替换 ASN/GLN 的转子构象 | 关 |

### 输出

**`patched_standard.pdb`**：被替换的残基保留 MCCE 构象中的全部原子（含氢），残基名保持标准名，CYD 写为 CYS。

**`patched_amber.pdb`**：残基名改为 Amber 命名，侧链只保留重原子，骨架只额外保留 H，适合直接交给 tleap 补氢。命名规则如下：

| 原残基 | 判断依据 | Amber 名 |
| --- | --- | --- |
| ASP | 有 HD1 或 HD2 | ASH，否则 ASP |
| GLU | 有 HE1 或 HE2 | GLH，否则 GLU |
| LYS | 有 HZ1、HZ2 但无 HZ3 | LYN，否则 LYS |
| HIS | HD1 + HE2 / 仅 HD1 / 仅 HE2 | HIP / HID / HIE |
| CYS | 无 HG | CYM，否则 CYS |
| CYD | — | CYX |

**`population_report.txt`**：列出每个残基选中的状态 ID 和占比，结果同时打印到屏幕。占比落在 40–60% 的残基会标记 `WARNING`，说明其状态不确定，建议人工检查。

两个 PDB 中的原子都会重新编号，`CONECT` 记录也会按新编号重新映射。

### 注意

- 原始 PDB 中的残基按链号和残基号匹配，**带插入码（icode）的残基不会被替换**。
- 如果 `step2_out.pdb` 里找不到某个残基的骨架构象 `_000` 或被选中的侧链构象，脚本会直接报错退出。

---

## export_residue_microstates.py

把一个残基在 `step2_out.pdb` 中的所有侧链构象分别与骨架 `_000` 拼接，每个构象生成一个 PyMOL 对象，沿指定轴等间距排开，并在每个对象上方标注状态名，最后保存为 `.pse` 会话。

如果提供了 fort.38，或者当前目录下存在 `fort.38`，脚本会用其中的标签（如 `HIS01B0071_002`）命名对象，并按 fort.38 中的顺序排列。没有 fort.38 时，标签取自 `step2_out.pdb` 中的构象代码。

### 用法

```bash
python export_residue_microstates.py --step2-pdb step2_out.pdb --residue B0071
```

| 参数 | 说明 | 默认值 |
| --- | --- | --- |
| `--step2-pdb` | MCCE `step2_out.pdb`（必填） | |
| `--residue` | 目标残基，链号 + 残基号，如 `B0071` 或 `C91`（必填） | |
| `--fort38` | fort.38 路径 | 当前目录的 `fort.38`（如存在） |
| `--out-prefix` | 输出文件前缀 | `{RES}_{链号}{残基号}_microstates` |
| `--spacing` | 相邻微观态之间的间距（Å） | `6.0` |
| `--axis` | 排列方向：`x` / `y` / `z` | `x` |
| `--pymol` | PyMOL 可执行文件 | `pymol` |
| `--keep-temp` | 在当前目录保留中间生成的 `.pdb` 和 `.pml` | 关 |

输出为 `{out-prefix}.pse`，用 PyMOL 直接打开即可。

---

## 推荐流程

1. 先运行 `patch_mcce_pdb.py`，得到质子化后的结构和占比报告。
2. 在 `population_report.txt` 中找出标记了 `WARNING` 的残基。
3. 对这些残基运行 `export_residue_microstates.py`，在 PyMOL 中对比各个构象，再决定是否需要手动调整。
