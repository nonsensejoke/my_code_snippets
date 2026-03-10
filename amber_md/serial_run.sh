#!/bin/bash
# =============================================================================
# AMBER MD Pipeline Script
# For PPI Dimer with grafted loop
# Usage: bash run_md.sh [options]
# =============================================================================

# ------------------------------------------------------------
# 用户可配置参数（User-configurable parameters）
# ------------------------------------------------------------
DRYRUN=false                            # dry-run 模式（仅打印命令，不实际运行）
AMBER_EXE_EM="${AMBER_EXE_EM:-}"       # EM 阶段可执行程序（优先）
AMBER_EXE_MD="${AMBER_EXE_MD:-}"       # NVT/NPT 阶段可执行程序（优先）
AMBER_EXE="${AMBER_EXE:-pmemd}"        # 通用 fallback（EM/MD 均未指定时使用）
TOPOLOGY=""                            # 拓扑文件 (.parm7 / .prmtop)
INIT_COORDS=""                         # 初始坐标文件 (.rst7 / .inpcrd)
REF_COORDS=""                          # 参考坐标（用于 restraint，通常与初始坐标相同）
INPUT_DIR="."                          # 输入文件所在目录
OUTPUT_DIR="."                         # 输出文件目录

# ------------------------------------------------------------
# 颜色输出
# ------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# ------------------------------------------------------------
# 解析命令行参数
# ------------------------------------------------------------
usage() {
    echo ""
    echo -e "${BOLD}Usage:${NC} bash run_md.sh [options]"
    echo ""
    echo "Options:"
    echo "  -e <executable>   AMBER executable for ALL stages (default: pmemd)"
    echo "                    Used as fallback when -E or -M are not specified"
    echo "  -E <executable>   AMBER executable for EM  stages only"
    echo "                    e.g. sander (EM 通常用 CPU 版本)"
    echo "  -M <executable>   AMBER executable for NVT/NPT stages only"
    echo "                    e.g. pmemd.cuda (MD 通常用 GPU 版本)"
    echo "  -p <topology>     Topology file (.parm7 / .prmtop)  [required]"
    echo "  -c <coordinates>  Initial coordinates (.rst7 / .inpcrd)  [required]"
    echo "  -r <ref_coords>   Reference coordinates for restraints"
    echo "                    (default: same as -c)"
    echo "  -i <input_dir>    Directory containing .txt input files (default: .)"
    echo "  -o <output_dir>   Directory for output files (default: .)"
    echo "  -n                Dry-run mode: print commands without executing
  -h                Show this help message"
    echo ""
    echo "Examples:"
    echo "  # EM 用 CPU sander，NVT/NPT 用 GPU pmemd.cuda"
    echo "  bash run_md.sh -E sander -M pmemd.cuda -p system.parm7 -c system.rst7"
    echo ""
    echo "  # 所有阶段统一使用同一个程序"
    echo "  bash run_md.sh -e pmemd.cuda -p system.parm7 -c system.rst7"
    echo ""
    echo "  # 也可通过环境变量指定"
    echo "  AMBER_EXE_EM=sander AMBER_EXE_MD=pmemd.cuda bash run_md.sh -p system.parm7 -c system.rst7

  # Dry-run：仅打印每步命令，不实际运行（无需真实文件存在）
  bash run_md.sh -n -e pmemd.cuda -p system.parm7 -c system.rst7"
    echo ""
    exit 1
}

while getopts "e:E:M:p:c:r:i:o:nh" opt; do
    case $opt in
        e) AMBER_EXE="$OPTARG" ;;
        E) AMBER_EXE_EM="$OPTARG" ;;
        M) AMBER_EXE_MD="$OPTARG" ;;
        p) TOPOLOGY="$OPTARG" ;;
        c) INIT_COORDS="$OPTARG" ;;
        r) REF_COORDS="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        n) DRYRUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

# -E / -M 未指定时，fallback 到 -e
[[ -z "$AMBER_EXE_EM" ]] && AMBER_EXE_EM="$AMBER_EXE"
[[ -z "$AMBER_EXE_MD" ]] && AMBER_EXE_MD="$AMBER_EXE"

# ------------------------------------------------------------
# 参数检查
# ------------------------------------------------------------
if [[ -z "$TOPOLOGY" || -z "$INIT_COORDS" ]]; then
    echo -e "${RED}[ERROR]${NC} Topology (-p) and initial coordinates (-c) are required."
    usage
fi

if [[ "$DRYRUN" == false ]]; then
    if [[ ! -f "$TOPOLOGY" ]]; then
        echo -e "${RED}[ERROR]${NC} Topology file not found: $TOPOLOGY"
        exit 1
    fi
    if [[ ! -f "$INIT_COORDS" ]]; then
        echo -e "${RED}[ERROR]${NC} Initial coordinates not found: $INIT_COORDS"
        exit 1
    fi
fi

# 默认 ref = 初始坐标
if [[ -z "$REF_COORDS" ]]; then
    REF_COORDS="$INIT_COORDS"
fi

# 检查 AMBER 可执行程序（EM 和 MD 分别检查，dry-run 时跳过）
if [[ "$DRYRUN" == false ]]; then
    if ! command -v "$AMBER_EXE_EM" &> /dev/null; then
        echo -e "${RED}[ERROR]${NC} AMBER EM executable not found: $AMBER_EXE_EM"
        echo "  Please check your AMBER installation or specify with -E or -e"
        exit 1
    fi
    if ! command -v "$AMBER_EXE_MD" &> /dev/null; then
        echo -e "${RED}[ERROR]${NC} AMBER MD executable not found: $AMBER_EXE_MD"
        echo "  Please check your AMBER installation or specify with -M or -e"
        exit 1
    fi
fi

# 检查 ambpdb
if ! command -v ambpdb &> /dev/null; then
    echo -e "${YELLOW}[WARNING]${NC} 'ambpdb' not found in PATH."
    echo "  PDB conversion will be skipped."
    AMBPDB_AVAILABLE=false
else
    AMBPDB_AVAILABLE=true
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------
# 辅助函数
# ------------------------------------------------------------

# 格式化秒数 -> 时分秒
format_time() {
    local seconds=$1
    local hours=$((seconds / 3600))
    local minutes=$(( (seconds % 3600) / 60 ))
    local secs=$((seconds % 60))
    printf "%02dh %02dm %02ds" $hours $minutes $secs
}

# 将 rst7 转换为 pdb
convert_to_pdb() {
    local rst7_file="$1"
    local pdb_file="$2"
    if [[ "$AMBPDB_AVAILABLE" == true ]]; then
        ambpdb -p "$TOPOLOGY" -c "$rst7_file" > "$pdb_file" 2>/dev/null
        if [[ $? -eq 0 ]]; then
            echo -e "  ${GREEN}[PDB]${NC} Saved: $pdb_file"
        else
            echo -e "  ${YELLOW}[WARNING]${NC} ambpdb conversion failed for $rst7_file"
        fi
    fi
}

# 运行单步 MD/EM
# 参数: step_name  input_file  prev_rst7  exe  step_ref
run_step() {
    local step="$1"
    local input_file="${INPUT_DIR}/$2"
    local prev_rst7="$3"
    local exe="$4"
    local step_ref="$5"

    local out_file="${OUTPUT_DIR}/${step}.out"
    local rst7_file="${OUTPUT_DIR}/${step}.rst7"
    local nc_file="${OUTPUT_DIR}/${step}.nc"
    local pdb_file="${OUTPUT_DIR}/${step}.pdb"
    local log_file="${OUTPUT_DIR}/${step}.log"

    echo ""
    echo -e "${CYAN}════════════════════════════════════════════════${NC}"
    echo -e "${BOLD}  Step: ${step}${NC}"
    echo -e "${CYAN}════════════════════════════════════════════════${NC}"

    # ── Dry-run 模式：最先检查，仅打印命令后立即返回 ────────────
    if [[ "$DRYRUN" == true ]]; then
        echo -e "  ${YELLOW}[DRY-RUN]${NC} Command that would be executed:"
        echo ""
        printf "  %s -O \\\\\n" "$exe"
        printf "      -i  %s \\\\\n" "$input_file"
        printf "      -o  %s \\\\\n" "$out_file"
        printf "      -p  %s \\\\\n" "$TOPOLOGY"
        printf "      -c  %s \\\\\n" "$prev_rst7"
        printf "      -r  %s \\\\\n" "$rst7_file"
        printf "      -x  %s \\\\\n" "$nc_file"
        printf "      -inf %s \\\\\n" "${OUTPUT_DIR}/${step}.info"
        printf "      -ref %s\n"     "$step_ref"
        echo ""
        echo -e "  ${BOLD}# PDB conversion (after run)${NC}"
        printf "  ambpdb -p %s -c %s > %s\n" "$TOPOLOGY" "$rst7_file" "$pdb_file"
        echo ""
        echo "$rst7_file"
        return 0
    fi

    # ── 正常模式 ─────────────────────────────────────────────────

    # 检查 PDB 是否已存在（跳过判断）
    if [[ -f "$pdb_file" ]]; then
        echo -e "  ${YELLOW}[SKIP]${NC} $pdb_file already exists. Skipping this step."
        echo -e "         (Remove $pdb_file to re-run this step)"
        echo "$rst7_file"
        return 0
    fi

    # 检查输入文件
    if [[ ! -f "$input_file" ]]; then
        echo -e "  ${RED}[ERROR]${NC} Input file not found: $input_file"
        exit 1
    fi

    # 检查前一步的 rst7
    if [[ ! -f "$prev_rst7" ]]; then
        echo -e "  ${RED}[ERROR]${NC} Previous restart file not found: $prev_rst7"
        exit 1
    fi

    echo -e "  Executable : ${exe}"
    echo -e "  Input      : ${input_file}"
    echo -e "  Topology   : ${TOPOLOGY}"
    echo -e "  Coords in  : ${prev_rst7}"
    echo -e "  Ref coords : ${step_ref}"
    echo -e "  Coords out : ${rst7_file}"
    echo -e "  Trajectory : ${nc_file}"
    echo ""

    # 计时开始
    local start_time=$(date +%s)
    echo -e "  ${GREEN}[RUN]${NC} Started at: $(date '+%Y-%m-%d %H:%M:%S')"

    # 运行 AMBER
    "$exe" -O \
        -i  "$input_file" \
        -o  "$out_file" \
        -p  "$TOPOLOGY" \
        -c  "$prev_rst7" \
        -r  "$rst7_file" \
        -x  "$nc_file" \
        -inf "${OUTPUT_DIR}/${step}.info" \
        -ref "$step_ref" \
        2>&1 | tee "$log_file"

    local exit_code=${PIPESTATUS[0]}

    # 计时结束
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))

    if [[ $exit_code -ne 0 ]]; then
        echo ""
        echo -e "  ${RED}[ERROR]${NC} Step ${step} FAILED (exit code: $exit_code)"
        echo -e "  Check log: $log_file"
        exit $exit_code
    fi

    echo ""
    echo -e "  ${GREEN}[DONE]${NC} Step ${step} completed."
    echo -e "  ${BOLD}  Elapsed time: $(format_time $elapsed)${NC}"

    # 转换为 PDB
    convert_to_pdb "$rst7_file" "$pdb_file"

    # 返回 rst7 路径（供下一步使用）
    echo "$rst7_file"
}

# ------------------------------------------------------------
# 打印配置摘要
# ------------------------------------------------------------
echo ""
echo -e "${BOLD}╔══════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}║         AMBER MD Pipeline - Configuration        ║${NC}"
echo -e "${BOLD}╚══════════════════════════════════════════════════╝${NC}"
echo -e "  AMBER EM  Exec   : ${CYAN}${AMBER_EXE_EM}${NC}
  AMBER MD  Exec   : ${CYAN}${AMBER_EXE_MD}${NC}"
echo -e "  Topology         : ${TOPOLOGY}"
echo -e "  Initial Coords   : ${INIT_COORDS}"
echo -e "  Initial Ref      : ${REF_COORDS}"
echo -e "  Input  Directory : ${INPUT_DIR}"
echo -e "  Output Directory : ${OUTPUT_DIR}"
echo -e "  ambpdb available : ${AMBPDB_AVAILABLE}"
echo -e "  Dry-run mode     : ${BOLD}${DRYRUN}${NC}"
echo ""

# 总计时开始
PIPELINE_START=$(date +%s)

# ------------------------------------------------------------
# 定义步骤序列
# 格式: "步骤名称|输入文件名|阶段类型(em/md)|完成后更新全局REF(yes/-)"
#
#   第4段 = yes  → 该步完成后，将全局 CURRENT_REF 更新为该步输出的 rst7
#   第4段 = -    → 不更新，后续步骤继续沿用当前的 CURRENT_REF
#
# 全局 CURRENT_REF 初始值 = 用户通过 -r 指定的参考结构（默认为初始坐标）
#
# 当前配置：
#   nvt-7 标记为 yes → nvt-7 完成后 CURRENT_REF 更新为 nvt-7.rst7
#   → npt-1/2 自动以弛豫后的构象为 restraint 锚点
#
STEPS=(
    "em-1|em-1_in.txt|em|-"
    "em-2|em-2_in.txt|md|-"
    "em-3|em-3_in.txt|md|-"
    "em-4|em-4_in.txt|md|-"
    "em-5|em-5_in.txt|md|-"
    "em-6|em-6_in.txt|md|yes"
    "nvt-1|nvt-1_in.txt|md|-"
    "nvt-2|nvt-2_in.txt|md|-"
    "nvt-3|nvt-3_in.txt|md|-"
    "nvt-4|nvt-4_in.txt|md|-"
    "nvt-5|nvt-5_in.txt|md|-"
    "nvt-6|nvt-6_in.txt|md|-"
    "npt-1|npt-1_in.txt|md|-"
    "npt-2|npt-2_in.txt|md|-"
    "npt-3|npt-3_in.txt|md|-"
    "npt-4|npt-4_in.txt|md|-"
)

# ------------------------------------------------------------
# 依次运行每一步
# ------------------------------------------------------------
CURRENT_RST7="$INIT_COORDS"
CURRENT_REF="$REF_COORDS"     # 全局 REF，初始值 = 用户指定的参考结构

for entry in "${STEPS[@]}"; do
    STEP_NAME=$(echo "$entry" | cut -d'|' -f1)
    INPUT_FNAME=$(echo "$entry" | cut -d'|' -f2)
    STAGE_TYPE=$(echo "$entry" | cut -d'|' -f3)
    UPDATE_REF=$(echo "$entry" | cut -d'|' -f4)

    # 根据阶段类型选择对应的可执行程序
    if [[ "$STAGE_TYPE" == "em" ]]; then
        STEP_EXE="$AMBER_EXE_EM"
    else
        STEP_EXE="$AMBER_EXE_MD"
    fi

    # 本步统一使用全局 CURRENT_REF
    STEP_REF="$CURRENT_REF"

    # run_step 会打印很多内容，最后一行是 rst7 路径
    # 用临时文件捕获 rst7 路径
    TMPFILE=$(mktemp)

    run_step "$STEP_NAME" "$INPUT_FNAME" "$CURRENT_RST7" "$STEP_EXE" "$STEP_REF" | tee /dev/tty | tail -1 > "$TMPFILE"
    NEW_RST7=$(cat "$TMPFILE")
    rm -f "$TMPFILE"

    # 更新 CURRENT_RST7
    # dry-run 下文件不存在，直接使用预期输出路径；正常模式下验证文件存在
    if [[ "$DRYRUN" == true ]]; then
        CURRENT_RST7="${OUTPUT_DIR}/${STEP_NAME}.rst7"
    elif [[ -n "$NEW_RST7" && -f "$NEW_RST7" ]]; then
        CURRENT_RST7="$NEW_RST7"
    elif [[ -f "${OUTPUT_DIR}/${STEP_NAME}.rst7" ]]; then
        CURRENT_RST7="${OUTPUT_DIR}/${STEP_NAME}.rst7"
    fi

    # 若该步标记了 update_ref=yes，则将全局 CURRENT_REF 更新为本步输出的 rst7
    if [[ "$UPDATE_REF" == "yes" ]]; then
        if [[ "$DRYRUN" == true ]]; then
            # dry-run 下 rst7 不存在，使用预期路径更新
            CURRENT_REF="${OUTPUT_DIR}/${STEP_NAME}.rst7"
            echo -e "  ${CYAN}[REF UPDATE]${NC} (dry-run) Global reference would be updated to: ${CURRENT_REF}"
        elif [[ -f "$CURRENT_RST7" ]]; then
            CURRENT_REF="$CURRENT_RST7"
            echo -e "  ${CYAN}[REF UPDATE]${NC} Global reference updated to: ${CURRENT_REF}"
        else
            echo -e "  ${YELLOW}[WARNING]${NC} update_ref=yes but rst7 not found: $CURRENT_RST7"
            echo -e "             Global reference NOT updated."
        fi
    fi
done

# ------------------------------------------------------------
# 总计时
# ------------------------------------------------------------
PIPELINE_END=$(date +%s)
PIPELINE_ELAPSED=$((PIPELINE_END - PIPELINE_START))

echo ""
echo -e "${BOLD}╔══════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}║              Pipeline Completed!                 ║${NC}"
echo -e "${BOLD}╚══════════════════════════════════════════════════╝${NC}"
echo -e "  Total elapsed time: ${BOLD}$(format_time $PIPELINE_ELAPSED)${NC}"
echo -e "  Final restart file: ${CURRENT_RST7}"
echo ""
echo -e "  Generated PDB checkpoints:"
for entry in "${STEPS[@]}"; do
    STEP_NAME="${entry%%|*}"
    PDB="${OUTPUT_DIR}/${STEP_NAME}.pdb"
    if [[ -f "$PDB" ]]; then
        echo -e "    ${GREEN}✓${NC}  $PDB"
    else
        echo -e "    ${YELLOW}✗${NC}  $PDB  (not found)"
    fi
done
echo ""
