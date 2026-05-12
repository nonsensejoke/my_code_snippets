#!/usr/bin/env bash

set -euo pipefail

# ── 默认参数 ────────────────────────────────────────────────────────────────
SIZE_LIMIT="200M"
ZIP_DIR="."
ZIP_NAME="big-files.zip"
TARGET_DIR="."
COMPRESS_LEVEL=0   # 默认不压缩（-0），打包速度快，适合视频/图片等已压缩文件

# ── 依赖检查 ────────────────────────────────────────────────────────────────
check_deps() {
  local missing=()
  for cmd in zip unzip find; do
    command -v "$cmd" >/dev/null 2>&1 || missing+=("$cmd")
  done
  if [[ ${#missing[@]} -gt 0 ]]; then
    echo "错误：缺少以下命令：${missing[*]}" >&2
    exit 1
  fi
}

# ── 帮助 ────────────────────────────────────────────────────────────────────
show_help() {
  cat <<EOF
用法:
  $0 [选项] [目标目录]

示例:
  $0 ~/Downloads
  $0 --size 500M ~/Movies
  $0 --compress --output-dir /Volumes/MyDrive --name movies-large.zip ~/Movies

选项:
  -s, --size        文件大小阈值，默认: 200M
                    单位: c(字节) k(KiB) M(MiB) G(GiB)
                    例如: 200M, 1G, 500k

  -o, --output-dir  zip 放置目录，默认: ~/Desktop

  -n, --name        zip 文件名，默认: big-files.zip
                    （自动补全 .zip 后缀）

  -c, --compress    启用压缩（默认为仅打包不压缩，适合视频/图片）
                    启用后使用压缩级别 1（速度与体积的平衡）

  -h, --help        显示帮助
EOF
}

# ── 参数解析 ────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--size)
      SIZE_LIMIT="$2"
      shift 2
      ;;
    -o|--output-dir)
      ZIP_DIR="$2"
      shift 2
      ;;
    -n|--name)
      ZIP_NAME="$2"
      shift 2
      ;;
    -c|--compress)
      COMPRESS_LEVEL=1
      shift
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    -*)
      echo "未知选项: $1" >&2
      show_help
      exit 1
      ;;
    *)
      TARGET_DIR="$1"
      shift
      ;;
  esac
done

# ── 参数校验 ────────────────────────────────────────────────────────────────

# SIZE_LIMIT 格式：数字 + find 支持的单位
if [[ ! "$SIZE_LIMIT" =~ ^[0-9]+[ckMG]$ ]]; then
  echo "错误：无效的大小格式 '$SIZE_LIMIT'" >&2
  echo "      应为数字 + 单位，例如: 200M、1G、500k、102400c" >&2
  exit 1
fi

# ZIP_NAME 自动补全 .zip 后缀
[[ "$ZIP_NAME" != *.zip ]] && ZIP_NAME="${ZIP_NAME}.zip"

# ── 路径处理 ────────────────────────────────────────────────────────────────
if [[ ! -d "$TARGET_DIR" ]]; then
  echo "错误：目标目录不存在：$TARGET_DIR" >&2
  exit 1
fi

TARGET_DIR="$(cd "$TARGET_DIR" && pwd)"
ZIP_DIR="$(mkdir -p "$ZIP_DIR" && cd "$ZIP_DIR" && pwd)"
ZIP_PATH="$ZIP_DIR/$ZIP_NAME"

BASE_NAME="${ZIP_NAME%.zip}"
ORIGINAL_PATH_FILE="$ZIP_DIR/${BASE_NAME}_ORIGINAL_DIRECTORY.txt"
FILE_LIST="$(mktemp)"

# ── 清理 ────────────────────────────────────────────────────────────────────
cleanup() {
  rm -f "$FILE_LIST"
}
trap cleanup EXIT

# ── 信息输出 ────────────────────────────────────────────────────────────────
echo "目标目录:   $TARGET_DIR"
echo "大小阈值:   $SIZE_LIMIT"
echo "zip 路径:   $ZIP_PATH"
echo "路径记录:   $ORIGINAL_PATH_FILE"
echo "压缩级别:   -${COMPRESS_LEVEL}$([ "$COMPRESS_LEVEL" -eq 0 ] && echo "（仅打包，不压缩）" || echo "（已启用压缩）")"
echo

# ── 防覆盖检查 ───────────────────────────────────────────────────────────────
if [[ -e "$ZIP_PATH" ]]; then
  echo "错误：zip 文件已存在，避免覆盖: $ZIP_PATH" >&2
  exit 1
fi
if [[ -e "$ORIGINAL_PATH_FILE" ]]; then
  echo "错误：路径记录文件已存在，避免覆盖: $ORIGINAL_PATH_FILE" >&2
  exit 1
fi

# ── 查找大文件 ───────────────────────────────────────────────────────────────
cd "$TARGET_DIR"

echo "正在查找超过 $SIZE_LIMIT 的文件..."
find . -type f -size +"$SIZE_LIMIT" -print > "$FILE_LIST"

if [[ ! -s "$FILE_LIST" ]]; then
  echo "没有找到超过 $SIZE_LIMIT 的文件。"
  exit 0
fi

echo
echo "将被打包的文件:"
cat "$FILE_LIST"
echo

FILE_COUNT=$(grep -c '' "$FILE_LIST")
echo "共找到 $FILE_COUNT 个文件。"
echo

read -r -p "继续创建 zip 吗？[y/N] " CONFIRM
if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
  echo "已取消。"
  exit 0
fi

# ── 打包 ────────────────────────────────────────────────────────────────────
printf '%s\n' "$TARGET_DIR" > "$ORIGINAL_PATH_FILE"

echo
echo "正在创建 zip（压缩级别 -${COMPRESS_LEVEL}）..."
zip "-${COMPRESS_LEVEL}" "$ZIP_PATH" -@ < "$FILE_LIST"

echo
echo "正在验证 zip 完整性..."
unzip -t "$ZIP_PATH"

echo
echo "zip 创建并验证完成。"
echo "zip:      $ZIP_PATH"
echo "路径记录: $ORIGINAL_PATH_FILE"
echo

# ── 删除原文件 ───────────────────────────────────────────────────────────────
read -r -p "确认删除这些原始大文件吗？此操作不可撤销。[y/N] " DELETE_CONFIRM
if [[ ! "$DELETE_CONFIRM" =~ ^[Yy]$ ]]; then
  echo "已保留原文件，没有删除。"
  exit 0
fi

echo
echo "正在删除原始大文件..."
while IFS= read -r file; do
  rm -f -- "$file"
  echo "  已删除: $file"
done < "$FILE_LIST"

echo
echo "完成。请将以下两个文件一起复制到外置硬盘："
echo "  $ZIP_PATH"
echo "  $ORIGINAL_PATH_FILE"
