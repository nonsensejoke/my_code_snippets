#!/usr/bin/env bash
#
# cluster_load.sh — 并行采集一组 Linux 节点的 CPU 负载与内存占用，
#                   按"空闲程度"（1 分钟负载率）从低到高排序输出，
#                   方便在提交计算任务前挑选最闲的机器。
#
# 特点:
#   * 所有节点并行 SSH，整体耗时约等于最慢的一台
#   * 远端只读 /proc/loadavg 与 /proc/meminfo，不依赖 free/uptime/top 的输出格式，
#     跨发行版稳定（这些命令在不同版本里列数会变）
#   * 连接失败自动重试一次（集群网络偶发抖动很常见，单次失败不代表节点宕机）
#   * 负载率按核数归一化，不同核数的机器可以横向对比
#   * 综合 CPU 与内存给出推荐节点，内存已用 >=85% 的直接剔除
#
# 依赖: bash、ssh、awk —— 无需在远端安装任何东西
#
# 用法示例:
#   ./cluster_load.sh -n "node01 node02 node03"
#   ./cluster_load.sh -n "$(echo node{01..10})"        # 借用 shell 的花括号展开
#   ./cluster_load.sh -f ~/.cluster_nodes              # 从文件读取节点清单
#   CLUSTER_NODES="gpu01 gpu02" ./cluster_load.sh      # 用环境变量
#   ./cluster_load.sh -n "node01" -d example.com -u alice -k ~/.ssh/cluster_key
#
# 建议把常用参数写进 ~/.bashrc 省得每次敲:
#   export CLUSTER_NODES="$(echo node{01..10})"
#   export CLUSTER_DOMAIN="example.com"
#   export CLUSTER_USER="alice"
#   export CLUSTER_KEY="$HOME/.ssh/cluster_key"
#
set -uo pipefail

KEY="${CLUSTER_KEY:-$HOME/.ssh/id_rsa}"
SSH_USER="${CLUSTER_USER:-$(id -un)}"
DOMAIN="${CLUSTER_DOMAIN:-}"          # 留空则直接使用主机名，不追加域名后缀
CONNECT_TIMEOUT="${CLUSTER_TIMEOUT:-15}"
NODES="${CLUSTER_NODES:-}"
NODEFILE=""

usage() {
  cat <<USAGE
用法: $(basename "$0") [选项]

选项:
  -n "<节点列表>"  主机名，空格分隔，例如 "node01 node02"
  -f <文件>        从文件读取节点清单（每行一个，# 开头为注释）
  -k <私钥路径>    SSH 私钥       (默认: \$HOME/.ssh/id_rsa)
  -u <用户名>      SSH 用户       (默认: 当前本地用户名)
  -d <域名后缀>    主机域名后缀   (默认: 空，即直接用主机名连接)
  -t <秒>          SSH 连接超时   (默认: 15)
  -h               显示本帮助

环境变量 CLUSTER_NODES / CLUSTER_KEY / CLUSTER_USER / CLUSTER_DOMAIN /
CLUSTER_TIMEOUT 同样生效，命令行参数优先级更高。

输出列说明:
  负载率 = load1 / 核数 x 100%，已按核数归一化，可跨机型直接比较
  趋势   = 对比 load1 与 load15，判断负载在上升还是消退
  状态   = 空闲 (<25% 且内存 <70%) / 中等 / 繁忙；内存 >=85% 额外标注
USAGE
}

while getopts "n:f:k:u:d:t:h" opt; do
  case "$opt" in
    n) NODES="$OPTARG" ;;
    f) NODEFILE="$OPTARG" ;;
    k) KEY="$OPTARG" ;;
    u) SSH_USER="$OPTARG" ;;
    d) DOMAIN="$OPTARG" ;;
    t) CONNECT_TIMEOUT="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 2 ;;
  esac
done

if [ -n "$NODEFILE" ]; then
  if [ ! -r "$NODEFILE" ]; then
    echo "错误: 无法读取节点清单文件: $NODEFILE" >&2
    exit 1
  fi
  NODES="$NODES $(sed -e 's/#.*//' "$NODEFILE" | tr '\n' ' ')"
fi

if [ -z "${NODES// /}" ]; then
  echo "错误: 未指定任何节点。请用 -n / -f 或环境变量 CLUSTER_NODES 指定。" >&2
  echo >&2
  usage >&2
  exit 2
fi

if [ ! -r "$KEY" ]; then
  echo "错误: 无法读取 SSH 私钥: $KEY  (可用 -k 指定)" >&2
  exit 1
fi

WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/cluster-load.XXXXXX")
trap 'rm -rf "$WORKDIR"' EXIT INT TERM

# 远端命令: 只读 /proc，输出一行
#   <核数> <load1> <load5> <load15> <MemTotal_kB> <MemAvailable_kB>
# 用 MemAvailable 而不是 MemFree —— 后者不含可回收的 cache，会严重高估内存压力。
REMOTE_CMD='
c=$(nproc 2>/dev/null || echo 0)
read -r l1 l5 l15 _rest < /proc/loadavg
mt=$(grep -m1 MemTotal /proc/meminfo | tr -dc 0-9)
ma=$(grep -m1 MemAvailable /proc/meminfo | tr -dc 0-9)
echo "$c $l1 $l5 $l15 ${mt:-0} ${ma:-0}"
'

target_of() {
  if [ -n "$DOMAIN" ]; then printf '%s@%s.%s' "$SSH_USER" "$1" "$DOMAIN"
  else                      printf '%s@%s'    "$SSH_USER" "$1"
  fi
}

probe_once() {
  local h="$1" out
  out=$(ssh -i "$KEY" \
            -o StrictHostKeyChecking=no \
            -o BatchMode=yes \
            -o LogLevel=ERROR \
            -o ConnectTimeout="$CONNECT_TIMEOUT" \
            -o ServerAliveInterval=5 \
            -o ServerAliveCountMax=3 \
            "$(target_of "$h")" "$REMOTE_CMD" 2>"$WORKDIR/$h.err" | tail -n1)
  printf '%s' "$out"
}

probe() {
  local h="$1" out attempt
  # 集群偶发连接抖动，失败后重试一次再判定为不可达
  for attempt in 1 2; do
    out=$(probe_once "$h")
    [ -n "$out" ] && break
    [ "$attempt" = 1 ] && sleep 2
  done
  if [ -n "$out" ]; then
    printf '%s %s\n' "$h" "$out" > "$WORKDIR/$h.ok"
  else
    local msg
    msg=$(tr -d '\r' < "$WORKDIR/$h.err" | grep -v '^$' | tail -n1)
    printf '%s\t%s\n' "$h" "${msg:-SSH 两次均无响应（连接超时或会话被中断）}" > "$WORKDIR/$h.bad"
  fi
}

count=$(printf '%s\n' $NODES | wc -l | tr -d ' ')
printf '正在并行采集 %s 个节点 (超时 %ss) ...\n' "$count" "$CONNECT_TIMEOUT" >&2

for h in $NODES; do probe "$h" & done
wait

# ---------------- 结果汇总 ----------------
printf '\n采集时间: %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S')"

if ! ls "$WORKDIR"/*.ok >/dev/null 2>&1; then
  echo "所有节点均无法访问。" >&2
else
  printf '%-12s %5s %7s %7s %7s %8s %6s %14s  %s\n' \
         节点 核数 load1 load5 load15 负载率 趋势 内存 状态
  printf '%s\n' "---------------------------------------------------------------------------------------"

  cat "$WORKDIR"/*.ok | awk '
    {
      host=$1; cpu=$2+0; l1=$3+0; l5=$4+0; l15=$5+0; mt=$6+0; ma=$7+0
      if (cpu <= 0) cpu = 1
      ratio    = l1 / cpu * 100
      memtot   = mt / 1048576                      # kB -> GiB
      memused  = (mt - ma) / 1048576
      memratio = (mt > 0) ? (mt - ma) / mt * 100 : 0

      # 趋势: 对比 1 分钟与 15 分钟负载
      if      (l15 > 0 && l1 < l15 * 0.8) trend = "↓ 降"
      else if (l15 > 0 && l1 > l15 * 1.2) trend = "↑ 升"
      else                                trend = "→ 平"

      if      (ratio < 25 && memratio < 70) st = "✅ 空闲"
      else if (ratio < 60 && memratio < 85) st = "🟡 中等"
      else                                  st = "🔴 繁忙"
      if (memratio >= 85) st = st "  ⚠ 内存紧张"

      # 首列为排序键（负载率），打印时用 cut 去掉
      printf "%012.4f\t%-12s %5d %7.2f %7.2f %7.2f %7.1f%% %6s %5.0f/%-4.0fG  %s\n", \
             ratio, host, cpu, l1, l5, l15, ratio, trend, memused, memtot, st
    }' | sort -n | cut -f2-
fi

# ---------------- 无法访问的节点 ----------------
if ls "$WORKDIR"/*.bad >/dev/null 2>&1; then
  printf '\n无法访问的节点:\n'
  cat "$WORKDIR"/*.bad | while IFS=$'\t' read -r h msg; do
    printf '  %-12s %s\n' "$h" "${msg:-未知错误}"
  done
fi

# ---------------- 推荐 ----------------
if ls "$WORKDIR"/*.ok >/dev/null 2>&1; then
  # 综合评分: CPU 负载率 + 内存占用率的一半; 内存 >=85% 的节点直接剔除
  best=$(cat "$WORKDIR"/*.ok | awk '
    { cpu=($2+0>0)?$2+0:1
      cr=$3/cpu*100
      mr=($6+0>0)?($6-$7)/$6*100:0
      if (mr >= 85) next
      printf "%012.4f %s %.0f %.0f\n", cr + mr/2, $1, cr, mr }' \
    | sort -n | head -3 \
    | awk '{ printf "%s(CPU %s%%, 内存 %s%%)  ", $2, $3, $4 }')
  if [ -n "$best" ]; then
    printf '\n推荐使用: %s\n' "$best"
  else
    printf '\n推荐使用: 无 — 所有可达节点内存均已超过 85%%\n'
  fi
fi
