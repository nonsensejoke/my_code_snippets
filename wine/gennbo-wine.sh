#!/usr/bin/env bash
set -Eeuo pipefail

program_name=${0##*/}

usage() {
  cat <<EOF
Usage:
  $program_name [options] INPUT.47 [-- GENNBO_ARG ...]

Run the Windows GenNBO/NBO7 executables under Wine and save clean text
output. A pseudo-terminal is used because the Windows Fortran runtime fails
when its console output is redirected directly.

Options:
  -o, --output FILE       Output file (default: INPUT.nbo)
  -f, --force             Overwrite an existing output file
  -i, --int i4|i8         Integer build (default: i8)
  -m, --memory SIZE       GenNBO memory, e.g. 50mw, 512mb, 1gb
                          (default: NBOMEM or 50mw)
  -d, --detail VALUE      Set NBODTL; use "verbose" for detailed output
      --verbose           Equivalent to --detail verbose
      --nbobin DIR        Directory containing gennbo.*.exe and nbo7.*.exe
                          (default: NBOBIN or the bin directory beside script)
      --wine COMMAND      Wine command/path (default: WINE_CMD or wine)
      --wine-prefix DIR   Set WINEPREFIX for this run
      --timeout SECONDS   Stop after this many seconds; 0 disables timeout
                          (default: NBO_TIMEOUT or 0)
      --stderr FILE       Save Wine diagnostics to FILE, even when empty
  -h, --help              Show this help

Environment variables:
  NBOBIN, NBOMEM, NBODTL, NBO_INT, NBO_TIMEOUT, WINE_CMD, WINEPREFIX,
  WINEDEBUG

Examples:
  $program_name test.47
  $program_name --verbose --memory 1gb -o test-verbose.nbo test.47
  $program_name --int i4 --timeout 600 calculation.47
EOF
}

die() {
  printf '%s: error: %s\n' "$program_name" "$*" >&2
  exit 2
}

note() {
  printf '%s: %s\n' "$program_name" "$*" >&2
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
input=
output=
integer_kind=${NBO_INT:-i8}
memory=${NBOMEM:-50mw}
detail=${NBODTL:-}
nbobin=${NBOBIN:-"$script_dir/bin"}
wine_command=${WINE_CMD:-wine}
wine_prefix=${WINEPREFIX:-}
timeout_seconds=${NBO_TIMEOUT:-0}
stderr_file=
force=0
extra_args=()

while (($#)); do
  case $1 in
    -o|--output)
      (($# >= 2)) || die "$1 requires a file"
      output=$2
      shift 2
      ;;
    -f|--force)
      force=1
      shift
      ;;
    -i|--int)
      (($# >= 2)) || die "$1 requires i4 or i8"
      integer_kind=$2
      shift 2
      ;;
    -m|--memory)
      (($# >= 2)) || die "$1 requires a size"
      memory=$2
      shift 2
      ;;
    -d|--detail)
      (($# >= 2)) || die "$1 requires a value"
      detail=$2
      shift 2
      ;;
    --verbose)
      detail=verbose
      shift
      ;;
    --nbobin)
      (($# >= 2)) || die "$1 requires a directory"
      nbobin=$2
      shift 2
      ;;
    --wine)
      (($# >= 2)) || die "$1 requires a command or path"
      wine_command=$2
      shift 2
      ;;
    --wine-prefix)
      (($# >= 2)) || die "$1 requires a directory"
      wine_prefix=$2
      shift 2
      ;;
    --timeout)
      (($# >= 2)) || die "$1 requires a number of seconds"
      timeout_seconds=$2
      shift 2
      ;;
    --stderr)
      (($# >= 2)) || die "$1 requires a file"
      stderr_file=$2
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      extra_args=("$@")
      break
      ;;
    -*)
      die "unknown option: $1"
      ;;
    *)
      [[ -z $input ]] || die "only one input file is allowed (use -- before extra GenNBO arguments)"
      input=$1
      shift
      ;;
  esac
done

[[ -n $input ]] || { usage >&2; exit 2; }
[[ $integer_kind == i4 || $integer_kind == i8 ]] || die "--int must be i4 or i8"
[[ $timeout_seconds =~ ^[0-9]+$ ]] || die "--timeout must be a non-negative integer"
memory=${memory,,}
[[ $memory =~ ^[1-9][0-9]*(mw|mb|gb)$ ]] || die "--memory must look like 50mw, 512mb, or 1gb"

command -v "$wine_command" >/dev/null 2>&1 || die "Wine command not found: $wine_command"
command -v script >/dev/null 2>&1 || die "required command not found: script"
command -v sed >/dev/null 2>&1 || die "required command not found: sed"
command -v tr >/dev/null 2>&1 || die "required command not found: tr"

input=$(realpath -e -- "$input") || die "input file not found: $input"
[[ -f $input && -r $input ]] || die "input is not a readable regular file: $input"
[[ ${input,,} == *.47 ]] || note "warning: input filename does not end in .47"

nbobin=$(realpath -e -- "$nbobin") || die "NBO binary directory not found: $nbobin"
genexe=$nbobin/gennbo.$integer_kind.exe
nboexe=$nbobin/nbo7.$integer_kind.exe
[[ -f $genexe && -r $genexe ]] || die "GenNBO executable not found: $genexe"
[[ -f $nboexe && -r $nboexe ]] || die "NBO7 executable not found: $nboexe"

if [[ -z $output ]]; then
  output=${input%.*}.nbo
fi
output_dir=$(dirname -- "$output")
output_base=$(basename -- "$output")
mkdir -p -- "$output_dir"
output_dir=$(cd -- "$output_dir" && pwd -P)
output=$output_dir/$output_base

if [[ -e $output && $force -ne 1 ]]; then
  die "output already exists: $output (use --force to overwrite)"
fi

if [[ -n $stderr_file ]]; then
  stderr_dir=$(dirname -- "$stderr_file")
  stderr_base=$(basename -- "$stderr_file")
  mkdir -p -- "$stderr_dir"
  stderr_dir=$(cd -- "$stderr_dir" && pwd -P)
  stderr_file=$stderr_dir/$stderr_base
  if [[ -e $stderr_file && $force -ne 1 ]]; then
    die "stderr file already exists: $stderr_file (use --force to overwrite)"
  fi
fi

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/gennbo-wine.XXXXXX")
tmp_output=$(mktemp "$output_dir/.${output_base}.XXXXXX")
tmp_stderr=$tmpdir/wine.stderr
cleanup() {
  rm -f -- "$tmp_output"
  rm -rf -- "$tmpdir"
}
trap cleanup EXIT INT TERM

# Wine exposes the Unix root directory as drive Z: by default.
nboexe_windows="Z:${nboexe//\//\\}"

wine_env=(
  "WINEDEBUG=${WINEDEBUG:--all}"
  "NBOEXE=$nboexe_windows"
  "NBOMEM=$memory"
)
[[ -z $detail ]] || wine_env+=("NBODTL=$detail")
[[ -z $wine_prefix ]] || wine_env+=("WINEPREFIX=$wine_prefix")

input_dir=$(dirname -- "$input")
input_name=$(basename -- "$input")
child=(env "${wine_env[@]}" "$wine_command" "$genexe" "$input_name" "${extra_args[@]}")
printf -v child_command '%q ' "${child[@]}"
printf -v quoted_stderr '%q' "$tmp_stderr"
child_command+="2>$quoted_stderr"

note "input:  $input"
note "output: $output"
note "build:  $integer_kind; memory: $memory${detail:+; detail: $detail}"

cd -- "$input_dir"
set +e
if ((timeout_seconds > 0)); then
  timeout "${timeout_seconds}s" script -qefc "$child_command" /dev/null |
    sed $'s/\033\\[[0-9;?]*[A-Za-z]//g' |
    tr -d '\r' >"$tmp_output"
else
  script -qefc "$child_command" /dev/null |
    sed $'s/\033\\[[0-9;?]*[A-Za-z]//g' |
    tr -d '\r' >"$tmp_output"
fi
pipeline_status=("${PIPESTATUS[@]}")
set -e

run_status=0
for status in "${pipeline_status[@]}"; do
  if ((status != 0)); then
    run_status=$status
    break
  fi
done

if ((run_status != 0)); then
  failure_stderr=${stderr_file:-"${output}.wine.stderr"}
  if [[ -s $tmp_stderr ]]; then
    cp -- "$tmp_stderr" "$failure_stderr"
    note "Wine diagnostics: $failure_stderr"
  fi
  if ((run_status == 124)); then
    die "timed out after ${timeout_seconds} seconds"
  fi
  die "GenNBO failed with exit status $run_status"
fi

if ! grep -aFq 'NBO analysis completed' "$tmp_output"; then
  failure_stderr=${stderr_file:-"${output}.wine.stderr"}
  if [[ -s $tmp_stderr ]]; then
    cp -- "$tmp_stderr" "$failure_stderr"
    note "Wine diagnostics: $failure_stderr"
  fi
  die "process exited successfully but completion marker was not found"
fi

if [[ -n $stderr_file ]]; then
  cp -- "$tmp_stderr" "$stderr_file"
elif [[ -s $tmp_stderr ]]; then
  cp -- "$tmp_stderr" "${output}.wine.stderr"
  note "Wine diagnostics: ${output}.wine.stderr"
fi

if ((force)); then
  mv -f -- "$tmp_output" "$output"
else
  mv -- "$tmp_output" "$output"
fi
note "completed successfully"
