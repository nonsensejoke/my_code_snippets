#!/bin/zsh
# macOS L2TP/IPSec VPN auto-connect & watchdog.
#
# Configure via environment variables (see --help) or edit the defaults below.
# No credentials are stored in this file: the password and IPSec shared secret
# are either prompted for at runtime or read from the macOS Keychain.

set -u

SERVICE_NAME="${SERVICE_NAME:-MyVPN}"        # Name of the VPN service in System Settings
VPN_USER="${VPN_USER:-vpn-user}"             # VPN account name
VPN_PEER="${VPN_PEER:-10.0.0.1}"             # Peer/gateway IP used for reachability checks
CHECK_INTERVAL_SECONDS="${CHECK_INTERVAL_SECONDS:-30}"
RETRY_INTERVAL_SECONDS="${RETRY_INTERVAL_SECONDS:-60}"

VPN_PASSWORD="${VPN_PASSWORD:-}"
VPN_SECRET="${VPN_SECRET:-}"
USE_KEYCHAIN=0
EXITING=0

usage() {
  cat <<EOF
Usage:
  ./macos-vpn-autoconnect.sh [--use-keychain]

Default behavior:
  Prompts once for the VPN password and IPSec shared secret, starts the
  "$SERVICE_NAME" L2TP VPN, then keeps monitoring it. If it disconnects, the
  script retries every $RETRY_INTERVAL_SECONDS seconds until the VPN is
  reachable again.

Options:
  --use-keychain   Do not prompt. Ask macOS to use credentials already saved
                   in System Settings / Keychain.

Environment overrides:
  SERVICE_NAME               Current: $SERVICE_NAME
  VPN_USER                   Current: $VPN_USER
  VPN_PEER                   Current: $VPN_PEER
  CHECK_INTERVAL_SECONDS     Current: $CHECK_INTERVAL_SECONDS
  RETRY_INTERVAL_SECONDS     Current: $RETRY_INTERVAL_SECONDS
  VPN_PASSWORD               Optional; avoids password prompt
  VPN_SECRET                 Optional; avoids shared secret prompt
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --use-keychain)
      USE_KEYCHAIN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

prompt_hidden() {
  local message="$1"
  /usr/bin/osascript \
    -e "set r to display dialog \"$message\" default answer \"\" with hidden answer buttons {\"Cancel\", \"OK\"} default button \"OK\" cancel button \"Cancel\"" \
    -e 'text returned of r'
}

ensure_credentials() {
  if [[ "$USE_KEYCHAIN" -eq 1 ]]; then
    return 0
  fi

  if [[ -z "$VPN_PASSWORD" ]]; then
    VPN_PASSWORD="$(prompt_hidden "Enter the VPN password for $SERVICE_NAME user $VPN_USER")" || exit 1
  fi

  if [[ -z "$VPN_SECRET" ]]; then
    VPN_SECRET="$(prompt_hidden "Enter the IPSec Shared Secret for $SERVICE_NAME")" || exit 1
  fi

  if [[ -z "$VPN_PASSWORD" || -z "$VPN_SECRET" ]]; then
    echo "VPN password and shared secret are required." >&2
    exit 1
  fi
}

start_vpn() {
  log "Starting VPN service: $SERVICE_NAME"

  /usr/sbin/scutil --nc select "$SERVICE_NAME" >/dev/null 2>&1 || true

  if [[ "$USE_KEYCHAIN" -eq 1 ]]; then
    /usr/sbin/scutil --nc start "$SERVICE_NAME" >/dev/null 2>&1
  else
    /usr/sbin/scutil --nc start "$SERVICE_NAME" \
      --user "$VPN_USER" \
      --password "$VPN_PASSWORD" \
      --secret "$VPN_SECRET" >/dev/null 2>&1
  fi
}

vpn_interface_up() {
  local iface

  for iface in $(/sbin/ifconfig -l | tr ' ' '\n' | grep '^ppp'); do
    /sbin/ifconfig "$iface" 2>/dev/null | grep -q -- "--> $VPN_PEER" && return 0
  done

  return 1
}

vpn_reachable() {
  vpn_interface_up || return 1
  /sbin/ping -c 1 -t 2 "$VPN_PEER" >/dev/null 2>&1
}

wait_until_connected() {
  while true; do
    if vpn_reachable; then
      log "VPN is connected and $VPN_PEER is reachable."
      return 0
    fi

    log "VPN is not reachable. Retrying in ${RETRY_INTERVAL_SECONDS}s."
    sleep "$RETRY_INTERVAL_SECONDS"
    start_vpn
  done
}

stop_vpn() {
  log "Stopping VPN service: $SERVICE_NAME"
  /usr/sbin/scutil --nc stop "$SERVICE_NAME" >/dev/null 2>&1 || true
}

cleanup_and_exit() {
  local exit_code="${1:-130}"

  if [[ "$EXITING" -eq 0 ]]; then
    EXITING=1
    stop_vpn
  fi

  exit "$exit_code"
}

trap 'cleanup_and_exit 130' INT
trap 'cleanup_and_exit 143' TERM
trap 'cleanup_and_exit 129' HUP

ensure_credentials
start_vpn
wait_until_connected

while true; do
  sleep "$CHECK_INTERVAL_SECONDS"

  if vpn_reachable; then
    continue
  fi

  log "VPN appears disconnected."
  start_vpn
  wait_until_connected
done
