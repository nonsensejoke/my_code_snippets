#!/bin/bash

TOPIC="here-you-define-ntfy-topic"

LAST_FILE="$HOME/.last_ip"

# Wi-Fi / LAN IP
LAN_IP=$(ipconfig getifaddr en0)

# VPN IP (utun5)
VPN_IP=$(ipconfig getifaddr utun5 2>/dev/null)
# ipconfig 对 utun 接口有时取不到，回退到 ifconfig
[ -z "$VPN_IP" ] && VPN_IP=$(ifconfig utun5 2>/dev/null | awk '/inet /{print $2}')

# Wi-Fi 未连接
[ -z "$LAN_IP" ] && exit 0

# 用一行文本记录两个 IP 的状态
CURRENT="LAN=${LAN_IP} VPN=${VPN_IP:-none}"

LAST=""
[ -f "$LAST_FILE" ] && LAST=$(cat "$LAST_FILE")

if [ "$CURRENT" != "$LAST" ]; then
    echo "$CURRENT" > "$LAST_FILE"

    curl \
      -H "Title: $(hostname) IP Changed" \
      -d "LAN: ${LAN_IP}
VPN: ${VPN_IP:-未连接}" \
      https://ntfy.sh/$TOPIC
fi
