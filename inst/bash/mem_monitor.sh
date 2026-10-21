#!/usr/bin/env bash
# Usage: mem_monitor.sh PID [OUTFILE] [INTERVAL_SEC] [PHASE_FILE]
#
# Samples the RAM usage of process PID every INTERVAL_SEC seconds and
# appends it to OUTFILE as CSV: time,phase,rss_mb,peak_mb
#
# Stops when:
#   - the process PID exits, or
#   - PHASE_FILE (if given) is deleted.

if [ -z "$1" ]; then
  echo "Usage: $0 PID [OUTFILE] [INTERVAL_SEC] [PHASE_FILE]" >&2
  exit 1
fi

pid=$1
out=${2:-mem_trace.csv}
dt=${3:-1}
phase_file=${4:-}

echo "time,phase,rss_mb,peak_mb" > "$out"

while kill -0 "$pid" 2>/dev/null; do
  phase=""
  if [ -n "$phase_file" ]; then
    [ -f "$phase_file" ] || break
    phase=$(head -n1 "$phase_file" 2>/dev/null | tr -d ',"\r')
  fi

  awk -v t="$(date +%s.%N)" -v ph="$phase" '
    /^VmRSS:/ {r=$2}
    /^VmHWM:/ {h=$2}
    END { if (r) printf "%s,%s,%.1f,%.1f\n", t, ph, r/1024, h/1024 }
  ' "/proc/$pid/status" 2>/dev/null >> "$out"

  sleep "$dt"
done
