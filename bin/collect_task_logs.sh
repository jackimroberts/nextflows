#!/bin/bash

# Nextflow Task Log Collector
# Collect all .command.out files, sort by timestamp, and deduplicate process blocks

LAUNCH_DIR="$1"
if [[ -z "$LAUNCH_DIR" ]]; then
    echo "Usage: $0 <nextflow_launch_directory> <output_directory>"
    exit 1
fi

OUT_DIR="$2"
if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR=$LAUNCH_DIR
fi

mkdir -p $OUT_DIR

# Find work directory
work_pattern="$LAUNCH_DIR"
[[ -d "$LAUNCH_DIR/work" ]] && work_pattern="$LAUNCH_DIR/work"
[[ -d "$LAUNCH_DIR/../../work" ]] && work_pattern="$LAUNCH_DIR/../../work"

# Collect all .command.out files, sort by timestamp, and deduplicate process blocks
{
    for out_file in "$work_pattern"/*/*/.command.out; do
        [[ ! -s "$out_file" ]] && continue
        done_time=$(stat -c %Y "$out_file")
        echo "$done_time|$out_file"
    done
} | sort -t'|' -k1,1n | \
while IFS='|' read -r done_time out_file; do
    timestamp=$(date -d "@$done_time" "+%Y-%m-%d %H:%M:%S")
    nextflow_id=$(echo "$out_file" | sed 's|.*/work/||' | cut -c1-9)
    echo "=== Nextflow ID: $nextflow_id $timestamp ==="
    cat "$out_file"
    echo
done | \
awk -v RS="====== PROCESS_SUMMARY" 'NR==1{print $0} NR>1 && !seen[$0]++{print $0}' > "$OUT_DIR/task_logs.txt"

# Collect all .command.err files, sort by timestamp, and deduplicate process blocks
{
    for err_file in "$work_pattern"/*/*/.command.err; do
        [[ ! -s "$err_file" ]] && continue
        done_time=$(stat -c %Y "$err_file")
        echo "$done_time|$err_file"
    done
} | sort -t'|' -k1,1n | \
while IFS='|' read -r done_time err_file; do
    timestamp=$(date -d "@$done_time" "+%Y-%m-%d %H:%M:%S")
    nextflow_id=$(echo "$err_file" | sed 's|.*/work/||' | cut -c1-9)
    echo "============================================"
    echo "=== Nextflow ID: $nextflow_id $timestamp ==="
    cat "$err_file"
    echo
done > "$OUT_DIR/error_logs.txt"
