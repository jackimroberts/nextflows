#!/bin/bash

# Nextflow Resource Usage Analyzer
# Compares requested vs actual resource usage for SLURM jobs

LAUNCH_DIR="$1"
if [[ -z "$LAUNCH_DIR" ]]; then
    echo "Usage: $0 <nextflow_launch_directory> <output_directory>"
    exit 1
fi

OUT_DIR="$2"
if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR=$LAUNCH_DIR
fi

echo "=== Nextflow Resource Usage Analysis ==="
echo "Launch directory: $LAUNCH_DIR"
echo "Analysis date: $(date)"
echo

# Convert memory strings to MB
convert_to_mb() {
    local mem="$1"
    if [[ "$mem" =~ ([0-9.]+)([KMGT]?B?) ]]; then
        local num="${BASH_REMATCH[1]}"
        local unit="${BASH_REMATCH[2]}"
        case "$unit" in
            KB|K) echo "$(echo "$num / 1024" | bc -l)" ;;
            MB|M|"") echo "$num" ;;
            GB|G) echo "$(echo "$num * 1024" | bc -l)" ;;
            TB|T) echo "$(echo "$num * 1024 * 1024" | bc -l)" ;;
            *) echo "$num" ;;
        esac
    else
        echo "0"
    fi
}


# Parse memory usage and return percentage
parse_memory_usage() {
    local max_rss="$1"
    local mem_req_mb="$2"
    
    if [[ -z "$max_rss" || -z "$mem_req_mb" || "$mem_req_mb" -le 0 ]]; then
        echo "?"
        return
    fi
    
    local mem_used_mb=""
    if [[ "$max_rss" =~ ([0-9]+)K ]]; then
        mem_used_mb=$(echo "${BASH_REMATCH[1]} / 1024" | bc -l)
    elif [[ "$max_rss" =~ ([0-9]+)M ]]; then
        mem_used_mb="${BASH_REMATCH[1]}"
    elif [[ "$max_rss" =~ ([0-9]+)G ]]; then
        mem_used_mb=$(echo "${BASH_REMATCH[1]} * 1024" | bc -l)
    fi
    
    if [[ -n "$mem_used_mb" ]]; then
        local mem_efficiency=$(echo "scale=1; $mem_used_mb * 100 / $mem_req_mb" | bc -l)
        printf "%.1f%%" "$mem_efficiency"
    else
        echo "?"
    fi
}

# Parse SLURM TotalCPU format to seconds
parse_total_cpu() {
    local total_cpu="$1"
    
    if [[ "$total_cpu" =~ ^([0-9]+):([0-9]+):([0-9]+)$ ]]; then
        # HH:MM:SS format
        local hours="${BASH_REMATCH[1]}"
        local mins="${BASH_REMATCH[2]}"
        local secs="${BASH_REMATCH[3]}"
        echo "$hours * 3600 + $mins * 60 + $secs" | bc -l
    elif [[ "$total_cpu" =~ ^([0-9]+):([0-9]+)\.([0-9]+)$ ]]; then
        # MM:SS.FFF format
        local mins="${BASH_REMATCH[1]}"
        local secs="${BASH_REMATCH[2]}"
        local frac="${BASH_REMATCH[3]}"
        echo "$mins * 60 + $secs.$frac" | bc -l
    elif [[ "$total_cpu" =~ ^([0-9]+)\.([0-9]+)$ ]]; then
        # SS.FFF format
        echo "$total_cpu"
    else
        echo "0"
    fi
}

# Calculate CPU efficiency from SLURM data
calculate_cpu_efficiency() {
    local total_cpu_time="$1"
    local allocated_cpu_time="$2"
    
    if [[ -z "$total_cpu_time" || -z "$allocated_cpu_time" || "$allocated_cpu_time" -le 0 ]]; then
        echo "unknown"
        return
    fi
    
    local cpu_efficiency=$(echo "scale=1; $total_cpu_time * 100 / $allocated_cpu_time" | bc -l)
    printf "%.1f%%" "$cpu_efficiency"
}

# Format time duration
format_duration() {
    local seconds="$1"
    local suffix="$2"
    
    if [[ -z "$seconds" || "$seconds" -le 0 ]]; then
        echo "?"
        return
    fi
    
    local hours=$((seconds / 3600))
    local mins=$(((seconds % 3600) / 60))
    local secs=$((seconds % 60))
    printf "%02d:%02d:%02d%s" $hours $mins $secs "${suffix:-}"
}

# Get job status from exit code
get_job_status() {
    local task_dir="$1"
    local job_id="$2"
    
    if [[ -f "$task_dir/.exitcode" ]]; then
        local exit_code=$(cat "$task_dir/.exitcode" 2>/dev/null | tr -d '\n\r' | head -c 10)
        case "$exit_code" in
            0) echo "COMPLETED" ;;
            127) echo "CMD_NOT_FOUND" ;;
            137) echo "OOM_KILLED" ;;
            143) echo "USER_CANCEL" ;;
            "") echo "COMPLETED" ;;
            *) echo "FAILED($exit_code)" ;;
        esac
    elif [[ "$job_id" =~ ^[0-9]+$ ]]; then
        echo "RUNNING"
    else
        echo "UNKNOWN"
    fi
}

# Header
printf "%-20s %-6s %-8s %-4s %-7s %-10s %-10s %-10s %-10s %-12s\n" \
    "Process" "Nf_ID" "SLURM_ID" "CPUs" "Mem(GB)" "%_Mem_Use" "%_CPU_Use" "Time_Use" "Queue_Wait" "Status"
echo "$(printf '%*s' 106 '' | tr ' ' '-')"

# Find work directory
work_pattern="$LAUNCH_DIR"
[[ -d "$LAUNCH_DIR/work" ]] && work_pattern="$LAUNCH_DIR/work"
[[ -d "$LAUNCH_DIR/../../work" ]] && work_pattern="$LAUNCH_DIR/../../work"

# Collect job data
declare -a job_data=()

for task_dir in "$work_pattern"/*/*; do
    [[ ! -d "$task_dir" || ! -f "$task_dir/.command.run" ]] && continue
    
    # Get submit time (for sorting)
    submit_time=$(stat -c %Y "$task_dir/.command.run" 2>/dev/null || echo "9999999999")
    
    # Extract basic info
    dir_path=$(echo "$task_dir" | sed 's|.*/work/||')
    nextflow_id="${dir_path:0:6}"
    
    # Get process name
    process_name="$nextflow_id"
    if process=$(grep -o "### name: '[^']*'" "$task_dir/.command.run" 2>/dev/null | sed "s/### name: '//" | sed "s/'//"); then
        [[ -n "$process" ]] && process_name="$process"
        # Remove subworkflow prefix (everything before and including the colon)
        process_name="${process_name##*:}"
    fi
    
    # Check if SLURM job
    is_slurm=false
    grep -q "#SBATCH\|sbatch" "$task_dir/.command.run" 2>/dev/null && is_slurm=true
    
    # Get SLURM job ID
    job_id="none"
    if slurm_file=$(find "$task_dir" -name "slurmjob-*.out-*" 2>/dev/null | head -1); then
        job_id=$(basename "$slurm_file" | grep -o 'slurmjob-[0-9]\+' | grep -o '[0-9]\+')
    fi
    
    # Parse resource requests
    cpus_req=$(grep -o '#SBATCH.*--cpus-per-task[= ]*[0-9]\+\|#SBATCH.*-c[[:space:]]*[0-9]\+' "$task_dir/.command.run" 2>/dev/null | grep -o '[0-9]\+$' | head -1 || echo "1")
    mem_req=$(grep -o '#SBATCH.*--mem[[:space:]]*[0-9]\+[KMGT]\?' "$task_dir/.command.run" 2>/dev/null | grep -o '[0-9]\+[KMGT]\?$' | head -1 || echo "1G")
    mem_req_mb=$(convert_to_mb "$mem_req")
    mem_req_gb=$(printf "%.1f" "$(echo "$mem_req_mb / 1024" | bc -l)")
    
    # Initialize defaults
    mem_usage="?"
    cpu_usage="?"
    elapsed="?"
    queue_wait="?"
    status="?"
    
    if [[ "$job_id" =~ ^[0-9]+$ ]]; then
        # Get job status
        status=$(get_job_status "$task_dir" "$job_id")
        
        # Get SLURM stats from batch step
        if job_info=$(sacct -j "$job_id" --format=JobID,State,Elapsed,MaxRSS,CPUTimeRAW,TotalCPU --parsable2 --noheader 2>/dev/null | grep "^$job_id\.batch" | head -1); then
            IFS='|' read -r _ slurm_state elapsed_sacct max_rss cpu_time_raw total_cpu <<< "$job_info"
            
            # Parse memory usage
            mem_usage=$(parse_memory_usage "$max_rss" "$mem_req_mb")
            
            # Calculate CPU efficiency
            if [[ -n "$total_cpu" && -n "$cpu_time_raw" && "$cpu_time_raw" -gt 0 ]]; then
                total_cpu_sec=$(parse_total_cpu "$total_cpu")
                [[ "$total_cpu_sec" != "0" ]] && cpu_usage=$(calculate_cpu_efficiency "$total_cpu_sec" "$cpu_time_raw")
            fi
            
            # Use SLURM elapsed time if available
            [[ -n "$elapsed_sacct" ]] && elapsed="$elapsed_sacct"
        fi
        
        # Calculate queue wait time
        if [[ -f "$task_dir/.command.begin" ]]; then
            start_time=$(stat -c %Y "$task_dir/.command.begin" 2>/dev/null)
            if [[ -n "$start_time" && "$start_time" -gt "$submit_time" ]]; then
                queue_wait=$(format_duration $((start_time - submit_time)))
            fi
        elif [[ "$is_slurm" == true ]]; then
            # Still in queue
            queue_wait=$(format_duration $(($(date +%s) - submit_time)) "*")
        fi
        
        # Calculate elapsed time from file timestamps if not from SLURM
        if [[ "$elapsed" == "unknown" && -f "$task_dir/.command.begin" ]]; then
            start_time=$(stat -c %Y "$task_dir/.command.begin" 2>/dev/null)
            if [[ -n "$start_time" ]]; then
                if [[ "$status" == "RUNNING" ]]; then
                    elapsed=$(format_duration $(($(date +%s) - start_time)) "*")
                elif [[ -f "$task_dir/.exitcode" ]]; then
                    end_time=$(stat -c %Y "$task_dir/.exitcode" 2>/dev/null)
                    [[ -n "$end_time" ]] && elapsed=$(format_duration $((end_time - start_time)))
                fi
            fi
        fi
    else
        # Handle non-SLURM jobs
        if [[ "$is_slurm" == false ]]; then
            status="local"
            job_id="local"
            queue_wait="N/A"
            
            # Calculate elapsed time for local jobs (submit to completion)
            if [[ -f "$task_dir/.exitcode" ]]; then
                end_time=$(stat -c %Y "$task_dir/.exitcode" 2>/dev/null)
                [[ -n "$end_time" ]] && elapsed=$(format_duration $((end_time - submit_time)))
            fi
        else
            status="CANCEL_QUEUED"
            job_id="none"
            # Calculate time spent in queue before cancellation
            if [[ -f "$task_dir/.command.log" ]]; then
                cancel_time=$(stat -c %Y "$task_dir/.command.log" 2>/dev/null)
            elif [[ -f "$task_dir/.command.err" ]]; then
                cancel_time=$(stat -c %Y "$task_dir/.command.err" 2>/dev/null)
            fi
            [[ -n "$cancel_time" && "$cancel_time" -gt "$submit_time" ]] && queue_wait=$(format_duration $((cancel_time - submit_time)))
        fi
    fi
    
    # Store job data with submit_time for sorting
    job_data+=("$submit_time|${process_name:0:20}|$nextflow_id|$job_id|$cpus_req|$mem_req_gb|$mem_usage|$cpu_usage|$elapsed|$queue_wait|$status")
done

# Sort by submit_time and display (skip submit_time in output)
IFS=$'\n' sorted_data=($(sort -t'|' -k1,1n <<< "${job_data[*]}"))
{
for line in "${sorted_data[@]}"; do
    IFS='|' read -r submit_time process nextflow_id job_id cpus_req mem_req_gb mem_usage cpu_usage elapsed queue_wait status <<< "$line"
    printf "%-20s %-6s %-8s %-4s %-7s %-10s %-10s %-10s %-10s %-12s\n" \
        "$process" "$nextflow_id" "$job_id" "$cpus_req" "$mem_req_gb" "$mem_usage" "$cpu_usage" "$elapsed" "$queue_wait" "$status"
done

echo
echo "=== Summary Notes ==="
echo "- CPUs_Req: Number of CPU cores requested"
echo "- Mem_Req: Memory in GB requested"
echo "- %_Mem_Use: Memory efficiency (used/requested * 100)"
echo "- %_CPU_Use: CPU efficiency (actual CPU time/allocated CPU time * 100)"
echo "- Time_Use: Actual runtime (* indicates job still running)"
echo "- Queue_Wait: Time spent waiting in SLURM queue (* indicates still queued)"
echo "- Status: Job completion status"
} > "$OUT_DIR/slurm_analysis.txt"