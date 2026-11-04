#!/bin/bash

# Recursively obtain child processes
recursive_pgrep() {
    local parent_pid="$1"
    local kid
    local kids

    kids=$(pgrep -P "$parent_pid" 2>/dev/null || true)

    for kid in $kids; do
        pid_to_check="$pid_to_check $kid"
        recursive_pgrep "$kid"
    done
}

# Validate arguments
if [[ $# -lt 1 ]]; then
    echo "ERROR: no command specified."
    exit 2
fi

# Last argument identifies the task
task_id="${@: -1}"

memory_file="memory.${task_id}.txt"
runtime_file="runtime.${task_id}.txt"
status_file="exit_status.${task_id}.txt"

start=$(date +%s)

echo "0" > "$memory_file"
echo "RUNNING" > "$status_file"

# Execute command in background
"$@" &
main_pid=$!

echo "task_wrapper.sh: monitoring PID $main_pid"

# Monitor memory while the main process is running
while true; do

    sleep 1

    # Check whether the main process still exists
    ps -p "$main_pid" -o rss= >/dev/null 2>&1 || break

    # Obtain all relevant PIDs
    pid_to_check="$main_pid"
    recursive_pgrep "$main_pid"

    # Calculate total RSS memory (KB)
    total_mem=0

    for pid in $pid_to_check; do

        rss=$(ps -p "$pid" -o rss= 2>/dev/null | tr -d ' ')

        if [[ "$rss" =~ ^[0-9]+$ ]]; then
            total_mem=$((total_mem + rss))
        fi

    done

    # Preserve maximum observed memory
    previous_mem=$(cat "$memory_file" 2>/dev/null || echo 0)

    if (( total_mem > previous_mem )); then
        echo "$total_mem" > "$memory_file"
    fi

    # Save memory history
    echo "$total_mem" >> "$main_pid.${task_id}.mem"

done

# Recover the actual exit code of the executed command
wait "$main_pid"
task_exit_code=$?

# Record execution status
echo "$task_exit_code" > "$status_file"

stop=$(date +%s)
echo $((stop - start)) > "$runtime_file"

echo "task_wrapper.sh: task completed"
echo "PID: $main_pid"
echo "Maximum RSS: $(cat "$memory_file") KB"
echo "Runtime: $(cat "$runtime_file") seconds"
echo "Exit status: $task_exit_code"

if (( task_exit_code > 0 )); then
    echo "task_wrapper.sh: ERROR - command failed."
else
    echo "task_wrapper.sh: SUCCESS - command completed."
fi

# Propagate the original exit code
exit "$task_exit_code"
