#!/bin/bash

recursive_pgrep ()
{
    kids=$(pgrep -P $1)
    for kid in $kids
    do
        pid_to_check=$pid_to_check" "$kid
        recursive_pgrep $kid
    done
}
echo "2+2" | bc
start=`date +%s`
"$@" & # run the command in background
main_pid=$! # get the corresponding pid
echo "0" > memory.${@: -1}.txt
while true; do
  sleep 1
  still_running="$(ps -o rss= $main_pid 2> /dev/null)" || break
  # obtain all relevant pid
  pid_to_check=$main_pid
  recursive_pgrep $main_pid
  # sum their rss values
  total_mem=0
  for pid in $pid_to_check
  do
      total_mem=$(echo $total_mem"+"$(ps -p $pid -o rss | grep -v RSS) | bc 2>/dev/null)
  done
  cat memory.${@: -1}.txt > $main_pid.${@: -1}.mem
  echo $total_mem >> $main_pid.${@: -1}.mem
  sort -nr $main_pid.${@: -1}.mem | head -n 1 > memory.${@: -1}.txt
done
echo "task_wrapper.sh:task completed"
stop=`date +%s`
echo $((stop-start)) > runtime.${@: -1}.txt
