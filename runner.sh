#!/bin/bash

# start and end values for the numbers
start=1
end=3

# step value for the numbers
step=0.1

# initialize the list
list=()

# loop through the numbers and convert them to strings
i=0
while [ $(bc <<< "$start + $i * $step <= $end") -eq 1 ]
do
  num=$(bc <<< "$start + $i * $step")
  str=$(printf "%.1f" "$num")
  list+=("$str")
  i=$((i+1))
done

# print the list
# echo "${list[@]}"

# name of the Python file
python_file="vis.py"

# arrays of arguments to pass to the Python file
arg1="50"
arg2=list
arg3="G"

# loop through the arrays of arguments and run the Python file with each set of arguments
for ((i=0; i<${#arg1[@]}; i++))
do
  python "$python_file" "${arg1}" "${arg2[i]}" "${arg3}"
done
