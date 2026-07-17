#!/bin/bash

# Get sorted unique indices
indices=$(ls slow_pullf*.xvg | sed -E 's/.*f([0-9]+)\.xvg/\1/' | sort -n)

new_index=0

for old_index in $indices; do

    # Rename to temporary names first
    mv "slow_pullf${old_index}.xvg" "tmpf_${new_index}.xvg"
    
    if [ -f "slow_pullx${old_index}.xvg" ]; then
        mv "slow_pullx${old_index}.xvg" "tmpx_${new_index}.xvg"
    fi

    ((new_index++))
done

# Rename temporary files back
for file in tmpf_*.xvg; do
    num=$(echo "$file" | sed -E 's/tmpf_([0-9]+)\.xvg/\1/')
    mv "$file" "slow_pullf${num}.xvg"
done

for file in tmpx_*.xvg; do
    num=$(echo "$file" | sed -E 's/tmpx_([0-9]+)\.xvg/\1/')
    mv "$file" "slow_pullx${num}.xvg"
done
