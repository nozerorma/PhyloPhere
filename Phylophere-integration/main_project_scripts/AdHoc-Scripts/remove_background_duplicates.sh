#!/bin/bash


for i in $(find . -name "background-*"); do
    echo $i
    cat $i >> global_background.txt
done

# Global background consists on a list of genes separated by \n
# Remove duplicates
sort global_background.txt | uniq > global_background_no_duplicates.txt
