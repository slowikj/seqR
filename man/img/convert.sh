#!/bin/bash

for filename in ./*.pdf; do
	without_path=${filename%.*}
	without_path=$(basename $without_path)
	pdftoppm $filename $without_path -png
	echo $filename
	echo $without_path
	echo "end"
done
