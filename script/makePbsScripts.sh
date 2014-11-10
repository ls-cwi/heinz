#!/bin/bash
#1 STP input dir
#2 python make*Pbs script
#3 result subdir
#4 PBS timelimit
#5 Heinz timelimit
for f in $1/*.stp
do
	echo $f
	$2 `basename $f .stp` $f `dirname $f | sed -e "s/heinz\/data\/DIMACS/result\/$3/g"` $4 $5 > `basename $f .stp`.pbs
done
