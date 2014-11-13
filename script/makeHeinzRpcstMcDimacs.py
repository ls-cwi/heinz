#!/usr/bin/python
import sys

if len(sys.argv) != 5:
    sys.stderr.write("Usage: " + sys.argv[0] + " <instances> <full_output_dir> <timelimit_pbs> <timelimit_heinz>\n")
    sys.exit(1)

lines = open(sys.argv[1]).readlines()
n = len(lines)
bins = n / 15
if n % 15 != 0:
  bins += 1

full_output_dir = sys.argv[2]
timelimit_pbs = sys.argv[3]
timelimit_heinz = sys.argv[4]

for i in range(bins):
  with open(str(i) + ".pbs", "w") as f:
    f.write("#PBS -N " + str(i)+"\n")
    f.write("#PBS -o " + str(i) + ".out\n")
    f.write("#PBS -e " + str(i) + ".err\n")
    f.write("#PBS -lwalltime=" + timelimit_pbs + "\n")
    f.write("#PBS -lnodes=1:cpu3\n")
    f.write("cd ~/DIMACS2014/\n")

    nn = 15
    if i == bins - 1:
      nn = n % 15
    for j in range(nn):
      idx = 15 * i + j
      s = lines[idx].rstrip("\n").split(" ")
      
      filename1 = full_output_dir + "/" + "4-" + s[0]
      f.write("( /usr/bin/time -o " + filename1 + ".time bin/heinz_rpcst_mc " + s[2] + " " + timelimit_heinz + " 2 " + filename1 + ".dimacs" + " > " + filename1 + ".out 2> " + filename1 + ".err ) &\n")
    f.write("wait\n")
