#!/usr/bin/python
import sys
import random
import os

# sys.argv[1] # samples
# sys.argv[2] tmp file
# sys.argv[3] nodes file
# sys.argv[4] heinz command

n = int(sys.argv[1])
tmp_file = sys.argv[2]
tmp_file_err = tmp_file + ".err"
nodes_file = sys.argv[3]
heinz_cmd = sys.argv[4]

# read nodes
scores = {}
f = open(nodes_file)
for line in f:
  if not line.startswith("#"):
    s = line.split()
    scores[s[0]] = s[1]

# shuffe
vals = scores.values()
key = scores.keys()[0]
UB_list = []

count = 0
for i in range(n):
  random.shuffle(vals)
  f = open(tmp_file, "w")
  for j, node in enumerate(scores):
    f.write(node + " " + vals[j] + "\n")
  f.close()

  # run heinz
  os.system(heinz_cmd + " -n " + tmp_file + " > /dev/null 2> " + tmp_file_err)
  #print 'grep "\[" ' + tmp_file_err + ' | sed -e "s/.*, \([0-9]*\.\?[0-9]*\)]/\\1/g" | sort -n -r | head -n 1'
  UB = float(os.popen('grep "\[" ' + tmp_file_err + ' | sed -e "s/.*, \([0-9]*\.\?[0-9]*\)]/\\1/g" | sort -n -r | head -n 1').read())
  UB_list.append(UB)
  print UB
#  if LB > UB:
#    count += 1

#print "p =", float(count) / float(n)
