#!/usr/bin/python
import sys

with open(sys.argv[1]) as f:
  sys.stdout.write('"'+sys.argv[1]+'"')
  sys.stdout.write(',"'+sys.argv[1].split("/")[2]+'"')
  for line in f:
    s = line.rstrip("\n").split()
    if len(s) == 0:
      continue
    if s[0] == "Name":
      sys.stdout.write(","+s[1])
    elif s[0] == "Problem":
      sys.stdout.write(","+s[1])
    elif s[0] == "Program":
      sys.stdout.write(","+s[1])
    elif s[0] == "Threads":
      sys.stdout.write(","+s[1])
    elif s[0] == "Time":
      sys.stdout.write(","+s[1])
    elif s[0] == "Dual":
      sys.stdout.write(","+s[1])
    elif s[0] == "Primal":
      sys.stdout.write(","+s[1])
sys.stdout.write("\n")
