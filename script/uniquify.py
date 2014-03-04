#!/usr/bin/python
import sys

edges = set()
self = 0
multiple = 0
with open(sys.argv[1]) as f:
  for line in f:
    if not line.startswith("#"):
      s = line.split()
      if s[0] < s[1]:
        if (s[0], s[1]) not in edges:
          edges.add((s[0], s[1]))
        else:
          multiple += 1
      elif s[0] > s[1]:
        if (s[1], s[0]) not in edges:
          edges.add((s[1], s[0]))
        else:
          multiple += 1
      else:
        self += 1

for (a,b) in edges:
  print a, b

if self != 0:
  sys.stderr.write("Self loops: " + str(self))
if multiple != 0:
  sys.stderr.write("Multiple edges: " + str(multiple))
