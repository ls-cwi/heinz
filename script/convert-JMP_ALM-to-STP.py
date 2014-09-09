#!/usr/bin/python

# convert dat format used by Alvarez et al. to STP format used in DIMACS challenge 2014
# GW Klau, CWI

import sys
import os

f = open(sys.argv[1])
i = 0
w = [] # weights of the nodes
E = [] # edges

for line in f:
   if line.startswith("node"): break
next (f) ## skip comment

#parse nodes
for line in f:
    if line.startswith("link"): break
    l = line.split()
    assert int(l[0]) == i, "node number mismatch"
    w.append(l[3])
    print i, w[i]
    i = i + 1    

next (f) ## skip comment

#parse edges
for line in f:
    l = line.split()
    print l[1], l[2]
    E.append([l[1], l[2]])
    

print "33D32945 STP File, STP Format Version 1.0\n"

print "SECTION Comment"
print "Name\t\"" + os.path.splitext(sys.argv[1])[0] + "\""
print "Creator\t\"CWI_LS_GWK\""
print "Comment\t\"converted from format used in Alvarez-Miranda et al., based on the random geometric instances of Johnson et al.\""
print "END\n"

print "SECTION Graph"
print "Nodes ", len(w)
print "Edges ", len(E)
for e in E:
    print "E ", e[0], e[1]
print "END"

print "SECTION Terminals"
print "Terminals ", len(w)
for i, v in enumerate(w):
    print "T ", i, v
print "END\n"

print "EOF"




