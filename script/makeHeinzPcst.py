#!/usr/bin/python
import sys

if len(sys.argv) != 5:
    sys.stderr.write("Usage: " + sys.argv[0] + " <stp_filename> <stp_ful_filename> <timelimit_pbs> <timelimit_heinz>\n")
    sys.exit(1)

filename = sys.argv[1]
filename1 = "1-" + filename
filename2 = "2-" + filename
filename3 = "3-" + filename
full_filename = sys.argv[2]
timelimit_pbs = sys.argv[3]
timelimit_heinz = sys.argv[4]

print "#PBS -N", filename
print "#PBS -o", filename + ".out"
print "#PBS -e", filename + ".err"
print "#PBS -l", "walltime=" + timelimit_pbs
print "cd ~/DIMACS2014/"
print "( /usr/bin/time -o " + filename1 + ".time bin/heinz -no-enum -p -m 2 -t " + timelimit_heinz + " -stp-pcst " + full_filename + " -o " + filename1 + ".hnz" + " > " + filename1 + ".out 2> " + filename1 + ".err ) &"
print "( /usr/bin/time -o " + filename2 + ".time bin/heinz -no-enum -m 2 -t " + timelimit_heinz + " -stp-pcst " + full_filename + " -o " + filename2 + ".hnz" + " > " + filename2 + ".out 2> " + filename2 + ".err ) &"
print "( /usr/bin/time -o " + filename3 + ".time bin/heinz -m 2 -t " + timelimit_heinz + " -stp-pcst " + full_filename + " -o " + filename3 + ".hnz" + " > " + filename3 + ".out 2> " + filename3 + ".err ) &"
print "wait"
