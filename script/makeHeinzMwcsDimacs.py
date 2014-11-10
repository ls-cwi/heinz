#!/usr/bin/python
import sys

if len(sys.argv) != 6:
    sys.stderr.write("Usage: " + sys.argv[0] + " <stp_filename> <stp_full_filename> <full_output_dir> <timelimit_pbs> <timelimit_heinz>\n")
    sys.exit(1)

filename = sys.argv[1]
full_filename = sys.argv[2]
full_output_dir = sys.argv[3]
filename1 = full_output_dir + "/" + "1-" + filename
filename2 = full_output_dir + "/" + "2-" + filename
filename3 = full_output_dir + "/" + "3-" + filename
timelimit_pbs = sys.argv[4]
timelimit_heinz = sys.argv[5]

print "#PBS -N", filename
print "#PBS -o", filename + ".out"
print "#PBS -e", filename + ".err"
print "#PBS -lwalltime=" + timelimit_pbs
print "#PBS -l:cpu3"
print "cd ~/DIMACS2014/"
print "( /usr/bin/time -o " + filename1 + ".time bin/heinz_mwcs_no_pre " + full_filename + " " + timelimit_heinz + " 2 " + filename1 + ".dimacs" + " > " + filename1 + ".out 2> " + filename1 + ".err ) &"
print "( /usr/bin/time -o " + filename2 + ".time bin/heinz_mwcs_no_dc " + full_filename + " " + timelimit_heinz + " 2 " + filename2 + ".dimacs" + " > " + filename2 + ".out 2> " + filename2 + ".err ) &"
print "( /usr/bin/time -o " + filename3 + ".time bin/heinz_mwcs_dc " + full_filename + " " + timelimit_heinz + " 2 " + filename3 + ".dimacs" + " > " + filename3 + ".out 2> " + filename3 + ".err ) &"
print "wait"
