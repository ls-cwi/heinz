#!/usr/bin/python
import sys
import subprocess

if len(sys.argv) != 5:
    sys.stderr.write("Usage: " + sys.argv[0] + " <executable> <check_executable> <input_file> <output_file>\n")
    sys.exit(1)

executable = sys.argv[1]
check_executable = sys.argv[2]
input_file = sys.argv[3]
output_file = sys.argv[4]
time_limit = 10
threads = 1

print "Running " + executable + " on " + input_file
print(executable + " " + input_file + " " + str(time_limit) + " " + str(threads) + " " + output_file)
print(executable + " " + input_file + " " + str(time_limit) + " " + str(threads) + " " + output_file)
status = subprocess.call(executable + " -stp " + input_file + " -t " + str(time_limit) + " -m " + str(threads) + " -o " + output_file, shell=True)
if status != 0:
    sys.exit(status)
else:
    print "Checking solution " + output_file
    print(check_executable + " -stp " + input_file + " -s " + output_file)
    status = subprocess.call(check_executable + " -stp " + input_file + " -s " + output_file, shell=True)
sys.exit(status)
