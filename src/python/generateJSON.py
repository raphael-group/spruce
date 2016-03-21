#!/usr/bin/python
import sys
import os
import subprocess

if len(sys.argv) != 3:
    sys.stderr.write("Usage " + sys.argv[0] + " <VISUALIZE_EXECUTABLE> <SOLUTION>\n")
    sys.exit(1)

viz_exec = sys.argv[1]
sol_file = sys.argv[2]

#subprocess.call([viz_exec, "-j", sol_file])
p = subprocess.Popen(viz_exec.split() + ["-s", sol_file], stdout=subprocess.PIPE)
q = subprocess.Popen(["dot", "-Tsvg"], stdin=p.stdout, stdout=subprocess.PIPE)

for i in range(6):
    q.stdout.readline()

r = subprocess.Popen(viz_exec.split() + ["-j", sol_file], stdout=subprocess.PIPE)
sys.stdout.write(r.stdout.readline())
sys.stdout.write('\t"svg": "')
for line in q.stdout:
    sys.stdout.write(line.rstrip("\n").replace('"', '\\"'))
print '",'

for line in r.stdout:
    sys.stdout.write(line)

