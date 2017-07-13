#!/usr/bin/python
import sys
import os
import subprocess

if len(sys.argv) != 3:
    sys.stderr.write("Usage " + sys.argv[0] + " <VISUALIZE_EXECUTABLE> <SOLUTION>\n")
    sys.exit(1)

viz_exec = sys.argv[1]
sol_file = sys.argv[2]

output_filename = os.path.basename(sol_file)

#subprocess.call([viz_exec, "-j", sol_file])
p = subprocess.Popen(viz_exec.split() + ["-s", sol_file], stdout=subprocess.PIPE)
q = subprocess.Popen(["dot", "-Tsvg"], stdin=p.stdout, stdout=subprocess.PIPE)

for i in range(6):
    q.stdout.readline()

r = subprocess.Popen(viz_exec.split() + ["-j", sol_file], stdout=subprocess.PIPE)

with open(output_filename + ".json", "w") as f:
    f.write(r.stdout.readline())
    f.write('\t"svg": "')
    for line in q.stdout:
        f.write(line.rstrip("\n").replace('"', '\\"'))
    f.write('",\n')

    for line in r.stdout:
        f.write(line)

with open(output_filename + ".html", "w") as f:
    header = ""
    with open("header.txt") as h:
        header = h.read()

    body = ""
    with open("body.txt") as b:
        body = b.read()

    footer = ""
    with open("footer.txt") as f2:
        footer = f2.read()

    f.write(header)
    f.write("<title>SPRUCE Visualization :: " + output_filename + "</title>\n")
    f.write(body)
    f.write("d3.json('%s', function(data){\n" % (output_filename + ".json"))
    f.write(footer)
