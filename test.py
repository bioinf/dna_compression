#!/usr/bin/python

from compress import compress
import sys

quick = False
if len(sys.argv) > 1:
    quick = (sys.argv[1].find('-q') != -1)


print 'Quick testing:', quick

for filename in open('test.lst').readlines():
    if filename[0] == '#': continue

    filename, quick2, cond_huffman = filename.split()

    if (not quick) or quick2 == 'True': 
        print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        compress(filename, [cond_huffman == 'True'])
