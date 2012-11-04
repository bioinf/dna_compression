#!/usr/bin/python

from compress import compress
from decompress import decompress

import colorama; colorama.init()
from colorama import Fore

import filecmp
import os
import sys

quick = False
if len(sys.argv) > 1:
    quick = (sys.argv[1].find('-q') != -1)


print 'Quick testing:', quick

for filename in open('test.lst').readlines():
    if filename[0] == '#': continue

    filename, quick2, cond_huffman = filename.split()

    if (not quick) or quick2 == 'True': 
        print '\n' + '%' * 90
        compress(filename, [cond_huffman == 'True'])
        print '-' * 60
        decompress(filename + '.z')
        path = os.path.dirname(filename) + '/'
        if filecmp.cmp(path + 'out3', path + 'out33'):
            print Fore.GREEN + 'Quality files are identical' + Fore.RESET
        else:
            print Fore.RED + 'Quality files are different!' + Fore.RESET

            
