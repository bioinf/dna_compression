#!/usr/bin/python

from compress import compress
from decompress import decompress

import colorama; colorama.init()
from colorama import Fore
from prettytable import PrettyTable

from time import time
import filecmp
import os
import sys



def testfile(filename, cond_huffman):

    print '\n' + '%' * 90

    report = {'name' : filename}
    start_time = time()
    report['before_size'] = str(os.path.getsize(filename)) + " bytes"

    # Compress
    compress(filename, [cond_huffman == 'True'])

    report['comp_time'] = '%.2f' % (time() - start_time) + 's'
    report['after_size'] = str(os.path.getsize(filename + '.z')) + " bytes"
    report['ratio'] = Fore.YELLOW + str(os.path.getsize(filename + '.z') * 100 /
                          os.path.getsize(filename)) + '%' + Fore.RESET

    start_time = time()

    print '-' * 60

    # Decompress
    decompress(filename + '.z', filename + '.new')
    report['decomp_time'] = '%.2f' % (time() - start_time) + 's'


    success = True

    path = os.path.dirname(filename) + '/'

    if filecmp.cmp(path + 'out2', path + 'out22'):
        print Fore.GREEN + 'Sequences files are identical' + Fore.RESET
    else:
        print Fore.RED + 'Sequences files are different!' + Fore.RESET
        success = False

    if filecmp.cmp(path + 'out3', path + 'out33'):
        print Fore.GREEN + 'Quality files are identical' + Fore.RESET
    else:
        print Fore.RED + 'Quality files are different!' + Fore.RESET
        success = False

    if filecmp.cmp(path + 'out1', path + 'out11'):
        print Fore.GREEN + 'Info files are identical' + Fore.RESET
    else:
        print Fore.RED + 'Info files are different!' + Fore.RESET
        success = False

    if filecmp.cmp(filename, filename + '.new'):
        print Fore.GREEN + 'Files are identical' + Fore.RESET
    else:
        print Fore.RED + 'Files are different!' + Fore.RESET
        success = False

    report['success'] = success

    return report



quick = False
if len(sys.argv) > 1:
    if sys.argv[1] == '-h':
        print 'Usage:'
        print 'test.py [-q]        -- test using test.lst table, -q -- only quick tests'
        print 'test.py file.fastq  -- test only this fastq'
        exit()
    else:
        if sys.argv[1] == '-q':
            lines = open('test.lst').readlines()
            quick = True
        else:
            lines = [sys.argv[1] + ' False True']
            quick = False
else:
    lines = open('test.lst').readlines()
    quick = False

print 'Quick testing:', quick



table = PrettyTable(["File name", "Size before", "Size after", "Ratio", "CTime", "DTime", "Result"])
for line in lines:
    if line[0] == '#': continue
    filename, quick2, cond_huffman = line.split()
    if (not quick) or quick2 == 'True': 
        report = testfile(filename, cond_huffman)
        table.add_row([report['name'], report['before_size'],
                   report['after_size'], report['ratio'], report['comp_time'],
                   report['decomp_time'], report['success']])


table.align["CTime"] = "r"
table.align["DTime"] = "r"
table.align["Ratio"] = "r"
table.align["Size before"] = "r"
table.align["Size after"] = "r"

print
print(table)
