#!/usr/bin/python
from Huffman3 import *

import sys
import os

###################### Parameters ###########################
cond_huffman = True

#############################################################

if len(sys.argv) < 2:
    print "Using: " + sys.argv[0] + " file.fastq"
    exit()

filename = sys.argv[1]
path, file = os.path.dirname(filename), os.path.basename(filename)
path += '/'


def analyze(filename, path):
    print "Processing: " + path + file
    print "File size: " + str(os.path.getsize(filename)) + " bytes"

    print "Generating out_1, out_2 and out_3"
    f = open(filename, 'r')
    out1 = open(path + 'out_1', 'w')
    out2 = open(path + 'out_2', 'w')
    out3 = open(path + 'out_3', 'w')
    line = f.readline()
    num_reads = 0
    while line:
        out1.write(line)
        line = f.readline()
        out2.write(line)
        line = f.readline()
        line = f.readline()
        out3.write(line)
        lng = len(line.replace('\n', ''))
        num_reads += 1
        line = f.readline()

    f.close()
    out1.close()
    out2.close()
    out3.close()



    print "Counting frequences..."
    print "Positions: ", 

    freq = [{} for i in range(lng)]
    cond_freq = [{} for i in range(lng)]

    out3 = open(path + 'out_3', 'r')
    line = out3.readline()
    while line:
        c = line[lng - 1]
        if not freq[lng - 1].has_key(c): freq[lng - 1][c] = 0
        freq[lng - 1][c] += 1

        c2 = line[0]
        for i in range(lng - 1):
            c1 = c2
            c2 = line[i+1]
            if not cond_freq[i].has_key(c1):
                cond_freq[i][c1] = {}
            if not cond_freq[i][c1].has_key(c2):
                cond_freq[i][c1][c2] = 0
            if not freq[i].has_key(c1):
                freq[i][c1] = 0
                
            freq[i][c1] += 1    
            cond_freq[i][c1][c2] += 1

        line = out3.readline()            
    out3.close()
    

    for c1 in cond_freq[i]:
        for c2 in cond_freq[i][c1]:
            cond_freq[i][c1][c2] = cond_freq[i][c1][c2] / float(freq[i][c1])
            
            
    for i in range(len(freq)):
        summ = float(sum(freq[i].values()))
        for sym in freq[i]:
            freq[i][sym] /= summ


    print "\nBuilding trees..."
    alph_card = 0
    symbols = [0 for i in range(lng)]
    for i in range(lng):
        pairs = [ (sym, freq[i][sym]) for sym in freq[i] ]
        alph_card = max(alph_card, len(pairs))
        symbols[i] = makenodes(pairs)
        iterate(symbols[i])

    cond_symbols = []
    for i in range(lng - 1):
        cond_symbols.append({})
        for c in cond_freq[i]:
            cf = cond_freq[i][c]
            pairs = [ (sym, cf[sym]) for sym in cf ]
            cond_symbols[i][c] = makenodes(pairs)
            iterate(cond_symbols[i][c])
    
    return num_reads, lng, alph_card, symbols, cond_symbols



def squeeze(path, num_reads, symbols, cond_symbols):
    print "\nSqueezing..."

    summ = 0; count = 0
    count_step = int(num_reads / 100.0) + 1
    f = open(path + 'out_3', 'r')
    line = f.readline()
    while line:
        if count % count_step == 0:
            print str(count / count_step) + "% ",
        count += 1

        cur = line[0]
        summ += len(encode(cur, symbols[0]))
        for i in range(lng - 1):
            prev = cur
            cur = line[i+1]
            if cond_huffman:
                summ += len(encode(cur, cond_symbols[i][prev]))
            else:
                summ += len(encode(cur, symbols[i + 1]))
            
        line = f.readline()

    f.close()
    
    return summ


num_reads, lng, alph_card, symbols, cond_symbols = analyze(filename, path)
summ = squeeze(path, num_reads, symbols, cond_symbols)



quality_bytes = summ / 8
header_bytes = 2*(lng * alph_card**2 if cond_huffman else lng * alph_card)
nucl_bytes = (num_reads * lng) / 4
info_bytes = 0    ## TBD
print "\nSqueezed to: " + str(quality_bytes + 
                              nucl_bytes + 
                              header_bytes +
                              info_bytes) + " bytes"

print "Caution: reads info lost!"
