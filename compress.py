#!/usr/bin/python

from Huffman3 import *
from progressbar import *

import sys
import os


def analyze(filename, path):
    print "Processing: " + filename
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
    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()
    count = 0

    freq = [{} for i in range(lng)]
    cond_freq = [{} for i in range(lng)]

    out3 = open(path + 'out_3', 'r')
    line = out3.readline()
    while line:
        count += 1

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
        pbar.update(count)

    pbar.finish()
    out3.close()
    

    for c1 in cond_freq[i]:
        for c2 in cond_freq[i][c1]:
            cond_freq[i][c1][c2] = cond_freq[i][c1][c2] / float(freq[i][c1])
            
            
    for i in range(len(freq)):
        summ = float(sum(freq[i].values()))
        for sym in freq[i]:
            freq[i][sym] /= summ


    print "Building trees..."
    alph_card = 0
    symbols = [0 for i in range(lng)]
    for i in range(lng):
        pairs = [ (sym, freq[i][sym]) for sym in freq[i] ]
        alph_card = max(alph_card, len(pairs))
        symbols[i] = makenodes(pairs)
        iterate(symbols[i])

    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = lng - 2).start()
    cond_symbols = []
    for i in range(lng - 1):
        pbar.update(i)
        cond_symbols.append({})
        for c in cond_freq[i]:
            cf = cond_freq[i][c]
            pairs = [ (sym, cf[sym]) for sym in cf ]
            cond_symbols[i][c] = makenodes(pairs)
            iterate(cond_symbols[i][c])
    pbar.finish()

    print "Building tables..."
    table = [{} for i in range(lng)]
    for i in range(lng):
        for c in freq[i]:
            table[i][c] = encode(c, symbols[i])

    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = lng - 2).start()
    cond_table = [{} for i in range(lng - 1)]
    for i in range(lng - 1):
        pbar.update(i)
        for c1 in cond_freq[i]:
            cond_table[i][c1] = {}
            for c2 in cond_freq[i][c1]:
                cond_table[i][c1][c2] = encode(c2, cond_symbols[i][c1])
    pbar.finish()

    return num_reads, lng, alph_card, symbols, cond_symbols



def squeeze(path, num_reads, lng, symbols, cond_symbols, cond_huffman):
    print "Squeezing..."

    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()

    summ = 0; count = 0
    f = open(path + 'out_3', 'r')
    line = f.readline()
    while line:
        count += 1
        pbar.update(count)

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

    pbar.finish()
    f.close()
    
    return summ


def compress(filename, parameters):
    cond_huffman = parameters[0]

    path, file = os.path.dirname(filename), os.path.basename(filename)
    path += '/'

    num_reads, lng, alph_card, symbols, cond_symbols = analyze(filename, path)
    summ = squeeze(path, num_reads, lng, symbols, cond_symbols, cond_huffman)

    quality_bytes = summ / 8
    header_bytes = 2*(lng * alph_card**2 if cond_huffman else lng * alph_card)
    nucl_bytes = (num_reads * lng) / 4
    info_bytes = 0
    print "Squeezed to: " + str(quality_bytes + 
                                  nucl_bytes + 
                                  header_bytes +
                                  info_bytes) + " bytes"

    print "Caution: reads info lost!"



if __name__ == '__main__':

    ###################### Parameters ###########################
    cond_huffman = True
    
    parameters = [cond_huffman]
    #############################################################

    if len(sys.argv) < 2:
        print "Using: " + sys.argv[0] + " file.fastq"
        exit()

    filename = sys.argv[1]


    compress(filename, parameters)
