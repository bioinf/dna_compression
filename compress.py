#!/usr/bin/python

from Huffman3 import *
from progressbar import *

import sys
import os
from struct import pack, calcsize

def analyze(filename, path):
    print "Generating out1, out2 and out3..."
    f = open(filename, 'r')
    out1 = open(path + 'out1', 'w')
    out2 = open(path + 'out2', 'w')
    out3 = open(path + 'out3', 'w')
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

    out3 = open(path + 'out3', 'r')
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

    return num_reads, lng, alph_card, table, cond_table



def write_tables(out, lng, table, cond_table):
    out.write(pack('H', lng))

    for i in range(lng):
        out.write(pack('B', len(table[i])))
        for c in table[i]:
            out.write(c)
            out.write(pack('B', len(table[i][c])))
            if len(table[i][c]) > 0:
                out.write(pack('H', int('0' + table[i][c], 2)))

    for i in range(lng - 1):
        out.write(pack('B', len(cond_table[i])))
        for prev in cond_table[i]:
            out.write(prev)
            out.write(pack('B', len(cond_table[i][prev])))
            for c in cond_table[i][prev]:
                out.write(c)
                out.write(pack('B', len(cond_table[i][prev][c])))
                if len(cond_table[i][prev][c]) > 0:
                    out.write(pack('H', int('0' + cond_table[i][prev][c], 2)))



def squeeze_quality(path, filename, num_reads, lng, table, cond_table, cond_huffman):

    out = open(path + filename, 'wb')
    out.write(pack('L', num_reads))
    
    # Write tables to file
    write_tables(out, lng, table, cond_table)

    # Squeezing
    print "Squeezing..."

    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()
    
    summ = 0; count = 0
    cache = ''; MaxNcache = 8 * calcsize('>Q')
    f = open(path + 'out3', 'r')
    line = f.readline()
    while line:
        count += 1
        pbar.update(count)

        cur = line[0]
        summ += len(table[0][cur])
        cache += table[0][cur]
        for i in range(lng - 1):
            prev = cur
            cur = line[i+1]
            if cond_huffman:
                summ += len(cond_table[i][prev][cur])
                cache += cond_table[i][prev][cur]
            else:
                summ += len(table[i + 1][cur])
                cache += table[i + 1][cur]

        line = f.readline()
        

        while len(cache) > MaxNcache:
            out.write(pack('>Q', int(cache[:MaxNcache], 2)))
            cache = cache[MaxNcache:]

    pbar.finish()
    f.close()
    
    # Cache tail processing
    bytesize = calcsize('B') * 8 # Byte size in bits
    while cache:
        if len(cache) < bytesize:
            cache += '0' * (bytesize - len(cache))
        out.write(pack('B', int(cache[:bytesize], 2)))
        cache = cache[bytesize:]
    out.write('\x00' * 8)
    out.close()

    return summ


def compress(filename, parameters):
    cond_huffman = parameters[0]

    print "Compressing: " + filename
    print "File size: " + str(os.path.getsize(filename)) + " bytes"

    path, _ = os.path.dirname(filename), os.path.basename(filename)
    path += '/'

    num_reads, lng, alph_card, table, cond_table = analyze(filename, path)
    qual_summ = squeeze_quality(path, filename + '.z', num_reads, lng, table, cond_table, cond_huffman)

    quality_bytes = qual_summ / 8


    header_bytes = 2*(lng * alph_card**2 if cond_huffman else lng * alph_card)
    nucl_bytes = (num_reads * lng) / 4
    info_bytes = 0
    print "Squeezed to: " + str(quality_bytes + 
                                  nucl_bytes + 
                                  header_bytes +
                                  info_bytes) + " bytes"
    print "Quality: " + str(quality_bytes) + " bytes"
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
