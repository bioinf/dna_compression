#!/usr/bin/python

from progressbar import *

import sys
import os
import struct



def read_tables(filein):
    lng = unpack('H', filein.read(2))[0]

    table = [{} for i in range(lng)]
    for i in range(lng):
        n = unpack('B', filein.read(1))
        for j in range(n):
            c = filein.read(1)
            lenc = unpack('B', filein.read(1))
            table[i][c] = bin(unpack('H', filein.read(2)))[2:].zfill(lenc)

    cond_table = [{} for i in range(lng-1)]
    for i in range(lng - 1):
        n = unpack('B', filein.read(1))
        for j in range(n):
            prev = filein.read(1)
            cond_table[i][prev] = {}
            m = unpack('B', filein.read(1))
            for k in range(m):
                c = filein.read(1)
                lenc = unpack('B', filein.read(1))
                if lenc > 0:
                    cond_table[i][prev][c] = bin(unpack('H', filein.read(2)))[2:].zfill(lenc)
                else:
                    cond_table[i][prev][c] = ''
                
    return lng, table, cond_table


def desqueeze_quality(path, filename, num_reads, lng, table, cond_table, cond_huffman):
    print "Squeezing..."

    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()
    
    summ = 0; count = 0
    cache = ''; MaxNcache = 8 * struct.calcsize('L')
    f = open(path + 'out_3', 'r')
    out = open(path + filename, 'wb')
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
            out.write(struct.pack('L', int(cache[:MaxNcache], 2)))
            cache = cache[MaxNcache:]

    pbar.finish()
    f.close()
    out.close()

    return summ


def decompress(filename, parameters):
    cond_huffman = parameters[0]

    path, file = os.path.dirname(filename), os.path.basename(filename)
    path += '/'

    num_reads, lng, alph_card, table, cond_table = analyze(filename, path)
    qual_summ = squeeze_quality(path, 'quality', num_reads, lng, table, cond_table, cond_huffman)

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
