#!/usr/bin/python

from progressbar import *
from struct import unpack, calcsize

import sys
import os



def read_tables(filein):
    lng = unpack('H', filein.read(2))[0]

    table = [{} for i in range(lng)]
    for i in range(lng):
        n = unpack('B', filein.read(1))[0]
        for j in range(n):
            c = filein.read(1)
            lenc = unpack('B', filein.read(1))[0]
            table[i][c] = bin(unpack('H', filein.read(2))[0])[2:].zfill(lenc)

    cond_table = [{} for i in range(lng-1)]
    for i in range(lng - 1):
        n = unpack('B', filein.read(1))[0]
        for j in range(n):
            prev = filein.read(1)
            cond_table[i][prev] = {}
            m = unpack('B', filein.read(1))[0]
            for k in range(m):
                c = filein.read(1)
                lenc = unpack('B', filein.read(1))[0]
                if lenc > 0:
                    cond_table[i][prev][c] = bin(unpack('H', filein.read(2))[0])[2:].zfill(lenc)
                else:
                    cond_table[i][prev][c] = ''
                
    return lng, table, cond_table


def desqueeze_quality(path, filename):

    filein = open(path + filename, 'rb')

    lng, table, cond_table = read_tables(filein)
    
    print lng, table, cond_table

    error

    print "Squeezing..."

    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()
    
    summ = 0; count = 0
    cache = ''; MaxNcache = 8 * struct.calcsize('L')
    f = open(path + 'out_3', 'r')
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
            filein.write(struct.pack('L', int(cache[:MaxNcache], 2)))
            cache = cache[MaxNcache:]

    pbar.finish()
    f.close()
    filein.close()

    return summ


def decompress(filename):

    path, filee = os.path.dirname(filename), os.path.basename(filename)
    path += '/'


    desqueeze_quality(path, filee)



if __name__ == '__main__':

    ###################### Parameters ###########################
    cond_huffman = True
    
    parameters = [cond_huffman]
    #############################################################

    if len(sys.argv) < 2:
        print "Using: " + sys.argv[0] + " file.fstq"
        exit()

    filename = sys.argv[1]


    decompress(filename)
