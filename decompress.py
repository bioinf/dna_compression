#!/usr/bin/python

from progressbar import *
from struct import unpack, calcsize

import sys
import os



def read_tables(filein, lng):

    table = [{} for i in range(lng)]
    for i in range(lng):
        n = unpack('B', filein.read(1))[0]
        for j in range(n):
            c = filein.read(1)
            lenc = unpack('B', filein.read(1))[0]
            if lenc > 0:
                table[i][c] = bin(unpack('H', filein.read(2))[0])[2:].zfill(lenc)
            else:
                table[i][c] = ''

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
                
    return table, cond_table

def prepare_table(lng, table, cond_table):
    
    print 'Preparing auxiliary tables...'
    
    table2 = [{'maxlen' : 0, 't' : []} for t in table]

    for i in range(len(table2)):
        maxlen = max([len(tt) for tt in table[i].values()])
        table2[i]['maxlen'] = maxlen
        table2[i]['t'] = [{'c': '?', 'len' : 0} for j in range(2**maxlen)]
        for c in table[i]:
            code = table[i][c]
            lng = len(code)
            diff = 2**(maxlen - lng)
            code = int('0' + code, 2) * diff
            table2[i]['t'][code : code + diff] = [{'c' : c, 'len' : lng}] * diff


    cond_table2 = [{c : {'maxlen' : 0, 't' : []} for c in t} for t in cond_table]

    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = len(cond_table2)).start()

    for i in range(len(cond_table2)):
        pbar.update(i)
        for cc in cond_table2[i]:
            maxlen = max([len(tt) for tt in cond_table[i][cc].values()])
            cond_table2[i][cc]['maxlen'] = maxlen
            cond_table2[i][cc]['t'] = [{'c': '?', 'len' : 0} for j in range(2**maxlen)]
            for c in cond_table[i][cc]:
                code = cond_table[i][cc][c]
                lng = len(code)
                diff = 2**(maxlen - lng)
                code = int('0' + code, 2) * diff
                cond_table2[i][cc]['t'][code:code+diff] = [{'c' : c, 'len' : lng}] * diff
        
    pbar.finish()

    return table2, cond_table2


def desqueeze_quality(path, filein, cond_huffman, num_reads, lng):

    table, cond_table = read_tables(filein, lng)
    table2, cond_table2 = prepare_table(lng, table, cond_table)

    print "Desqueezing..."

    def get_next(table, filein, cache):
        cachesize = calcsize('Q') * 8  # Cache size (64 bits)
        cachelen = cache['len']        # Substantial bits of cache
        cachev = cache['v']            # Cache value
        bytesize = calcsize('B') * 8   # Byte size in bits

        maxlen = table['maxlen']
        table = table['t']


        if cachelen < maxlen:
            cachev = cachev >> (cachesize - cachelen)
            while cachesize - cachelen >= bytesize:
                byte = unpack('B', filein.read(1))[0]
                cachev = (cachev << bytesize) | byte
                cachelen += bytesize
            cachev = cachev << (cachesize - cachelen)

        index = cachev >> (cachesize - maxlen)
        c = table[index]
        reallen = c['len']
        c = c['c']

        cache['v'] = (cachev & ~((~(1 << reallen)) << (cachesize - reallen))) << reallen
        cache['len'] = cachelen - reallen
        
        return c


    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()
    
    cache = {'v' : 0, 'len' : 0}
    out = open(path + 'out33', 'w')
    reads = 0
    while reads < num_reads:
        pbar.update(reads)

        cur = get_next(table2[0], filein, cache)
        line = cur
        for i in range(lng - 1):
            if cond_huffman:
                cur = get_next(cond_table2[i][cur], filein, cache)
            else:
                cur = get_next(table2[i + 1], filein, cache)

            line += cur
        
        out.write(line + '\n')
        reads += 1

    pbar.finish()
    out.close()
    filein.close()


def decompress(filename):

    print "Decompressing: " + filename
    print "File size: " + str(os.path.getsize(filename)) + " bytes"

    path, filename = os.path.dirname(filename) + '/', os.path.basename(filename)

    # Input file
    filein = open(path + filename, 'rb')

    # Read parameters
    cond_huffman = bool(unpack('B'*1, filein.read(1))[0])
    num_reads = unpack('L', filein.read(calcsize('L')))[0]
    lng = unpack('H', filein.read(2))[0]

    # Read quality
    desqueeze_quality(path, filein, cond_huffman, num_reads, lng)



if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "Using: " + sys.argv[0] + " file.fastq.z"
        exit()

    filename = sys.argv[1]


    decompress(filename)
