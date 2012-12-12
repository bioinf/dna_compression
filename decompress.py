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
                    if lenc <= 16:
                        cond_table[i][prev][c] = bin(unpack('H', filein.read(2))[0])[2:].zfill(lenc)
                    else:
                        cond_table[i][prev][c] = bin(unpack('I', filein.read(4))[0])[2:].zfill(lenc)
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


def desqueeze_quality(path, filein, cond_huffman, num_reads, lng, table, cond_table):

    table2, cond_table2 = prepare_table(lng, table, cond_table)

    print "Desqueezing quality..."

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

    cachesize = calcsize('Q') * 8     
    len, tail = (cache['len'] / 8) * 8, cache['len'] % 8
    rest = bin(cache['v'])[2:].zfill(cachesize)[tail:][:len]

    return rest

def desqueeze_seq(path, filein, num_reads, lng, rest):


    print "Desqueezing sequences..."

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
        #print cachelen, cachev, cachesize, maxlen, index,
        c = table[index]
        reallen = c['len']
        c = c['c']

        cache['v'] = (cachev & ~((~(1 << reallen)) << (cachesize - reallen))) << reallen
        cache['len'] = cachelen - reallen
        
        return c


    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()


    table = {'maxlen' : 3, 't' : [{'len' : 2, 'c' : 'A'}, {'len' : 2, 'c' : 'A'}, 
                                  {'len' : 2, 'c' : 'T'}, {'len' : 2, 'c' : 'T'}, 
                                  {'len' : 2, 'c' : 'G'}, {'len' : 2, 'c' : 'G'}, 
                                  {'len' : 3, 'c' : 'C'}, {'len' : 3, 'c' : 'N'}]}
    cache = {'v' : int(rest + '0' * (64 - len(rest)), 2), 'len' : len(rest)}
    out = open(path + 'out22', 'w')
    reads = 0
    while reads < num_reads:
        pbar.update(reads)

        line = ''
        for i in range(lng):
            line += get_next(table, filein, cache)

        out.write(line + '\n')
        reads += 1

    pbar.finish()
    out.close()
    filein.close()



def desqueeze_info(path, filein, num_reads, pattern):

    print "Desqueezing info..."
    
    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()

    d = pattern['d']
    pat = pattern['pat']
    pats = pat.split('(\d*)')
    mins = pattern['mins']
    dmins = pattern['dmins']
    use_diff = pattern['use_diff']
    prev = [0] * len(d)

    count = 0
    f = open(path + 'out11', 'w'); 
    while count < num_reads:
        count += 1
        pbar.update(count)

        out = pats[0]

        for i in range(len(d)):
            c = d[i]
            num = unpack(c, filein.read(calcsize(c)))[0]
            if use_diff[i]:
                prev[i], num = num + prev[i], num + dmins[i] + prev[i]
            else:
                num += mins[i]

            out += str(num) + pats[i+1]

        f.write(out)


    pbar.finish()
    f.close()

    return


def assemble(path, filename, num_reads):

    f1 = open(path + 'out11')
    f2 = open(path + 'out22')
    f3 = open(path + 'out33')

    fileout = open(filename, 'wb')
    
    
    for i in range(num_reads):
        fileout.write(f1.readline())
        fileout.write(f2.readline())
        fileout.write('+\n')
        fileout.write(f3.readline())

    f1.close()
    f2.close()
    f3.close()

    fileout.close()
    

    return


def read_parameters(filein):

    cond_huffman = bool(unpack('B'*1, filein.read(1))[0])
    num_reads = unpack('L', filein.read(calcsize('L')))[0]
    lng = unpack('H', filein.read(2))[0]
    pat = filein.read(unpack('B', filein.read(1))[0])
    lend = unpack('B', filein.read(1))[0]
    d = ''; mins = []; dmins = []; use_diff = []
    for i in range(lend):
        d += filein.read(1) #fileout.write(pattern['d'][i])
        mins.append(unpack('I', filein.read(4))[0])
        dmins.append(unpack('i', filein.read(4))[0])
        use_diff.append(bool(unpack('B', filein.read(1))[0]))

    pattern = {'d' : d, 'pat' : pat, 'mins' : mins, 
               'dmins' : dmins, 'use_diff' : use_diff} 

    return cond_huffman, num_reads, lng, pattern


# Decompressing: from filename1 to filename2 
def decompress(filename1, filename2):

    print "Decompressing: " + filename1
    print "File size: " + str(os.path.getsize(filename1)) + " bytes"

    path, filename1 = os.path.dirname(filename1) + '/', os.path.basename(filename1)

    # Input file
    filein = open(path + filename1, 'rb')

    

    # Read parameters
    cond_huffman, num_reads, lng, pattern = read_parameters(filein)

    # Read info headers
    desqueeze_info(path, filein, num_reads, pattern)
     
    # Read tables
    table, cond_table = read_tables(filein, lng)

    # Read quality
    rest = desqueeze_quality(path, filein, cond_huffman, num_reads, lng, table, cond_table)

    # Read sequence
    desqueeze_seq(path, filein, num_reads, lng, rest)

    filein.close()

    # Assemble out[11-33]
    assemble(path, filename2, num_reads)


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "Using: " + sys.argv[0] + " file.fastq.z"
        exit()

    filename1 = sys.argv[1]
    filename2 = sys.argv[2]

    decompress(filename1, filename2)
