#!/usr/bin/python

from Huffman3 import *
from progressbar import *

import sys
import os
import re
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


    # TBD: analyzing header pattern 
    print "Pattern extraction..."

    def get_min_common_pattern(pat, line):
        ipat = 0; iline = 0
        lngpat = len(pat); lngline = len(line)
        while ipat < lngpat and iline < lngline:
            if pat[ipat] == line[iline]:
                ipat += 1; iline += 1
                continue

            if ipat + 5 <= lngpat and pat[ipat:ipat+5] == '(\d*)':
                while iline < lngline and line[iline].isdigit(): 
                    iline += 1
                ipat += 5
                continue
            
            if pat[ipat].isdigit() and line[iline].isdigit():
                ipat2 = ipat; iline2 = iline
                while ipat2 < lngpat and pat[ipat2].isdigit(): 
                    ipat2 += 1
                while iline2 < lngline and line[iline2].isdigit(): 
                    iline2 += 1
                
                pat = pat[:ipat] + '(\d*)' + pat[ipat2:]
                ipat += 5
                iline = iline2
                lngpat = len(pat)
                continue

            else: error

        return pat


    out1 = open(path + 'out1')
    line = out1.readline()
    pat = line
    while line:
        patre = re.compile(pat)
        if not patre.match(line):
            pat = get_min_common_pattern(pat, line)
        line = out1.readline()

    out1.close()

    pat = "@ERR001268.(\d*) 080821_HWI-EAS301_0002_30ALBAAXX:1:(\d*):(\d*):(\d*)/(\d*)";
    pattern = {
        're' : re.compile(pat),
        'd' : [],
        }



    return num_reads, lng, alph_card, table, cond_table, pattern



def write_tables(out, lng, table, cond_table):

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



def squeeze_quality(path, fileout, num_reads, lng, table, cond_table, cond_huffman):


    print "Squeezing quality..."

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
            fileout.write(pack('>Q', int(cache[:MaxNcache], 2)))
            cache = cache[MaxNcache:]

    pbar.finish()
    f.close()
    
    # Cache tail processing
    bytesize = calcsize('B') * 8 # Byte size in bits
    while cache:
        if len(cache) < bytesize:
            cache += '0' * (bytesize - len(cache))
        fileout.write(pack('B', int(cache[:bytesize], 2)))
        cache = cache[bytesize:]

    return summ / 8 + 1


def squeeze_seq(path, fileout, num_reads, lng):

    print "Squeezing sequences..."


    def seq_to_bits(s):
        return s.replace('A', '00').replace('T', '01').replace('G', '10').replace('C', '110').replace('N', '111')


    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()

    summ = 0; count = 0
    cache = ''; MaxNcache = calcsize('>Q') * 8 # Max number of bits in cache
    f = open(path + 'out2', 'r'); line = f.readline()
    while line:
        count += 1
        pbar.update(count)

        summ += len(line.replace('\n', '').replace('\r', ''))
        cache += seq_to_bits(line.replace('\n', '').replace('\r', ''))

        while len(cache) > MaxNcache:
            fileout.write(pack('>Q', int(cache[:MaxNcache], 2)))
            cache = cache[MaxNcache:]

        line = f.readline()

    pbar.finish()
    f.close()

    # Cache tail processing
    summ += len(cache) 
    cache = seq_to_bits(cache)
    bytesize = calcsize('B') * 8 # Byte size in bits
    while cache:
        if len(cache) < bytesize:
            cache += '0' * (bytesize - len(cache))
        fileout.write(pack('B', int(cache[:bytesize], 2)))
        cache = cache[bytesize:]


    return summ / 4 + 1


def squeeze_info(path, fileout, num_reads, pattern):

    print "Squeezing info..."
    pat = pattern['re']

    widgets = [Bar('#'), ' ', ETA()]
    pbar = ProgressBar(widgets = widgets, maxval = num_reads).start()

    info_bytes = 0
    count = 0
    #cache = ''; MaxNcache = calcsize('>Q') * 8 # Max number of bits in cache
    f = open(path + 'out1', 'r'); line = f.readline()
    while line:
        count += 1
        pbar.update(count)

        m = pat.match(line).groups()

        fileout.write(pack('I', int(m[0])))
        fileout.write(pack('B', int(m[1])))
        fileout.write(pack('H', int(m[2])))
        fileout.write(pack('H', int(m[3])))
        fileout.write(pack('B', int(m[4])))

        info_bytes += 4
#        cache += seq_to_bits(line.replace('\n', '').replace('\r', ''))

#        while len(cache) > MaxNcache:
#            fileout.write(pack('>Q', int(cache[:MaxNcache], 2)))
#            cache = cache[MaxNcache:]

        line = f.readline()

    pbar.finish()
    f.close()

    # Cache tail processing
#    summ += len(cache) 

    return info_bytes


def compress(filename, parameters):

    print "Compressing: " + filename
    filesize = os.path.getsize(filename)
    print "File size: " + str(filesize) + " bytes"

    path, _ = os.path.dirname(filename) + '/', os.path.basename(filename)

    # Analyzing and splitting out to three files: out1, out2, out3
    num_reads, lng, alph_card, table, cond_table, pattern = analyze(filename, path)

    

    # Output file
    fileout = open(filename + '.z', 'wb')

    # Write parameters to file
    cond_huffman = parameters[0]
    fileout.write(pack('B', cond_huffman))
    fileout.write(pack('L', num_reads))
    fileout.write(pack('H', lng))


    # Write headers
    info_bytes = squeeze_info(path, fileout, num_reads, pattern)

    # Write Huffman tables to file
    write_tables(fileout, lng, table, cond_table)

    # Write quality
    qual_bytes = squeeze_quality(path, fileout, num_reads,
                                lng, table, cond_table, cond_huffman)

    # Write sequence
    seq_bytes = squeeze_seq(path, fileout, num_reads, lng)


 
    fileout.write('0' * 8)
    fileout.close()

    
    # Delete temporary out[1-3] files


    # Print results
    tabl_bytes = 2*(lng * alph_card**2 if cond_huffman else lng * alph_card)
    nucl_bytes = (num_reads * lng) / 4
    print "Squeezed from " + str(filesize) + " bytes to " + str(qual_bytes +
                                nucl_bytes +
                                tabl_bytes +
                                info_bytes) + " bytes"
    print "Quality: "   + str(qual_bytes) + " bytes"
    print "Sequences: " + str(nucl_bytes) + " bytes"
    print "Tables: " + str(tabl_bytes)  + " bytes"
    print "Headers info: " + str(info_bytes) + " bytes"
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
