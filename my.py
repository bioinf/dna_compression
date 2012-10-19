from Huffman3 import *
import sys
import os

if len(sys.argv) < 2:
    print "Using: " + sys.argv[0] + " file.fastq"
    exit()


filename = sys.argv[1]
path, file = os.path.dirname(filename), os.path.basename(filename)
path += '/'
print "Processing: " + path + file


print "Generating out and out_2"
f = open(filename, 'r')
out = open(path + 'out', 'w')
out2 = open(path + 'out_2', 'w')
line = f.readline()
line = f.readline()
line = f.readline()
num_reads = 0
while line:
    line = f.readline()
    n = len(line.replace('\n', ''))
    out.write(line)
    out2.write(line.replace('\n', ''))
    line = f.readline()
    line = f.readline()
    line = f.readline()
    num_reads += 1


f.close()
out.close()
out2.close()


print "Generating out[0.." + str(n - 1) + ']'
f = open(path + 'out', 'r')
out = []
for i in range(n):
    out.append(open(path + 'out' + str(i), 'w'))

for line in f.readlines():
    for i in range(n):
        out[i].write(line[i])
        
for ff in out:
    ff.close()

f.close()



print "Building trees..."

lng = n

freq = [{} for i in range(lng)]
cond_freq = [{} for i in range(lng)]


f1 = open(path + 'out0', 'r')
c1 = f1.read(1)
while c1:
    if not freq[0].has_key(c1):
        freq[0][c1] = 0
    freq[0][c1] += 1
    c1 = f1.read(1)

f1.close()

f1 = open(path + 'out' + str(lng - 1), 'r')
c1 = f1.read(1)
while c1:
    if not freq[lng - 1].has_key(c1):
        freq[lng - 1][c1] = 0
    freq[lng - 1][c1] += 1
    c1 = f1.read(1)

f1.close()
    
for i in range(lng - 1):
    print i
    f1 = open(path + 'out' + str(i), 'r')
    f2 = open(path + 'out' + str(i + 1), 'r')

    c1 = f1.read(1)
    c2 = f2.read(1)
    while c1:
        if not cond_freq[i].has_key(c1):
            cond_freq[i][c1] = {}
        if not cond_freq[i][c1].has_key(c2):
            cond_freq[i][c1][c2] = 0
        if not freq[i].has_key(c1):
            freq[i][c1] = 0
        
        freq[i][c1] += 1    
        cond_freq[i][c1][c2] += 1
        c1 = f1.read(1)
        c2 = f2.read(1)

    for c1 in cond_freq[i]:
        for c2 in cond_freq[i][c1]:
            cond_freq[i][c1][c2] = cond_freq[i][c1][c2] / float(freq[i][c1])
    
    f1.close()
    f2.close()
    
for i in range(len(freq)):
    summ = float(sum(freq[i].values()))
    for sym in freq[i]:
        freq[i][sym] /= summ

symbols = [0 for i in range(lng)]
for i in range(lng):
    pairs = [ (sym, freq[i][sym]) for sym in freq[i] ]
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
               

f = open(path + 'out_2', 'r')

summ = 0
count = 0
cur = f.read(1)
while cur:
    if count % 1000 == 0:
        print count/float(num_reads)
    count += 1
    
    summ += len(encode(cur, symbols[0]))
    for i in range(lng - 1):
        prev = cur
        cur = f.read(1)
        summ += len(encode(cur, cond_symbols[i][prev]))
        #summ += len(encode(cur, symbols[i + 1]))
        
    cur = f.read(1)

f.close()
print "Squeezed to: " + str(summ + num_reads * n / 4) + " bytes"