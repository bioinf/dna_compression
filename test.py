from compress import compress

quick = True

for filename in open('test.lst').readlines():
    if filename[0] == '#': continue

    filename, quick2, cond_huffman = filename.split()
    if (not quick) or bool(quick2): 
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        compress(filename, [bool(cond_huffman)])
