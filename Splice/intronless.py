import sys
import csv
#import pysam
#from collections import OrderedDict
#import pandas as pd

gtffile=sys.argv[1]
#fh = open(gtffile, 'r').readlines()

the_file=open('novel.tsv', 'w')

mixlist=list()
with open(gtffile , 'r') as fh:
    next(fh)
    next(fh)
    (fh)
    for line in fh:
        fields = line.split('\t')

        line2 = fh.next()
        line3 = fh.next()
        fields3 = line3.split('\t')
        if fields[2]=='transcript' and fields3[2]=='transcript':
            
            cub_this = [fields[0],fields[3],fields[4],fields[6]]
            cub_next = [fields3[0],fields3[3],fields3[4],fields3[6]]
            if cub_this not in mixlist:
                mixlist.append(cub_this)
                
            if cub_next not in mixlist:
                mixlist.append(cub_next)
print(mixlist)                
with open('novel.tsv', 'a') as the_file:
#                print(type(fields3))
    wr= csv.writer(the_file,delimiter='\t',escapechar=' ', quoting=csv.QUOTE_NONE)                 
    for line in mixlist:
        wr.writerow(line)
