import sys
import csv
#import pysam
#from collections import OrderedDict
#import pandas as pd

gtffile=sys.argv[1]
f=open(gtffile , 'r')

the_file=open('novel_%s.tsv'%gtffile, 'w')

fh=f.readlines()[2:(-1)]
i=1
print(len(fh))
for i in range (2,(len(fh)-2)):
    with open('novel_%s.tsv'%gtffile, 'a') as the_file:
        wr= csv.writer(the_file,delimiter='\t',escapechar=' ', quoting=csv.QUOTE_NONE)
        fields = fh[i].split('\t')
        ref=fields[8].split(';')[3].strip()
        if not ref.startswith("ref"):
            fields3 = fh[i+2].split('\t')
            if fields[2] == 'transcript' and fields3[2] == 'exon':
                line = [fields[0],fields[3],fields[4],fields[6],fields[8].split(';')[1]]
                wr.writerow(line)

