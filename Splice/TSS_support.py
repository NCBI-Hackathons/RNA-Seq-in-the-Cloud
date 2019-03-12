import sys
import csv
from collections import defaultdict

novel_file=sys.argv[1]
TSS_file=sys.argv[2]
tss_list=dict()
novel_list=dict()
the_file=open('test_tight_novel_%s.tsv'%novel_file, 'w')
dist=10
with open(TSS_file , 'r') as tss_f:
    for line in tss_f:
        fields = line.split('\t')
        fields2=fields[8].split(';')
        if fields[2]=='region':
            tss_list[(fields[0],fields[6],fields2[2])]=(fields[3],fields[4])
            #chr, strand, id:start,stop
print(tss_list.items())
with open(novel_file,'r') as novel_f:
    for line in novel_f:
        fields = line.split('\t')
        if fields[3].strip()=='-':
            novel_list[(fields[0],fields[3],fields[4])]=fields[2]
        if fields[3].strip()=='+':
            novel_list[(fields[0],fields[3],fields[4])]=fields[1]
#chr, strand, id:start
the_file=open('test_tight_novel_%s.tsv'%novel_file, 'a') 
wr= csv.writer(the_file,delimiter='\t',escapechar=' ', quoting=csv.QUOTE_NONE)
for key_t,loc_t in tss_list.items():
    for key_n,loc_n in novel_list.items():
        if key_t[0]==key_n[0] and key_t[1].strip()==key_n[1].strip():
            if int(loc_n)>=(int(loc_t[0])-dist) and int(loc_n)<=(int(loc_t[1])-dist):
                wr.writerow([key_n[0],key_n[1],key_n[2],key_t[1],key_t[2]])
