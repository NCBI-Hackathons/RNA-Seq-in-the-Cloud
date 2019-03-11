import sys
import csv

novel_file=sys.argv[1]
TSS_file=sys.argv[2]

tss_list=list()

#mix_list=list()
the_file=open('novel_tss.tsv', 'w')
with open(TSS_file , 'r') as tss_f:
    for line in tss_f:
        fields = line.split('\t')
        tss_list.append([fields[0],fields[2],fields[3],fields[4],fields[8]])

with open(novel_file,'r') as novel_f:
    for line in novel_f:
        fields = line.split('\t')
        for fields1 in tss_list:
            if fields1[1]=='transcription_start_site':
                if fields[0] == fields1[0]:
                    with open('novel_tss.tsv', 'a') as the_file:
                        wr= csv.writer(the_file,delimiter='\t',escapechar=' ', quoting=csv.QUOTE_NONE)
                        if fields[1]>=fields1[2]:
                            if fields[2]<=fields1[3]:
#                                mix_list.append([fields[0],fields[1],fields[2],fields1[2],fields1[3],fields1[4],'TSS inside'])
                                wr.writerow([fields[0],fields[1],fields[2],fields1[2],fields1[3],fields1[4],'TSS inside'])
                        elif fields[2]>=fields1[2]:
#                            mix_list.append([fields[0],fields[1],fields[2],fields1[2],fields1[3],fields1[4],'before the TSS'])
                            wr.writerow([fields[0],fields[1],fields[2],fields1[2],fields1[3],fields1[4],'before the TSS'])
                        elif fields[1]<=fields1[3]:
#                            mix_list.append([fields[0],fields[1],fields[2],fields1[2],fields1[3],fields1[4],'after the TSS'])
                            wr.writerow([fields[0],fields[1],fields[2],fields1[2],fields1[3],fields1[4],'after the TSS'])

                                

