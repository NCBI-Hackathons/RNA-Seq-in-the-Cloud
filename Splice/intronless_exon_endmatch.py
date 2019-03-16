import sys
import csv
from collections import defaultdict

gtffile=sys.argv[1]
f=open(gtffile , 'r')
ref_unfiltered_exon_dict=defaultdict(list)
ref_filtered_exon_dict=defaultdict(list)
final_exon_dict=defaultdict(list)

i=1
fh=f.readlines()[2:(-1)]       

for i in range (2,(len(fh)-2)):
    with open('exon_novel_%s.tsv'%gtffile, 'a') as the_file:
        fields = fh[i].split('\t')
        if fields[2] == 'transcript':
            ref1=fields[8].split('.')[1] #gene_id "STRG.48"; transcript_id "STRG.48
            ref2=fields[8].split(';')[2] #reference_id "XM_005244719.4" -- unique
            fullname=fields[8].split(';')[1] #transcript_id "STRG.48.6"       
            j=i+1
            if fh[i+2].split('\t')[2]=='exon':
                if fh[j].split('\t')[6]=='+':
                    if (fh[j+1].split('\t')[2] == 'exon' and fh[j-1].split('\t')[2] == 'transcript' and fh[j].split('\t')[2] == 'exon'):
                        fieldsx=fh[j].split('\t')
                        new_key=(fields[0],fields[6],ref1)
                        ref_filtered_exon_dict[new_key].append((fieldsx[8],int(fieldsx[3]),int(fieldsx[4])))
     
                if fh[j].split('\t')[6]=='-':
                    while fh[j].split('\t')[2] == 'exon':
                        j=j+1
                    if (fh[j-2].split('\t')[2] == 'exon' and fh[j].split('\t')[2] == 'transcript' and fh[j-1].split('\t')[2] == 'exon'):
                        fieldsx=fh[j].split('\t')
                        new_key=(fields[0],fields[6],ref1)
                        ref_filtered_exon_dict[new_key].append((fieldsx[8],int(fieldsx[3]),int(fieldsx[4])))
                        
referenced_exon_dict=defaultdict(list)                  
not_referenced_exon_dict=defaultdict(list)
common_exon_dict=defaultdict(list)
for k,v in ref_filtered_exon_dict.items():
    for vv in v:
        if "reference_id" in vv[0]:
            referenced_exon_dict[k].append(v)
        else:
            not_referenced_exon_dict[k].append(v)
            
#collapes all referenced exon ends into common_exon_dict
for k,v in referenced_exon_dict.items():
    if k not in not_referenced_exon_dict.keys():
        final_exon_dict[k].append(not_referenced_exon_dict[k])
    if k in not_referenced_exon_dict.keys():
        for vv in referenced_exon_dict[k]:
            if k[1]=='+':
                vv_end=str(vv).split(',')
                common_exon_dict[k].append(vv_end[2])
            elif k[1]=='-':
                vv_end=str(vv).split(',')
                common_exon_dict[k].append(vv_end[1])
                
            
for k, v in not_referenced_exon_dict.items():
    if k[1]=='+':
        for vv in v:
            vv_end=str(vv).split(',')
            if vv_end[2] not in common_exon_dict[k]:
                final_exon_dict[k].append(vv)
    if k[1]=='-':
        for vv in v:
            vv_end=str(vv).split(',')
            if vv_end[1] not in common_exon_dict[k]:
                final_exon_dict[k].append(vv)    
            

the_file=open('exon_novel_endmatch_%s.tsv'%gtffile, 'w')
with open('exon_novel_endmatch_%s.tsv'%gtffile, 'a') as the_file:
    wr= csv.writer(the_file,delimiter='\t',escapechar=' ', quoting=csv.QUOTE_NONE)
    for k,v in final_exon_dict.items():
        for vv in v:
            if len(vv) >0:
                vv_end=str(vv).split(',')
                wr.writerow([k[0],str(vv_end[1]).split(',')[0],str(vv_end[2]).split(',')[0].replace(')]', ''),k[1],str(vv_end).split(';')[1].strip()])

                    
            

         
























        
    
            
