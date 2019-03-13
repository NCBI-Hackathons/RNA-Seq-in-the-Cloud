# -*- coding: utf-8 -*-

####Splice filtering
####Nicole Lama
####RNA_Seq_In_Cloud
####2019/03/11

################# STEPS #######################################################
## input: gtf files

## filtering steps:
#i. remove known txs A
#ii. remove intronless tx
#iii. keep only internal exons (remove begining and terminal exons)
#iv. Filter retained introns

## output: tsv (chr,exon_start,exon_stop,strand)
#import os,sys,csv
#from glob import glob
#import numpy as np
#import pandas as pd



#define important variables and arrays/dictionaries

gtf = []
numOfReads = []
allExonCoords = {}
#rejectTx = [] ## store rejected transcripts with less than 2 introns
exonTxH = {}
exonKey = "NaN"
exonCount = 0
uniqueExonId = 0

#################### FUNCTIONS ###############################################
def examineRetainedIntrons(exonH,exonCoords):
    delItems = set()
    delId = []
    for exonInfo in exonH:
        ex1_id = exonInfo.split(",")[0]
        ex1_start=exonInfo.split(",")[2]
        ex1_end=exonInfo.split(",")[3]

        for s in exonCoords.values():
            if ex1_start == s[0]:
#                print(s[0],ex1_id)
                for e,exId in zip(exonCoords.values(),exonCoords.keys()):
#                    print(exId,ex1_id)
                    if ex1_end == e[1] and int(ex1_id) != exId:
                        if str(exId) not in delId:
#                            print(s[0],e[1],exId,ex1_id)
                            delItems.add(exonInfo)
                            delId.append(ex1_id)
                 
    for exon in delItems:  
        print("removing: " + exon)
        del exonH[exon]
    return exonH

def generateGenCodeHash(genEx):
    genCode = {}
    with open(genEx,'r') as g:
        for row in g:
            tx = row.split("\t")
            if row[0] == "#": #filter out header and gaps
                continue
            if tx[2] == "exon":
                exonKey = "{0},{1},{2},{3}".format(tx[0],tx[3],tx[4],tx[6])
                tx_id=tx[8].split(";")[1].split(" ")[2]
                if exonKey in exonTxH.keys():
                    genCode[exonKey].append(tx_id)
                else:
                    genCode[exonKey] = []
                    genCode[exonKey].append(tx_id)
    return genCode

    
def queryGenCode(genEx,eEx):
    genVerifiedExons=[]
    for exon in eEx.keys():
        exon = ",".join(exon.split(",")[1:5])
        if genEx.get(exon,0) != 0:
            genVerifiedExons.append(exon)
    return genVerifiedExons

            
            
    
def writeOutFile(exons,tsvOut):
    with open(tsvOut, "w+", newline = '') as tOut:
            tOut.write("chromosome\texon_start\texon_end\tstrand_+/-\n")
            for row in exons:
                row=row.split(",")
                print(row)
                chrm = row[0]
                exS = row[1]
                exE = row[2]
                strd = row[3]
                tOut.write("{0}\t{1}\t{2}\t{3}\n".format(chrm,exS,exE,strd))
            


########################## INIT ###############################################

humGenCode="gencode.v29.annotation.gtf"

with open("sample.500.gtf",'r') as g:
    for row in g:
        tx = row.split("\t")
        if row.split(" ")[0] == "#": #filter out header and gaps
            continue
        if tx[2] == "transcript":
            if exonTxH.get(exonKey,0) != 0 and exonCount > 3:
#                print(exonCount,exonKey)
                exonTxH.popitem()
                exonCount = 0
        if any("reference_id" in attr for attr in tx[8].split(";")): #filter known tx
            allExonCoords[uniqueExonId]=[tx[3],tx[4]]
            uniqueExonId += 1
            continue
        if tx[2] == "exon":
            if not any("1" in exNum for exNum in tx[8].split(";")[2]): #filter for less than 2 introns and remove first exon
                exonKey = "{0},{1},{2},{3},{4}".format(uniqueExonId,tx[0],tx[3],tx[4],tx[6])
                tx_id=tx[8].split(";")[1].split(" ")[2]
                if exonKey in exonTxH.keys():
                    exonTxH[exonKey].append(tx_id)
                    allExonCoords[uniqueExonId]=[tx[3],tx[4]]
                    uniqueExonId += 1
                else:
                    exonTxH[exonKey] = []
                    exonTxH[exonKey].append(tx_id)
                    allExonCoords[uniqueExonId]=[tx[3],tx[4]]
                    uniqueExonId += 1
                exonCount+=1
#exonTxH.popitem() #pop off last item of hash since it is external exon

print(len(exonTxH))   
filteredTx = examineRetainedIntrons(exonTxH,allExonCoords)
print(len(filteredTx))
numOfReads.append(len(filteredTx))
genH = generateGenCodeHash(humGenCode)
genVerExons = queryGenCode(genH,filteredTx)
        #rejectTx.append(row) 
print(genVerExons)
print(len(genVerExons))

writeOutFile(genVerExons,"test.tsv")





