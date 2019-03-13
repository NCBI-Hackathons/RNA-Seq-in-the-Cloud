# -*- coding: utf-8 -*-

NOTE: USE DETECTNOVELEXONS.PY
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
from glob import glob
import gzip,csv,os
from timeit import timeit

#define important variables and arrays/dictionaries

gtf = []
numOfReads = []
allExonCoords = {}
uniqueExonId = 0
exonTxH = {}
exonKey = "NaN"
exonCount = 0
uniqueExonId = 0

#################### FUNCTIONS ###############################################

def checkCoords(d,coord):
    if d.get(coord,"NaN") == "NaN": 
        d[coord] = 1
    else:
        d[coord] +=1
    return d

   
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
##                            print(s[0],e[1],exId,ex1_id)
                            delItems.add(exonInfo)
                            delId.append(ex1_id)
#                 
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
                if exonKey in genCode.keys():
                    genCode[exonKey].append(tx_id)
                else:
                    genCode[exonKey] = []
                    genCode[exonKey].append(tx_id)
    print("finished Gen Hash")
    return genCode

def genAnnotatedHash(a):
    aHash = {}
    with gzip.open(a,'r') as f:
        for row in f:
		tx = row.split("\t")
                tx[3]=tx[3].strip("\n")
		#print(tx[0])
		#chr1    +       155308807       155308965
		#print("start", tx[3].strip("\n"), "end")
		exonKey = "{0},{1},{2},{3}".format(tx[0],tx[2],tx[3],tx[1])
		#print(exonKey)
		aHash[exonKey]=1
    #print(aHash)
    return aHash
    
def queryRef(refH,eEx):
    genVerifiedExons = {}
    for exon in eEx.items():
	#print(exon[0])
        if refH.get(exon[0],0) != 0:
            genVerifiedExons[exon[0]] = exon[1]  
    return genVerifiedExons

            
            
    
def writeOutFile(exons,tsvOut):
    with open(tsvOut, "w+") as tOut:
            tOut.write("chromosome\texon_start\texon_end\tstrand\n")
            for row in exons.items():
                ex=row[0].split(",")
                chrm = ex[0]
                exS = ex[1]
                exE = ex[2]
                strd = ex[3]
                txId = row[1]
                tOut.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrm,exS,exE,strd,txId))           

def writeLogFile(f,txC,tEx,filEx,anEx):
    with open("{0}.logFile.out".format(f), "w+") as tOut:
            tOut.write("total_lines\ttotal_exons\tfiltered_exons\tannotated_exons\n")
            tOut.write("{0}\t{1}\t{2}\t{3}\n".format(txC,tEx,filEx,anEx))
########################## INIT ###############################################

humGenCode = "./GenCode/gencode.v29.annotation.gtf"
print("Generating Gencode Exon Information")
genH = generateGenCodeHash(humGenCode)
print("Generating Annotated Exon Information")
annotated = "./ourAlignments/all_tx_unique_exons.tsv.gz"
aH = genAnnotatedHash(annotated)

expGTFs = glob("/data/study_gtf/SRP*.gtf")
print("parsing through experimental gtfs")
for f in expGTFs:
    fileName = os.path.basename(f)
    numOfReads = []
    allExonCoords = {}
    uniqueExonId = 0
    exonTxH = {}
    exonKey = "NaN"
    exonCount = 0
    uniqueExonId = 0
    txCount = 0
    numOfExons = 0
    print("working on file: "+ fileName)
    with open(f,'r') as g:
        for row in g:
            txCount += 1
            tx = row.split("\t")
    
            if row.split(" ")[0] == "#": #filter out header and gaps
                continue
    
            if tx[2] == "transcript":
                if exonTxH.get(exonKey,0) != 0 and exonCount > 3:
    #                print(exonCount,exonKey)
                    exonTxH.popitem()
                    exonCount = 0
            #exonStartCoordsDict = checkCoords(exonStartCoordsDict,tx[3])
            #exonEndCoordsDict = checkCoords(exonEndCoordsDict,tx[4])                
            attributes = tx[8].split(";")
            if any("reference_id" in attr for attr in attributes): #filter known tx
                allExonCoords[uniqueExonId]=[tx[3],tx[4]]
                uniqueExonId += 1
                continue
            if tx[2] == "exon":
		numOfExons += 1
	        if int(attributes[2].split(" ")[2].strip('\"')) > 1: #filter 2 introns and remove first exon
                    exonKey = "{0},{1},{2},{3}".format(tx[0],tx[3],tx[4],tx[6])
                    #print(tx[6])
                    tx_id=attributes[1].split(" ")[2]
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
    if len(exonTxH) > 2 :        
        exonTxH.popitem() #pop off last item of hash since it is external exon
    
    print("# of total lines: {0}".format(txCount))
    print("# of total exons: {0}".format(numOfExons))
#    print("preprocessing found {} novel exons".format(len(exonTxH)))   
    filteredTx = examineRetainedIntrons(exonTxH,allExonCoords)
    print("preprocessing found {} novel exons".format(len(filteredTx)))
    print("querying genCode for novel exons")
    genVerExons = queryRef(genH,filteredTx)
    print("querying annotated file for novel exons")
    annVerExons = queryRef(aH,filteredTx)
    print("novel exons verified by genCode:{0} ".format(len(genVerExons)))
    print("novel exons verified by annotat:{0} ".format(len(annVerExons)))
    
    writeOutFile(filteredTx,"{0}.allExon.rmRetInts.G.tsv".format(fileName))
    writeOutFile(annVerExons,"{0}.Ann.Gen.rmRetInts.tsv".format(fileName))
    fileName += ".GenCode"
    writeLogFile(fileName,txCount,numOfExons,len(filteredTx),len(annVerExons))





