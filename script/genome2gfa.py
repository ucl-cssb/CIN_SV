#!/usr/bin/env python3

#########################################################################
# Originally written by Tass Abdessalem, revised by Bingxin Lu
# Description: converting the genome_* files to GFA format
#########################################################################


import sys
import pprint
from genome import GenomeFile

if (len(sys.argv) != 2):
    print("Usage: python gfa.py <inputfile>")
    sys.exit(1)

data = GenomeFile(sys.argv[1]).records

gfaFile = []
csvFile = ""

gfaFile.append(f"H\tVN:Z:1.0\n")
csvFile += (f"Ref,Type,Chr,Hap,Start,End\n")

# add suffix to allow duplicated nodes
nodeEnds = {}
nodeCount = {}
nodeCount_adj = {}

addedSegmentIDs = []

for record in data:
    # print(vars(record))
    for node in record.nodes:  
        # print(vars)      
        gfaLine = []
        if (node.type == "INTERVAL" ):           
            # used to find adjacent intervals of an adjacency edge based on positions
            # Due to duplications, the same position may have different ends
            key_start = f"{node.fromChr}/{node.fromHaplotype}/{node.fromPosition}"
            key_end = f"{node.toChr}/{node.toHaplotype}/{node.toPosition}"
            nodeCount[key_start] = nodeCount.get(key_start, 0) + 1
            nodeCount[key_end] = nodeCount.get(key_end, 0) + 1     
            
            nodeRef = f"{node.type}_Chr{node.fromChr}:Hap{node.fromHaplotype}:{nodeCount[key_start]}:{nodeCount[key_end]}:{node.fromPosition}-{node.toPosition}"
                              
            if nodeCount[key_start] == 1:
                nodeEnds[key_start] = nodeRef
            else:
                nodeEnds[f"{key_start}/{nodeCount[key_start]}"] = nodeRef
                
            if nodeCount[key_end] == 1:
                nodeEnds[key_end] = nodeRef
            else:
                nodeEnds[f"{key_end}/{nodeCount[key_end]}"] = nodeRef           

            # cannot simply remove duplicated segment due to duplications
            # if (nodeRef not in addedSegmentIDs):
            #     addedSegmentIDs.append(nodeRef)
            # else:
            #     continue
            
            # write segments or intervals first
            csvFile += (f"{nodeRef},{node.type},{node.fromChr},{node.fromHaplotype},{node.fromPosition},{node.toPosition}\n")
            gfaLine.append("S")
            gfaLine.append((nodeRef))
            gfaLine.append("*")
            size = abs(int(node.toPosition) - int(node.fromPosition))               
            gfaLine.append(f"LN:i:{size}")                     
            gfaFile.append("\t".join(gfaLine) + "\n")  
            
# pprint.pprint(nodeEnds)
            
for record in data:
    # print(vars(record))
    for node in record.nodes:  
        # print(vars)      
        gfaLine = []                      
        if (node.type == "VAR" or node.type == "REF"):
            key_start = f"{node.fromChr}/{node.fromHaplotype}/{node.fromPosition}"
            key_end = f"{node.toChr}/{node.toHaplotype}/{node.toPosition}"
            
            nodeCount_adj[key_start] = nodeCount_adj.get(key_start, 0) + 1
            nodeCount_adj[key_end] = nodeCount_adj.get(key_end, 0) + 1
            
            # print(nodeCount_adj[key_start])
            # print(nodeCount_adj[key_end])
            
            # find the adjacent intervals based on count of appearance
            if nodeCount_adj[key_start] == 1:
                fromRef = nodeEnds[key_start]
            else:
                fromRef = nodeEnds[f"{key_start}/{nodeCount_adj[key_start]}"]
                
            if nodeCount_adj[key_end] == 1:
                toRef = nodeEnds[key_end]
            else:
                toRef = nodeEnds[f"{key_end}/{nodeCount_adj[key_end]}"]
                                            
            gfaLine.append("L")
            gfaLine.append(fromRef)
            gfaLine.append("+")
            gfaLine.append(toRef)
            gfaLine.append("+")
            gfaLine.append("0M")
            gfaLine.append(f"ID:Z:{node.raw}")
            gfaFile.append("\t".join(gfaLine) + "\n")


with open(f"{sys.argv[1].replace('.tsv', '')}.gfa", "w+") as file:
    file.writelines(gfaFile)
    
with open(f"{sys.argv[1].replace('.tsv', '')}.csv", "w+") as file:
    file.write(csvFile)
    


