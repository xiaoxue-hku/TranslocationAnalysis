import os
import sys
import tokenize
import random
import gzip
import math

#Usage: makeCoverageProfile.py [coverageBed output] [outfile prefix] [bins]
#Makes profile matrix output from coverageBed depth command
#[bins] is optional. It is used to specify how many bins to place the profile into.
#If a value is selected, an extra file "outfile_prefix"_binned.txt will be generated
#with binned profiles for each bed entry and the average profile at the end.
#
#It is also possible to pipe the output of coverageBed directly to the command (e.g.):
# coverageBed -a 5utr.total.bed -b ~/Data/RLoop-Gquad/mutations/SCC/SCC_XPC_Mutations.bed -d | python ~/programming/2016scripts/makeCoverageProfile.py - 5utr.total.XPC_prof 20



if len(sys.argv) < 2:
   print "Usage: makeCoverageProfile.py [coverageBed output] [outfile prefix] [bins]"
   exit()

bins = 0
if len(sys.argv) > 3:
   bins = int(sys.argv[3])
if sys.argv[1] <> "-":
   f = open(sys.argv[1],'r')
else:
   f = sys.stdin
outM = open(sys.argv[2]+"_whole.txt",'w')
if bins <> 0:
   outB = open(sys.argv[2]+"_binned.txt",'w')

profiles = {} #{name:[]}
last = None
laststrand = None
curprof = []
sumprof = []
sumvals = []
for line in f:
    tmp = line[:-1].split("\t")
    cur = tmp[0]+":"+tmp[1]+"-"+tmp[2]
    if last <> None and last <> cur:
       if laststrand == "-":
          profiles[last+"_"+str(len(curprof))+"_"+laststrand] = curprof[::-1]
       else:
          profiles[last+"_"+str(len(curprof))+"_"+laststrand] = curprof[:]
       pos = 0
       for i in profiles[last+"_"+str(len(curprof))+"_"+laststrand]:
           try:
               sumprof[pos]+=i
               sumvals[pos]+=1
           except IndexError:
               sumprof.append(i)
               sumvals.append(1)
           pos+=1
       curprof = []
    last = cur
    laststrand = tmp[-3]                            #This is the strand column, assumes that it is the last column in original bed annotation file. Might need to change if different
    curprof.append(int(tmp[-1]))

for i in profiles.keys():
    outM.write(i+"\t")
    for j in profiles[i]:
        outM.write(str(j)+"\t")
    outM.write("\n")
outM.write("Average\t")
print len(sumprof)
for i in range(0,len(sumprof)):
    outM.write(str(float(sumprof[i])/float(sumvals[i]))+"\t")


if bins <> 0:
   sumbinprof = [0]*bins
   sumbinval = [0]*bins
   for i in profiles.keys():
       if len(profiles[i]) >= bins:
         outB.write(i+"\t")
         binprof = [0]*bins
         binvals = [0]*bins
         c = 0
         for j in profiles[i]:
             pos = int(math.floor((float(c)/float(len(profiles[i])))*float(bins)+0.00000001))        #0.000000001 offset is to account of floating point errors
             binprof[pos]+=j
             binvals[pos]+=1
             sumbinprof[pos]+=j
             sumbinval[pos]+=1
             c+=1
         for j in range(0,bins):
             outB.write(str(float(binprof[j])/float(binvals[j]))+"\t")
         outB.write("\n")
   outB.write("Average\t")
   for j in range(0,bins):
       try:
           outB.write(str(float(sumbinprof[j])/float(sumbinval[j]))+"\t")
       except ZeroDivisionError:
           outB.write("0\t")
   outB.write("\n")

f.close()
outM.close()
if bins <> 0:
   outB.close()


