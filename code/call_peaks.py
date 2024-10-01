#!/usr/bin/env python

import re
import sys

if len(sys.argv) < 3:
	print("Usage: <keep> <lrt table (including global significance)>")
	sys.exit(0)

keep = open(sys.argv[1])
tab = open(sys.argv[2])
tab.readline()

sigbits = []

tab_tag = tab.readline()
ttag = tab_tag.split()[0]

last = ""

for l in keep.readlines():
	tag = l.split()[0]
	sig = 0
	logFC = 0
	pVal = 1
	if tag == ttag:
		sig = int(float(tab_tag.split()[5]))
		logFC = float(tab_tag.split()[1])
		pVal = float(tab_tag.split()[4])
		tab_tag = tab.readline()
		if not tab_tag:
			tab_tag = "XX XX"
		ttag = tab_tag.split()[0]
	sigbits.append((tag, sig, logFC, pVal))


start = -1
end = -1
logFC = 0
size = 0
sig = 0



lapseOneTag = False
for i in range(0, len(sigbits), 2):
	if i+3 < len(sigbits) and (sigbits[i+3][3]<0.05 and  sigbits[i][3]<0.05 and sigbits[i+3][2]>0 and sigbits[i][2]>0):
		if sigbits[i][1] > 0:
			sig += 1
		if sigbits[i+3][1] > 0:
			sig += 1
		size += 2
		logFC += sigbits[i][2] + sigbits[i+3][2]
                #peak+=sigbits[i][0]+","+sigbits[i+3][0]+","
		if start < 0:
			start = i
		end = i+3
		lapseOneTag = True
	elif lapseOneTag:
		lapseOneTag = False
	elif start >= 0:
		s=int(re.split("[\_:-]",sigbits[start][0])[1])
		e=int(re.split("[\_:-]",sigbits[end][0])[1])
		m=(s+e)/2
		unmeth = 0
		for j in range(end+3, len(sigbits),1):
			x=int(float(re.split("[\_:-]",sigbits[j][0])[1]))
			if x-m >500:
				break
			if not sigbits[j][1] and not sigbits[j-1][1]:
				unmeth +=1
		for j in range(start -1, 0, -1):
			x=int(re.split("[\_:-]",sigbits[j][0])[1])
			if m-x >500:
				break
			if not sigbits[j][1] and not sigbits[j+1][1]:
				unmeth +=1

		ll = logFC/size
		print ("%s\t%d\t%d\t%d\t%f\t%f\t%d" %(re.split("[\_:-]",sigbits[start][0])[1], s, e, size, (float(size))/(unmeth+size), ll, sig))
		last = 0
		skip = 0
		logFC = 0
		start = -1
		end = -1
		size = 0
		sig = 0
              
