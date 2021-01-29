import pandas as pd
import tskit
import msprime
import os
import sys
import numpy as np


i = int(sys.argv[1])
dirName = str(sys.argv[2])
#Filename of tools used for simulate data
#Either msprime or slim
data_sim = sys.argv[3]


datalistNewick = []
tslistNewick = []
relistNewick = []
mslistNewick2 = []
datafile = str(dirName) + "/" + str(data_sim) + str(i) + ".trees"
tsfile = str(dirName) + "/ts" + str(i) + ".trees"
refile = str(dirName) + "/relate" + str(i) + ".trees"
msfile2 = str(dirName) + "/msprime" + str(i) + "_2.trees"
dataNewick = str(dirName) + "/" + str(data_sim) + str(i) + ".newick"
tsNewick = str(dirName) + "/ts" + str(i) + ".newick"
reNewick = str(dirName) + "/relate" + str(i) + ".newick"
msNewick2 = str(dirName) + "/msprime" + str(i) + "_2.newick"
ms = tskit.load(datafile)
ts = tskit.load(tsfile)
re = tskit.load(refile)
ms2 = tskit.load(msfile2)

for tree in ms.trees():
	datalistNewick.append(tree.newick())

for tree in ts.trees():
	tslistNewick.append(tree.newick())

for tree in re.trees():
	relistNewick.append(tree.newick())

for tree in ms2.trees():
	mslistNewick2.append(tree.newick())

writer = open(dataNewick, "w")
tmpNewick = ""
for tree in range(len(datalistNewick)):
	tmpNewick = tmpNewick + str(datalistNewick[tree]) + "\n"
writer.write(tmpNewick)
writer.close

writer = open(tsNewick, "w")
tmpNewick = ""
for tree in range(len(tslistNewick)):
	tmpNewick = tmpNewick + str(tslistNewick[tree]) + "\n"
writer.write(tmpNewick)
writer.close

writer = open(reNewick, "w")
tmpNewick = ""
for tree in range(len(relistNewick)):
	tmpNewick = tmpNewick + str(relistNewick[tree]) + "\n"
writer.write(tmpNewick)
writer.close

writer = open(msNewick2, "w")
tmpNewick = ""
for tree in range(len(mslistNewick2)):
	tmpNewick = tmpNewick + str(mslistNewick2[tree]) + "\n"
writer.write(tmpNewick)
writer.close

