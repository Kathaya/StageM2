import argparse
import tskit
import numpy as np
import cyvcf2
import tsinfer
import tqdm
import json
import pandas as pd


def chromosome_length(vcf):
    """
    Get the chromosome_length
    """
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]

def add_diploid_sites(vcf, samples):
    progressbar = tqdm.tqdm(total=samples.sequence_length, desc="Read VCF", unit='bp')
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
    pos = 0
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        progressbar.update(variant.POS - pos)
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF] + variant.ALT
        ancestral = variant.INFO.get('AA', variant.REF)
        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        allele_index = {old_index: ordered_alleles.index(allele)
            for old_index, allele in enumerate(alleles)}
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [allele_index[old_index]
            for row in variant.genotypes for old_index in row[0:2]]
        samples.add_site(pos, genotypes=genotypes, alleles=alleles)


def kc(file1, file2):
    """
    Compute the Kendall Coljin distance between 2 set of tree sequence with ponderation between surface recrouvrement of 2 trees. 
    Input :
        - file1 : treesequence file
        - file2 : treesequence file
    Output : 
        - None, used to write via bash in a results file as vcf format
    """
    ts1 = tskit.load(file1).simplify()
    ts2 = tskit.load(file2).simplify()
    #For each tree in real data we cumpute distance with each tree in tsinfer simulated data that share at least one SNP in covering range of both trees

    distance=[]
    chevauchement=[]
    i=j=0
    while(True):
        t1=ts1.at_index(i, sample_lists=True)
        t2=ts2.at_index(j, sample_lists=True)
        
        #Si jamais la borne sup de arbre 2 inférieur borne inf de l'arbre 1, augmente l'indic de 2
        if(t1.interval[0]>t2.interval[1]):
            j += 1
        #Pareil pour 1
        elif(t2.interval[0]>t1.interval[1]):
            i += 1

        #Maintenant les arbres ont au moins un SNP en commun
        else :
        #The value 0 indicate that we don't consider branch length but only topology
        #Calcul de la distance
            distance.append(t2.kc_distance(t1,0))
        #Calcul du chevauchement

            #if t1 inside t2
            if (t1.interval[0]<=t2.interval[1] and t1.interval[0]>=t2.interval[0]) and (t1.interval[1]<=t2.interval[1] and t1.interval[1]>=t2.interval[0]):
                chevauchement.append(t1.interval[1]-t1.interval[0])
                i += 1
            #if t2 inside t1
            elif (t2.interval[0]<=t1.interval[1] and t2.interval[0]>=t1.interval[0]) and (t2.interval[1]<=t1.interval[1] and t2.interval[1]>=t1.interval[0]):
                chevauchement.append(t2.interval[1]-t2.interval[0])
                j += 1
            #if t1 part inside 2
            elif (t1.interval[1]<=t2.interval[1] and t1.interval[1]>= t2.interval[0]):
                chevauchement.append(t1.interval[1]-t2.interval[0])
                i += 1
            #if t2 part inside t1    
            else:
                chevauchement.append(t2.interval[1]-t1.interval[0])
                j += 1
        #Condition de sortie de boucle
        if(i>= ts1.last().index or j >= ts2.last().index):
            break
    distList=np.array(distance)
    chevList=np.array(chevauchement)
    chevxdistList=distList*chevList
    
    mean=sum(chevxdistList)/sum(chevList)
    l1bis = np.repeat(mean, len(distList))
    l1bis = distList-l1bis
    l1bis = np.square(l1bis)
    l3bis = l1bis*chevList
    var = sum(l3bis)/sum(chevList)
    #Repeter la valeur de l3 sur une len(l1)
    #soustrait cette l1 avec cette liste
    #Mise au carré
    #division par sum(l2)
    
    print("{}\\t{}".format(mean, var))


#function used to do the inference with tsinfer from a vcf file
def ts_sim(file1, file2, file3):
    """
    Inference via tsinfer of coalescente trees via a vcf file as summary of a population genetics

    Input:
        -file1 : vcf file summarizing genetical state of a population
        -file2: sample files used to summarize the diploid_site
        -file3: trees file where to write the infered treesequence
    """
    vcf = cyvcf2.VCF(file1)
    with tsinfer.SampleData(path=file2,
                        sequence_length=chromosome_length(vcf)) as samples:
        add_diploid_sites(vcf, samples)

    print("Sample file created for {} samples ".format(samples.num_samples) +
          "({} individuals) ".format(samples.num_individuals) +
          "with {} variable sites.".format(samples.num_sites), flush=True)

    # Do the inference
    ts = tsinfer.infer(samples, num_threads=1)
    print("Inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        ts.num_trees, ts.sequence_length/1e6, ts.num_edges))
    ts.dump(file3)

    #Writing of simulation.samples file
    with tsinfer.SampleData(
            path=file2, sequence_length=ts.sequence_length,
            num_flush_threads=1) as sample_data:
        for var in ts.variants():
            sample_data.add_site(var.site.position, var.genotypes, var.alleles)


#Function used to create .map file needed for relate execution
def map(file1, file2):
    with open(file1, 'r') as fillin, open(file2, 'w') as fillout:
        lignes = fillin.readlines()
        for ligne in lignes:
            snp = float(ligne.split()[2])
            fillout.write("{}\t2e-8\t{}\n".format(int(snp), round(snp/10000*0.1, 5)))

def ntree(file1):
    """
    Function calculating the number of tree in a TreeSequence file
    
    Input : 
        -file1 : TreeSquence file
    """
    ts = tskit.load(file1)
    print(ts.last().index)

if __name__ == "__main__":
    function_map = {'kc' : kc, 'ts_sim' : ts_sim, 'map' : map, 'ntree' : ntree}

    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=function_map.keys(), help="Script to use, choices are kc; map; ts_sim")
    parser.add_argument("--file1", help="First file")
    parser.add_argument("--file2", help="Second file")
    parser.add_argument("--file3", help="Thrid file")
    parser.add_argument("--threads", help="Number of threads use for simulations", type=int)

    args = parser.parse_args()



    func = function_map[args.command]

    if args.command == 'kc':
        if args.file1 == None or args.file2 == None:
            print("Kc distance function usage : python function_map.py ts_sim --file1 firstfile --file2 secondfile.")
        else:
            func(args.file1, args.file2)

    if args.command == 'ts_sim':
        if args.file1 == None or args.file2 == None or args.file3 == None:
            print("ts_sim function usage : python function_map.py kc --file1 vcffile  --file2 temporarysamplefilepath --file3 treefile.")
        else:
            func(args.file1, args.file2, args.file3)

    if args.command == 'map':
        if args.file1 == None or args.file2 == None:
            print("Map function usage : python function_map.py map --file1 hapsfile --file2 outputfile")
        else:
            func(args.file1, args.file2)
    if args.command == 'ntree':
        if args.file1 == None:
            print("Ntree usage v: python function_map.py ntree --file1 treesequencefile")
        else:
            func(args.file1)






