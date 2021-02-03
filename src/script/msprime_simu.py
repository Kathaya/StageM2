import msprime
import numpy as np
import sys
import tsinfer
import tskit
import pandas as pd
import argparse

def get_arguments():
    """
    Lecture des arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-vp', "--vcf_path", help="Where to write the .vcf file")
    parser.add_argument('-tp', "--tree_path", help="Where to write the .trees file")
    parser.add_argument('-p', type=int, default=100, help="Population size")
    parser.add_argument('-s', type=int, default=30, help="Sample size")
    parser.add_argument('-f', type=int, default=1, help="Factor of demegraphic changes")
    parser.add_argument('-g', type=int, default=5, help="Number of generation before change in demogaphy")
    parser.add_argument('-rho', type=float, default=2e-8, help="Recombination rate")
    parser.add_argument('-l', type=float, default=1e5, help="Chromosome length")
    parser.add_argument('-mu', type=float, default=1e-8, help="Mutation rate")
    args = parser.parse_args()
    return args
    
def verification(tree_sequence):
    """
    Verification of mutation at the first and last position of the genome. 
    In that case, we will regenerate the information.
    """
    flag=0
    valeur = []
    for var in tree_sequence.variants():
        index = int(round(var.site.position))
        if index == 0 or index == L:
            flag=1
            break
        else:
            valeur.append(index)
    if flag == 0 and len(valeur)>=10:
        return 0
    else:
        return 1


if __name__ == "__main__":
    args = get_arguments()
    
    tree_path = args.tree_path
    vcf_path = args.vcf_path
    popSize = args.p # present day population size
    sampleSize = args.s # sample size
    factorChange = args.f # factor for demographic changes, see use below
    genChange = args.g # number of generations before present the change in demography occured

    rho = args.rho # number of recombination events per base pair = recombination rate
    L = args.l # length of the chromosome
    mu = args.mu # rate of mutations

    # Population size before the change
    # Thus : dividing by the factor : sudden expansion, multiplying by the factor : sudden contraction
    # In the example here : the population size is 10 today, but was 100*10 just 5 generations ago
    popSizeBefore = popSize * factorChange



    # List of events : here we have one event only - you can add some
    events = [msprime.PopulationParametersChange(time=genChange, initial_size=popSizeBefore)]


    x=0
    while 1:
            x += 1
            tree_sequence = msprime.simulate(sample_size=sampleSize, Ne=popSize,
                                             length=L, recombination_rate=rho, mutation_rate=mu,
                                             demographic_events = events)
            if verification(tree_sequence)==0:
                    break
            if x==100:
                    print("PROBLEMEEEEEEEEE")
                    break

    tree_sequence.dump(tree_path+".trees")
    with open(vcf_path+".vcf", "w") as vcf_file: 
            tree_sequence.write_vcf(vcf_file, ploidy=2, position_transform="legacy")
    if x != 1:
        print("Nombre d'essaie : {}\n".format(x))



"""
    tree_sequence = msprime.simulate(sample_size=sampleSize, Ne=popSize,length=L,
            recombination_rate=rho, mutation_rate=mu,demographic_events = events)


    tree_sequence.dump(tree_path+".trees")

    with open(vcf_path+".vcf", "w") as vcf_file: 
            tree_sequence.write_vcf(vcf_file, ploidy=2)

    data = pd.read_csv(vcf_path+".vcf", sep="\t", header=5)
    data.sort_values("POS", inplace = True)
    data.drop_duplicates(subset="POS", keep=False, inplace=True)
    data2=data[data.POS != 100000]
    data3=data2[data.POS != 0]



    with open(vcf_path+".vcf","r") as fi:
        header = []
        for ln in fi:
            if ln.startswith("##"):
                header.append(ln)

    with open(vcf_path+".vcf", "w") as fillout:
        for ligne in header:
            fillout.write(ligne)
    data3.to_csv(vcf_path+".vcf", sep="\t", mode='a', index=False)
        #fillout.write(data.to_string(index=False))
"""

