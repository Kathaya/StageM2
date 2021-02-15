#!/usr/bin/python3

import sys, os, subprocess, argparse, random
import tskit, msprime, pyslim, numpy as np, pandas as pd
"""
Manage the SLiM simulation pipeline
- generate slim simulation file with correct parameters
- Do slim simulation and generate vcf
- Recapitate slim simulation and write vcf 
"""

def get_arguments():
    """
    Lecture des arguments
    """
    parser = argparse.ArgumentParser(description='')

    #Mandatory arguments :
    parser.add_argument('--path', required=True ,help="Generic name of output files")
    #parser.add_argument('-n', '--name', type=str, required=True,  help="<Required> Simulation_name")

    #Optionnal arguments :
    parser.add_argument('-p', '--popsize', type=int, default=1000, help="Population size")
    parser.add_argument('-s', '--sampsize', type=int, default=30, help="Sample size")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads. Default: 4")
    parser.add_argument('-r', '--rho', type=float, default=2e-8, help="Recombination rate")
    parser.add_argument('-m', '--mu', type=float, default=1e-8, help="Mutation rate")
    parser.add_argument('-l', '--klen', type=float, default=1e5, help="Chromosome length")
    parser.add_argument('-o', '--output', type=int, default=10, help="Number of generation to calculate summary statistics ?")
    parser.add_argument('--tsr_start', type=int, default=0, help="Generation number for starting TSR")
    parser.add_argument('--tsr_end', type=int, default=10, help="Generation number for stop TSR")
    parser.add_argument('-a', '--alpha', type=float, default=0, help="TSR power, between 0 (no TSR) and 2 (Max TSR)")
    parser.add_argument('--nb_hist', type=int, default=21, help="Je sais pas")
    parser.add_argument('-d', type=int, default=1, help="Ditribution law for reproduction, either Poisson like (1) or gama like (2), default 1: ")
    parser.add_argument('-tp', '--trans_parent', default=2, help="Parental transmission either 0 for women, 1 for men or 2 for biparental")
    parser.add_argument('-gv', '--generations_vec', nargs='+', default=10, type=int, help="Condition d'arret de sortie d'arbres passer en vector sous la forme -gv x y z")
    args = parser.parse_args()
    return args

def simplification(sampleSize):
    """
    Reduce the numbers of individuals to write
    """
    keep_indivs = np.random.choice(ts.individuals_alive_at(0), sampleSize, replace=False)
    keep_nodes = []
    for i in keep_indivs:
       keep_nodes.extend(ts.individual(i).nodes)
    sts = ts.simplify(keep_nodes)
    return sts

def recapitation(simu_name, mu, Ne, sampleSize):
    """
    Function used to recapitate the SLiM simulation and write appropriate vcf file
    """
    ts = pyslim.load(simu_name + ".trees")
    ts = ts.recapitate(recombination_rate = mu, Ne = popSize)
    ts = pyslim.SlimTreeSequence(msprime.mutate(ts, mu))
    halfsampleSize = int(sampleSize/2)
    sts=simplification(halfsampleSize)
    while(sts.num_individuals != halfsampleSize):
        sts=simplification(halfsampleSize)
    ts.dump(simu_name+"_full.trees")
    sts.dump(simu_name+".trees")
    with open(simu_name+".vcf", "w") as vcf_file: 
        sts.write_vcf(vcf_file)

    data = pd.read_csv(simu_name+".vcf", sep="\t", header=5)
    data.sort_values("POS", inplace = True)
    data.drop_duplicates(subset="POS", keep='first', inplace=True)
    data2=data[data.POS < 100000]
    print(data2)
    data3=data2[data.POS > 0]
    print(data3)
    with open(simu_name+".vcf","r") as fi:
        header = []
        for ln in fi:
            if ln.startswith("##"):
                header.append(ln)

    with open(simu_name+".vcf", "w") as fillout:
        for ligne in header:
            fillout.write(ligne)
    data3.to_csv(simu_name+".vcf", sep="\t", mode='a', index=False)


if __name__ == "__main__":
    wdpath = sys.argv[0].replace("slimEngine.py", "")
    args = get_arguments()
    path=args.path
    popSize = args.popsize
    #name = args.name
    sampleSize = args.sampsize
    threads = args.threads
    rho = args.rho
    mu = args.mu
    L = int(args.klen)
    output = args.output
    tsr_start = args.tsr_start
    tsr_end = args.tsr_end
    alpha = args.alpha
    nb_hist = args.nb_hist
    distrib = args.d
    trans_parent = args.trans_parent
    generations_vec = args.generations_vec
    if isinstance(generations_vec, list) or len(generations_vec) == 1:
        generations_vec = int(generations_vec[0])

    #print("slimEngine simulation : {}\n".format(i))

    # Dictionnary creation for managing the arguments
    dict_replace = { "enr_Bash" : output, "deb_Bash" : tsr_start,
    "fin_Bash" : tsr_end, "trans_parent" : trans_parent, "L_Bash" : L,
    "rec_Bash" : rho, "tailleBash" : popSize, "nb_Bash" : nb_hist,
    "trans_Bash" : alpha, "distrib_Bash" : distrib }


    # Write into the SLiM code the choosen parameters
    with open(wdpath+"simupop.c", "r") as fillin:
        slim_string = fillin.read()
    
    for key in dict_replace:    
        slim_string = slim_string.replace(str(key), str(dict_replace[key]))

    simu_name = path
    # Adding in the end of the SLiM code the choosen output_trees
    
    if isinstance(generations_vec, int):
        slim_string =  slim_string + str(int(generations_vec)) + " late() {sim.treeSeqOutput('./" + str(simu_name) + ".trees');}" + "\n"
    else:
        for generation in generations_vec:
            slim_string =  slim_string + str(generation) + " late() {sim.treeSeqOutput('./" + str(simu_name) + "-" +  str(generation) + ".trees');}" + "\n"



    # Writing the new SLiM code

    slim_new_name = simu_name + ".c"
    with open(slim_new_name, 'w') as f:
        f.write(slim_string)

    # Launching SLiM

    slim_output = os.popen("~/slim/build/slim " + slim_new_name).read()


    # Writing output in csv format

    csv_header = "Generation Alpha Mere_enfant Pere_enfant Parent_enfant" # English
    for hist in range(0, int(nb_hist)):
        csv_header = csv_header + " Sibship" + str(hist)

    # -> make it a csv !
    with open(str(simu_name)+'.txt', 'w') as f:
        f.write(csv_header + "\n" + "\n".join(slim_output.splitlines()[17:]))


    # Slim recapitation
    # Passer en fonction et boucler dessus
    if isinstance(generations_vec, int):
        ts = pyslim.load(simu_name + ".trees")
        ts = ts.recapitate(recombination_rate = mu, Ne = popSize)
        ts = pyslim.SlimTreeSequence(msprime.mutate(ts, mu))

        halfsampleSize = int(sampleSize/2)
        sts=simplification(halfsampleSize)
        while(sts.num_individuals != halfsampleSize):
            sts=simplification(halfsampleSize)
        ts.dump(simu_name+"_full.trees")
        sts.dump(simu_name+".trees")
        with open(simu_name+".vcf", "w") as vcf_file: 
            sts.write_vcf(vcf_file, position_transform="legacy")

        data = pd.read_csv(simu_name+".vcf", sep="\t", header=5)
        data.sort_values("POS", inplace = True)
        data.drop_duplicates(subset="POS", keep='first', inplace=True)
        data2=data[data.POS < 100000]
        print(data2)

        with open(simu_name+".vcf","r") as fi:
            header = []
            for ln in fi:
                if ln.startswith("##"):
                    header.append(ln)

        with open(simu_name+".vcf", "w") as fillout:
            for ligne in header:
                fillout.write(ligne)
        data2.to_csv(simu_name+".vcf", sep="\t", mode='a', index=False)




"""
# Write into the SLiM code the choosen parameters
with open('../simupop.c') as f:
    slim_string = f.read()
    
for key in dict_replace:    
    slim_string = slim_string.replace(str(key), str(dict_replace[key]))


# Adding in the end of the SLiM code the choosen output_trees

for generation in param_dict["generations_vec"]:
    slim_string =  slim_string + str(generation) + " late() {sim.treeSeqOutput('./" + str(simu_name) + "-" +  str(generation) + ".trees');}" + "\n"



# Writing the new SLiM code

slim_new_name = simu_name + ".c"

with open(slim_new_name, 'w') as f:
    f.write(slim_string)


# Launching SLiM

slim_output = os.popen("../../build/slim " + slim_new_name).read()

# Writing output in csv format

csv_header = "Generation Alpha Mere_enfant Pere_enfant Parent_enfant" # English
for hist in range(0, int(param_dict["nb_hist"])):
    csv_header = csv_header + " Sibship" + str(hist)

# -> make it a csv !
with open(str(simu_name)+'.txt', 'w') as f:
    f.write(csv_header + "\n" + "\n".join(slim_output.splitlines()[17:]))
"""





