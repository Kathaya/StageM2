import argparse
import os
import time
import shutil
"""
pop=$1 #Population size
nb_sim=$2 #Number of simulation
nb_ech=$3 #Number of echantillion : population size/2
threads=$4 #Number of threads for computation
dirName=$5 #Main diretory name
variable=$6 #secondary directory name
rho=$7 #Recombinasion rate
mu=$8 #Muation rate
gen_length=$9 #Genome length
"""

def get_arguments():
    """
    Lecture des arguments
    """
    parser = argparse.ArgumentParser(description='')

    #Mandatory arguments :
    parser.add_argument('d', help="Path to the working directory")
    parser.add_argument('sd', help="Sub_directory name")

    #Optionnal arguments :
    parser.add_argument('-p', type=int, default=10000, help="Population size")
    parser.add_argument('-i', type=int, default=1000, help="Number of simulation")
    parser.add_argument('-s', type=int, default=30, help="Sample size")
    parser.add_argument('-t', type=int, default=4, help="Number of threads. Default: 4")
    parser.add_argument('-r', type=float, default=2e-8, help="Recombination rate")
    parser.add_argument('-m', type=float, default=1e-8, help="Mutation rate")
    parser.add_argument('-l', type=float, default=1e5, help="Chromosome length")
    #parser.add_argument('-f', type=int, default=1, help="Factor of demegraphic changes")
    #parser.add_argument('-g', type=int, default=5, help="Number of generation before change in demogaphy")
    args = parser.parse_args()
    return args

def computation_time(start_time):
    """
    Time calculation in second between 2 dates
    ===========================
    Input : 
        - start_time : starting time of a task
    ===========================
    Output :
        - Time difference between start_time and current time
    """
    return time.time() - start_time

if __name__ == "__main__":
    args = get_arguments()
    path=args.d+"/"+args.sd

    popSize = args.p 
    nsim = args.i
    sampleSize = args.s
    threads = args.t
    rho = args.r
    mu = args.m
    L = args.l

    start_time=time.time()
    os.environ['threads'] = str(threads)

    #cleaning of directory
    if os.path.exists(path) and os.path.isdir(path):
        shutil.rmtree(path)
    os.mkdir(path)

    os.system("rm -rf launch*")

    #Preparing output file
    with open(path+"/kc_distance.csv", "w") as fillout:
        fillout.write("Simulation\tMean ms/ts\tVar ms/ts\tMean ms/re\tVar ms/re\tMean ts/re\tVar ts/re\tMean ms/out\tVar ms/out\tMean re/out\tVar re/out\tMean ts/out\tVar ts/out\n")
    with open(path+"/ntree.csv", "w") as fillout:
        fillout.write("msprime\ttsinfer\trelate\tout\n")

    #Preparing msprime simulation
    with open("launch_ms.txt", "w") as fillout:
        for i in range(1, args.i+1):
            vcf_path = path+"/simulation_"+str(i)
            tree_path = path+"/msprime"+str(i)

            fillout.write("echo msprime {}\n".format(i))
            fillout.write('python script/msprime_simu.py -p {} -s {} -vp {} -tp {} -rho {} -mu {} -l {}\n'.format(popSize, sampleSize, vcf_path, tree_path, rho, mu, L))

    ms_time = time.time()
    os.system("bash -c 'parallel -a launch_ms.txt -j $threads'")
    ms_time = computation_time(ms_time)
    
    #Preparing outgroup via msprime
    with open("launch_ms_2.txt", "w") as fillout:
        for i in range(1, args.i+1):
            vcf_path = path+"/simulation_"+str(i)+"_2"
            tree_path = path+"/msprime"+str(i)+"_2"
            #fillout.write(python script/msprime_simu.py -p popSize -s sampleSize -vp vcf_path -tp tree_path -rho rho -mu mu)
            fillout.write("echo msprime2 {}\n".format(i))
            fillout.write('python script/msprime_simu.py -p {} -s {} -vp {} -tp {} -rho {} -mu {}\n'.format(popSize, sampleSize, vcf_path, tree_path, rho, mu))
    os.system("bash -c 'parallel -a launch_ms_2.txt -j $threads'")


    #Preparing tsinfer simulation
    with open("launch_ts.txt", "w") as fillout:
        for i in range(1, args.i+1):
            vcf_path = path+"/simulation_"+str(i)+".vcf"
            sample_path = path+"/sample"+str(i)+".samples"
            ts_path = path+"/ts"+str(i)+".trees"

            fillout.write("python script/function_map2.py ts_sim --file1 {} --file2 {} --file3 {}\n".format(vcf_path, sample_path, ts_path))
    
    ts_time = time.time()
    os.system("bash -c 'parallel -a launch_ts.txt -j $threads'")
    ts_time = computation_time(ts_time)
    
    #Preparing file for relate
    with open("launch_pre_relate.txt", "w") as pre, open(path+"/launch_map_relate.txt", "w") as map_re, open(path+"/launch_resim.txt", "w") as re, open(path+"/relate_conversion.txt", "w") as conv:
        for i in range(1, args.i+1):
            haps_path = path+"/relate"+str(i)+".haps"
            sample_path = path+"/relate"+str(i)+".sample"
            simu= path+"/simulation_"+str(i)
            map_path = path+"/relate"+str(i)+".map"
            relate_path = path+"/relate"+str(i)
            relate_curr = "relate"+str(i)

            pre.write("~/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps {} --sample {} -i {}\n".format(haps_path, sample_path, simu))
            map_re.write("python script/function_map2.py map --file1 {} --file2 {}\n".format(haps_path, map_path))
            re.write("~/relate/scripts/RelateParallel/RelateParallel.sh --mode All -m {} -N {} --haps {} --sample {} --map {} -o {} --threads {}\n".format(mu, popSize, haps_path, sample_path, map_path, relate_curr, threads))
            conv.write("~/relate/bin/RelateFileFormats --mode ConvertToTreeSequence -i {} -o {}\n".format(relate_curr, relate_path))

    re_time = time.time()
    os.system("bash -c 'parallel -a launch_pre_relate.txt -j {}'".format(threads))
    os.system("bash -c 'parallel -a {}/launch_map_relate.txt -j {}'".format(path, threads))
    os.system("bash -c 'parallel -a {}/launch_resim.txt -j {}'".format(path, threads))
    os.system("bash -c 'parallel -a {}/relate_conversion.txt -j {}'".format(path, threads))
    re_time = computation_time(re_time)
    
    #Summary stats cuomputation
    with open(path+"/launch_kc.txt", "w") as kc, open(path+"/launch_ntree.txt", "w") as ntree:
        for i in range(1, args.i+1):
            mstree = path+"/msprime"+str(i)+".trees"
            tstree = path+"/ts"+str(i)+".trees"
            retree = path+"/relate"+str(i)+".trees"
            m2tree = path+"/msprime"+str(i)+"_2.trees"

            kc.write("echo {}'\t'$(python script/function_map2.py kc --file1 {} --file2 {})'\t'".format(i, mstree, tstree))
            kc.write("$(python script/function_map2.py kc --file1 {} --file2 {})'\t'".format(mstree,retree))
            kc.write("$(python script/function_map2.py kc --file1 {} --file2 {})'\t'".format(tstree,retree))
            kc.write("$(python script/function_map2.py kc --file1 {} --file2 {})'\t'".format(mstree,m2tree))
            kc.write("$(python script/function_map2.py kc --file1 {} --file2 {})'\t'".format(retree,m2tree))
            kc.write("$(python script/function_map2.py kc --file1 {} --file2 {})\n".format(tstree,m2tree))
            #kc.write("\n")

            ntree.write("echo $(python script/function_map2.py ntree --file1 {})'\t'".format(mstree))
            ntree.write(str("$(python script/function_map2.py ntree --file1 {})'\t'".format(tstree)))
            ntree.write(str("$(python script/function_map2.py ntree --file1 {})'\t'".format(retree)))
            ntree.write(str("$(python script/function_map2.py ntree --file1 {})\n".format(m2tree)))

    ss_time = time.time()
    os.system("bash -c 'parallel -a {} -j {}' >> {}".format(path+"/launch_kc.txt", threads, path+"/kc_distance.csv"))
    os.system("bash -c 'parallel -a {} -j {}' >> {}".format(path+"/launch_ntree.txt", threads, path+"/ntree.csv"))
    os.system("python script/imbalance2.py {} -i {} -t {} -l {}".format(path, nsim, threads, int(L)))
    ss_time = computation_time(ss_time)

    #Time calculation
    start_time = computation_time(start_time)
    #Clearing tempory files
    os.system("rm -rf relate*")

    with open(path+"/time.txt", "w") as fillout:
        fillout.write("Total = {}s\nmsprime = {}s\ntsinfer = {}s\nRelate = {}s\nSS = {}s\n".format(start_time, ms_time, ts_time, re_time, ss_time))
    