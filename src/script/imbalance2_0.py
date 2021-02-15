import tskit, msprime, os, sys, argparse
import pandas as pd
import numpy as np

def get_arguments():
    """
    Lecture des arguments
    """
    parser = argparse.ArgumentParser(description='')
    #Mandatory arguments :
    parser.add_argument('-p', required=True, help="Path to the working directory")
    
    #Optionnal arguments :
    parser.add_argument('-t', type=int, default=4, help="Number of threads. Default: 4")
    parser.add_argument('-i', type=int, default=1000, help="Number of simulation")
    parser.add_argument('-l', type=float, default=1e5, help="Chromosome length")
    parser.add_argument('--sim_name', type=str, help="Generique name of the simulation")
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
	args = get_arguments()
	dirName = args.p
	par = args.t
	nb_sim = args.i
	k_len = args.l
	#Filename of tools used for simulate data
	#Either msprime or slim
	data_sim = args.sim_name
	wdpath = sys.argv[0].replace("imbalance2_0.py", "")
	os.environ['curr'] = wdpath



	# conversion en Newick pour index Imbalance

	os.environ['path'] = dirName
	os.environ['par'] = str(par)
	os.environ['data_sim'] = data_sim

	os.system("bash -c 'rm launch_conv.txt'")
	for i in range(1,nb_sim+1):
		os.environ['i']=str(i)
		os.system("echo python ${curr}newick_conv.py $i $path $data_sim >> launch_conv.txt")

	os.system("bash -c 'parallel -a launch_conv.txt -j $par'")
	os.system("bash -c 'parallel -a    launchms.txt -j $par'")


	# conversion en format adapté au programme de JT pour l'index Imbalance

	os.system ("bash -c 'rm launchms.txt'")
	os.system ("bash -c 'rm launchts.txt'")
	os.system ("bash -c 'rm launchre.txt'")
	os.system ("bash -c 'rm launchms2.txt'")
	for i in range(1,nb_sim+1):
		dataNewick = str(dirName) + "/" + str(data_sim) + str(i) + ".newick"
		tsNewick = str(dirName) + "/ts" + str(i) + ".newick"
		reNewick = str(dirName) + "/relate" + str(i) + ".newick"
		msNewick2 = str(dirName) + "/msprime" + str(i) + "_2.newick"
		print(dataNewick)
		os.environ['dataNewick'] = dataNewick
		os.environ['tsNewick'] = tsNewick
		os.environ['reNewick'] = reNewick
		os.environ['msNewick2'] = msNewick2
		os.system("echo python ${curr}format_netwick.py $dataNewick $dataNewick'.out' >> launchms.txt")
		os.system("echo python ${curr}format_netwick.py $tsNewick $tsNewick'.out' >> launchts.txt")
		os.system("echo python ${curr}format_netwick.py $reNewick $reNewick'.out' >> launchre.txt")
		os.system("echo python ${curr}format_netwick.py $msNewick2 $msNewick2'.out' >> launchms2.txt")
		#os.system("bash -c './format_netwick.py '$msNewick' '$msNewick'.out' >> lauchms.txt")
		#os.system("bash -c './format_netwick.py '$tsNewick' '$tsNewick'.out' >> lauchts.txt")
		#os.system("bash -c './format_netwick.py '$reNewick' '$reNewick'.out' >> lauchre.txt")


	os.system("bash -c 'parallel -a launchms.txt -j $par'")
	os.system("bash -c 'parallel -a launchts.txt -j $par'")
	os.system("bash -c 'parallel -a launchre.txt -j $par'")
	os.system("bash -c 'parallel -a launchms2.txt -j $par'")
				
	# calcul index Imbalance

	os.system ("bash -c 'rm launch2ms.txt'")
	os.system ("bash -c 'rm launch2ts.txt'")
	os.system ("bash -c 'rm launch2re.txt'")
	os.system ("bash -c 'rm launch2ms2.txt'")

	for i in range(1,nb_sim+1):
		dataNewick = str(dirName) + "/" + str(data_sim) + str(i) + ".newick"
		tsNewick = str(dirName) + "/ts" + str(i) + ".newick"
		reNewick = str(dirName) + "/relate" + str(i) + ".newick"
		msNewick2 = str(dirName) + "/msprime" + str(i) + "_2.newick"
		os.environ['dataNewick'] = dataNewick
		os.environ['tsNewick'] = tsNewick
		os.environ['reNewick'] = reNewick
		os.environ['msNewick2'] = msNewick2
		os.system("echo './${curr}CalcIPrimeV2_linux64.exe -f '$dataNewick'.out > '$dataNewick'.results' >> launch2ms.txt")
		os.system("echo './${curr}CalcIPrimeV2_linux64.exe -f '$tsNewick'.out > '$tsNewick'.results' >> launch2ts.txt")
		os.system("echo './${curr}CalcIPrimeV2_linux64.exe -f '$reNewick'.out > '$reNewick'.results' >> launch2re.txt")
		os.system("echo './${curr}CalcIPrimeV2_linux64.exe -f '$msNewick2'.out > '$msNewick2'.results' >> launch2ms2.txt")
		#os.system("bash -c './script/CalclPrimeV2_linux64.exe -f '$msNewick'.out' > $msNewick'.results' >>launch2ms.txt")
		#os.system("bash -c './script/CalclPrimeV2_linux64.exe -f '$tsNewick'.out' > $tsNewick'.results' >>launch2ts.txt")
		#os.system("bash -c './script/CalclPrimeV2_linux64.exe -f '$reNewick'.out' > $reNewick'.results' >>launch2re.txt")
	os.environ['par'] = str(par)
	os.system("bash -c 'parallel -a launch2ms.txt -j $par'")
	os.system("bash -c 'parallel -a launch2ts.txt -j $par'")
	os.system("bash -c 'parallel -a launch2re.txt -j $par'")
	os.system("bash -c 'parallel -a launch2ms2.txt -j $par'")

		

	def calc_mean_var_imb(file1, k_len):
		recouvList=[]
		meanList=[]
		Res=file1+".newick.results"
		Tfile=file1+".trees"
		tsfile=tskit.load(Tfile)
		data=pd.read_table(Res)
		imbList=data['IprimeNonBin2N']

		for tree in tsfile.trees():
			recouvList.append(tree.interval[1]-tree.interval[0])
		imbList=np.array(imbList)
		recouvList=np.array(recouvList)
		l3=recouvList*imbList
		mean=sum(l3)/k_len
		meanList = np.repeat(mean, len(imbList))

		imbListbis = imbList-meanList
		imbListbis = np.square(imbListbis)
		varList = imbListbis*recouvList
		var = sum(varList)/k_len

		return mean,var



	msimbMeanList = []
	msimbVarList = []
	tsimbMeanList = []
	tsimbVarList = []
	reimbMeanList = []
	reimbVarList = []
	outimbMeanList = []
	outimbVarList = []
	for i in range(1,nb_sim+1):
		msRes = str(dirName) + "/" + str(data_sim) + str(i)
		tsRes = str(dirName) + "/ts" + str(i)
		reRes = str(dirName) + "/relate" + str(i)
		msRes2 = str(dirName) + "/msprime" + str(i) +"_2"
		mstmp = calc_mean_var_imb(msRes, k_len)
		msimbMeanList.append(mstmp[0])
		msimbVarList.append(mstmp[1])
		tstmp = calc_mean_var_imb(tsRes, k_len)
		tsimbMeanList.append(tstmp[0])
		tsimbVarList.append(tstmp[1])
		retmp = calc_mean_var_imb(reRes, k_len)
		reimbMeanList.append(retmp[0])
		reimbVarList.append(retmp[1])
		outtmp = calc_mean_var_imb(msRes2, k_len)
		outimbMeanList.append(outtmp[0])
		outimbVarList.append(outtmp[1])


	JTmsimbMeanList = []
	JTtsimbMeanList = []
	JTreimbMeanList = []
	JToutimbMeanList = []
	for i in range(1,nb_sim+1):
		msRes = str(dirName) + "/" + str(data_sim) + str(i) + ".newick.results"
		tsRes = str(dirName) + "/ts" + str(i) + ".newick.results"
		reRes = str(dirName) + "/relate" + str(i) + ".newick.results"
		msRes2 = str(dirName) + "/msprime" + str(i) + "_2.newick.results"
		print(msRes2)
		msdata = pd.read_table(msRes)
		tsdata = pd.read_table(tsRes)
		redata = pd.read_table(reRes)
		outdata = pd.read_table(msRes2)

		JTmsimbMeanList.append(np.mean(msdata['IprimeNonBin2N']))
		#msimbVarList.append(np.var(msdata['IprimeNonBin2N']))
		JTtsimbMeanList.append(np.mean(tsdata['IprimeNonBin2N']))  
		#tsimbVarList.append(np.var(tsdata['IprimeNonBin2N']))
		JTreimbMeanList.append(np.mean(redata['IprimeNonBin2N']))
		#reimbVarList.append(np.var(redata['IprimeNonBin2N']))
		JToutimbMeanList.append(np.mean(outdata['IprimeNonBin2N']))


	# Creation du data frame
	d = {'msimbMean': msimbMeanList, 'msimbVar': msimbVarList, 'tsimbMean': tsimbMeanList, 'tsimbVar': tsimbVarList, 'reimbMean': reimbMeanList, 'reimbVar': reimbVarList, 'outimbMean': outimbMeanList, 'outimbVarList': outimbVarList, 
	'JTmsimbMeanList': JTmsimbMeanList, 'JTtsimbMeanList': JTtsimbMeanList, 'JTreimbMeanList': JTreimbMeanList, 'JToutimbMeanList': JToutimbMeanList}
	df = pd.DataFrame(data=d)
	nameDF = str(dirName) + "/Imbalance.csv"
	df.to_csv(nameDF, index=False, sep="\t")



