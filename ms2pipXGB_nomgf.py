import sys
import numpy as np
import pandas as pd
import ms2pipfeatures_cython
import argparse
import multiprocessing

#Cython compiled modules that contain the predictive models
import vectors_b_pkl
import vectors_y_pkl

def main():

	parser = argparse.ArgumentParser(description='MS2PIP on XGBoost')    
	parser.add_argument('pep_file', metavar='<.PEPREC file>',
					 help='file containing peptide identifications')
	parser.add_argument('-c','--num_cpu', metavar='INT',action="store", dest='num_cpu',default='23',
					 help='number of cores')

	args = parser.parse_args()

	num_cpu = int(args.num_cpu)
		
	#read peptide file into Pandas DataFrame
	data = pd.read_csv(args.pep_file,sep=' ',index_col=False)
	data = data.fillna('-')
	
	#make sure this column is read as a string
	data['spec_id'] = [str(x) for x in data['spec_id']]

	#distribute work over workers
	num_psms_per_cpu = len(data)/num_cpu		
	sys.stderr.write('starting workers...\n')	
	myPool = multiprocessing.Pool(num_cpu)
	results = []
	for i in range(num_cpu):
		#process_file(data.iloc[i*num_psms_per_cpu:(i+1)*num_psms_per_cpu])
		results.append(myPool.apply_async(process_file,args=(data.iloc[i*num_psms_per_cpu:(i+1)*num_psms_per_cpu])))
	myPool.close()
	myPool.join()
	
	#workers done...merging results
	sys.stderr.write('\nmerging results and writing files...\n')
	for r in results:
		(r_pep,r_mz,r_peaks) = r.get()
		for (pep,mz,peaks) in zip(r_pep,r_mz,r_peaks):
			print "%s %s %i" % (pep[0],pep[1],pep[2])
			tmp = np.argsort(mz)
			for j in tmp:
				print "%f %f" % (mz[j],peaks[j])
		
#main function that predicts the MS2 peak intensities for the peptides in 'data'
def process_file(data):
	
	#a_map converts amio acids to integers
	aminos = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	a_map = {}
	for i,a in enumerate(aminos):
		a_map[a] = i
	#a_masses maps a amino acid integer to mz
	a_masses = [71.037114,103.009185,115.026943,129.042593,147.068414,57.021464,137.058912,113.084064,
				128.094963,113.084064,131.040485,114.042927,97.052764,128.058578,156.101111,87.032028,
				101.047679,99.068414,186.079313,163.063329]
	
	result_pep = [] #stores peptide info
	result_mz = [] #stores fragment ion masses
	result_peaks = [] #stores predicted peak intensities
	
	for i,r in data.iterrows():
		peptide = r['peptide']
		mods = r['modifications']
		charge = int(r['charge'])
		result_pep.append([peptide,mods,charge])

		#prepare modification numpy array for indexing
		modified = np.zeros(len(peptide),dtype=np.uint16)
		if mods != '-':
			l = mods.split('|')
			for i in range(0,len(l),2):
				modified[int(l[i])] = float(l[i+1])

		#prepare peptide numpy array for indexing
		peptide_i = np.array([a_map[x] for x in peptide],dtype=np.uint16)
		
		#compute fragment ion masses and store
		tmp_mz = []
		mz = 0
		for j in range(len(peptide_i)-1):
			mz += a_masses[peptide_i[j]] + modified[j]
			tmp_mz.append(mz+1.007324)
		mz = 18.010526
		tmpp = peptide_i[::-1]
		tmpm = modified[::-1]
		for j in range(len(tmpp)-1):
			mz += a_masses[tmpp[j]] + tmpm[j]
			tmp_mz.append(mz+1.007324)
		result_mz.append(tmp_mz)
		
		#compute feature vectors, predict and store the MS2 peak intensities
		X = []
		v = ms2pipfeatures_cython.compute_vector_nopos(peptide_i,modified,charge,0)
		for j,vv in enumerate(v):
			X.append(vectors_b_pkl.myscore(vv))
		v = ms2pipfeatures_cython.compute_vector_nopos(peptide_i,modified,charge,1)
		for j,vv in enumerate(v):
			X.append(vectors_y_pkl.myscore(vv))
		result_peaks.append(X)
		
	return (result_pep,result_mz,result_peaks)


def print_logo():
	logo = """
 _____ _____ ___ _____ _____ _____ 
|     |   __|_  |  _  |     |  _  |
| | | |__   |  _|   __|-   -|   __|
|_|_|_|_____|___|__|  |_____|__|   
                                   
           """
	print logo
	print "Version 2.1\nby Sven Degroeve\n"

if __name__ == "__main__":
	print_logo()
	main()        

