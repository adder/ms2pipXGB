import sys
import numpy as np
import pandas as pd
import utils_cython as uc
import ms2pipfeatures_cython
import pickle
import argparse
import multiprocessing

from scipy.stats import pearsonr

import vectors_b_pkl
import vectors_y_pkl

def main():

	parser = argparse.ArgumentParser(description='MS2PIP on XGBoost')    
	parser.add_argument('pep_file', metavar='<.PEPREC file>',
					 help='file containing peptide identifications')
	parser.add_argument('spec_file', metavar='<.PEPREC.mgf file>',
					 help='file containing ms2 spectra')
	parser.add_argument('-c','--num_cpu', metavar='INT',action="store", dest='num_cpu',default='23',
					 help='number of cores')

	args = parser.parse_args()

	num_cpu = int(args.num_cpu)

	#ptm_map = {}
	#ptm_map['CAM'] = 57.02146
	#ptm_map['Oxidation'] = 15.994915
	
	#load tables with chemical properties of amino acids
	chemdata = pickle.load(open("amino_acid_properties.pkl","rb" ))

	#a_map converts the peptide amio acids to integers
	aminos = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	a_map = {}
	for i,a in enumerate(aminos):
		a_map[a] = i
		
	data = pd.read_csv(args.pep_file,sep=' ',index_col=False)
	data['peplen'] = data['peptide'].apply(return_length)
	data = data.fillna('-')
	#make sure this column is a string
	data['spec_id'] = [str(x) for x in list(data['spec_id'])]
		
	sys.stderr.write('scanning spectrum file...\n')
	titles = scan_spectrum_file(args.spec_file)
	num_spectra = len(titles)
	num_spectra_per_cpu = int(num_spectra/(num_cpu-1))
	sys.stderr.write("%i spectra (%i per cpu)\n"%(num_spectra,num_spectra_per_cpu))
						
	sys.stderr.write('starting workers...\n')
	
	myPool = multiprocessing.Pool(num_cpu)
	
	results = []
	i = 0
	for i in range(num_cpu-1):
		tmp = titles[i*num_spectra_per_cpu:(i+1)*num_spectra_per_cpu]
		"""
		process_file(
										i,
										args.spec_file,
										tmp,
										data[data.spec_id.isin(tmp)],
										a_map
										)
		"""
		results.append(myPool.apply_async(process_file,args=(
										i,
										args.spec_file,
										tmp,
										data[data.spec_id.isin(tmp)],
										a_map
										)))		
	i+=1
	tmp = titles[i*num_spectra_per_cpu:]
	results.append(myPool.apply_async(process_file,args=(
									i,
									args.spec_file,
									tmp,
									data[data.spec_id.isin(tmp)],
									a_map
									)))		
	myPool.close()
	myPool.join()
	
	#workers done...merging results
	sys.stderr.write('\nmerging results and writing files...\n')
	p_b = []
	p_y = []
	for r in results:
		(sp_b,sp_y) = r.get()
		p_b.extend(sp_b)
		p_y.extend(sp_y)
	
	tmp = pd.DataFrame()
	tmp['b'] = p_b
	tmp['y'] = p_y
	tmp.to_csv("pearson.csv",index=False)
		
def process_file(pid,spec_file,titles,data,a_map):
	
	specdict = data[['spec_id','peptide','modifications']].set_index('spec_id').to_dict()
	peptides = specdict['peptide']
	modifications = specdict['modifications']

	p_b = []
	p_y = []

	title = ""
	parent_mz = 0.
	charge = 0
	msms = []
	peaks = []
	f = open(spec_file)
	skip = False
	vectors_b = []
	vectors_y = []
	targets_b_1 = []
	#targets_b_2 = []
	targets_y_1 = []
	#targets_y_2 = []
	psmids = []
	psmids_count = 0
	while (1):
		rows = f.readlines(2000000)
		sys.stderr.write('.')
		if not rows: break
		for row in rows:
			row = row.rstrip()
			if skip:
				if row[:10] == "BEGIN IONS":
					skip = False
				else:
					continue
			if row == "": continue
			if row[:5] == "TITLE":
				title = row[6:].replace(' ','')
				if not title in peptides:
					skip = True
					continue
			elif row[0].isdigit():
				tmp = row.split()
				msms.append(float(tmp[0]))
				peaks.append(float(tmp[1]))
			elif row[:10] == "BEGIN IONS":
				msms = []
				peaks = []
			elif row[:6] == "CHARGE":
				charge = int(row[7:9].replace("+",""))
			elif row[:7] == "PEPMASS":
				parent_mz = float(row[8:].split()[0])
			elif row[:8] == "END IONS":
				#process
				if not title in peptides: continue
				peptide = peptides[title]
				mods = modifications[title]
				peplen = len(peptide)
				modified = [0.]*(peplen+2)
				#k = False
				#if mods != '-':
				#	l = mods.split('|')
				#	for i in range(0,len(l),2):
				#		if not l[i+1] in ptm_map:
				#			print "PTM not found: %s" % (l[i+1])
				#			k = True
				#			break
				#		modified[int(l[i])] = ptm_map[l[i+1]]
				#if k: 
				#	continue
				modified = np.array(modified)
				peptide = np.array([a_map[x] for x in peptide],dtype=np.uint16) #change datatype?
				msms = np.array(msms)
				peaks = np.array(peaks)
				peaks = peaks / np.sum(peaks)
				peaks *= 64000
				peaks = peaks.astype(np.uint16)
				targets = np.array(uc.find_ion_peaks(peptide,modified,msms,peaks),dtype=np.uint16)
				l = len(targets)/2
				targets_b_1.append(targets[:l])
				#targets_b_2.append(targets[l:2*l])
				targets_y_1.append(targets[l:])
				#targets_y_2.append(targets[3*l:])
				modified = modified.astype(np.uint16)
				#vectors_b.append(ms2pipfeatures_cython.compute_vector_nopos(peptide,modified,charge,0))
				#vectors_y.append(ms2pipfeatures_cython.compute_vector_nopos(peptide,modified,charge,1))
				v = ms2pipfeatures_cython.compute_vector_nopos(peptide,modified,charge,0)
				tmpp = []
				for j,vv in enumerate(v):
					tmpp.append(vectors_b_pkl.myscore(vv))
				p_b.append(pearsonr(tmpp,targets[:l])[0])
				v = ms2pipfeatures_cython.compute_vector_nopos(peptide,modified,charge,1)
				tmpp = []
				for j,vv in enumerate(v):
					tmpp.append(vectors_y_pkl.myscore(vv))
				p_y.append(pearsonr(tmpp,targets[l:])[0])

	return (p_b,p_y)

def get_feature_names():
	names = ['peplen','ionnumber','ionnumber_rel','pmz']
	for c in ['bas','hydro','heli','pI']:
		names.append('mean_'+c)
	names.append("mz_ion")
	names.append("mz_ion_other")
	names.append("charge")
	names.append("iontype")
	for k in ['0','1','-2','-1']:
		#for a in aminos:
		#	names.append(a+'_'+k)
		for c in ['bas','hydro','heli','pI']:
			names.append(c+'_'+k)
		#names.append('mod_'+k)
	for k in ['i','i-1','i+1','i+2']:
		#for a in aminos:
		#	names.append(a+'_'+k)	
		for c in ['bas','hydro','heli','pI']:
			names.append(c+'_'+k)
		#names.append('mod_'+k)
	for c in ['bas','hydro','heli','pI']:
		names.append('mean_ion_'+c)
		names.append('mean_ion_other_'+c)
		names.append('max_ion_'+c)
		names.append('max_ion_other_'+c)
		names.append('min_ion_'+c)
		names.append('min_ion_other_'+c)
	return names

def scan_spectrum_file(filename):
	titles = []
	f = open(filename)
	while (1):
		rows = f.readlines(100000)
		if not rows: break
		for row in rows:
			if row[:5] == "TITLE":
				titles.append(row.rstrip()[6:])
	f.close()
	return titles

def return_length(x):
	return len(x)	

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

