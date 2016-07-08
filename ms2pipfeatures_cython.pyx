import sys
import numpy as np
cimport numpy as np

masses = np.array([71,103,115,129,147,57,137,113,128,113,131,114,97,128,156,87,101,99,186,163],dtype=np.uint16)
chem = np.array([[10, 23, 10, 17, 37, 27, 0, 61, 23, 55, 20, 30, 29, 34, 33, 100, 14, 26, 17, 39, 21, 30, 35],[51, 18, 75, 25, 35, 100, 16, 3, 94, 0, 97, 82, 12, 0, 22, 22, 21, 39, 80, 98, 95, 70, 28],[93, 49, 31, 45, 39, 95, 79, 56, 100, 43, 98, 90, 52, 0, 54, 53, 60, 72, 97, 69, 100, 75, 47],[40, 100, 28, 0, 5, 33, 40, 60, 40, 87, 40, 37, 33, 44, 36, 100, 36, 35, 39, 39, 40, 36, 40]],dtype=np.uint16)

def compute_vector_nopos(peptide, modifications, int charge, int ionpart):
	cdef int i,j,k
	cdef unsigned short tmp,s,maxs1,mins1,maxs2,mins2
	cdef unsigned short mzb
	
	cdef unsigned short[:] pep = peptide
	cdef unsigned short[:] modis = modifications

	cdef int peplen = len(pep)
	
	cdef unsigned short mz = 0
	
	cdef np.ndarray[np.uint16_t,ndim=1] chem_total = np.array([0,0,0,0],dtype=np.uint16)
	for i from 0 <= i < peplen:
		mz = mz + masses[pep[i]]
		for j from 0 <= j < 4:
			chem_total[j] = chem_total[j] + chem[j][pep[i]]

	for i from 0 <= i < peplen+2:
		mz = mz + modis[i]
		
	vector = []
		
	mzb = 0
	for i from 0 <= i < peplen-1:
		vector.append([])
		vector[i].append(peplen)
		if ionpart == 0:
			vector[i].append(i)
			vector[i].append(int(float(100*i)/peplen))
		else:
			vector[i].append(peplen-1-i)
			vector[i].append(int(float(peplen-1-i)/peplen))
		vector[i].append(mz)
		for j from 0 <= j < 4:
			vector[i].append(chem_total[j]/peplen)
		mzb = mzb + masses[pep[i]] + modis[i]
		if ionpart == 0:
			vector[i].append(mz)
			vector[i].append(mz-mzb)		
		else:
			vector[i].append(mz-mzb)		
			vector[i].append(mz)
		vector[i].append(charge)
		vector[i].append(ionpart)
		for k in [0,1,-2,-1]:
			#for j from 0 <= j < 20:
			#	if j == pep[k]:
			#		vector[i].append(1)
			#	else:
			#		vector[i].append(0)
			for j from 0 <= j < 4:
					vector[i].append(chem[j][pep[k]])
			#vector[i].append(modis[k])
		for k in [i,i-1,i+1,i+2]:
			if k < 0: k = 0
			if k >= peplen: k = peplen - 1
			#for j from 0 <= j < 20:
			#	if j == pep[k]:
			#		vector[i].append(1)
			#	else:
			#		vector[i].append(0)
			for j from 0 <= j < 4:
					vector[i].append(chem[j][pep[k]])
			#vector[i].append(modis[k])
		for j from 0 <= j < 4:
			s = 0
			maxs1 = 0
			mins1 = 9999			
			for k from 0 <= k < i:
				tmp = chem[j][pep[k]]
				s = s + tmp
				if tmp > maxs1:
					maxs1 = tmp
				if tmp < mins1:
					mins1 = tmp
			s=0
			maxs2 = 0
			mins2 = 9999			
			for k from i <= k < peplen:
				tmp = chem[j][pep[k]]
				s = s + tmp
				if tmp > maxs2:
					maxs2 = tmp
				if tmp < mins2:
					mins2 = tmp
			if ionpart == 0:			
				vector[i].append(s/(i+1))
				vector[i].append((chem_total[j]-s)/(peplen-i))
				vector[i].append(maxs1)
				vector[i].append(mins1)			
				vector[i].append(maxs2)
				vector[i].append(mins2)			
			else:
				vector[i].append((chem_total[j]-s)/(peplen-i))
				vector[i].append(s/(i+1))
				vector[i].append(maxs2)
				vector[i].append(mins2)			
				vector[i].append(maxs1)
				vector[i].append(mins1)			
	

	return np.array(vector,dtype=np.uint16)
	
