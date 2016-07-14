import sys
import numpy as np
cimport numpy as np

#amino acid masses as integers
masses = np.array([71,103,115,129,147,57,137,113,128,113,131,114,97,128,156,87,101,99,186,163],dtype=np.uint16)
#amino acid chemical properties as positive integers
chem = np.array([[10, 23, 10, 17, 37, 27, 0, 61, 23, 55, 20, 30, 29, 34, 33, 100, 14, 26, 17, 39, 21, 30, 35],[51, 18, 75, 25, 35, 100, 16, 3, 94, 0, 97, 82, 12, 0, 22, 22, 21, 39, 80, 98, 95, 70, 28],[93, 49, 31, 45, 39, 95, 79, 56, 100, 43, 98, 90, 52, 0, 54, 53, 60, 72, 97, 69, 100, 75, 47],[40, 100, 28, 0, 5, 33, 40, 60, 40, 87, 40, 37, 33, 44, 36, 100, 36, 35, 39, 39, 40, 36, 40]],dtype=np.uint16)

#computes feature vectors for all nterm-ions
#
# peptide: peptide sequence as numpy np.uint16 array
# modifications: masses to add as numpy np.uint16 array
# charge: spectrum charge state
def compute_vectors(peptide, modifications, int charge):
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

	for i from 0 <= i < peplen:
		mz = mz + modis[i]
		
	vector = []
		
	mzb = 0
	for i from 0 <= i < peplen-1:
		buf = np.zeros(164,dtype=np.uint16)
		buf[0] = peplen
		buf[1] = i
		buf[2] = int(float(100*i)/peplen)
		buf[3] = mz
		for j from 0 <= j < 4:
			buf[4+j] = int(chem_total[j]/peplen)
		mzb = mzb + masses[pep[i]] + modis[i]
		buf[8] = mzb
		buf[9] = mz-mzb
		buf[10] = charge
		buf[11] = modis[i]
		for j from 0 <= j < 4:
			buf[12+j] = chem[j][pep[0]]
		for j from 0 <= j < 4:
			buf[16+j] = chem[j][pep[1]]
		for j from 0 <= j < 4:
			buf[20+j] = chem[j][pep[-2]]
		for j from 0 <= j < 4:
			buf[24+j] = chem[j][pep[-1]]
		for j from 0 <= j < 4:
			buf[28+j] = chem[j][pep[i]]
		if (i-1) < 0:
			for j from 0 <= j < 4:
				buf[32+j] = chem[j][pep[0]]
		else:
			for j from 0 <= j < 4:
				buf[32+j] = chem[j][pep[i-1]]
		for j from 0 <= j < 4:
			buf[36+j] = chem[j][pep[i+1]]
		if (i+2) >= peplen:
			for j from 0 <= j < 4:
				buf[40+j] = chem[j][pep[peplen-1]]
		else:
			for j from 0 <= j < 4:
				buf[40+j] = chem[j][pep[i+2]]

		for j from 0 <= j < 4:
			s = 0
			maxs = 0
			mins = 9999			
			for k from 0 <= k <= i:
				tmp = chem[j][pep[k]]
				s = s + tmp
				if tmp > maxs:
					maxs = tmp
				if tmp < mins:
					mins = tmp
			buf[44+(j*8)] = s
			buf[45+(j*8)] = int(s/(i+1))
			buf[46+(j*8)] = maxs
			buf[47+(j*8)] = mins
			
			s = 0
			maxs = 0
			mins = 9999			
			for k from (i+1) <= k < peplen:
				tmp = chem[j][pep[k]]
				s = s + tmp
				if tmp > maxs:
					maxs = tmp
				if tmp < mins:
					mins = tmp
			buf[48+(j*8)] = s
			buf[49+(j*8)] = int(s/(peplen-(i+1)))
			buf[50+(j*8)] = maxs
			buf[51+(j*8)] = mins

		s = 0
		maxs = 0
		mins = 9999			
		for k from 0 <= k <= i:
			tmp = masses[pep[k]]
			s = s + tmp
			if tmp > maxs:
				maxs = tmp
			if tmp < mins:
				mins = tmp
		buf[76] = s
		buf[77] = int(s/(i+1))
		buf[78] = maxs
		buf[79] = mins
		
		s = 0
		maxs = 0
		mins = 9999			
		for k from (i+1) <= k < peplen:
			tmp = masses[pep[k]]
			s = s + tmp
			if tmp > maxs:
				maxs = tmp
			if tmp < mins:
				mins = tmp
		buf[80] = s
		buf[81] = int(s/(peplen-(i+1)))
		buf[82] = maxs
		buf[83] = mins

		buf[84+pep[0]] = 1
		buf[104+pep[-1]] = 1
		buf[124+pep[i]] = 1
		buf[144+pep[i+1]] = 1
				
		vector.append(buf)

	return np.array(vector,dtype=np.uint16)
	
