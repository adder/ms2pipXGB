from cpython cimport array
import array
import numpy as np
cimport numpy as np

cdef double* masses = [71.037114,103.009185,115.026943,129.042593,147.068414,57.021464,137.058912,113.084064,128.094963,113.084064,131.040485,114.042927,97.052764,128.058578,156.101111,87.032028,101.047679,99.068414,186.079313,163.063329]

def find_ion_peaks(speptide,smodified,smsms,speaks):
	cdef unsigned short[:] peptide = speptide
	cdef double[:] modified = smodified
	cdef double[:] msms = smsms
	cdef unsigned short[:] peaks = speaks

	cdef int ip
	
	cdef int msms_pos = 0
	cdef int mzlookup_pos = 0
	cdef double mz
	cdef double tol = 0.8
	cdef double tolmz = 0.8
	cdef unsigned short max,tmp2
	cdef int tmp
	
	cdef int lenmsms = len(msms)
	cdef int lenmzlookup_sorted

	cdef int peplen = len(peptide)
	cdef int num_ions = (peplen-1)*2

	cdef np.ndarray[np.float32_t,ndim=1] mzlookup = np.empty(num_ions,dtype=np.float32)
	cdef np.ndarray[np.float32_t,ndim=1] mzlookup_sorted = np.empty(num_ions,dtype=np.float32)
	cdef np.ndarray[np.uint16_t,ndim=1] result = np.zeros(num_ions,dtype=np.uint16)
	cdef np.ndarray[np.int_t,ndim=1] mzlookup_sorted_index

	cdef double mass

	cdef int mpos = 0
	mass = 0.
	for ip from 0 <= ip < peplen-1:
		mass += masses[peptide[ip]] + modified[ip+1]
		mzlookup[mpos] = mass+1.
		mpos += 1
	#mass = 0.
	#for ip from 0 <= ip < peplen-1:
	#	mass += masses[peptide[ip]] + modified[ip+1]
	#	mzlookup[mpos] = (mass+2.)/2.
	#	mpos += 1
	peptide = peptide[::-1]
	mass = 18.010526
	for ip from 0 <= ip < peplen-1:
		mass += masses[peptide[ip]] + modified[ip+1]
		mzlookup[mpos] = (mass+1.)
		mpos += 1
	#mass = 18.010526
	#for ip from 0 <= ip < peplen-1:
	#	mass += masses[peptide[ip]] + modified[ip+1]
	#	mzlookup[mpos] = (mass+2.)/2.
	#	mpos += 1

	mzlookup_sorted_index = np.argsort(mzlookup)
	lenmzlookup_sorted = mzlookup_sorted_index.shape[0]
	for il from 0 <= il < lenmzlookup_sorted:
		mzlookup_sorted[il] = mzlookup[mzlookup_sorted_index[il]]		
	msms_pos = 0
	mzlookup_pos = 0

	data=[]
	indices = []
	
	while True:
		if msms_pos == lenmsms: break
		if mzlookup_pos == lenmzlookup_sorted: break
		mz = mzlookup_sorted[mzlookup_pos]
		#tolmz = (tol*mz)/1000000.
		if msms[msms_pos] > (mz+tolmz):
			mzlookup_pos += 1
		elif msms[msms_pos] < (mz-tolmz):
			msms_pos += 1
		else:
			max = peaks[msms_pos]
			tmp = msms_pos + 1
			if tmp < lenmsms:
				while (msms[tmp] <= (mz+tolmz)):
					tmp2 = peaks[tmp]
					if max < tmp2:
						max = tmp2
					tmp += 1
					if tmp == lenmsms: 
						break
			#data.append(max)
			#indices.append(mzlookup_sorted_index[mzlookup_pos])
			result[mzlookup_sorted_index[mzlookup_pos]] = max
			#result[mzlookup_sorted_index[mzlookup_pos]] = int(100*msms[msms_pos])
			mzlookup_pos += 1
		
	return result
