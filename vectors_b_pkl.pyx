cdef extern from "vectors_b_pkl_c.c":
	float score_b(short unsigned short[164] v)

def myscore(sv):
	cdef unsigned short[164] v = sv
	return score_b(v)
