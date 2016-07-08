cdef extern from "vectors_y_pkl_c.c":
	float score_y(short unsigned short[68] v)

def myscore(sv):
	cdef unsigned short[68] v = sv
	return score_y(v)
