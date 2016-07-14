import os
import sys
import argparse
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
import random
import operator
import pickle
import matplotlib.pyplot as plt

def main():

	parser = argparse.ArgumentParser(description='XGBoost training')    
	parser.add_argument('vectors',metavar='<vectors.pkl>',
					 help='feature vector file')
	parser.add_argument('targets',metavar='<targets.pkl>',
					 help='target file')
	parser.add_argument('psmids',metavar='<psmids.pkl>',
					 help='PSM groups file')
	parser.add_argument('type',metavar='<model type>',
					 help='{b,y}')
	parser.add_argument('-p',metavar='INT', action="store", dest='num_cpu', default=23, 
					 help='number of cpu\'s to use')
	args = parser.parse_args()
	
	sys.stderr.write('loading data\n')
	
	vectors = pd.read_pickle(args.vectors)
	targets = pd.read_pickle(args.targets)
	psmids = pd.read_pickle(args.psmids)
	#psmids = psmids.PSMid
	
	
	#vectors = vectors[targets.target>0]
	#psmids = psmids[targets.target>0].PSMid
	#targets = targets[targets.target>0]

	psmids = psmids.PSMid

	
	targets = np.log2(targets+10)
	
	#selecting charge +2 peptides only!!
	np.random.seed(1)
	upeps = psmids[vectors.charge==2].unique()
	num_psms = len(upeps)
	np.random.shuffle(upeps)

	#creating train/test split
	test_psms = upeps[:20000]
	test_vectors = vectors[psmids.isin(test_psms)]
	test_targets = targets[psmids.isin(test_psms)]
	train_vectors = vectors[~psmids.isin(test_psms)]
	train_targets = targets[~psmids.isin(test_psms)]
	
	sys.stderr.write('loading data done\n')
		
	#rename features to understand decision tree dump
	train_vectors.columns = ['Feature'+str(i) for i in range(len(train_vectors.columns))]
	test_vectors.columns = ['Feature'+str(i) for i in range(len(train_vectors.columns))]
	numf = len(train_vectors.columns.values)
	
	#create XGBoost datastructure
	sys.stderr.write('creating DMatrix\n')
	xtrain = xgb.DMatrix(train_vectors, label=train_targets)
	xeval = xgb.DMatrix(test_vectors, label=test_targets)
	sys.stderr.write('creating DMatrix done\n')
	
	evallist  = [(xeval,'eval')]
	
	#set XGBoost parameters; make sure to tune well!
	param = {"objective":"reg:linear",
	         "nthread":int(args.num_cpu),
	         "silent":1,
	         "eta":0.3,
	         #"max_delta_step":1,
	         "max_depth":10,
			 "gamma":0.1,
			 "min_child_weight":10,
			 "subsample":1,
			 "colsample_bytree":1,
			 #"scale_pos_weight":num_neg/num_pos
			 #"scale_pos_weight":2
	         }
	plst = param.items()          
	plst += [('eval_metric', 'mae')]
	
	#train XGBoost
	#bst = xgb.cv( plst, xtrain, 200,nfold=5,callbacks=[xgb.callback.print_evaluation(show_stdv=False),xgb.callback.early_stop(3)])
	bst = xgb.train( plst, xtrain, 50, evallist,early_stopping_rounds=10)
	#bst = xgb.train( plst, xtrain, 50, evallist)
	
	#save model
	bst.save_model(args.vectors+'.xgboost')
	
	#get feature importances
	importance = bst.get_fscore()
	importance = sorted(importance.items(), key=operator.itemgetter(1))
	print importance
	
	predictions = bst.predict(xeval)
	
	tmp = pd.DataFrame()
	tmp['target'] = test_targets.target.values
	tmp['predictions'] = predictions
	tmp.to_csv('predictions.csv',index=False)
	
	#plt.scatter(x=test_targets,y=predictions)
	#plt.show()

	#dump model to .c code
	convert_model_to_c(bst,args,numf)
	
def convert_model_to_c(bst,args,numf):
	#dump model and write .c file
	bst.dump_model('dump.raw.txt') 
	num_nodes = []
	mmax = 0
	with open('dump.raw.txt') as f:
		for row in f:
			if row.startswith('booster'):
				if row.startswith('booster[0]'):
					mmax = 0
				else:
					num_nodes.append(mmax+1)
					mmax = 0
				continue			
			l=int(row.rstrip().replace(' ','').split(':')[0])
			if l > mmax:
				mmax = l
	num_nodes.append(mmax+1)
	forest = []
	tree = None
	b = 0
	with open('dump.raw.txt') as f:
		for row in f:
			if row.startswith('booster'):
				if row.startswith('booster[0]'):
					tree = [0]*num_nodes[b]
					b += 1
				else:
					forest.append(tree)
					tree = [0]*num_nodes[b]
					b+=1
				continue		
			#if b == len(num_nodes)-10: break	
			l=row.rstrip().replace(' ','').split(':')
			if l[1][:4] == "leaf":
				tmp = l[1].split('=')
				tree[int(l[0])] = [-1,float(tmp[1]),-1,-1] #!!!!
			else:	
				tmp = l[1].split('yes=')
				tmp[0]=tmp[0].replace('[Features','')
				tmp[0]=tmp[0].replace('[Feature','')
				tmp[0]=tmp[0].replace(']','')
				tmp2 = tmp[0].split('<')
				if float(tmp2[1]) < 0: tmp2[1] = 1
				tmp3 = tmp[1].split(",no=")
				tmp4 = tmp3[1].split(',')
				tree[int(l[0])] = [int(tmp2[0]),float(tmp2[1]),int(tmp3[0]),int(tmp4[0])]
		forest.append(tree)
	
		tmp = args.vectors.replace('.','_')
		tmp2 = tmp.split('/')
		with open(tmp+'_c.c','w') as fout:
			fout.write("static float score_"+args.type+"(unsigned short* v){\n")
			fout.write("float s = 0.;\n")
			#for tt in [0]:
			#for tt in range(len(forest)-10):
			for tt in range(len(forest)):
				fout.write(tree_to_code(forest[tt],0,1))
			fout.write("\nreturn s;}\n")
	
		with open(tmp+'.pyx','w') as fout:
			fout.write("cdef extern from \"" + tmp2[-1] + "_c.c\":\n")
			fout.write("\tfloat score_%s(short unsigned short[%i] v)\n\n"%(args.type,numf))
			fout.write("def myscore(sv):\n")
			fout.write("\tcdef unsigned short[%i] v = sv\n"%numf)
			fout.write("\treturn score_%s(v)\n"%args.type)
	
	os.remove('dump.raw.txt')


def tree_to_code(tree,pos,padding):
	p = "\t"*padding
	if tree[pos][0] == -1:
		if tree[pos][1] < 0:
			return p+"s = s %f;\n"%tree[pos][1]
		else:
			return p+"s = s + %f;\n"%tree[pos][1]			
	return p+"if (v[%i]<%i){\n%s}\n%selse{\n%s}"%(tree[pos][0],tree[pos][1],tree_to_code(tree,tree[pos][2],padding+1),p,tree_to_code(tree,tree[pos][3],padding+1))
	
def print_logo():
	logo = """                                   
 _____ _____ _____ _____ _____ 
|_   _| __  |  _  |     |   | |
  | | |    -|     |-   -| | | |
  |_| |__|__|__|__|_____|_|___|
                               
           """
	print logo
	print "Version 2.1\nby Sven Degroeve\n"

if __name__ == "__main__":
	print_logo()
	main()        

	
	
