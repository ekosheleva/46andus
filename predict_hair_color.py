import numpy as np
import pandas as pd
import numpy.random as rand
from sklearn.preprocessing import OneHotEncoder

def predict_phenotype(allsnps, forest):
	#fill child snps

	nsnps = len(allsnps.index)

	parentGenotype = np.array(allsnps[['altXX', 'altXY']])
	parentGenotype = np.concatenate([parentGenotype[:, 0], parentGenotype[:, 1]])

	#store heterozygous sites	
	hets = parentGenotype == 1
	nhets = np.sum(hets)

	childGenotype = np.zeros((min(2**nhets, 256), 2*len(allsnps)))

	#childGenotype[parentGenotype == 0, :] = 0
	childGenotype[:, parentGenotype == 2] = 1

	print 'Number of heterozygous sites is', nhets

	totalProbabilities = np.zeros(3)

	normalize = 256.0

	if nhets <= 8: #fewer than 256 possible genotypes, enumerate all possibilities
		normalize = 2**nhets.astype(float)
           	childGenotype[:, hets] = np.array([[int(x) for x in bin(g)[2:].zfill(nhets)] for g in range(0, 2**nhets)])
    		fullChildGenotype = childGenotype[:, 0:nsnps] + childGenotype[:, nsnps:]
		fullChildGenotype = OneHotEncoder(n_values=3).fit(fullChildGenotype).transform(fullChildGenotype).toarray()
    		marginal = forest.predict_proba(fullChildGenotype)  
    		totalProbabilities = np.sum(marginal, 0)/normalize 
	else: #otherwise, randomly sample 256 genotypes with replacement
    		genotypes = rand.randint(0, 2**nhets, 256)
    		childGenotype[:, hets] = np.array([[int(x) for x in bin(g)[2:].zfill(nhets)] for g in genotypes])
    		fullChildGenotype = childGenotype[:, 0:nsnps] + childGenotype[:, nsnps:]
		fullChildGenotype = OneHotEncoder(n_values=3).fit(fullChildGenotype).transform(fullChildGenotype).toarray()    		
		marginal = forest.predict_proba(fullChildGenotype)  
    		totalProbabilities = np.sum(marginal, 0)/normalize

	#print totalProbabilities	
	return totalProbabilities

