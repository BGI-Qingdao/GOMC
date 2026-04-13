import numpy as np
# from sklearn.metrics.pairwise import paired_distances
from sklearn.metrics import pairwise_distances
from sys import argv

bgc_input = argv[1]
gcf_input = argv[2]
cos_out = argv[3]
bgc = np.loadtxt(bgc_input)
gcf = np.loadtxt(gcf_input)

# matrix = paired_distances(bgc, gcf, metric="cosine", n_jobs=1)
matrix = pairwise_distances(bgc, gcf, metric="cosine", n_jobs=1)
# print(matrix)

metrixs = open(cos_out,'w+')
for i in range(len(matrix)):
	for j in range(len(matrix[i])):
		metrixs.write(str(matrix[i][j]))
		metrixs.write('\t')   
	metrixs.write('\n')      
metrixs.close()