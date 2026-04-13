from Bio.Cluster import treecluster
import numpy
from numpy import random
from sys import argv

# # 随机生成矩阵并保存
# data = random.random(size=(10,12))
# print(data)
# metrixs = open('data.txt','w+')
# for i in range(len(data)):
# 	for j in range(len(data[i])):
# 		metrixs.write(str(data[i][j]))
# 		metrixs.write('\t')   
# 	metrixs.write('\n')      
# metrixs.close()

# 从文件中读取数据
bgc_input = argv[1]
tree_output = argv[2]


data = numpy.loadtxt(bgc_input)
# print(data)

tree = treecluster(data, method='a', dist='u')
# tree.scale()
# print(tree)
# print(type(tree))

trees = open(tree_output,'w+')
trees.write(str(tree))
trees.close()