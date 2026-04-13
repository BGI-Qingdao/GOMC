# if your antismash results were obtained by antismash v6.x, use GCF_Analysis_with_Bigslice1.pl
#!! if you want to calculate the cosine distance between the BGC gbks and the BiG-FAM database, due to limitations of the BiG-FAM database, you can only use the BGC gbk from antismash v6.x or v5.x.
perl GCF_pipline/GCF_Analysis_with_Bigslice1.pl ./GCF_clustering_parameter.txt

#if your antismash results were obtained by antismash v7.x, use GCF_Analysis_with_Bigslice2.pl
perl GCF_pipline/GCF_Analysis_with_Bigslice2.pl ./GCF_clustering_parameter.txt
