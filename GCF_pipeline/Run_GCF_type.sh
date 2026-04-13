bin=/dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/GCF_type
r_dir=/dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/envs/gcf/bin/Rscript
#手动修改
run_gcf_result_dir=
output_dir=

if [ ! -d ${output_dir}/tmp ];then mkdir ${output_dir}/tmp;fi

#1. 获取每个BGC的GCF, MAG和type信息
echo "extract BGC type info..."
${r_dir} ${bin}/0_get_bgc_type.R ${run_gcf_result_dir}/2_Bigslice_Output/result/data.db ${output_dir}/tmp/BGC_type.txt
echo "done"
echo "merge BGC info..."
perl ${bin}/1_merge_bgc_info.pl ${run_gcf_result_dir}/2_Bigslice_Output/GCF_info.txt ${output_dir}/tmp/BGC_type.txt ${output_dir}/1_BGC_type_info.txt
echo "done"

#2. 按phylum分别统计不同type下BGC和GCF的数目
echo "count number of phylum..."
perl ${bin}/2_get_mags_type_info.pl ${output_dir}/1_BGC_type_info.txt ${output_dir}/2_Phylum_type_info.txt
echo "done"
