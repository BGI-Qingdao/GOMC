#修改gcf结果路径
gcf_result=
output=

plot=yes
gcf_fun=
fun_order=
phylum_order=

python=/dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/envs/gcf/bin/python

#1.检查BigFam输入矩阵是否存在
if [ ! -e /dellfsqd2/ST_OCEAN/USER/zhouchanghao/Database/GCF/Big_Fam/Big_fam_BGC_feature_input.txt ]
then
$python /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/extract_bgc_features_matrix.py /dellfsqd2/ST_OCEAN/USER/zhouchanghao/metagenome/get_gcf/Big_Fam/full_run_result /dellfsqd2/ST_OCEAN/USER/zhouchanghao/Database/GCF/Big_Fam/Big_fam_BGC_feature.txt
perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/prepare_bgc_for_cal_cos.pl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/Database/GCF/Big_Fam/Big_fam_BGC_feature.txt /dellfsqd2/ST_OCEAN/USER/zhouchanghao/Database/GCF/Big_Fam/Big_fam_BGC_feature_input.txt
fi

#2.提取并处理BGC矩阵
if [ ! -d ${output}/tmp/ ];then mkdir ${output}/tmp;fi
$python /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/extract_bgc_features_matrix.py ${gcf_result}/2_Bigslice_Output ${output}/tmp/All_BGC_feature.txt
perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/Novelty/0_sort_hmm_order.pl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/Database/GCF/Big_Fam/Big_fam_BGC_feature.txt ${output}/tmp/All_BGC_feature.txt ${output}/tmp/All_BGC_feature_sorted.txt
perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/prepare_bgc_for_cal_cos.pl ${output}/tmp/All_BGC_feature_sorted.txt ${output}/0_All_BGC_feature_sorted_input.txt

#3.计算BGC, GCF与BigFam的距离
$python /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/Novelty/1_cal_dis.py ${output}/0_All_BGC_feature_sorted_input.txt /dellfsqd2/ST_OCEAN/USER/zhouchanghao/Database/GCF/Big_Fam/Big_fam_BGC_feature_input.txt ${output}/1_All_BGC_Dist.txt
perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/Novelty/2_get_min_cos_dis.pl ${gcf_result}/2_Bigslice_Output/GCF_info.txt ${output}/1_All_BGC_Dist.txt ${output}/2_All_BGC_Min_Dist.txt ${output}/3_All_GCF_Novelty.txt 0.2

#4.绘制function和phylum的柱形图
if [ $plot == yes ]
then
if [ ! -d ${output}/plot ];then mkdir ${output}/plot;fi
perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/Novelty/3_get_gcf_function_info.pl $fun_order $gcf_fun ${output}/3_All_GCF_Novelty.txt ${output}/plot/plot_fun_input.txt
perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/Novelty/4_get_gcf_phylum_info.pl $phylum_order ${gcf_result}/2_Bigslice_Output/GCF_info.txt ${output}/3_All_GCF_Novelty.txt ${output}/plot/plot_phylum_input.txt
perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/Novelty/5_convert2matrix.pl $fun_order ${output}/plot/plot_fun_input.txt ${output}/plot/Function_matrix.txt
perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/GCF_pipline/bin/Novelty/5_convert2matrix.pl $phylum_order ${output}/plot/plot_phylum_input.txt ${output}/plot/Phylum_matrix.txt
fi

