#use white matter analysis soft to cluser 800 fibers
##############################################
############################################
#write by wjq-2020-3-30
#input address
Input=/media/wjq/brain1/3paper_diabetes/dti_data_remain_1122
Output=/media/zeng/0003D399000FC463/wjq
wma=/home/zeng/whitematteranalysis-master/bin
atlas=/home/wjq/soft/ORG-800FC-100HCP
hcp1_path=/home/wjq/soft/ORG-800FC-100HCP1
hardtrans_path=/home/wjq/anaconda2/bin
slicer_path=/home/wjq/soft/Slicer-4.10.2-linux-amd64/Slicer

###########step1注册到图谱空间#######################
for dir in `ls $Input`
do
echo "tckconvert..............." $dir
tckconvert    $Input/$dir/csd_WMA/csd_40_thouand_fibers.tck  $Input/$dir/csd_WMA/csd_40_thouand_fibers.vtk  -nthreads 16   -force

done
#################################################

###########step2注册到图谱空间#######################
for dir in `ls $Input`
do
echo "wm_register_to_atlas_new.py..............." $dir
 wm_register_to_atlas_new.py -l 40 -mode affine  $Input/$dir/csd_WMA/csd_40_thouand_fibers.vtk   $atlas/atlas.vtp     $Input/$dir/csd_WMA
done
##########


###########step3聚类（图谱空间）#######################
for dir in `ls $Input`
do
echo "wm_cluster_from_atlas.py..............." $dir
wm_cluster_from_atlas.py -l 40    -norender       $Input/$dir/csd_WMA/csd_40_thouand_fibers/output_tractography/csd_40_thouand_fibers_reg.vtk    $atlas    $Input/$dir/csd_WMA  -j 16

done
##########


###########step4移除噪点（图谱空间）#######################
for dir in `ls $Input`
do
echo "wm_cluster_remove_outliers.py..............." $dir
wm_cluster_remove_outliers.py  $Input/$dir/csd_WMA/csd_40_thouand_fibers_reg   $atlas    $Input/$dir/csd_WMA/20-2/
done
##########

###########step5分l、r、c（图谱空间）#######################
for dir in `ls $Input`
do
echo "wm_cluster_from_atlas.py..............." $dir
wm_separate_clusters_by_hemisphere.py   -clusterLocationFile   $hcp1_path/cluster_hemisphere_location.txt -atlasMRML  $atlas/clustered_tracts_display_100_percent.mrml  $Input/$dir/csd_WMA/20-2/csd_40_thouand_fibers_reg_outlier_removed   $Input/$dir/csd_WMA/hemisphere_clusters
done
######################################
 
###########step6转换到原始individual空间#######################
for dir in `ls $Input`
do
echo "wm_harden_transform.py..............." $dir

python  $hardtrans_path/wm_harden_transform.py -i  -t   $Input/$dir/csd_WMA/csd_40_thouand_fibers/output_tractography/itk_txform_csd_40_thouand_fibers.tfm $Input/$dir/csd_WMA/hemisphere_clusters/tracts_commissural $Input/$dir/csd_WMA/indispace_hemisphere_clusters/tracts_commissural  $slicer_path

 python  $hardtrans_path/wm_harden_transform.py -i  -t   $Input/$dir/csd_WMA/csd_40_thouand_fibers/output_tractography/itk_txform_csd_40_thouand_fibers.tfm $Input/$dir/csd_WMA/hemisphere_clusters/tracts_left_hemisphere $Input/$dir/csd_WMA/indispace_hemisphere_clusters/tracts_left_hemisphere  $slicer_path

 python  $hardtrans_path/wm_harden_transform.py -i  -t   $Input/$dir/csd_WMA/csd_40_thouand_fibers/output_tractography/itk_txform_csd_40_thouand_fibers.tfm $Input/$dir/csd_WMA/hemisphere_clusters/tracts_right_hemisphere $Input/$dir/csd_WMA/indispace_hemisphere_clusters/tracts_right_hemisphere  $slicer_path


done
######################################
