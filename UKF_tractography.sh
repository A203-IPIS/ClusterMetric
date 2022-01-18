################## - ReadME - ######################
## wjq use proprocess AD  data.
#2019-10-27
################ - Step.0 - input and output data tempr -  ######################
input=/media/wjq/brain1/3paper_diabetes/DTI_data_75
output=/media/zeng/0003D399000FC463/wjq
slicer_path=/home/wjq/soft/Slicer-4.10.2-linux-amd64/lib/Slicer-4.10
ukf_path=/home/wjq/soft/ukf_build_run/UKFTractography-build/UKFTractography/bin
######## - Step.1 - transform .grad file use bval and bvec file ###########
for case in `ls $input`
do
echo "transform .grad file use bval and bvec file....................."  $case
python  /media/zeng/0003D399000FC463/wjq/sh/grad2bval.py   $input/$case/grad    $input/$case     $input/$case
done

########### - Step.2 - transform .nrrd file use predti.nii and grad file ########


###################################################
for case in `ls $input`
do
echo "transform .nrrd file use predti.nii and grad file.................." $case
 cd  $slicer_path/cli-modules
./DWIConvert    --inputBVectors  $input/$case/dti.bvec  --inputBValues  $input/$case/dti.bval   -o  $input/$case/predti.nhdr   --conversionMode   FSLToNrrd     --inputVolume  $input/$case/EPI_eddy_dti.nii.gz  --allowLossyConversion
done
###################################################


########### - Step3. - extract b0 image and brain_mask  ########
for case in `ls $input`
do
echo "extract b0 image and brain_mask................" $case
#fslroi  $input/$case/Eddy_EPI_dti.nii.gz   $input/$case/b0   0  1
#bet  $input/$case/b0.nii.gz  $input/$case/Brain  -m -f 0.15
done
#################################################

########### - Step4. - transform the brainmask.nii to brainmask.nrrd  ########
for case in `ls $input`
do
echo "transform the brainmask.nii to brainmask.nrrd................" $case
 cp    /media/wjq/brain1/2paper/mask_dti.bval  $input/$case
 cp    /media/wjq/brain1/2paper/mask_dti.bvec  $input/$case
 cd  $slicer_path/cli-modules
 ./DWIConvert    --inputBVectors  $input/$case/mask_dti.bvec  --inputBValues  $input/$case/mask_dti.bval   -o  $input/$case/mask_dti.nrrd   --conversionMode   FSLToNrrd     --inputVolume  $input/$case/Brain_mask.nii.gz  --allowLossyConversion
done
#################################################


########### - Step5. - transform the brainmask.nii to brainmask.nrrd  ########
for case in `ls $input`
do
echo "transform the brainmask.nii to brainmask.nrrd................" $case

cd  $slicer_path/cli-modules
./DWIConvert    --inputBVectors  /mnt/hgfs/E/3paper/mask_dti.bvec  --inputBValues  /mnt/hgfs/E/3paper/mask_dti.bval   -o  $input/$case/Brain_mask_mask.nrrd   --conversionMode   FSLToNrrd     --inputVolume  $input/$case/Brain_mask_mask.nii.gz  --allowLossyConversion
done
#################################################

#     --stoppingFA 0.15   --seedsFile    $input/$case/mask_dti.nrrd

########### - Step6. - ukftractography  ########
for case in `ls $input`
do
echo "ukftractography.........going......." $case
cd  $ukf_path
$ukf_path/UKFTractography   --dwiFile   $input/$case/predti.nhdr  --maskFile $input/$case/mask_dti.nrrd    --tracts    $input/$case/ukf_DB_wholebrain.vtk    --seedsPerVoxel  5          --numThreads  16   --numTensor  2   --stepLength  0.3
done
#################################################



##########################over######################################3
