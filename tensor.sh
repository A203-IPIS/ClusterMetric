Input=/media/wjq/brain2/4paper_fiber_segmentation_according_location/AD_data/AD


for dir in `ls $Input`
do
dwi2tensor  $Input/$dir/dti.nii.gz   -fslgrad  $Input/$dir/dti.bvec   $Input/$dir/dti.bval  $Input/$dir/tensor_6.nii  -mask $Input/$dir/brain_mask.nii.gz      -force

tensor2metric  $Input/$dir/tensor_6.nii  -fa   $Input/$dir/FA_image.nii    -ad  $Input/$dir/AD_image.nii   -rd  $Input/$dir/RD_image.nii   -adc   $Input/$dir/MD_image.nii  -mask $Input/$dir/brain_mask.nii.gz -force
done
