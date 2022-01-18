#  ClusterMetric
![image](https://github.com/A203-IPIS/ClusterMetric/blob/main/bundle%20profile.png)

(a) displays the process of whole-brain finer parcellation and node identification.The code is 'bundle_profile_example.py'.
(b) displays the atlas nodes. (c) displays the statistical results and the significant regions within  a cluster. The code is 'bunbdle_profile.m'  and 'Mark_significant_regions within_bundle.py'. 

#  Dependance 
This code depends on whitematteranalysis, Dipy, Fixel, BUAN, Mrtrix3, UKF Tracrography.

#  Steps 
 (1) Whole-brain tractography using UKF['UKF_tractography.sh' in this code].
 
 
 (2) Fiber parcellation['WMA_Parcellation.sh'] node identification using ClusterMetric[bundle_profile_example.py].
 
 
 (3)DTI and Fixel computation using Mrtrix3.
 
 
 (4)building bundle profiles and statistical analysis['bunbdle_profile.m'].
 
