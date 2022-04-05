#  ClusterMetric
![image](https://github.com/A203-IPIS/ClusterMetric/blob/main/node_clusterometric.jpg)


(A) displays the process of whole-brain finer parcellation and node identification. 


(B) displays the atlas nodes and the statistical results and the significant regions within  a cluster. 

(c) displays the node identification in clusterometric analysis. 
#  Dependance 
This code depends on [whitematteranalysis](https://github.com/SlicerDMRI/whitematteranalysis), [Dipy](https://dipy.org/), [Mrtrix3](https://mrtrix.readthedocs.io/en/latest/index.html), [UKF Tracrography](https://github.com/pnlbwh/ukftractography).
 
#  Steps 
 (1) Whole-brain tractography using UKF['UKF_tractography.sh' ].
 
 
 (2) Fiber parcellation['WMA_Parcellation.sh'] node identification using ClusterMetric[bundle_profile_example.py].
 
 
 (3)DTI parameter and Fixel computation using Mrtrix3.
 
 
 (4)Building bundle profiles and statistical analysis['bunbdle_profile.m'].
 
 (5) Marking the results ['Mark_significant_regions within_bundle.py']. 

 
