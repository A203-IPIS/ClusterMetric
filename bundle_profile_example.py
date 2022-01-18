import numpy as np
from dipy.stats.analysis import assignment_map
from dipy.io.streamline import load_trk
import os
import nibabel as nib
from dipy.stats.analysis import afq_profile
from dipy.tracking.streamline import values_from_volume
from dipy.tracking.streamlinespeed import set_number_of_points
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline
import os.path as op
import AFQ.models.dti as dti
from fury import actor, window
from dipy.io.streamline_wjq import save_tractogram, load_tractogram_wjq
from dipy.io.stateful_tractogram import Space
from math import acos
from scipy import stats
import random

import scipy.io as scio
from ClusterMetric import test_clsutermetric
from dipy.tracking.streamline import transform_streamlines

def  get_FA_Value_func(fa_img,streams,indexs,affine_w):
     values = values_from_volume(fa_img, streams, affine_w)
     value_one_dim = []
     fa_vale_100 = []
     for num_line in range(len(values)):
         for line_point_num in range(len(values[num_line])):
             value_one_dim.append(values[num_line][line_point_num])
     for node_num in range(0,120):
         inde_num = np.where(indexs==node_num)
         inde_num = np.array(inde_num)
         value_one_dim = np.array(value_one_dim)
         value_fa = value_one_dim[inde_num]
         value_fa[np.isnan(value_fa)] = 0
         idx = np.argwhere(np.all(value_fa[..., :] < 0.1, axis=0))
         value_fa = np.delete(value_fa, idx, axis=1)
         divide_d = len(value_fa[0,:])
         if divide_d==0:
            divide_d =1
         mean_value_fa = sum(value_fa[0,:])/divide_d
         fa_vale_100.append(mean_value_fa)
     fa_vale_100 = np.array(fa_vale_100)
     return fa_vale_100

################
def get_MD_Value_func(fa_img, streams, indexs, affine_w):
    values = values_from_volume(fa_img, streams, affine_w)
    value_one_dim = []
    fa_vale_100 = []
    for num_line in range(len(values)):
        for line_point_num in range(len(values[num_line])):
            value_one_dim.append(values[num_line][line_point_num])
    for node_num in range(0, 120):
        inde_num = np.where(indexs == node_num)
        inde_num = np.array(inde_num)
        value_one_dim = np.array(value_one_dim)
        value_fa = value_one_dim[inde_num]
        value_fa[np.isnan(value_fa)] = 0
        idx = np.argwhere(np.all(value_fa[..., :] < 0.0000001, axis=0))
        value_fa = np.delete(value_fa, idx, axis=1)
        divide_d = len(value_fa[0, :])
        if divide_d == 0:
            divide_d = 1
        mean_value_fa = sum(value_fa[0, :]) / divide_d
        fa_vale_100.append(mean_value_fa)
    fa_vale_100 = np.array(fa_vale_100)
    return fa_vale_100


def get_fc_value_func(FC_index_data,FC_direction_data,FC_value_data,bundle_line_individual_space,indx):
    length_streams = len(bundle_line_individual_space)
    current_line_fc_value = []
    for curen_line_num in range(length_streams):
        currline_coor = bundle_line_individual_space[curen_line_num]
        currleng_line = len(currline_coor)
        currline_direct = np.zeros((currleng_line,3))
        currline_direct[0:currleng_line-2,:] = currline_coor[0:currleng_line-2,:] - currline_coor[1:currleng_line-1,:]
        currline_direct[currleng_line-1] = currline_direct[currleng_line-2]
        for currpoint in range(len(currline_coor)):
            curr_point_coord = np.round(currline_coor[currpoint])
            curr_point_dirct = currline_direct[currpoint]
            x_coor = int(curr_point_coord[0])
            y_coor = int(curr_point_coord[1])
            z_coor = int(curr_point_coord[2])
            if z_coor >=60:
               z_coor = 59
            if y_coor >= 60:
               y_coor = 59
            currondexnum = FC_index_data[x_coor,y_coor,z_coor,:]
            ind_from = int(currondexnum[1])
            dst = int(currondexnum[1] + currondexnum[0])
            fc_direction = np.squeeze(FC_direction_data[ind_from:dst])
            icon_empty = []
            conr_val_multip = np.abs(np.dot(fc_direction,curr_point_dirct))
            gg= conr_val_multip.size
            if gg > 0:
               max_inde = np.argmax(conr_val_multip)
               fc_value_point = FC_value_data[ind_from+max_inde]
               current_line_fc_value.append(fc_value_point)
            else:
                current_line_fc_value.append(0)
    array_fcvalue =  np.array(current_line_fc_value)
    ###############################
    fa_vale_100 = []
    for node_num in range(0, 120):
        inde_num = np.where(indx == node_num)
        inde_num = np.array(inde_num)
        value_fa = array_fcvalue[inde_num]
        dividee = len(value_fa[0, :])
        if dividee == 0:
           dividee = 1
        mean_value_fa = sum(value_fa[0, :]) /dividee
        fa_vale_100.append(mean_value_fa)
    fa_vale_100 = np.squeeze(np.array(fa_vale_100))

    return fa_vale_100
############################
def remove_zrros_negtive_val(input_mat):
    copy_mat = input_mat
    lengthy_in = len(input_mat)
    for curr_pointsd in range(lengthy_in):
        cur_vals = input_mat[curr_pointsd]
        if  cur_vals<0:
            copy_mat[curr_pointsd]= 0
        elif cur_vals > 0.004:
            copy_mat[curr_pointsd] = 0
    return copy_mat

#################
def read_FA_FC_date_func(cluster_name_numb,sphere_icon):
    all_five_property_mat = np.zeros((1, 5, 120))
    five_property_mat = np.zeros((5,120))
    value_45_5_100 = {}
    path_dir_AD = '/media/wjq/brain2/4paper_fiber_segmentation_according_location/code_4_paper/github/AD_data/'
    if sphere_icon =='L':
        BUNDLE_PATH_AD = '/wma/hemisphere-clusters/tracts_left_hemisphere_vtk/'
        BUNDLE_PATH_ATLAS = '/tracts_left_hemisphere/'
        individual_BUNDLE_PATH_AD = '/wma/indispace-hemisphere-clusters_vtk/tracts_left_hemisphere/'
    elif sphere_icon =='R':
        BUNDLE_PATH_AD = '/wma/hemisphere-clusters/tracts_right_hemisphere_vtk/'
        BUNDLE_PATH_ATLAS = '/tracts_right_hemisphere/'
        individual_BUNDLE_PATH_AD = '/wma/indispace-hemisphere-clusters_vtk/tracts_right_hemisphere/'
    else:
        BUNDLE_PATH_AD = '/wma/hemisphere-clusters/tracts_commissural_vtk/'
        BUNDLE_PATH_ATLAS = '/tracts_commissural/'
        individual_BUNDLE_PATH_AD = '/wma/indispace-hemisphere-clusters_vtk/tracts_commissural/'
    AD_subjects_all = os.listdir(path_dir_AD)
    AD_subjects_all.sort(key=lambda x: x[:1])
    subjects_num = len(AD_subjects_all)
    for current_subj in range(0,subjects_num):
        curre_subje_name = AD_subjects_all[current_subj]
        #######################
        famdrdadfc_all_path = path_dir_AD + curre_subje_name + '/property/'
        save_property_path = famdrdadfc_all_path + cluster_name_numb + '_' + sphere_icon
        isExists = os.path.exists(famdrdadfc_all_path)
        if not isExists:
            os.makedirs(famdrdadfc_all_path)
        ########################
        cluster_na = 'cluster_00' + cluster_name_numb +'.vtk'
        cluster_na_atlas = 'cluster_00' + cluster_name_numb + '.vtk'
        full_vtk_path_subject = path_dir_AD + curre_subje_name + BUNDLE_PATH_AD + cluster_na
        full_vtk_path_atlas = '/media/wjq/brain2/4paper_fiber_segmentation_according_location/code_4_paper/github/AD_data/subject/atlas' + BUNDLE_PATH_ATLAS + cluster_na_atlas
        data_nii_data_path =path_dir_AD +  curre_subje_name + "/data.nii.gz"
        img = nib.load(data_nii_data_path)
        affine_data = img.affine
        ######################
        bundle_individual = load_tractogram_wjq(full_vtk_path_subject, img, to_space=Space.VOX, bbox_valid_check=False)
        bundle_individual_line = bundle_individual.streamlines
        length_fiber = len(bundle_individual_line)
        if length_fiber ==0:
           f = open(save_property_path, 'a')
           np.save(save_property_path,five_property_mat)
           print(curre_subje_name)
           print(current_subj)
           continue
        bundle_atlas = load_tractogram_wjq(full_vtk_path_atlas, img, to_space=Space.VOX, bbox_valid_check=False)
        bundle_atlas_line = bundle_atlas.streamlines
        number_srg  = 101
        # indx ,_= assignment_map(bundle_individual_line, bundle_atlas_line, n)
        # indx = np.array(indx)
        indx = test_clsutermetric(bundle_individual_line, bundle_atlas_line,number_srg)
        indx = np.array(indx)
        ##########################################
        full_vtk_path_subject = path_dir_AD + curre_subje_name + individual_BUNDLE_PATH_AD + cluster_na
        bundle_individual_space = load_tractogram_wjq(full_vtk_path_subject, img, to_space=Space.VOX, bbox_valid_check=False)
        bundle_line_individual_space = bundle_individual_space.streamlines
        ###################
        data_FA_data_path = path_dir_AD + curre_subje_name + "/FA_image.nii"
        img = nib.load(data_FA_data_path)
        FA_data = img.get_fdata()
        FA_data_100 = get_FA_Value_func(FA_data,bundle_line_individual_space,indx,affine_data)
        ####################
        data_FA_data_path = path_dir_AD + curre_subje_name + "/MD_image.nii"
        img = nib.load(data_FA_data_path)
        FA_data = img.get_fdata()
        MD_data_100 = get_MD_Value_func(FA_data,bundle_line_individual_space,indx,affine_data)
        # ##################
        # data_FA_data_path = path_dir_AD + curre_subje_name + "/RD_image.nii"
        # img = nib.load(data_FA_data_path)
        # FA_data = img.get_fdata()
        # RD_data_100 = get_FA_Value_func(FA_data,bundle_line_individual_space,indx,affine_data)
        ####################
        # data_FA_data_path = path_dir_AD + curre_subje_name + "/AD_image.nii"
        # img = nib.load(data_FA_data_path)
        # FA_data = img.get_fdata()
        # AD_data_100 = get_FA_Value_func(FA_data,bundle_line_individual_space,indx,affine_data)

        ######### compute FC #####################
        fc_all_path = path_dir_AD + curre_subje_name + '/FC_ONLY_SUBJECT/'
        FC_index_path = fc_all_path + 'index.nii'
        FC_DIRECTION_path = fc_all_path + 'directions.nii'
        FC_VALUE_path = fc_all_path + 'fc_individual.nii'
        ##############################
        FC_index = nib.load(FC_index_path)
        FC_index_data = FC_index.get_fdata()
        FC_index_data =np.array(FC_index_data)
        ############################
        FC_direction = nib.load(FC_DIRECTION_path)
        FC_direction_data = FC_direction.get_fdata()
        FC_direction_data = np.array(FC_direction_data)
        ############################
        FC_values = nib.load(FC_VALUE_path)
        FC_value_data = FC_values.get_fdata()
        FC_value_data = np.array(FC_value_data)
        ############################
        bundle_line_individual_space = transform_streamlines(bundle_line_individual_space,
                                            np.linalg.inv(affine_data))
        FC_NODE_values = get_fc_value_func(FC_index_data,FC_direction_data,FC_value_data,bundle_line_individual_space,indx)

        five_property_mat[0, :] = FA_data_100
        five_property_mat[1, :] = remove_zrros_negtive_val(MD_data_100)
        # five_property_mat[2, :] = remove_zrros_negtive_val(RD_data_100)
        # five_property_mat[3, :] = remove_zrros_negtive_val(AD_data_100)
        five_property_mat[4, :] = FC_NODE_values

        print(current_subj)
        print(cluster_name_numb)
        np.save(save_property_path, five_property_mat)
        all_five_property_mat[current_subj,:,:] = five_property_mat

    all_path_resul = '/media/wjq/brain2/4paper_fiber_segmentation_according_location/code_4_paper/github/AD_data/subject/mat/'
    final_p = all_path_resul   + cluster_name_numb +  '_' + sphere_icon + '.mat'
    array_proprr = np.array(all_five_property_mat)
    value_45_5_100['vale_name'] = array_proprr
    scio.savemat(final_p, {'vale_name': value_45_5_100['vale_name']})
    okkk=0

    return 1


##################################################
def remove_zeross(x_data):
    new_data = []
   # x_data =np.array([1,0,0,9])
    idx = np.array(np.where(x_data != 0))
    length_d = len(idx)
    for curr_d in range(length_d):
        indee = idx[curr_d]
        cur_ddd = x_data[indee]
        new_data.append(cur_ddd)
    return np.array(new_data)
####################################################
def remove_nan(x_data):
    new_data = []
   # x_data =np.array([1,0,0,9])
    idx = np.array(np.where(x_data != np.nan))
    length_d = len(idx)
    for curr_d in range(length_d):
        indee = idx[curr_d]
        cur_ddd = x_data[indee]
        new_data.append(cur_ddd)
    return np.array(new_data)

hemisphere_icon = 'C'
cluster_name_numb = '364'
yyL = read_FA_FC_date_func(cluster_name_numb, hemisphere_icon)
















