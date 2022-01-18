
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
from mpl_toolkits.mplot3d import Axes3D
import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage.interpolation import map_coordinates
from scipy.spatial.distance import mahalanobis

from dipy.utils.optpkg import optional_package
from dipy.io.utils import save_buan_profiles_hdf5
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import AveragePointwiseEuclideanMetric
from dipy.tracking.streamline import (set_number_of_points,
                                      values_from_volume,
                                      orient_by_streamline,
                                      transform_streamlines,
                                      Streamlines)



def generate_3D_circle(sample_num , r , center , n,directions_line_center,curr_R0,curr_R):

    r = curr_R
    begin_pai  = (np.pi/2)- (curr_R0/curr_R)*(np.pi/6)
    end_pai =  (np.pi/2)+ (curr_R0/curr_R)*(np.pi/6)

    sample = np.linspace( begin_pai, end_pai , sample_num)

    #a = np.cross(n, [1, 0, 0])  # np.cross(), 向量叉积
    a = np.cross(n, directions_line_center)
    # if np.all(a == 0):  # a 是否为0向量
    #     a = np.cross(n, [0, 1, 0])
    b = np.cross(n, a)
    # 归一化a，b（圆面两个互垂直向量）
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    # 利用空间圆的参数方程生成圆
    c1 = center[0] * np.ones((sample_num, 1))
    c2 = center[1] * np.ones((sample_num, 1))
    c3 = center[2] * np.ones((sample_num, 1))
    [c_x, c_y, c_z] = c1 + r * a[0] * np.cos(sample) + r * b[0] * np.sin(sample), \
                      c2 + r * a[1] * np.cos(sample) + r * b[1] * np.sin(sample), \
                      c3 + r * a[2] * np.cos(sample) + r * b[2] * np.sin(sample)
    return c_x, c_y, c_z

def     draw_circle(sample_num,R,center,directions,directions_line,curr_R0,curr_R):
        x = []
        x.append(generate_3D_circle(sample_num,R,center,directions,directions_line,curr_R0,curr_R))

        for i in x:
            points_x = i[0]
            points_y = i[1]
            points_z = i[2]
        return points_x,points_y, points_z

def bundle_tractometric_new_w(segment_number,curr_R,curr_R0,ever_add_circle,bundle_atlas_line):



    #segment_number = 60
    ########################
    n_fea=segment_number*2
    indx_ss = assignment_map(bundle_atlas_line, bundle_atlas_line, n_fea)
    indx = np.array(indx_ss[0])
    center_line = np.array(indx_ss[1])
    colors = [np.random.rand(3) for si in range(n_fea)]
    disks_color = []
    for j in range(len(indx)):
        disks_color.append(tuple(colors[indx[j]]))
    # scene = window.Scene()
    # scene.SetBackground(1, 1, 1)
    # scene.add(actor.streamtube(bundle_atlas_line, colors=disks_color))
    # scene.set_camera(focal_point=(-18.17281532, -19.55606842, 6.92485857),
    #                  position=(-360.11, -340.46, -40.44),
    #                  view_up=(-0.03, 0.028, 0.89))
    #window.show(scene, size=(600, 600))
    ###########  I can  play #############
    bundle_num = len(bundle_atlas_line)
    first_points_lines_all = np.zeros((bundle_num,3))

    first_line = center_line
    first_points_lines_all[0]=first_line[0,:]
    first_line_length = len(first_line)
    first_line_start = first_line[0:4]
    first_line_end = first_line[first_line_length-4:first_line_length]
    bundle_lines = []
    for line_id in range(1,bundle_num):
        cuuren_line_data = bundle_atlas_line[line_id]
        length_cu_l = len(cuuren_line_data)
        start_three_point_in_a_line = cuuren_line_data[0:4]
        end_three_point_in_a_line = cuuren_line_data[length_cu_l-4:length_cu_l]
        ##############
        start_err = first_line_start-start_three_point_in_a_line
        sum_two_line_start = sum(np.sqrt(sum(np.dot(start_err, np.transpose(start_err)))))
        end_err = first_line_start-end_three_point_in_a_line
        sum_two_line_end = sum(np.sqrt(sum(np.dot(end_err, np.transpose(end_err)))))
        if sum_two_line_start >sum_two_line_end:
            new_line_curr = cuuren_line_data[::-1]
        else:
            new_line_curr = cuuren_line_data

        bundle_lines.append(new_line_curr)
        first_points_lines_all[line_id] = new_line_curr[0,:]

    #plt.show()
    R_radius = 65
    lamba = 0.8
    directions = center_line[1]-center_line[0]

    x2 = center_line[:,0]
    y2 = center_line[:,1]
    z2 = center_line[:,2]
    # 创建系数矩阵A
    A = np.zeros((3, 3))
    for i in range(0, n_fea):
        A[0, 0] = A[0, 0] + x2[i] ** 2
        A[0, 1] = A[0, 1] + x2[i] * y2[i]
        A[0, 2] = A[0, 2] + x2[i]
        A[1, 0] = A[0, 1]
        A[1, 1] = A[1, 1] + y2[i] ** 2
        A[1, 2] = A[1, 2] + y2[i]
        A[2, 0] = A[0, 2]
        A[2, 1] = A[1, 2]
        A[2, 2] = n_fea    #  100
    # print(A)

    b = np.zeros((3, 1))
    for i in range(0, n_fea):
        b[0, 0] = b[0, 0] + x2[i] * z2[i]
        b[1, 0] = b[1, 0] + y2[i] * z2[i]
        b[2, 0] = b[2, 0] + z2[i]

    A_inv = np.linalg.inv(A)
    X = np.dot(A_inv, b)
    #print('平面拟合结果为：z = %.3f * x + %.3f * y + %.3f' % (X[0, 0], X[1, 0], X[2, 0]))

    R = 0
    for i in range(0, n_fea):
        R = R + (X[0, 0] * x2[i] + X[1, 0] * y2[i] + X[2, 0] - z2[i]) ** 2
    #print('方差为：%.*f' % (3, R))

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("z")
    #ax1.scatter(x2, y2, z2, c='r', marker='o')

    x_p = np.linspace(-50, 200, 100)
    y_p = np.linspace(-50, 80, 50)
    x_p, y_p = np.meshgrid(x_p, y_p)
    z_p = X[0, 0] * x_p + X[1, 0] * y_p + X[2, 0]
    plan_direction = X
    A_d = -plan_direction[0]
    B_d = -plan_direction[1]
    C_d = 1
    divide_low = np.sum(A_d*A_d+B_d*B_d+C_d*C_d)
    D_plane = -plan_direction[2]
    new_line_in_plane = np.zeros((n_fea,3))
    for  num_point_id in range(0,n_fea):
        curr_points = center_line[num_point_id]
        c_x = curr_points[0]
        c_y = curr_points[1]
        c_z = curr_points[2]
        new_X = (B_d*B_d+C_d*C_d)*c_x-A_d*(B_d*c_y+c_z*C_d+D_plane)
        new_Y = (A_d * A_d + C_d * C_d) * c_y - B_d * (A_d * c_x + c_z * C_d+D_plane)
        new_Z = (A_d * A_d + B_d*B_d) * c_z - C_d * (A_d * c_x + c_y * B_d+D_plane)
        d_dist = A_d*new_X+B_d*new_Y+C_d*new_Z
        new_line_in_plane[num_point_id, 0] = new_X/divide_low
        new_line_in_plane[num_point_id, 1] = new_Y/divide_low
        new_line_in_plane[num_point_id, 2] = new_Z/divide_low

    ax1.plot3D(new_line_in_plane[:,0], new_line_in_plane[:,1], new_line_in_plane[:,2], 'blue')
    ###########################3
    length_lines = len(new_line_in_plane)
    line_filter = new_line_in_plane

    for filter_numb in range(0,30):
       for  p_id in range(1,length_lines-1):
            line_filter[p_id,:] = (new_line_in_plane[p_id-1] + new_line_in_plane[p_id+1])/2
       new_line_in_plane = line_filter

    ################## EXTANT ##################
    start_direc = new_line_in_plane[0,:] - new_line_in_plane[1,:]
    start_direc = start_direc / np.sqrt(np.sum(np.dot(start_direc, np.transpose(start_direc))))
    first_ADD_one_points = new_line_in_plane[0,:]+ 5*np.sqrt(np.sum(np.dot(start_direc, np.transpose(start_direc))))*start_direc

    start_direc = new_line_in_plane[n_fea-1,:] - new_line_in_plane[n_fea-2,:]
    start_direc = start_direc / np.sqrt(np.sum(np.dot(start_direc, np.transpose(start_direc))))
    end_ADD_one_points = new_line_in_plane[n_fea-1,:]+ 5*np.sqrt(np.sum(np.dot(start_direc, np.transpose(start_direc))))*start_direc

    all_line_cent = np.zeros((n_fea+2,3))
    all_line_cent[0,:] = first_ADD_one_points
    all_line_cent[n_fea+1,:] = end_ADD_one_points
    all_line_cent[1:n_fea+1,:] = new_line_in_plane

    final_ppp = np.zeros((n_fea,3))
    for  fina_p_id in range(0,n_fea):
        index_ids = int(fina_p_id*1.01)
        final_ppp[fina_p_id,:] = all_line_cent[index_ids,:]


    ax1.scatter(final_ppp[:, 0], final_ppp[:, 1], final_ppp[:, 2], 'green')

    new_line_in_plane = final_ppp
    directions_cir =np.zeros((3))
    directions_cir[0] = A_d
    directions_cir[1] = B_d
    directions_cir[2] = C_d
    directions_cir = directions_cir / np.sqrt(np.sum(np.dot(directions_cir, np.transpose(directions_cir))))
    #############################3
    points_and_color = []
    ########### parameters ###################
    # segment_number = 49       #   x2
    # curr_R = 30
    # curr_R0 = 50
    # ever_add_circle = 9
    ################################
    all_line_curr_point_cir_line = np.zeros((segment_number*2+1,2500, 3))
    for points_lines_ID in range(0,segment_number):
        points_lines = points_lines_ID
        curr_R = curr_R+ever_add_circle
        points_lines = int(points_lines*(segment_number/segment_number))
        directions = new_line_in_plane[points_lines+1,:]-new_line_in_plane[points_lines,:]
        directions = directions / np.sqrt(np.sum(np.dot(directions, np.transpose(directions))))
        next_points = new_line_in_plane[points_lines]+ curr_R*directions
        uuu = np.dot(directions_cir,np.transpose(directions))
        points_x,points_y, points_z = draw_circle(50, curr_R, next_points, directions_cir,directions,curr_R0,curr_R)
        ax1.plot3D(points_x[0, :], points_y[0, :], points_z[0, :], 'blue')
        lambbaa = 0.4
        curr_point_cir_line = np.zeros((2500, 3))
        for circle_vertical_id in range(0, 25):
            new_circle_center = next_points + lambbaa*circle_vertical_id * directions_cir
            points_x_NEXT,points_y_NEXT, points_z_NEXT = draw_circle(50, curr_R, new_circle_center, directions_cir,directions,curr_R0,curr_R)
           # ax1.plot3D(points_x_NEXT[0, :], points_y_NEXT[0, :], points_z_NEXT[0, :], 'green')

            begin_in = circle_vertical_id * 50
            end_in = (circle_vertical_id+1) * 50

            curr_point_cir_line[begin_in:end_in, 0] = points_x_NEXT[0, :]
            curr_point_cir_line[begin_in:end_in, 1] = points_y_NEXT[0, :]
            curr_point_cir_line[begin_in:end_in, 2] = points_z_NEXT[0, :]
            hhh=0
        for circle_vertical_id in range(0, 25):
            new_circle_center = next_points + lambbaa*circle_vertical_id * (-directions_cir)
            points_x_NEXT,points_y_NEXT, points_z_NEXT = draw_circle(50, curr_R, new_circle_center, directions_cir,directions,curr_R0,curr_R)
           # ax1.plot3D(points_x_NEXT[0, :], points_y_NEXT[0, :], points_z_NEXT[0, :], 'green')
            begin_in = (circle_vertical_id+25) * 50
            end_in   = (circle_vertical_id+25+1)  * 50
            curr_point_cir_line[begin_in:end_in, 0] = points_x_NEXT[0, :]
            curr_point_cir_line[begin_in:end_in, 1] = points_y_NEXT[0, :]
            curr_point_cir_line[begin_in:end_in, 2] = points_z_NEXT[0, :]
        all_line_curr_point_cir_line[points_lines_ID,:,:] = curr_point_cir_line
        ax1.scatter(curr_point_cir_line[:,0], curr_point_cir_line[:,1], curr_point_cir_line[:,2], 'green')

        hhh=0
    #plt.show()
    curr_R = 30
    #curr_R0 = 60
    #############################################
    new_line_in_plane = new_line_in_plane[::-1]
    directions_cir =np.zeros((3))
    directions_cir[0] = A_d
    directions_cir[1] = B_d
    directions_cir[2] = C_d
    directions_cir = directions_cir / np.sqrt(np.sum(np.dot(directions_cir, np.transpose(directions_cir))))
    curr_point_cir_line = np.zeros((2500, 3))
    for points_lines_ID in range(0,segment_number+1):   #segment_number+1
        points_lines = points_lines_ID
        curr_R = curr_R+ever_add_circle
        points_lines = int(points_lines*(segment_number/segment_number))
        directions = new_line_in_plane[points_lines+1,:]-new_line_in_plane[points_lines,:]
        directions = directions / np.sqrt(np.sum(np.dot(directions, np.transpose(directions))))
        next_points = new_line_in_plane[points_lines]+ curr_R*directions
        uuu = np.dot(directions_cir,np.transpose(directions))
        points_x,points_y, points_z = draw_circle(50, curr_R, next_points, directions_cir,directions,curr_R0,curr_R)
        ax1.plot3D(points_x[0, :], points_y[0, :], points_z[0, :], 'blue')
        lambbaa = 0.4
        for circle_vertical_id in range(0, 25):
            new_circle_center = next_points + lambbaa*circle_vertical_id * directions_cir
            points_x_NEXT,points_y_NEXT, points_z_NEXT = draw_circle(50, curr_R, new_circle_center, directions_cir,directions,curr_R0,curr_R)
           # ax1.plot3D(points_x_NEXT[0, :], points_y_NEXT[0, :], points_z_NEXT[0, :], 'green')
            begin_in = circle_vertical_id * 50
            end_in = (circle_vertical_id + 1) * 50
            curr_point_cir_line[begin_in:end_in, 0] = points_x_NEXT[0, :]
            curr_point_cir_line[begin_in:end_in, 1] = points_y_NEXT[0, :]
            curr_point_cir_line[begin_in:end_in, 2] = points_z_NEXT[0, :]
            hhh=0
        for circle_vertical_id in range(0, 25):
            new_circle_center = next_points + lambbaa*circle_vertical_id * (-directions_cir)
            points_x_NEXT,points_y_NEXT, points_z_NEXT = draw_circle(50, curr_R, new_circle_center, directions_cir,directions,curr_R0,curr_R)
            #ax1.plot3D(points_x_NEXT[0, :], points_y_NEXT[0, :], points_z_NEXT[0, :], 'green')
            begin_in = (circle_vertical_id+25) * 50
            end_in   = (circle_vertical_id+25+1)  * 50
            curr_point_cir_line[begin_in:end_in, 0] = points_x_NEXT[0, :]
            curr_point_cir_line[begin_in:end_in, 1] = points_y_NEXT[0, :]
            curr_point_cir_line[begin_in:end_in, 2] = points_z_NEXT[0, :]
        curr_line_right = segment_number*2-points_lines_ID   # segment_number*2+1
        all_line_curr_point_cir_line[curr_line_right, :, :] = curr_point_cir_line
        hhh=0

    # for  circle_plane_ID in range(0,20):
    #     circle_plane_ID = circle_plane_ID*2
    #     curr_pla = all_line_curr_point_cir_line[circle_plane_ID,:,:]
       # ax1.scatter(curr_pla[:, 0], curr_pla[:, 1], curr_pla[:, 2], 'red')
    #plt.show()
    jjj = np.reshape(all_line_curr_point_cir_line, (-1, 3))
    for  dis_id in range(1,bundle_num-1):
        display_line = bundle_lines[dis_id]
        #ax1.plot3D(display_line[:,0], display_line[:,1], display_line[:,2], 'blue')
    ax1.set_xlim(-50,100)
    ax1.set_ylim(-50,100)
    ax1.set_zlim(0,100)
    #plt.show()
    _, indx_50000 = cKDTree(jjj, 1,
                      copy_data=True).query(bundle_atlas_line.get_data(), k=1)
    colors = [np.random.rand(3) for si in range(n_fea+1)]
    disks_color = []
    ID_X = []
    for j in range(len(indx_50000)):
        disks_color.append(tuple(colors[int(indx_50000[j]/2500)]))
        ID_X.append(int(indx_50000[j]/2500))

    scene = window.Scene()
    scene.SetBackground(1, 1, 1)
    scene.add(actor.streamtube(bundle_atlas_line, colors=disks_color))
    scene.set_camera(focal_point=(-18.17281532, -19.55606842, 6.92485857),
                     position=(-360.11, -340.46, -40.44),
                     view_up=(-0.03, 0.028, 0.89))
    #window.show(scene, size=(600, 600))
    return  ID_X,  jjj
################################################
def  test_clsutermetric(bundle_bundle_in,bundle_atlas_line,segment_number):
    segment_number = int(segment_number/2)
    curr_R = 40
    curr_R0 = 60
    ever_add_circle = 3
    ############################################

    ############################################
    #bundle_bundle_in = bundle_atlas_line
    index_ID,segment_plane = bundle_tractometric_new_w(segment_number,curr_R,curr_R0,ever_add_circle,bundle_atlas_line)

    _, indx_50000 = cKDTree(segment_plane, 1,
                            copy_data=True).query(bundle_bundle_in.get_data(), k=1)
    colors = [np.random.rand(3) for si in range(segment_number*2+1)]
    disks_color = []
    ID_X = []
    for j in range(len(indx_50000)):
        disks_color.append(tuple(colors[int(indx_50000[j] / 2500)]))
        ID_X.append(int(indx_50000[j] / 2500))

    scene = window.Scene()
    scene.SetBackground(1, 1, 1)
    scene.add(actor.streamtube(bundle_bundle_in, colors=disks_color))
    scene.set_camera(focal_point=(-18.17281532, -19.55606842, 6.92485857),
                     position=(-360.11, -340.46, -40.44),
                     view_up=(-0.03, 0.028, 0.89))
    #window.show(scene, size=(600, 600))

    kkk=0
    return ID_X



def  AFQ_buan(bundle_bundle_indivi,bundle_atlas):
    tractogram_lines = bundle_bundle_indivi

    orient_by=tractogram_lines[0]
    N=50
    bundle = orient_by_streamline(tractogram_lines, orient_by)
    fgarray = set_number_of_points(bundle, N)
    length_bund = len(bundle)
    color_inde = np.zeros((N*length_bund))
    in_numb = np.zeros((N))
    for gen_i in range(0,N):
        in_numb[gen_i] = int(gen_i)

    for curr_id in range(0,length_bund):

        start_id = curr_id*N
        end_id = (curr_id+1) * N
        color_inde[start_id:end_id] =in_numb

    colors = [np.random.rand(3) for si in range(100)]
    disks_color = []
    ID_X = []
    for j in range(len(color_inde)):
        disks_color.append(tuple(colors[int(color_inde[j])]))

    scene = window.Scene()
    scene.SetBackground(1, 1, 1)
    scene.add(actor.streamtube(fgarray, colors=disks_color))
    scene.set_camera(focal_point=(-18.17281532, -19.55606842, 6.92485857),
                     position=(-360.11, -340.46, -40.44),
                     view_up=(-0.03, 0.028, 0.89))
    window.show(scene, size=(600, 600))


    n = N
    indx =  assignment_map(bundle_bundle_indivi, bundle_atlas, n)
    indx = np.array(indx[0])
    colors = [np.random.rand(3) for si in range(100)]
    disks_color = []
    for j in range(len(indx)):
        disks_color.append(tuple(colors[int(indx[j])]))
    scene = window.Scene()
    scene.SetBackground(1, 1, 1)
    scene.add(actor.streamtube(bundle_bundle_indivi, colors=disks_color))
    scene.set_camera(focal_point=(-18.17281532, -19.55606842, 6.92485857),
                     position=(-360.11, -340.46, -40.44),
                     view_up=(-0.03, 0.028, 0.89))
    window.show(scene, size=(600, 600))


def  test_364():
    cluater_name = '621'
    path_dir_DM = '/media/wjq/brain1/3paper_diabetes/DTI_data_75'
    AD_subjects_all = os.listdir(path_dir_DM)
    length_subjs = len(AD_subjects_all)
    for cur_subj in range(1,2):
        print(cur_subj)
        curr_subject_name = AD_subjects_all[cur_subj]
        path_dir_DM_ms = path_dir_DM + '/' + curr_subject_name
        BUNDLE_PATH_ATLAS =path_dir_DM_ms + '/wma/atlas_Native_space_vtk/tracts_commissural/'  # tracts_left_hemisphere_trk   tracts_commissural_trk
        cluster_na = 'cluster_00' + cluater_name + '.vtk'
        full_vtk_path_atlas = BUNDLE_PATH_ATLAS  + cluster_na
        data_nii_data_path = path_dir_DM_ms + "/data.nii.gz"
        img = nib.load(data_nii_data_path)
        bundle_atlas = load_tractogram_wjq(full_vtk_path_atlas, img, to_space=Space.VOX, bbox_valid_check=False)
        bundle_atlas_line = bundle_atlas.streamlines
        #########################
        indivodual_vtk_path_ = path_dir_DM_ms +'/wma/Native_space_vtk/tracts_commissural/' + cluster_na
        bundle_indivi = load_tractogram_wjq(indivodual_vtk_path_, img, to_space=Space.VOX, bbox_valid_check=False)
        bundle_indivi_line = bundle_indivi.streamlines
        seg_numb = 110
        test_clsutermetric(bundle_atlas_line, bundle_atlas_line,seg_numb)
        #AFQ_buan(bundle_atlas_line, bundle_atlas_line)
        hhh=0

#test_364()


