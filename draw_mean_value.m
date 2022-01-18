
 
%  [dir_name_ad  ,num_ad] = GetFiles('/media/wjq/brain2/4paper_fiber_segmentation_according_location/finally_results_2021_12_29/tettt/'); %test_draw    results_paper4
 [dir_name_ad  ,num_ad] = GetFiles('/media/wjq/brain2/4paper_fiber_segmentation_according_location/code_4_paper/github/results/MAT/'); %test_draw    results_paper4



for temp_num = 1: num_ad  
    current_name = dir_name_ad{1,temp_num};
    str_pro = ['FA';'MD';'RD';'AD';'FC']
    %data_m=load(current_name, '-ascii')  
    data= load(current_name)  
   data = data.vale_name
    lenff =  length(data)
    if lenff == 0
        continue
    end
   length_name =  length(current_name)
   hemi_str = current_name(length_name-4)
   cluster_num = current_name(length_name-8:length_name-6)
   save_direc = [ '/media/wjq/brain2/4paper_fiber_segmentation_according_location/code_4_paper/github/results/PIC/']
   mkdir(save_direc)
   save_direc_p_value = [ '/media/wjq/brain2/4paper_fiber_segmentation_according_location/code_4_paper/github/results/PIC/']
 
   
 for currrent_proper = 1:5

     if    currrent_proper==3|| currrent_proper==4
         continue;
     end
     
     
    diffusion = str_pro(currrent_proper,:)
    DIFFU_NAME = str_pro(currrent_proper, :)
    db_all_csd_peak_val = data( :, currrent_proper,:)
    db_all_csd_peak_val = squeeze(db_all_csd_peak_val) 
    hold off;
    close(gcf);
    befor_p_valeus =  draw_curve_mean(db_all_csd_peak_val(24:45,10:110), db_all_csd_peak_val(1:23,10:110),save_direc,hemi_str,cluster_num,diffusion,mmse_all_subject_values)

    [FDR_p_valeus ] = mafdr(befor_p_valeus ,'BHFDR', true);
    x_ind = find(FDR_p_valeus<=0.05);
    y_pvss = FDR_p_valeus(x_ind);
    len_icon = size(x_ind,2)

      if  len_icon<1
          continue
      end
    
    if  len_icon>1 &&  len_icon < 30
 
       %close(gcf);
      % befor_p_valeus =  draw_curve_mean(db_all_csd_peak_val(24:45,:), db_all_csd_peak_val(1:23,:),save_direc,hemi_str,cluster_num,diffusion,mmse_all_subject_values)
       cluster_belong_tract = tract_seg(cluster_num);
      % ok = display_his(befor_p_valeus, FDR_p_valeus,cluster_belong_tract,save_direc_p_value,hemi_str,number_cluser,diffusion,currrent_proper)
        save_name_str = ['/media/wjq/brain2/4paper_fiber_segmentation_according_location/code_4_paper/github/results/PIC/',hemi_str,'_',cluster_num,'_',diffusion,'.mat']
      save(save_name_str,'FDR_p_valeus')     
    end
 
    
 end

end


