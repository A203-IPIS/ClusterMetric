function [ok] = draw_curve_mean(all_fa_mat_co_L, all_fa_mat_HFS_L,dir_save,hemi,cluster_num,diffusion,mmse_all_subject_values)

weight_size = size(all_fa_mat_co_L,2);
x=1:1:weight_size;
all_fa_mat_co_L(find(isnan(all_fa_mat_co_L)==1)) = 0;
all_fa_mat_HFS_L(find(isnan(all_fa_mat_HFS_L)==1)) = 0;
 %%%%%%%%%%%%%%%%%%%%
 all_fa_mat_SD_co_l  = zeros(1,weight_size)
all_fa_mat_SD_hfs_l = zeros(1,weight_size)
%%%%%%%%%%%%%%%
co_mean  =   zeros(1,weight_size)
hfs_mean =   zeros(1,weight_size)
%%%%%%%%%%%%%%%%%%%%%
p_values_mat = zeros(1,weight_size)
%%%%%%%%%%%%%%
for  current_conumb = 1:weight_size
        co_curr_nod = all_fa_mat_co_L(:,current_conumb)
         x_ind = find(co_curr_nod<=0)
         co_curr_nod(x_ind)=[]         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%      
         hfs_curr_nod = all_fa_mat_HFS_L(:,current_conumb)
         x_ind = find(hfs_curr_nod<=0)
         hfs_curr_nod(x_ind)=[]
         curr_co_mean =  mean(co_curr_nod);
         curr_hfs_mean = mean(hfs_curr_nod);
         
        if   isnan(curr_co_mean )
             curr_co_mean=0      
        end      
         if   isnan(curr_hfs_mean)
             curr_hfs_mean=0      
         end      
 
         curr_all_fa_mat_SD_co_l = std(co_curr_nod);
         curr_all_fa_mat_SD_hfs_l = std(hfs_curr_nod);
         if   isnan(curr_all_fa_mat_SD_co_l )
             curr_all_fa_mat_SD_co_l=0      
         end
         
         if   isnan(curr_all_fa_mat_SD_hfs_l)
             curr_all_fa_mat_SD_hfs_l=0      
         end  

        %%%%%%%%%%%%%%%%%%%%%%%%%%
         co_mean(current_conumb) =  curr_co_mean
         hfs_mean(current_conumb) = curr_hfs_mean
         %%%%%%%%%%%%%%%%%%%%%%%%%%
         all_fa_mat_SD_co_l(current_conumb) = curr_all_fa_mat_SD_co_l
         all_fa_mat_SD_hfs_l(current_conumb) = curr_all_fa_mat_SD_hfs_l        
         x_inut_test =  co_curr_nod(:,1)
         y_inut_test =  hfs_curr_nod(:,1)
         if  size(x_inut_test)<5
             continue;
         end
             
             
         [H,P,CI]= ttest2(x_inut_test,y_inut_test)
         p_values_mat(current_conumb) =  P
end
 
scale = 0.05
all_fa_mat_SD_co_l  = all_fa_mat_SD_co_l*scale
all_fa_mat_SD_hfs_l = all_fa_mat_SD_hfs_l*scale

for  loop_i = 1:1
        all_fa_mat_SD_co_l = media_filter(all_fa_mat_SD_co_l)
        all_fa_mat_SD_hfs_l = media_filter(all_fa_mat_SD_hfs_l)
end
 
for  loop_i = 1:1
        co_mean = media_filter(co_mean)
        hfs_mean = media_filter(hfs_mean)
end





%%% co
pro_low_co_l = media_filter(co_mean) - all_fa_mat_SD_co_l;
pro_high_co_l = media_filter(co_mean) + all_fa_mat_SD_co_l;

%%%%hfs
pro_low_hfs_l = hfs_mean - all_fa_mat_SD_hfs_l;
pro_high_hfs_l = hfs_mean + all_fa_mat_SD_hfs_l;
%%%%%%%%%%% co 
plot_sd_x = [(1:weight_size),fliplr(1:weight_size)];
plot_sd_y = [pro_low_co_l,fliplr(pro_high_co_l)];
h=fill(plot_sd_x ,plot_sd_y ,[255/255 192/255 193/255]);
set(h,'edgealpha',1,'facealpha',1) 
 hold on;
 %plot(x,co_mean,'LineWidth',4.5,'Color',[1 0 0]) ;
%hold on;
 %%%% hfs123,104,238
plot_sd_x = [(1:weight_size),fliplr(1:weight_size)];
plot_sd_y = [pro_low_hfs_l,fliplr(pro_high_hfs_l)];
h=fill(plot_sd_x,plot_sd_y,[123/255 104/255  238/255]);
set(h,'edgealpha',1,'facealpha',0.89) 
 hold on;
 
 
 
save_path = [dir_save,'mean_',hemi,'-',diffusion,'-',cluster_num,'.bmp']
cluster_nam = tract_seg(cluster_num);
[FDR_p_valeus ] = mafdr(p_values_mat ,'BHFDR', true);
ok = display_his(p_values_mat(1:100), FDR_p_valeus(1:100),cluster_nam,save_path,hemi,cluster_num,diffusion,3);

hold on;
plot(x,hfs_mean,'LineWidth',2.9,'Color',[0,0,1]) ;
hold on;
plot(x,co_mean,'LineWidth',2.9,'Color',[1 0 0]) ;
hold on;
 %%%%%%%%%%save  

 %set(gca,'YTick',[0:0.2:1]);
 STR1 = strcat(hemi,'-');
STR2 = strcat(STR1,cluster_num);
STR3 = [STR2,'-',diffusion,'/Node'];
xlabel(STR3,'Fontsize',25,'Fontname','Times New Roman');
%ylabel('Mean and SD','Fontsize',25,'Fontname','Times New Roman');
%set(gca, 'Fontname', 'bold','FontSize',38);
 %set(gca,'fontweight','bold', 'FontSize',50)
 set(get(gca,'xlabel'),'fontsize',30)
 set(gca,'FontName','Times New Roman','FontSize',30);%  'FontWeight','bold'
  xlim([0,100]); 
 set(gca,'XTick',[0:20:100]);

  title([cluster_nam,' Tract'],'fontname','Times New Roman','FontSize',25);
 %margin = get(gca, 'TightInset');
%set(gca, 'Position', [0+margin(1) 0+margin(2) 1-margin(1)-margin(3) 1-margin(2)-margin(4)]);
% h = legend('Healthy Controls (SD)' ,'Alzheimers Disease (SD)','Alzheimers Disease (Mean)','Healthy Controls (Mean)' );
% set(h,'Fontname','Times New Roman','Fontsize',6);
%ylim([0,2]); 
%   backColor = [0 0 0];
%   set(gca, 'color', backColor);


len_icon=0;
 if  diffusion== 'FA'
        hold on 
        [FDR_p_valeus ] = mafdr(p_values_mat ,'BHFDR', true);
        x_ind = find(FDR_p_valeus<=0.05);
        y_pvss_size = size(x_ind,2); 
        y_pvss = 0.05*ones(y_pvss_size);
       plot(x_ind,y_pvss,'.','Color','g','MarkerSize',15); 
     ylim([0.0,1]); 
      set(gca,'YTick',[0: 0.2:  1]);
 end

 if  diffusion== 'MD'
        ylim([0,2.5e-3]); 
        hold on 
        [FDR_p_valeus ] = mafdr(p_values_mat ,'BHFDR', true);
        x_ind = find(FDR_p_valeus<=0.05);
        y_pvss_size = size(x_ind,2); 
        y_pvss = (1e-4)*ones(y_pvss_size);
       plot(x_ind,y_pvss,'.','Color','g','MarkerSize',15); 
      set(gca,'YTick',[0: 10e-4:  4e-3]);
 
 end

  if  diffusion== 'FC'
     ylim([0,1.8]); 
       hold on 
        [FDR_p_valeus ] = mafdr(p_values_mat ,'BHFDR', true);
        x_ind = find(FDR_p_valeus<=0.05);
        y_pvss_size = size(x_ind,2); 
        y_pvss = (0.3)*ones(y_pvss_size);
       plot(x_ind,y_pvss,'.','Color','g','MarkerSize',15); 
         set(gca,'YTick',[0: 0.4:  1.6]);
 
 end
 %%%%%%%%%%%%%%%%%%%%5
%grid on;
     len_icon = size(x_ind,2)
      if  len_icon>1  &&  len_icon < 30
         save_path = [dir_save,'mean_',hemi,'-',diffusion,'-',cluster_num,'.bmp']
         saveas(gcf,save_path);
      end
 % hold off;
%close(gcf);
 ok = p_values_mat;

end