 
function [ok] = display_his(p_values_before, FDR_VALUE,cluster_nam,save_direc,hemi_str,number_cluser,DIFFU_NAME,CURRENT_SEQU)


if     hemi_str=='L'
    CURRENT_SEQU=5;
elseif   hemi_str=='R'
     CURRENT_SEQU=3;
else
     CURRENT_SEQU=6;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
if  CURRENT_SEQU==1
    dd_colour = [ 245/255    222/255 179/255];
    c_colour =    [ 216/255   191/255  216/255  ];
end

if  CURRENT_SEQU==2
    c_colour = [ 255/255    182/255 193/255];
    %c_colour =    [ 72/255   209/255  204/255  ];
end

if  CURRENT_SEQU==3
    dd_colour = [ 135/255    206/255 235/255];
    c_colour =    [ 147/255   112/255  219/255  ];
end

if  CURRENT_SEQU==4
    dd_colour = [ 144/255    188/255 143/255];
    c_colour =    [ 173/255   255/255  47/255  ];
end

if  CURRENT_SEQU==5
    dd_colour = [ 173/255    255/255 47/255];
    c_colour =    [ 244/255   233/255  10/255  ];
end

if  CURRENT_SEQU==6
    dd_colour = [ 240/255    230/255 140/255];
    c_colour =    [ 255/255   165/255  5/255  ];
end

if  CURRENT_SEQU==7
    dd_colour = [ 210/255    180/255 140/255];
    c_colour =    [ 255/255   127/255  80/255  ];
end

if  CURRENT_SEQU==8
    dd_colour = [ 244/255    164/255 96/255];
    c_colour =    [ 255/255   99/255  71/255  ];
end

if  CURRENT_SEQU==9
    dd_colour = [ 135/255    206/255  250/255];
    c_colour =    [ 250/255   215/255  5/255  ];
end

if  CURRENT_SEQU==10
    dd_colour = [ 238/255    130/255  220/255];
    c_colour =    [ 65/255   105/255  225/255  ];
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  DIFFU_NAME == 'FA'
    p_values_before= -log10(p_values_before)*0.15;
    FDR_VALUE= -log10(FDR_VALUE)*0.15;
    x=[1:1:100];
    line([0,100],[-log10(0.05)*0.15,-log10(0.05)*0.15],'linestyle','--','color','g')
end

if  DIFFU_NAME == 'MD'
    p_values_before= -log10(p_values_before)*0.001*0.5*0.7;
    FDR_VALUE= -log10(FDR_VALUE)*0.001*0.5*0.7;
    x=[1:1:100];
    line([0,100],[-log10(0.05)*0.001*0.5*0.7,-log10(0.05)*0.001*0.5*0.7],'linestyle','--','color','g')
end

if  DIFFU_NAME == 'FC'
    p_values_before= -log10(p_values_before)*0.28;
    FDR_VALUE= -log10(FDR_VALUE)*0.28;
    x=[1:1:100];
    line([0,100],[-log10(0.05)*0.28,-log10(0.05)*0.28],'linestyle','--','color','g')
end






% dd=  bar(x,p_values_before,0.9)
% dd.FaceColor =dd_colour;
% hold on ;
c=  bar(x,FDR_VALUE,0.9)
c.FaceColor = c_colour;
set(c,'edgecolor','none');
% set(dd,'edgecolor','none');
% set(h, 'edgecolor', [0.8,0.8,0.8]);
% x_label_str = [hemi_str,'-',number_cluser,'-',DIFFU_NAME]
% xlabel([x_label_str,'/Node'],'Fontsize',25,'Fontname','Times New Roman');
% ylabel('-log10(Pvalues)','Fontsize',25,'Fontname','Times New Roman');

%  hold on
%   line([0,100],[-log10(0.01),-log10(0.01)],'linestyle','--','color','b')
%  
 
% h = legend('p-value' ,'p-value (MCC)','p-value = 0.05','p-value = 0.01');
% 
% 
% set(h,'Fontname','Times New Roman','Fontsize',10);
% hold on;
% % legend('TF Bundle','Location','SouthEast');  
%  %text(2,4.5,'TF Bundle')
% set(gca,'FontName','Times New Roman','FontSize',30);%  'FontWeight','bold'
%  
% set(gca,'XTick',[0:20:100]); 
% 
%   cluster_nam = tract_seg(number_cluser);

% title([cluster_nam,' Tract'],'fontname','Times New Roman','FontSize',25);
% %   backColor = [0 0 0];
% %   set(gca, 'color', backColor);
%  str1 = [cluster_nam,'_',hemi_str,'_',number_cluser,'_histogram_',DIFFU_NAME,'.bmp'];
%  mat_sav = [save_direc,str1]
% saveas(gcf,mat_sav);

ok = 1;
hold off
end