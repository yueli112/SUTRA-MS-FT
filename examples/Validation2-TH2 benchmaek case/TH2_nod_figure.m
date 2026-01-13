clear
number=6;
for i=1:number
    buchang(1)=0;  
    buchang(2)=1;  
    buchang(3)=10;  
    buchang(4)=20;  
    buchang(5)=30;  
    buchang(6)=40; 
end
for j=1:number
    name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.09';
    name1=[name00,'\fr.ele'];
    name2=[name00,'\fr.nod'];
    jj=buchang(j);
    weizhi1= [17,13,10,9.5];
    weizhi2=[0.12 0.1 0.65 0.8];
    xnods=181;ynods=61;
    dpi = 300;jiedian=xnods*ynods;yuansu=(xnods-1)*(ynods-1);
    tt=15;
    fid = fopen(name1,'r');
    q=157+(yuansu+5)*jj;
    FormatString=repmat('%s ',1,7);
    ppp=textscan(fid,FormatString,jiedian,'HeaderLines',q);
    for i = 1:yuansu
        for w = 1:7
            ele(i, w) = str2double(ppp{w}{i});
        end
    end
    fclose(fid);
    qq1=1080;
    elev = zeros(qq1,7); 
    for i = 1:qq1
        elev(i, :) = ele(i * 10-5, :);
    end
    sudu=120;
    matrixv= zeros(sudu,7);
    for i = 1:sudu/6-1
        matrixv(6*i-5:6*i, :) = elev(54*(i+1)-53:54*(i+1)-48, :);
    end
    xe = matrixv(:, 2);
    ye = matrixv(:, 3);
    u = matrixv(:, 4);
    v = matrixv(:, 5);
    ue=u*4e2;ve=v*4e2;
    fid = fopen(name2,'r');
    q=151+(jiedian+5)*jj;
    FormatString=repmat('%s ',1,12);
    qqq=textscan(fid,FormatString,jiedian,'HeaderLines',q);
    for i = 1:jiedian
        for w = 1:12
            matrix(i, w) = str2double(qqq{w}{i});
        end
    end
    fclose(fid);
    x=zeros(ynods,xnods);
    y=zeros(ynods,xnods);
    wd=zeros(ynods,xnods);yl=zeros(ynods,xnods);
    si=zeros(ynods,xnods);siout=zeros(ynods,xnods);
    sl=zeros(ynods,xnods);
    msalt=zeros(ynods,xnods);
    xx=matrix(:,2);yy=matrix(:,3);yl1=matrix(:,4);
    yd1=matrix(:,6);wd1=matrix(:,5);si1=matrix(:,9);sl1=matrix(:,8);siout1=matrix(:,10);msalt1=matrix(:,12);
    for i=1:xnods
        x(:,i)=xx(ynods*i-ynods+1:ynods*i);
        y(1:ynods,i)=yy(ynods*i-ynods+1:ynods*i);
        wd(1:ynods,i)=wd1(ynods*i-ynods+1:ynods*i);
        yl(1:ynods,i)=yl1(ynods*i-ynods+1:ynods*i)/9810;
        si(1:ynods,i)=si1(ynods*i-ynods+1:ynods*i);
        siout(1:ynods,i)=siout1(ynods*i-ynods+1:ynods*i);
        sl(1:ynods,i)=sl1(ynods*i-ynods+1:ynods*i);
    end
    totwater=sl+si+siout;totice=si+siout;
    wd_orig = wd;  
    wd_lower = wd_orig(end:-1:2, :);
    wd = [wd_lower; wd_orig];   
    x_lower = x(end:-1:2, :);
    x = [x_lower; x];   
    y_lower = -y(end:-1:2,:);
    y= [y_lower; y];  
    xe_orig = xe;  
    ye_orig = ye;  
    ue_orig = ue;  
    ve_orig = ve;  
    num_ve_lower = length(ye_orig) - 1;
    ye_lower = -ye_orig;
    ve_lower = -ve_orig;
    xe= [xe_orig; xe_orig];
    ye = [ye_lower; ye_orig];
    ue = [ue_orig; ue_orig];
    ve= [ve_lower; ve_orig];
    xx1=0;xx2=3;yy1=-0.5;yy2=0.5;
    kuan=12;gao=3;size1=18;size2=16;zuo=2;ap=1;
    if j==1
        set (gcf,'unit','centimeters','position',[10 0.5 15 25]);
        axes('unit','centimeters','position',  [zuo 20+ap kuan gao]);
        wd=5*ones(121,181);wd(41:81,51:71)=-5;
        [c2,h2]=contourf(x,y,wd,'linestyle','none','LevelList', -5:0.1:5);hold on
        clabel(c2,h2)
        [c3 ,h3]=contour(x,y,wd,[0 0],'w-', 'LineWidth', 1.5);
        set(gca, 'FontSize', 15); 
        colormap(gca,jet(20));
        quiver(xe, ye, ue, ve, 'AutoScale', 'off','Color', 'k');
        x0 = 2.7; y0 = 0.3;
        arrow_len = 0.0002 ;
        arrow_len_scaled = arrow_len * 4e2;
        quiver(x0, y0, arrow_len_scaled, 0, 'k', 'MaxHeadSize', 2, 'AutoScale', 'off');
        text(x0-0.7 , y0, '0.0002 m/s', 'FontSize', 14, 'Color', 'k');
        caxis([-5 5]);
        title('Temperature distribution','FontSize', size1);
        ylabel('{\itz} (m)','FontSize', size2);
        xlim([xx1,xx2]); ylim([yy1,yy2]);
        set(gca,'XTick',[], 'FontSize', size2) ;
        text(gca, 0,0.5, ' {\bfa)} {\itt} = 0 h', 'BackgroundColor', 'white',  'EdgeColor', 'none',  'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 15);
    end
    if j==2
        axes('unit','centimeters','position',  [zuo 16.5+ap kuan gao]);
        [c2,h2]=contourf(x,y,wd,'linestyle','none','LevelList', -5:0.1:5);hold on
        clabel(c2,h2)
        [c3 ,h3]=contour(x,y,wd,[0 0],'w-', 'LineWidth', 1.5);
        set(gca, 'FontSize', 15); 
        colormap(gca,jet(20));
        quiver(xe, ye, ue, ve, 'AutoScale', 'off','Color', 'k');
        caxis([-5 5]);
        ylabel('{\itz} (m)','FontSize', size2);
        xlim([xx1,xx2]); ylim([yy1,yy2]);
        set(gca,'XTick',[], 'FontSize', size2) ;
        text(gca, 0,0.5, ' {\bfb)} {\itt} = 1 h', 'BackgroundColor', 'white',  'EdgeColor', 'none',  'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 15);
    end
    if j==3
        axes('unit','centimeters','position',  [zuo 13+ap kuan gao]);
        [c2,h2]=contourf(x,y,wd,'linestyle','none','LevelList', -5:0.1:5);hold on
        clabel(c2,h2)
        [c3 ,h3]=contour(x,y,wd,[0 0],'w-', 'LineWidth', 1.5);
        set(gca, 'FontSize', 15); 
        colormap(gca,jet(20));
        quiver(xe, ye, ue, ve, 'AutoScale', 'off','Color', 'k');
        caxis([-5 5]);
        ylabel('{\itz} (m)','FontSize', size2);
        xlim([xx1,xx2]); ylim([yy1,yy2]);
        set(gca,'XTick',[], 'FontSize', size2) ;
        text(gca, 0,0.5, ' {\bfc)} {\itt} = 5 h', 'BackgroundColor', 'white',  'EdgeColor', 'none',  'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 15);
    end
    if j==4
        axes('unit','centimeters','position',  [zuo 9.5+ap kuan gao]);
        [c2,h2]=contourf(x,y,wd,'linestyle','none','LevelList', -5:0.1:5);hold on
        clabel(c2,h2)
        [c3 ,h3]=contour(x,y,wd,[0 0],'w-', 'LineWidth', 1.5);
        set(gca, 'FontSize', 15); 
        colormap(gca,jet(20));
        quiver(xe, ye, ue, ve, 'AutoScale', 'off','Color', 'k');
        caxis([-5 5]);
        ylabel('{\itz} (m)','FontSize', size2);
        xlim([xx1,xx2]); ylim([yy1,yy2]);
        set(gca,'XTick',[], 'FontSize', size2) ;
        text(gca, 0,0.5, ' {\bfd)} {\itt} = 10 h', 'BackgroundColor', 'white',  'EdgeColor', 'none',  'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 15);
    end
    if j==5
        axes('unit','centimeters','position',  [zuo 6+ap kuan gao]);
        [c2,h2]=contourf(x,y,wd,'linestyle','none','LevelList', -5:0.1:5);hold on
        clabel(c2,h2)
        [c3 ,h3]=contour(x,y,wd,[0 0],'w-', 'LineWidth', 1.5);
        set(gca, 'FontSize', 15); 
        colormap(gca,jet(20));
        quiver(xe, ye, ue, ve, 'AutoScale', 'off','Color', 'k');
        caxis([-5 5]);
        ylabel('{\itz} (m)','FontSize', size2);
        xlim([xx1,xx2]); ylim([yy1,yy2]);
        set(gca,'XTick',[], 'FontSize', size2) ;
        text(gca, 0,0.5, ' {\bfe)} {\itt} = 15 h', 'BackgroundColor', 'white',  'EdgeColor', 'none',  'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 15);
    end
    if j==6
        axes('unit','centimeters','position',  [zuo 2.5+ap kuan gao]);
        [c2,h2]=contourf(x,y,wd,'linestyle','none','LevelList', -5:0.1:5);hold on
        clabel(c2,h2)
        [c3 ,h3]=contour(x,y,wd,[0 0],'w-', 'LineWidth', 1.5);
        set(gca, 'FontSize', 15); 
        colormap(gca,jet(20)); 
        quiver(xe, ye, ue, ve, 'AutoScale', 'off','Color', 'k');
        caxis([-5 5]);
        xlabel('{\itx} (m)','FontSize', size2);
        ylabel('{\itz} (m)','FontSize', size2);
        xlim([xx1,xx2]); ylim([yy1,yy2]);
        text(gca, 0,0.5, ' {\bff)} {\itt} = 20 h', 'BackgroundColor', 'white',  'EdgeColor', 'none',  'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 15);
        axes('unit','centimeters','position',  [zuo 0.4 kuan 6]);
        axis off;
        colorbarHandle=colorbar('horiz','YTick', -5:1:5, 'FontSize', 15);
        colormap(gca,jet(20));
        title(colorbarHandle, '\circC','position',  [350 10 ]);
        caxis([-5 5]);
    end
    j
end