clc; clear all;
filePath = 'E:\SUTRA\Validation4-Solute variations\experiment data for experiment 1.xlsx';
data = readmatrix(filePath, 'Range', 'B2');
num_points_T = 6481;
y_temperature1_85 = data(1:end, 1);
y_temperature1_75 = data(1:end, 2);
y_temperature1_65 = data(1:end, 3);
y_temperature1_55 = data(1:end, 4);
y_moisture1_85 = data(1:end, 5);
y_moisture1_75 = data(1:end, 6);
y_moisture1_65 = data(1:end, 7);
y_moisture1_55 = data(1:end, 8);
y_conductivity1_85= data(1:end, 9);
y_conductivity1_75= data(1:end, 10);
y_conductivity1_65= data(1:end, 11);
y_conductivity1_55= data(1:end, 12);
time=6480;
name1='E:\SUTRA\Validation4-Solute variations\Column experiment 1-H=60cm';
path1=[name1,'\y_wd.FIL'];
path2=[name1,'\y_sl.FIL'];
path3=[name1,'\y_yd.FIL'];
fid = fopen(path1,'r');
p=0;
FormatString=repmat('%s ',1,102);
ppp=textscan(fid,FormatString,time,'HeaderLines',p);
for i = 1:time
    for w = 1:102
        wd(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
fid = fopen(path2,'r');
p=0;
FormatString=repmat('%s ',1,102);
ppp=textscan(fid,FormatString,time,'HeaderLines',p);
for i = 1:time
    for w = 1:102
        sl(i, w) = str2double(ppp{w}{i})*0.39;
    end
end
fclose(fid);
fid = fopen(path3,'r');
p=0;
FormatString=repmat('%s ',1,102);
ppp=textscan(fid,FormatString,time,'HeaderLines',p);
for i = 1:time
    for w = 1:102
        yd(i, w) = str2double(ppp{w}{i})*1000;
    end
end
fclose(fid);
tt=wd(:,1)/60;
wd_90=wd(:,102);wd_85=wd(:,97);wd_80=wd(:,92);
wd_75=wd(:,87);wd_70=wd(:,81);wd_65=wd(:,76);wd_60=wd(:,68);
wd_62=wd(:,70);wd_55=wd(:,64);wd_50=wd(:,58);wd_45=wd(:,52);
wd_40=wd(:,47);wd_30=wd(:,35);wd_10=wd(:,11);
sl_85=sl(:,97);sl_75=sl(:,87);sl_65=sl(:,76);sl_55=sl(:,64);
sl_62=sl(:,70);
yd_85=yd(:,97);yd_75=yd(:,87);yd_65=yd(:,76);yd_55=yd(:,64);
yd_62=yd(:,70);yd_60=yd(:,68);
sl = readmatrix('E:\SUTRA\Validation4-Solute variations\SHAW for experiment1.xlsx', 'Sheet', 'moisture');
wd = readmatrix('E:\SUTRA\Validation4-Solute variations\SHAW for experiment1.xlsx', 'Sheet', 'temp');
yd = readmatrix('E:\SUTRA\Validation4-Solute variations\SHAW for experiment1.xlsx', 'Sheet', 'conc');
shawt=linspace(1,108,108);
shawwd_85=wd(14:end,6);
shawwd_75=wd(14:end,11);
shawwd_65=wd(14:end,16);
shawwd_55=wd(14:end,21);
shawsl_85=sl(14:end,6);
shawsl_75=sl(14:end,11);
shawsl_65=sl(14:end,16);
shawsl_55=sl(14:end,21);
shawyd_85=yd(14:end,7);
shawyd_75=yd(14:end,12);
shawyd_65=yd(14:end,17);
shawyd_55=yd(14:end,22);
axes('unit','centimeters','position',[2 12 8 7]);
size1=8;
t2=linspace(0,108,num_points_T);
set (gcf,'unit','centimeters','position',[13,5,25,20]);
idx = [1:60:600, 601:150:3600, 3601:60:4200, 4201:150:6481];
plot(t2(idx),y_temperature1_85(idx),'ko', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx),y_temperature1_75(idx),'ro', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx),y_temperature1_65(idx),'bo', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx),y_temperature1_55(idx),'co', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(tt,wd_85,'k-','linewidth',1.5);hold on
plot(tt,wd_75,'r-','linewidth',1.5);hold on
plot(tt,wd_65,'b-','linewidth',1.5);hold on
plot(tt,wd_55,'c-','linewidth',1.5);hold on
plot(shawt,shawwd_85,'k--','linewidth',1.5);hold on
plot(shawt,shawwd_75,'r--','linewidth',1.5);hold on
plot(shawt,shawwd_65,'b--','linewidth',1.5);hold on
plot(shawt,shawwd_55,'c--','linewidth',1.5);hold on
xlabel('Elapsed time (h)');
ylabel('Temperature (\circC)');
title('{\bfa)} Temperature');
xlim([0, 108]);
set(gca,'XTick',[0 20 40 60 80 100], 'FontSize', 14) ;
axes('unit','centimeters','position',[13 12 8 7]);
bu2=0.012;
idx85=[1:10:7.5*60, 7.5*60+1:150:60*60, 60*60+1:150:70*60, 70*60+1:150:6481];
idx75=[1:150:12.5*60, 12.5*60+1:20:17.5*60, 17.5*60+1:150:70*60, 70*60+1:150:6481];
idx65=[1:150:30*60, 30*60+1:20:37.5*60, 37.5*60+1:150:70*60, 70*60+1:150:6481];
idx55=[1:150:30*60, 30*60+1:150:37.5*60, 37.5*60+1:150:70*60, 70*60+1:150:6481];
set (gcf,'unit','centimeters','position',[13,5,25,20]);
plot(t2(idx85),y_moisture1_85(idx85),'ko', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx75),y_moisture1_75(idx75),'ro', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx65),y_moisture1_65(idx65),'bo', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx55),y_moisture1_55(idx55),'co', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(tt,sl_85,'k-','linewidth',1.5);hold on
plot(tt,sl_75,'r-','linewidth',1.5);hold on
plot(tt,sl_65,'b-','linewidth',1.5);hold on
plot(tt,sl_55,'c-','linewidth',1.5);hold on
plot(shawt,shawsl_85,'k--','linewidth',1.5);hold on
plot(shawt,shawsl_75,'r--','linewidth',1.5);hold on
plot(shawt,shawsl_65,'b--','linewidth',1.5);hold on
plot(shawt,shawsl_55,'c--','linewidth',1.5);hold on
xlabel('Elapsed time (h)');
ylabel('Unfrozen water content (m^3\cdot{}m^{-3})');
title('{\bfb)} Unfrozen water content');
xlim([0, 108]);ylim([0.05, 0.45]);
set(gca,'XTick',[0 20 40 60 80 100], 'FontSize', 14) ;
set(gca,'yTick',[0.05 0.15 0.25 0.35 0.45], 'FontSize', 14) ;
axes('unit','centimeters','position',[2 2 8 7]);
bu3=1;
idx85=[1:150:6481];
idx75=[1:150:6481];
idx65=[1:150:6481];
idx55=[1:150:6481];
plot(t2(idx85),y_conductivity1_85(idx85),'ko', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx75),y_conductivity1_75(idx75),'ro', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx65),y_conductivity1_65(idx65),'bo', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(t2(idx55),y_conductivity1_55(idx55),'co', 'MarkerSize', size1,'linewidth',1.5);hold on
plot(tt,yd_85,'k-','linewidth',1.5);hold on
plot(tt,yd_75,'r-','linewidth',1.5);hold on
plot(tt,yd_65,'b-','linewidth',1.5);hold on
plot(tt,yd_55,'c-','linewidth',1.5);hold on
plot(shawt,shawyd_85,'k--','linewidth',1.5);hold on
plot(shawt,shawyd_75,'r--','linewidth',1.5);hold on
plot(shawt,shawyd_65,'b--','linewidth',1.5);hold on
plot(shawt,shawyd_55,'c--','linewidth',1.5);hold on
xlabel('Elapsed time (h)');
ylabel('Porewater salinity (ppt)');
title('{\bfc)} Porewater salinity');
xlim([0, 108]);
ylim([0, 50]);
set(gca,'XTick',[0 20 40 60 80 100], 'FontSize', 14) ;
legend_labels = {'85 cm (obs.)', '75 cm (obs.)', '65 cm (obs.)', '55 cm (obs.)', ...
                 '85 cm (SUTRA-MS-FT)', '75 cm (SUTRA-MS-FT)', '65 cm (SUTRA-MS-FT)', '55 cm (SUTRA-MS-FT)','85 cm (SHAW)', '75 cm (SHAW)', '65 cm (SHAW)', '55 cm (SHAW)'};
legend_colors = {'ko', 'ro', 'bo', 'co', 'k-', 'r-', 'b-', 'c-', 'k--', 'r--', 'b--', 'c--'};
ax4 = gca;
hold on;
for i = 1:length(legend_labels)
    plot(nan, nan, legend_colors{i}, 'MarkerSize', size1, 'linewidth', 1.5);
end
set (gcf,'unit','centimeters','position',[15 3 22 20]);
hLegend = legend(legend_labels, 'NumColumns', 1, 'FontSize', 12, 'Location', 'northwest');
legendPos = get(hLegend, 'Position');
set(hLegend, 'Position', [legendPos(1) + 0.48, legendPos(2)+0.01, legendPos(3), legendPos(4)]);