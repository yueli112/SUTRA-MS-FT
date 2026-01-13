clear
set (gcf,'unit','centimeters','position',[15 3 22 20]);
time=401;
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0';
name1=[name00,'\TH2_MINITEP_liquid.DAT'];
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,7);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);
for i = 1:time
    for w = 1:7
        mintem0(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.03';
name1=[name00,'\TH2_MINITEP_liquid.DAT'];
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,7);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);
for i = 1:time
    for w = 1:7
        mintem03(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.09';
name1=[name00,'\TH2_MINITEP_liquid.DAT'];
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,7);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);
for i = 1:time
    for w = 1:7
        mintem09(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.15';
name1=[name00,'\TH2_MINITEP_liquid.DAT'];
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,7);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);
for i = 1:time
    for w = 1:7
        mintem15(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
th2data = 'E:\SUTRA\Validation2-TH2 benchmaek case\TH2data.xlsx';
gr0 = readmatrix(th2data, 'Sheet', 2, 'Range','A2:B353'); 
t00=gr0(:,1);min00=gr0(:,2);
gr03 = readmatrix(th2data, 'Sheet', 3, 'Range', 'A2:B353'); 
t03=gr03(:,1);min03=gr03(:,2);
gr09 = readmatrix(th2data, 'Sheet', 4, 'Range','A2:B353'); 
t09=gr09(:,1);min09=gr09(:,2);
gr15 = readmatrix(th2data, 'Sheet', 5, 'Range','A2:B353'); 
t15=gr15(:,1);min15=gr15(:,2);
axes('unit','centimeters','position',[2 12 8 7]);
t=linspace(1,(time-1)*20*30,time);
width=1.5;
plot(t,mintem0(:,2),'k-', 'LineWidth', width);hold on
plot(t,mintem03(:,2),'r-', 'LineWidth', width);hold on
plot(t,mintem09(:,2),'b-', 'LineWidth', width);hold on
plot(t,mintem15(:,2),'g-', 'LineWidth', width);hold on
plot(t00,min00,'k--', 'LineWidth', width)
plot(t03,min03,'r--', 'LineWidth', width)
plot(t09,min09,'b--', 'LineWidth', width)
plot(t15,min15,'g--', 'LineWidth', width)
xlim([0,210000]);
ylim([-5,5]); 
set(gca,'YTick',[-5 -3 -1 1 3 5]) ;
set(gca, 'FontSize', 14);
xlabel('Time (s)','FontSize', 16);
ylabel('Temperature minimum (\circC)','FontSize', 16);
title('{\bfa) TH2-PM1}');
time=401;
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0';
name1=[name00,'\ALL_ENERGY_CHANGE_AND_BOUNDARY_FLUX.DAT'];
fid = fopen(name1,'r');
FormatString=repmat('%s ',1,16);
ppp=textscan(fid,FormatString,time,'HeaderLines',1);
for i = 1:time
    for w = 1:16
        flux0(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.03';
name1=[name00,'\ALL_ENERGY_CHANGE_AND_BOUNDARY_FLUX.DAT'];
fid = fopen(name1,'r');
FormatString=repmat('%s ',1,16);
ppp=textscan(fid,FormatString,time,'HeaderLines',1);
for i = 1:time
    for w = 1:16
        flux03(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.09';
name1=[name00,'\ALL_ENERGY_CHANGE_AND_BOUNDARY_FLUX.DAT'];
fid = fopen(name1,'r');
FormatString=repmat('%s ',1,16);
ppp=textscan(fid,FormatString,time,'HeaderLines',1);
for i = 1:time
    for w = 1:16
        flux09(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.15';
name1=[name00,'\ALL_ENERGY_CHANGE_AND_BOUNDARY_FLUX.DAT'];
fid = fopen(name1,'r');
FormatString=repmat('%s ',1,16);
ppp=textscan(fid,FormatString,time,'HeaderLines',1);
for i = 1:time
    for w = 1:16
        flux15(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
gr0 = readmatrix(th2data, 'Sheet', 7, 'Range','A2:B352'); 
t00=gr0(:,1);min00=gr0(:,2);
gr03 = readmatrix(th2data, 'Sheet', 8, 'Range', 'A2:B352'); 
t03=gr03(:,1);min03=gr03(:,2);
gr09 = readmatrix(th2data, 'Sheet', 9, 'Range','A2:B702'); 
t09=gr09(:,1);min09=gr09(:,2);
gr15 = readmatrix(th2data, 'Sheet', 10, 'Range','A2:B1502'); 
t15=gr15(:,1);min15=gr15(:,2);
axes('unit','centimeters','position',[13 12 8 7]);
t=linspace(1,(time-1)*20*30,time-1);
width=1.5;
plot(t,-flux0(2:end,16)*2,'k-', 'LineWidth', width);hold on
plot(t,-flux03(2:end,16)*2,'r-', 'LineWidth', width);hold on
plot(t,-flux09(2:end,16)*2,'b-', 'LineWidth', width);hold on
plot(t,-flux15(2:end,16)*2,'g-', 'LineWidth', width);hold on
plot(t00,min00,'k--', 'LineWidth', width)
plot(t03,min03,'r--', 'LineWidth', width)
plot(t09,min09,'b--', 'LineWidth', width)
plot(t15,min15,'g--', 'LineWidth', width)
xlim([0,210000]);
set(gca, 'FontSize', 14);
xlabel('Time (s)','FontSize', 16);
ylabel('Net flux (W/m^2)','FontSize', 16);
title('{\bfb) TH2-PM2}');
time=401;
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0';
name1=[name00,'\TH2_MINITEP_liquid.DAT'];
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,7);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);
for i = 1:time
    for w = 1:7
        water0(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.03';
name1=[name00,'\TH2_MINITEP_liquid.DAT'];
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,7);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);
for i = 1:time
    for w = 1:7
        water03(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.09';
name1=[name00,'\TH2_MINITEP_liquid.DAT'];
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,7);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);
for i = 1:time
    for w = 1:7
        water09(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
name00='E:\SUTRA\Validation2-TH2 benchmaek case\CaseTH2Frozen Inclusion gradients=0.15';
name1=[name00,'\TH2_MINITEP_liquid.DAT'];
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,7);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);
for i = 1:time
    for w = 1:7
        water15(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
gr0 = readmatrix(th2data, 'Sheet', 12, 'Range','A2:B353'); 
t00=gr0(:,1);min00=gr0(:,2);
gr03 = readmatrix(th2data, 'Sheet', 13, 'Range', 'A2:B353'); 
t03=gr03(:,1);min03=gr03(:,2);
gr09 = readmatrix(th2data, 'Sheet', 14, 'Range','A2:B703'); 
t09=gr09(:,1);min09=gr09(:,2);
gr15 = readmatrix(th2data, 'Sheet', 15, 'Range','A2:B1503'); 
t15=gr15(:,1);min15=gr15(:,2);
axes('unit','centimeters','position',[2 2 8 7]);
t=linspace(1,(time-1)*20*30,time);
width=1.5;
plot(t,water0(:,3)*0.37/11041*1.5*2,'k-', 'LineWidth', width);hold on
plot(t,water03(:,3)*0.37/11041*1.5*2,'r-', 'LineWidth', width);hold on
plot(t,water09(:,3)*0.37/11041*1.5*2,'b-', 'LineWidth', width);hold on
plot(t,water15(:,3)*0.37/11041*1.5*2,'g-', 'LineWidth', width);hold on
plot(t00,min00,'k--', 'LineWidth', width)
plot(t03,min03,'r--', 'LineWidth', width)
plot(t09,min09,'b--', 'LineWidth', width)
plot(t15,min15,'g--', 'LineWidth', width)
xlim([0,210000]);
ylim([1.07,1.12]); 
set(gca, 'FontSize', 14);
xlabel('Time (s)','FontSize', 16);
ylabel('Total water volume (m^3)','FontSize', 16);
title('{\bfc) TH2-PM3}');
h=legend('0% (Simulation)','3% (Simulation)','9% (Simulation)','15% (Simulation)','0% (Grenier et al.2018)','3% (Grenier et al.2018)','9% (Grenier et al.2018)','15% (Grenier et al.2018)','FontSize', 15);
lgd = legend('Location', 'southeast');