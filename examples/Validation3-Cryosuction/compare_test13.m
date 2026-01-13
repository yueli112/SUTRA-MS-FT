clear
name00='E:\SUTRA\Validation3-Cryosuction\Test 13';
name1=[name00,'\cryosuction_all2.DAT'];
time=72;
nods=6002;
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,nods+2);
ppp=textscan(fid,FormatString,time*5,'HeaderLines',q);
for i = 1:time*5
    for w = 1:nods+2
        all(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);
q=3;p=4;
tem=(all(1:5:end,q:2:end)+all(1:5:end,p:2:end))/2;
pre=(all(2:5:end,q:2:end)+all(2:5:end,p:2:end))/2;
sii=(all(3:5:end,q:2:end)+all(3:5:end,p:2:end))/2;
sll=(all(4:5:end,q:2:end)+all(4:5:end,p:2:end))/2;
siout=(all(5:5:end,q:2:end)+all(5:5:end,p:2:end))/2;
sw=sii+sll;totwater=sw+siout;
h=linspace(30,0,3001);
totwater(:,1)=totwater(:,5);totwater(:,2)=totwater(:,5);totwater(:,3)=totwater(:,5);
t6 = [-5.251, -1.391, 1.458, 4.262, 6.823, 8.806, 11.192, 13.262, 15.159, 16.826, 18.580]';
th6 = [2.513, 5.057, 7.601, 10.145, 12.575, 15.083, 17.585, 20.093, 22.403, 25.060, 27.490]';
t12 = [-7.320, -5.409, -3.267, -1.314, 1.153, 3.986, 6.316, 8.806, 11.436, 14.455, 17.229]';
th12 = [2.513, 5.100, 7.523, 9.989, 12.341, 15.041, 17.507, 20.051, 22.481, 25.025, 27.490]';
t72 = [-8.197, -6.760, -4.618, -2.707, -1.356, -0.005, 3.488, 6.578, 9.755, 13.334, 16.280]';
th72 = [2.669, 5.135, 7.488, 10.109, 12.497, 15.083, 17.663, 20.093, 22.516, 25.145, 27.568]';
w6 = [20.842, 20.121, 19.766, 20.121, 19.949, 21.132, 18.345, 12.514, 12.933, 13.590, 13.407, 13.342, 13.880, 14.300, 14.418, 13.644, 14.117, 14.655, 14.956, 14.655, 14.838, 14.892, 14.838, 14.956, 14.773, 14.773, 14.300, 14.773, 15.193]';
wh6 = [-0.266, 0.746, 1.807, 2.628, 3.874, 4.815, 5.834, 6.924, 7.943, 8.997, 9.896, 10.951, 11.970, 12.989, 13.966, 14.900, 16.039, 16.895, 18.035, 18.969, 19.946, 21.120, 21.984, 23.194, 24.050, 25.190, 26.046, 28.119, 27.178]';
w12 = [21.498, 20.842, 20.960, 20.487, 20.960, 21.670, 21.132, 21.907, 23.155, 22.090, 12.460, 12.277, 12.815, 13.052, 12.815, 12.815, 13.224, 13.342, 13.826, 13.407, 13.590, 13.826, 13.407, 13.826, 13.471, 13.289, 13.471, 13.170, 13.471]';
wh12 = [0.123, 0.944, 2.076, 3.095, 3.952, 4.971, 5.990, 7.080, 8.177, 8.997, 9.939, 11.071, 11.857, 12.989, 14.001, 14.900, 16.039, 16.973, 17.914, 18.856, 19.910, 21.120, 22.019, 22.918, 24.128, 24.949, 26.088, 27.065, 28.084]';
w72 = [20.842, 19.830, 19.647, 20.487, 20.002, 20.777, 21.132, 20.960, 21.788, 22.326, 23.811, 27.620, 31.536, 19.109, 9.006, 9.071, 9.071, 9.598, 9.189, 8.705, 9.663, 9.308, 9.899, 9.598, 9.243, 9.071, 9.243, 8.651]';
wh72 = [0.080, 0.944, 1.998, 2.862, 4.383, 5.013, 6.301, 7.157, 8.063, 9.196, 10.130, 12.437, 13.222, 14.079, 14.942, 15.806, 17.016, 17.992, 18.891, 20.023, 21.000, 21.906, 22.996, 24.291, 24.914, 25.968, 28.162, 27.100]';
size1=8;
set (gcf,'unit','centimeters','position',[13,2,21,22]);
axes('unit','centimeters','position',[2 12 18 9]);
plot(th6, t6, 'ko', 'MarkerSize', size1, 'linewidth', 1.5); hold on;
plot(th12, t12, 'bo', 'MarkerSize', size1, 'linewidth', 1.5); hold on;
plot(th72, t72, 'go', 'MarkerSize', size1, 'linewidth', 1.5); hold on;
plot(h, tem(6,:), 'k-', 'linewidth', 1.5); hold on;
plot(h, tem(12,:), 'b-', 'linewidth', 1.5); hold on;
plot(h, tem(72,:), 'g-', 'linewidth', 1.5); hold on;
ylim([-15 20]);
set(gca, 'YTick', [-15 -10 -5 0 5 10 15 20], 'FontSize', 14);
set(gca, 'xTick', [ ], 'FontSize', 14);
ylabel('Temperature (\circC)');
title('{\bfa)} Temperature');
legend('t=6h (Jame 1977)', 't=24h (Jame 1977)', 't=72h (Jame 1977)', ...
       't=6h (Simulation)', 't=24h (Simulation)', 't=72h (Simulation)');
axes('unit','centimeters','position',[2 2 18 9]);
plot(wh6, w6/100*2.71, 'ko', 'MarkerSize', size1, 'linewidth', 1.5); hold on;
plot(wh12, w12/100*2.71, 'bo', 'MarkerSize', size1, 'linewidth', 1.5); hold on;
plot(wh72, w72/100*2.71, 'go', 'MarkerSize', size1, 'linewidth', 1.5); hold on;
plot(h, totwater(6,:), 'k-', 'linewidth', 1.5); hold on;
plot(h, totwater(12,:), 'b-', 'linewidth', 1.5); hold on;
plot(h, totwater(72,:), 'g-', 'linewidth', 1.5); hold on;
plot(h, 0.415*ones(3001,1), 'k--', 'linewidth', 1); hold on;
xlim([0 30]); ylim([0 1]);
set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1], 'FontSize', 14);
xlabel('Distance from cold end (cm)');
ylabel('Total water saturation (-)');
title('{\bfb)} Total water saturation');
legend('t=6h (Jame 1977)','t=24h (Jame 1977)','t=72h (Jame 1977)',...
       't=6h (Simulation)','t=24h (Simulation)','t=72h (Simulation)');
%print(fullfile('E:\SUTRA\论文freeze-thaw processes induce convective fingering\指流一审\模型验证\验证吸力\','test13_new2.png'), '-dpng',['-r', num2str(300)]);
