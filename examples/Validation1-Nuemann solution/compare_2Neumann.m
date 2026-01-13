clear
name00='E:\SUTRA\Validation1-Nuemann solution\T_i=-0.001℃';
save_path = [name00,'\']; % Replace with your save path
name1=[name00,'\ICENOD.DAT'];
time=480;  % 240
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,2002);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);

for i = 1:time
    for w = 1:2002
        % Convert string to double using str2double; assign 0 to the position if conversion fails
        si0(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);

si=si0(:,2:2002);
t0=si0(:,1)/2/86400;
result_cols = zeros(time, 1);

for i = 1:time
    % Iterate through each column of the current row
    for j = 1:2001
        % Use eps to handle floating-point comparisons to avoid precision issues
        not_0 = si(i, j)>0;
        if not_0
            result_cols(i) = j;  % Record the column number
            break;  % Break the inner loop once found
        end
        if j == 2001
            result_cols(i) = 0;  % Use 0 to indicate that the entire row is 0.9999
            fprintf('Warning: All elements in row %d are 0.9999\n', i);
        end
    end
end
sim0=(2001*ones(time,1)-result_cols)/1000;

name00='E:\SUTRA\Validation1-Nuemann solution\T_i=-5℃';
save_path = [name00,'\']; % Replace with your save path
name1=[name00,'\ICENOD.DAT'];
time=480;  % 240
fid = fopen(name1,'r');
q=0;
FormatString=repmat('%s ',1,2002);
ppp=textscan(fid,FormatString,time,'HeaderLines',q);

for i = 1:time
    for w = 1:2002
        % Convert string to double using str2double; assign 0 to the position if conversion fails
        si0(i, w) = str2double(ppp{w}{i});
    end
end
fclose(fid);

si=si0(:,2:2002);
t0=si0(:,1)/2/86400;
result_cols = zeros(time, 1);

for i = 1:time
    % Iterate through each column of the current row
    for j = 1:2001
        % Use eps to handle floating-point comparisons to avoid precision issues
        not_0 = si(i, j)>0;
        if not_0
            result_cols(i) = j;  % Record the column number
            break;  % Break the inner loop once found
        end
        if j == 2001
            result_cols(i) = 0;  % Use 0 to indicate that the entire row is 0.9999
            fprintf('Warning: All elements in row %d are 0.9999\n', i);
        end
    end
end
sim5=(2001*ones(time,1)-result_cols)/1000;

% Analytical solutions calculation
t=linspace(0,time/24/2,time);
sol0=3.2671e-4*(t*3600*24).^0.5;
sol5=2.85e-4*(t*3600*24).^0.5;

% Plotting
plot(t0,sim0,'k--', 'LineWidth', 1.5); hold on
plot(t,sol0,'r-', 'LineWidth', 1.5);
plot(t0,sim5,'k--', 'LineWidth', 1.5); hold on
plot(t,sol5,'r-', 'LineWidth', 1.5);

xlim([0,10]);
ylim([0,0.35]); 
%set(gca,'YTick',[-5 -3 -1 1 3 5]) ;
set(gca, 'FontSize', 14);
xlabel('Time (day)','FontSize', 16);
ylabel('Depths to thawing front (m)','FontSize', 16);
legend('Simulation','Neumann solution','FontSize', 15);
text(7, 0.2, ' {\itT_i} = -5\circC', 'FontSize', 15);
text(4, 0.27, ' {\itT_i} = -0.001\circC', 'FontSize', 15);

%print(fullfile('F:\Neumann solution','Neumann_2compare.png'), '-dpng',['-r', num2str(300)]);