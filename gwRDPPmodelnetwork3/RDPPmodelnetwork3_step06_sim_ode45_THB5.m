clear all
close all
clc

addpath('./gw_NetPPscript/');

tic
start_time = datestr(now)

global L laplacian_matrix a b c d1 d2


%%
network_name = 'ws_n=10k=4p=0.05number=69';
datanetwork_path = ['./network3 ' network_name '/'];
results_path = strrep(datanetwork_path,'network3','RDPPmodel results network3');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

net_number = 69;
mat_name = [datanetwork_path 'adj_matrix_wsnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);

cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); 
disp(degree')

laplacian_matrix = adjacent_matrix - diag(degree);

laplacian_matrix = full(laplacian_matrix);

[eigenvectors, matrix_eigenvalues] = eig(laplacian_matrix);
eigenvalues = diag(matrix_eigenvalues);
eigenvectors = fliplr(eigenvectors);
eigenvalues = flipud(eigenvalues);
eigenvalues(1) = 0;
eigenvectors(:,1) = - eigenvectors(:,1);

load([results_path 'THBNormalForm'],'n')
k = n;
special_eigenvector = eigenvectors(:,k);
%%

load([results_path 'THBNormalForm'],'b_star','a_star')

domain_n = 5;
mu1 = -0.05; mu2 = -0.04795; % 5

b = b_star + mu1; 
a = a_star + mu2; 
c = 0.2;


%%
d1 = 0.001; d2 = 10*d1;
Times = 0.01;
d1 = d1 / Times;
d2 = d2 / Times;

%% 
u_star = b * (a + (c-1) * b) / a; 
v_star = b * (1 - c) / c * u_star;


%% 
L = size(laplacian_matrix,1);
x = [1:L];



%% 

% % Initial condition 
init_n = 1;
Initial_U = u_star*ones(size(x)) + 0.01*ones(size(x)) + 0.05 * special_eigenvector';
Initial_V = v_star*ones(size(x)) + 0.01*ones(size(x)) + 0.05 * special_eigenvector';

% % Initial condition 
init_n = 2;
Initial_U = u_star*ones(size(x)) - 0.01*ones(size(x)) - 0.05 * special_eigenvector';
Initial_V = v_star*ones(size(x)) - 0.01*ones(size(x)) - 0.05 * special_eigenvector';



Y0 = [Initial_U'; Initial_V'];
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(2*L,1));
T = 20e4; Time_interval = 1;
[time,solY1] = ode45(@NetworkRDsystem,[0:Time_interval:T],Y0,options);

U_series = solY1(:,[1:L])'; V_series = solY1(:,[1+L:L+L])';

Final_U = U_series(:,end)'; Final_V = V_series(:,end)';

% % plot and save 
% plot_onesim(time, x, mu1, mu2, U_series, V_series)

simulated_u1 = U_series';
simulated_v1 = V_series';

figure
subplot(2,1,1)
hold on
for i=1:size(laplacian_matrix,1)
    plot(time,simulated_u1(:,i));
end
xlabel('time t', 'FontSize',18, 'Interpreter','latex')
ylabel('$u_{i}$', 'FontSize',18, 'Interpreter','latex')

ax = gca;
% % xticks = xticks(ax1);
% xticks = [1.595:0.005:1.61];
% xticklabels = arrayfun(@(x) sprintf('%.3f', x), xticks, 'UniformOutput', false);
% % xticklabels(1) = {'0'};
% set(gca, 'xtick', xticks);
% set(gca, 'xticklabel', xticklabels);
% 
% % yticks = yticks(ax);
% yticks = [0.7:0.1:1.1];
% yticklabels = arrayfun(@(x) sprintf('%.1f', x), yticks, 'UniformOutput', false);
% % yticklabels(1) = {'0'};
% set(gca, 'ytick', yticks);
% set(gca, 'yticklabel', yticklabels);

ax_FontSize = 12;
ax.XAxis.FontSize = ax_FontSize;  
ax.XAxis.FontName = 'Times New Roman';  
ax.YAxis.FontSize = ax_FontSize;  
ax.YAxis.FontName = 'Times New Roman';  
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'in';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

subplot(2,1,2)
hold on
for i=1:size(laplacian_matrix,1)
    plot(time,simulated_v1(:,i));
end
xlabel('time t', 'FontSize',18, 'Interpreter','latex')
ylabel('$v_{i}$', 'FontSize',18, 'Interpreter','latex')

ax = gca;
% % xticks = xticks(ax1);
% xticks = [1.595:0.005:1.61];
% xticklabels = arrayfun(@(x) sprintf('%.3f', x), xticks, 'UniformOutput', false);
% % xticklabels(1) = {'0'};
% set(gca, 'xtick', xticks);
% set(gca, 'xticklabel', xticklabels);
% 
% % yticks = yticks(ax);
% yticks = [0.7:0.1:1.1];
% yticklabels = arrayfun(@(x) sprintf('%.1f', x), yticks, 'UniformOutput', false);
% % yticklabels(1) = {'0'};
% set(gca, 'ytick', yticks);
% set(gca, 'yticklabel', yticklabels);

% ax_FontSize = 18;
ax.XAxis.FontSize = ax_FontSize;  
ax.XAxis.FontName = 'Times New Roman';  
ax.YAxis.FontSize = ax_FontSize;  
ax.YAxis.FontName = 'Times New Roman';  
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'in';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

save([results_path 'data_ODE45net_THB_microdomain_nis' num2str(domain_n) '_Init' num2str(init_n)'],'domain_n','init_n','x','mu1','mu2','time','U_series','V_series','Initial_U','Initial_V','Final_U','Final_V')

start_time
end_time1 = datestr(now)



toc

