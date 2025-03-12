clear all
close all
clc

addpath('./gw_NetFHNscripts/');

tic
start_time = datestr(now)

global L laplacian_matrix a b c d1 d2


network_name = 'brain';
datanetwork_path = ['./network0 ' network_name '/'];
results_path = strrep(datanetwork_path,'network0','FHN results network0');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

mat_name = [datanetwork_path 'adj_matrix_brainnetwork'];
load(mat_name)
adjacent_matrix = full(adj_matrix);

%%
load([results_path 'SSBNormalForm'],'k')


%%
cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); 
disp(degree')

laplacian_matrix = adjacent_matrix - diag(degree);

N = size(laplacian_matrix,1);
[eigenvectors, eigenvalues] = eig(laplacian_matrix);
eigenvalues = diag(eigenvalues);
eigenvectors = fliplr(eigenvectors);
eigenvalues = flipud(eigenvalues);
eigenvalues(1) = 0;
eigenvectors(:,1) = - eigenvectors(:,1);

special_eigenvector = eigenvectors(:,k);

% figure
% hold on
% plot(special_eigenvector)
%% 

load([results_path 'SSBNormalForm'],'a_star','c')
para_n = 1; 
mu_n = 1; mu=0.002;
a = a_star + mu;

%% 
load([results_path 'SSBNormalForm'],'du','dv')
d1 = du; d2 = dv;



%% 
u_star = 0; v_star = 0;

%% 
L = size(laplacian_matrix,1);
x = [1:L];


%% 
% T = 30; Time_interval = 1;
T = 1e4; Time_interval = 0.1;

% Initial condition 1
init_n = 1;
Initial_U = u_star*ones(1,L) + 0.05 * special_eigenvector';
Initial_V = v_star*ones(1,L) + 0.05 * special_eigenvector';

% % Initial condition 1
% init_n = 1;
% Initial_U = u_star*ones(1,L) + 0.01 * special_eigenvector';
% Initial_V = v_star*ones(1,L) + 0.01 * 10 * special_eigenvector';

% load datastage1_ODE45delta_T=0dot001_THB_pis0dot05_microdomain_nis3_Init2 Final_U Final_V
% Initial_U = Final_U;
% Initial_V = Final_V;

Y0 = [Initial_U'; Initial_V'];
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(2*L,1));

[time,solY1] = ode45(@NetworkFHNsystem,[0:Time_interval:T],Y0,options);

U_series = solY1(:,[1:L])'; V_series = solY1(:,[1+L:L+L])';

Final_U = U_series(:,end)'; Final_V = V_series(:,end)';

% plot and save 
% plot_oneSSBsim(time, x, mu, U_series, V_series)
plot_oneSSBsim2(time, x, mu, U_series, V_series)

mat_pathname = [results_path 'data_SSBFHN_para' num2str(para_n) '_mu_nis' num2str(mu_n) '_Init' num2str(init_n) '.mat'];
save(mat_pathname,'para_n', 'mu_n', 'init_n', 'x', 'mu', 'time', 'U_series', 'V_series', 'Initial_U', 'Initial_V', 'Final_U', 'Final_V')

% pause(0.5)
% close all
start_time
end_time1 = datestr(now)

