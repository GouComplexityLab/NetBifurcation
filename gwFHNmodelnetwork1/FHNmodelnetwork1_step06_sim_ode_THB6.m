clear all
close all
clc

addpath('./gw_NetFHNscripts/');

tic
start_time = datestr(now)

global L laplacian_matrix a b c d1 d2


%%
network_name = 'npbcws_n=10k=2p=0number=1';
datanetwork_path = ['./network1 ' network_name '/'];
results_path = strrep(datanetwork_path,'network1','FHN results network1');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

net_number = 1;
mat_name = [datanetwork_path 'adj_matrix_nwnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);
%%
load([results_path 'THBNormalForm'],'n')
k = n;
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

load([results_path 'THBNormalForm'],'a_star','c_star')
region_i = 6;
mu1=0.005; mu2=0.002;
a = a_star + mu1;
c = c_star + mu2;

%%
load([results_path 'THBNormalForm'],'du','dv')
d1 = du; d2 = dv;



%%
u_star = 0; v_star = 0;

%%
L = size(laplacian_matrix,1);
x = [1:L];


%% 
% T = 30; Time_interval = 1;
T = 2e4; Time_interval = 0.1;

% Initial condition 1
init_n = 1;
Initial_U = u_star*ones(1,L) + 0.01 + 0.01 * special_eigenvector';
Initial_V = v_star*ones(1,L) + 0.01 + 0.01 * special_eigenvector';

% % Initial condition 2
% init_n = 2; 
% Initial_U = u_star*ones(1,L) - 0.01 - 0.01 * special_eigenvector';
% Initial_V = v_star*ones(1,L) - 0.01 - 0.01 * special_eigenvector';

% % Initial condition 3
% init_n = 3;
% Initial_U = u_star*ones(1,L) + 0.12 + 0.2 * special_eigenvector';
% Initial_V = v_star*ones(1,L) + 0.12 + 0.2 * special_eigenvector';

% % Initial condition 4
% init_n = 4;
% Initial_U = u_star*ones(1,L) - 0.12 - 0.2 * special_eigenvector';
% Initial_V = v_star*ones(1,L) - 0.12 - 0.2 * special_eigenvector';

Y0 = [Initial_U'; Initial_V'];
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(2*L,1));

[time,solY1] = ode45(@NetworkFHNsystem,[0:Time_interval:T],Y0,options);

U_series = solY1(:,[1:L])'; V_series = solY1(:,[1+L:L+L])';

Final_U = U_series(:,end)'; Final_V = V_series(:,end)';

% plot and save 
plot_onesim2(time, x, mu1, mu2, U_series, V_series)

mat_pathname = [results_path 'data_THBFHN_region_i' num2str(region_i) '_Init' num2str(init_n) '.mat'];
save(mat_pathname,'region_i','init_n','x','mu1','mu2','time','U_series','V_series','Initial_U','Initial_V','Final_U','Final_V')

% pause(0.5)
% close all
start_time
end_time6 = datestr(now)

toc

