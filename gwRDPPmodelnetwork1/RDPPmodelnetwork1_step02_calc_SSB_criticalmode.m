% clear all
close all
clc

tic
digits = 100;
format long

%%
perturbating_mu =  0.0002;
maxmu = 0.0004;
max_y = 0.05;
qu_n = 6;

a = 1.0; 
c = 0.2;
d1 = 0.001; d2 = 10*d1;
Times = 0.01;

network_name = 'npbcws_n=10k=2p=0number=1';
datanetwork_path = ['./network1 ' network_name '/'];
results_path = strrep(datanetwork_path,'network1','RDPPmodel results network1');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

net_number = 1;
mat_name = [datanetwork_path 'adj_matrix_nwnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);

cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); 
disp(degree')

laplacian_matrix = adjacent_matrix - diag(degree);

%%
[eigenvectors, matrix_eigenvalues] = eig(laplacian_matrix);
eigenvalues = diag(matrix_eigenvalues);
eigenvectors = fliplr(eigenvectors);
eigenvalues = flipud(eigenvalues);
eigenvalues(1) = 0;
eigenvectors(:,1) = - eigenvectors(:,1);


L = size(laplacian_matrix,1);
d1 = d1 / Times;
d2 = d2 / Times;

%%
network_bk = [];
network_ak = [];
N = size(laplacian_matrix,1);
for k=1:N-1
    one_bk = d1/d2 - (eigenvalues(k)+eigenvalues(k+1))/c/(1-c)*d1 + eigenvalues(k)*eigenvalues(k+1)/c/c/(1-c)/(1-c)*d1*d2;
    network_bk = [network_bk,one_bk];
    
    s_1k = (-d2*(1-c^2)*eigenvalues(k)+c*(1-c)^2) / (-d2*eigenvalues(k)+c*(1-c));
    s_2k = - (d1*d2*eigenvalues(k)^2-d1*c*(1-c)*eigenvalues(k)) / (-d2*eigenvalues(k)+c*(1-c));
    one_ak = s_1k*one_bk + s_2k;
    network_ak = [network_ak,one_ak];

end

postion = find(network_ak>a);
k = postion(1);

%%
% k=19;
b_star = (-d2*eigenvalues(k)+c*(1-c)) / (-d2*(1-c^2)*eigenvalues(k)+c*(1-c)^2) *a ...
    + (d1*d2*eigenvalues(k)^2-d1*c*(1-c)*eigenvalues(k))/(-d2*(1-c^2)*eigenvalues(k)+c*(1-c)^2);
b_star_s=sprintf('%20.10f',b_star);
b_star_s=b_star;

%%
K11 = b_star*(1-c^2) - a;
K12 = -c^2;
K21 = b_star*(1-c)^2;
K22 = c*(c-1);
K_local = [K11,K12;K21,K22];
D_diffusion = [d1,0;0,d2];
M_large = kron(diag(ones(1,10)),K_local) + kron(laplacian_matrix, D_diffusion);
MK = K_local + eigenvalues(k)*D_diffusion;

%%
[MK_eigenvectors, MK_eigenvalues] = eig(MK);
MK_eigenvalues = diag(MK_eigenvalues);
[v index] = sort(real(MK_eigenvalues));
MK_eigenvalues = MK_eigenvalues(index);
MK_eigenvectors = MK_eigenvectors(:,index);

disp(MK_eigenvalues')
disp(MK_eigenvalues(end))
MK_q = -MK_eigenvectors(:,end);

[transposed_MK_eigenvectors, transposed_MK_eigenvalues] = eig(MK');
transposed_MK_eigenvalues = diag(transposed_MK_eigenvalues);
[v index] = sort(real(transposed_MK_eigenvalues));
transposed_MK_eigenvalues = transposed_MK_eigenvalues(index);
transposed_MK_eigenvectors = transposed_MK_eigenvectors(:,index);
% disp(transposed_J_at_point_eigenvalues')
disp(transposed_MK_eigenvalues(end))
MK_q_star = transposed_MK_eigenvectors(:,end);


MK_q = MK_q ./ (MK_q(1));
MK_q_star = MK_q_star ./ sum(MK_q.*MK_q_star);
% k=10;
phi = eigenvectors(:,k);
psi = eigenvectors(:,k);

figure
set(gcf,"Position",[300 300 500 320])
set(gca,'Position',[0.146 0.129047619047619 0.814 0.745]);
% axes1 = axes('Position',[0.16625 0.21375 0.75 0.75]);
% set(gcf,'Position',get(0,'ScreenSize'))
hold on
box on

plot([1:10],phi,'k--')
plot([1:10],zeros(1,10),'k:')

% gwcolors = lines(10);
gwcolors = hsv(11);
gwcolors(3,:) = [];

for i=1:10
    plot(i,phi(i),'MarkerSize',18,'Marker','.','LineStyle','none',...
    'Color',gwcolors(i,:));
end

xlim([1 10])
ylim([-0.8 0.8])

ax = gca;
% xticks = xticks(ax);
xticks = [1,3,5,8,10];
xticklabels = arrayfun(@(x) sprintf('%.0f' , x), xticks, 'UniformOutput', false);
% xticklabels(1) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-0.5,0, 0.5];
yticklabels = arrayfun(@(x) sprintf('%.2f', x), yticks, 'UniformOutput', false);
yticklabels(find(yticks==0)) = {'0'};
% yticklabels(end) = {'$a$'};
set(gca, 'ytick', yticks);
set(gca, 'yticklabel', yticklabels);

ax_FontSize = 24;
ax.XAxis.FontSize = ax_FontSize;  
ax.XAxis.FontName = 'Times New Roman';  
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickLabelInterpreter = 'latex';


ax.YAxis.FontSize = ax_FontSize;  
ax.YAxis.FontName = 'Times New Roman';  
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickLabelInterpreter = 'latex';


set(gca,'Box','on');
% set(gca,'Box','on','BoxStyle','full');
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'FontSize',24,'linewidth',1.0,'layer','top');


figure_name = [results_path 'Figure Step02 SSBNF critical mode npbcwsnetworks_n=10k=2p=0number=1.eps'];
saveas(gcf, figure_name, 'epsc');
