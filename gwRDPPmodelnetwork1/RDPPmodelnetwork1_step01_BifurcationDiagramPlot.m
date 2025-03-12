clear all
close all
clc

format long
digits = 50;

%%

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

%%
cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); 
disp(degree')

laplacian_matrix = adjacent_matrix - diag(degree);

[eigenvectors, matrix_eigenvalues] = eig(laplacian_matrix);
eigenvalues = diag(matrix_eigenvalues);
eigenvectors = fliplr(eigenvectors);
eigenvalues = flipud(eigenvalues);
eigenvalues(1) = 0;
eigenvectors(:,1) = - eigenvectors(:,1);

%%
c = 0.2;
d1 = 0.001; d2 = 10*d1;
Times = 0.01;

N = length(eigenvalues);
d1 = d1 / Times;
d2 = d2 / Times;

%%
if d1 < c*(1-c)*d2/(c*(1-c) - eigenvalues(2)*d2)
    disp('Network case: Turing instability occurs.')
else
    disp('Network case: Turing instability dose not occur.')
end

%%
[network_v network_p]=min(abs((eigenvalues - ((d1-d2)*c*(1-c))/(2*d1*d2))));
% [network_v network_p]

n = network_p;
network_b_star = - d1*d2/(c^2*(1-c)^2)*eigenvalues(n)^2 + (d1-d2)/c/(1-c)*eigenvalues(n) + 1;
network_a_star = (1-c^2)*network_b_star + c*(c-1);

format rat
network_b_star;
network_a_star;
network_J_n = d1*d2*eigenvalues(n)^2 + (d1*c*(c-1) + d2*(network_b_star*(1-c^2) - network_a_star))*eigenvalues(n) + c*(1-c)*(network_a_star + network_b_star*(c-1));
network_T_1 = -(d1+d2)*eigenvalues(1) - network_b_star * (1-c^2) + network_a_star + c*(1-c);
fprintf('Network Turing Hopf condition: network_J_n=%20.20f, network_T_1=%20.20f.\n', network_J_n, network_T_1)
format short

if network_a_star > network_b_star*(1-c)
    disp('there exists the positive equilibrium in Network.')
else
    disp('there exists no positive equilibrium in Network. Error for Network!')
end

% connecting points
network_bk = 0;
for k=1:N-1
    one_bk = d1/d2 - (eigenvalues(k)+eigenvalues(k+1))/c/(1-c)*d1 + eigenvalues(k)*eigenvalues(k+1)/c/c/(1-c)/(1-c)*d1*d2;
    network_bk = [network_bk,one_bk];
end

%%
Set_M = 50;
intersect_n = [];
for k=1:N
    
    b_HT1 = [];
    for j=1:N
    %     one_b_HT = - d1*d2/c^2/(1-c)^2 * eigenvalues(j)^2 + (d1 - d2)/c/(1-c) * eigenvalues(j) + 1;
        one_b_HT1 = (-d1*d2*eigenvalues(j)^2 + (d1-d2)*c*(1-c)*eigenvalues(j) - (d1+d2)*c*(1-c)*eigenvalues(k) ...
                      + d2*(d1+d2)*eigenvalues(j)*eigenvalues(k) + c^2*(1-c)^2) / c^2/(1-c)^2;
        b_HT1 = [b_HT1,one_b_HT1];
    end
    [max_b_HT intersect_n] = max(b_HT1);
    intersect_n = [intersect_n,intersect_n];

end

%%
choosed_k0 = 1;
b_HT0 = [];
for j=1:N
%     one_b_HT0 = - d1*d2/c^2/(1-c)^2 * eigenvalues(j)^2 + (d1 - d2)/c/(1-c) * eigenvalues(j) + 1;
    one_b_HT0 = (-d1*d2*eigenvalues(j)^2 + (d1-d2)*c*(1-c)*eigenvalues(j) - (d1+d2)*c*(1-c)*eigenvalues(choosed_k0) ...
                  + d2*(d1+d2)*eigenvalues(j)*eigenvalues(choosed_k0) + c^2*(1-c)^2) / c^2/(1-c)^2;
    b_HT0 = [b_HT0,one_b_HT0];
end
[intersect_kn_b0 intersect_kn0] = max(b_HT0);
intersect_kn_a0 = (1 - c^2) * intersect_kn_b0 - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k0);

%%
choosed_k1 = 2;
b_HT1 = [];
for j=1:N
%     one_b_HT1 = - d1*d2/c^2/(1-c)^2 * eigenvalues(j)^2 + (d1 - d2)/c/(1-c) * eigenvalues(j) + 1;
    one_b_HT1 = (-d1*d2*eigenvalues(j)^2 + (d1-d2)*c*(1-c)*eigenvalues(j) - (d1+d2)*c*(1-c)*eigenvalues(choosed_k1) ...
                  + d2*(d1+d2)*eigenvalues(j)*eigenvalues(choosed_k1) + c^2*(1-c)^2) / c^2/(1-c)^2;
    b_HT1 = [b_HT1,one_b_HT1];
end
[intersect_kn_b1 intersect_kn1] = max(b_HT1);
intersect_kn_a1 = (1 - c^2) * intersect_kn_b1 - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k1);

%%
choosed_k2 = 3;
b_HT2 = [];
for j=1:N
%     one_b_HT2 = - d1*d2/c^2/(1-c)^2 * eigenvalues(j)^2 + (d1 - d2)/c/(1-c) * eigenvalues(j) + 1;
    one_b_HT2 = (-d1*d2*eigenvalues(j)^2 + (d1-d2)*c*(1-c)*eigenvalues(j) - (d1+d2)*c*(1-c)*eigenvalues(choosed_k2) ...
                  + d2*(d1+d2)*eigenvalues(j)*eigenvalues(choosed_k2) + c^2*(1-c)^2) / c^2/(1-c)^2;
    b_HT2 = [b_HT2,one_b_HT2];
end
[intersect_kn_b2 intersect_kn2] = max(b_HT2);
intersect_kn_a2 = (1 - c^2) * intersect_kn_b2 - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k2);




%%
% M=n;
intersect_kn1 = 8;
M=intersect_kn1+1;
% M=8;

if M == N
    max_b = network_bk(M)*2;
else
    max_b = network_bk(M+1);
end
b_value = linspace(0,max_b,100);
T_a_values = zeros(M, length(b_value));

H_a_values = zeros(M, length(b_value));

points_number = 100;
segment_b_values = zeros(M, points_number);
segment_a_values = zeros(M, points_number);

s_1k = [];
s_2k = [];
connecting_b = 0;
connecting_a = 0;

for k=1:M-1
%     i
    one_s_1k = (-d2*(1-c^2)*eigenvalues(k)+c*(1-c)^2) / (-d2*eigenvalues(k)+c*(1-c));
    one_s_2k = - (d1*d2*eigenvalues(k)^2-d1*c*(1-c)*eigenvalues(k)) / (-d2*eigenvalues(k)+c*(1-c));
    
    T_a_values(k,:) = one_s_1k * b_value + one_s_2k;
    H_a_values(k,:) = (1 - c^2) * b_value - c*(1-c) + (d1 + d2)*eigenvalues(k);
    
    connecting_b = [connecting_b,network_bk(k+1)];
    one_connecting_a = one_s_1k * network_bk(k+1) + one_s_2k;
    connecting_a = [connecting_a,one_connecting_a];
    
    segment_b_values(k,:) = linspace(network_bk(k),network_bk(k+1),points_number);
    segment_a_values(k,:) = one_s_1k * segment_b_values(k,:) + one_s_2k;
    
    s_1k = [s_1k,one_s_1k];
    s_2k = [s_2k,one_s_2k];
end

intersect_1n_b = - d1*d2/c^2/(1-c)^2 * eigenvalues(n)^2 + (d1 - d2)/c/(1-c) * eigenvalues(n) + 1;
intersect_1n_a = (1 - c^2) * intersect_1n_b - c*(1-c) + (d1 + d2)*eigenvalues(1);

% linspace(network_bk(k),network_bk(k+1),points_number);
segment_H1_a_values = (1 - c^2) * segment_b_values(n,:) - c*(1-c) + (d1 + d2)*eigenvalues(1);
segment_b_values_to_end_1n = linspace(intersect_1n_b,max_b,points_number);
segment_H1_a_values_to_end = (1 - c^2) * segment_b_values_to_end_1n - c*(1-c) + (d1 + d2)*eigenvalues(1);

%%
intersect_11_b = - d1*d2/c^2/(1-c)^2 * eigenvalues(1)^2 + (d1 - d2)/c/(1-c) * eigenvalues(1) + 1;
intersect_11_a = (1 - c^2) * intersect_11_b - c*(1-c) + (d1 + d2)*eigenvalues(1);
intersect_11_b_to_intersect_1n_b = linspace(intersect_11_b,intersect_1n_b,points_number);
intersect_11_a_to_intersect_1n_a = (1 - c^2) * intersect_11_b_to_intersect_1n_b - c*(1-c) + (d1 + d2)*eigenvalues(1);

%%
segment_Hk_a_values1 = (1 - c^2) * segment_b_values(intersect_kn1,:) - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k1);
segment_b_values_to_end_kn1 = linspace(intersect_kn_b1,max_b,points_number);
segment_Hk_a_values1_to_end = (1 - c^2) * segment_b_values_to_end_kn1 - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k1);


%%
gwLineWidth = 0.5;
figure
set(gcf,"Position",[300 300 400 400])
axes1 = axes('Position',[0.16625 0.21375 0.75 0.75]);
hold on
dq=lines(M); 
for i=1:M
    plot(segment_b_values(i,:),segment_a_values(i,:),'LineStyle','-','LineWidth',gwLineWidth,'color',dq(i,:))
end
plot(b_value,T_a_values(1,:),'LineStyle','--','color','r','LineWidth',gwLineWidth)

plot(segment_b_values_to_end_1n,segment_H1_a_values_to_end,'LineStyle','-','color','b','LineWidth',gwLineWidth)
plot(segment_b_values_to_end_kn1,segment_Hk_a_values1_to_end,'LineStyle','-.','color','r','LineWidth',gwLineWidth)

s=100;scatter(connecting_b,connecting_a,s,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','.')

s=100;scatter(intersect_1n_b,intersect_1n_a,s,'MarkerFaceColor','r','MarkerEdgeColor','r','Marker','pentagram')

% s=100;scatter(intersect_kn_b1,intersect_kn_a1,s,'MarkerFaceColor','r','MarkerEdgeColor','r','Marker','pentagram')

% s=100;scatter(intersect_11_b,intersect_11_a,s,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','.')
plot(intersect_11_b_to_intersect_1n_b,intersect_11_a_to_intersect_1n_a,'LineStyle','--','color','b','LineWidth',gwLineWidth)
 


% plot(1.16,0.5,'MarkerSize',3,'Marker','.','Color','k')
% plot(4.3,2.5,'MarkerSize',3,'Marker','.','Color','k')

plot([0:0.5:25],3.0*ones(size([0:0.5:25])),'LineStyle','-.','Color','k','LineWidth',gwLineWidth)
plot([0:0.5:25],1.0*ones(size([0:0.5:25])),'LineStyle','-.','Color','k','LineWidth',gwLineWidth)


ax = gca;
% xticks = xticks(ax1);
xticks = [0:1.0:5];
xticklabels = arrayfun(@(x) sprintf('%.0f', x), xticks, 'UniformOutput', false);
xticklabels(1) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [0:1.0:5];
yticklabels = arrayfun(@(x) sprintf('%.0f', x), yticks, 'UniformOutput', false);
yticklabels(1) = {'0'};
% yticklabels(end) = {'$a$'};
set(gca, 'ytick', yticks);
set(gca, 'yticklabel', yticklabels);

ax_FontSize = 24;
ax.XAxis.FontSize = ax_FontSize;  
ax.XAxis.FontName = 'Times New Roman';  
ax.YAxis.FontSize = ax_FontSize;  
ax.YAxis.FontName = 'Times New Roman';  
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'in';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

xlim([0 4.5])
ylim([0 4.5])

xlabel('$b$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
ylabel('$a$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');


box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'linewidth',1.0,'layer','top');


figure_name = [results_path 'Figure Step01 BifurcationDiagram npbcwsnetworks_n=10k=2p=0number=1.eps'];
saveas(gcf, figure_name, 'epsc');
