clear all
close all
clc

%%

network_name = 'ws_n=10k=2p=0.05number=796';
datanetwork_path = ['./network2 ' network_name '/'];
results_path = strrep(datanetwork_path,'network2','RDPPmodel results network2');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

net_number = 796;
mat_name = [datanetwork_path 'adj_matrix_wsnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);

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

digits = 50;

weimu = 0.01;

a = 3.0;
c = 0.2;
d1 = 0.001; d2 = 10*d1;
Times = 0.01;

N = length(eigenvalues);
d1 = d1 / Times;
d2 = d2 / Times;


k = 1;

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
k = 1;
b_star = a/(1-c^2) + c/(1+c) - (d1+d2)*eigenvalues(k)/(1-c^2);
T_k = -(d1+d2)*eigenvalues(k) - b_star * (1-c^2) + a + c*(1-c);



format long
b_star;

omega_k = sqrt(d1*d2*eigenvalues(k)^2 + (d1*c*(c-1)+d2*(b_star*(1-c^2)-a))*eigenvalues(k) + c*(1-c)*(a+b_star*(c-1)));

pk1 = 1; 
pk2 = (d1*eigenvalues(k) + b_star*(1-c^2) - a - i*omega_k)/c^2;
pk = [pk1; pk2];

qk1 = (-d2*eigenvalues(k) + c*(1-c) + i*omega_k)/(2*i*omega_k);
qk2 = - c^2/(2*i*omega_k);
qk = [qk1; qk2];

f200 =[-(2*a*(- b_star*c^3 + 2*b_star*c^2 + a - b_star))/(b_star*(a - b_star + b_star*c)); 
       -(2*a*c*(c - 1)^2)/(a - b_star + b_star*c)];
f020 =[ (2*a*c^3)/(b_star^2*(a - b_star + b_star*c));
       -(2*a*c^3)/(b_star^2*(a - b_star + b_star*c))];  
f002 =[ 0; 0];  

f110 =[ (2*a*c^2*(c - 1))/(b_star*(a - b_star + b_star*c));
       -(2*a*c^2*(c - 1))/(b_star*(a - b_star + b_star*c))];
f101 = [1 - c^2; (c - 1)^2];
f011 =[ 0; 0];  

%  B_k1 = ((1-c^2)*qk1 + (1-c)^2*qk2) * pk1
% B_k1 = pk1*(qk.'*f101) + pk2*(qk.'*f011)
B_k1 = qk.'*(f101*pk1 + f011*pk2);
nv_k1 = real(B_k1);
 
f300 = [-(6*a^2*c^2*(c - 1)^2)/(b_star*(a - b_star + b_star*c)^2);
          (6*a^2*c^2*(c - 1)^2)/(b_star*(a - b_star + b_star*c)^2)];
f030 = [-(6*a^2*c^4)/(b_star^4*(a - b_star + b_star*c)^2);
         (6*a^2*c^4)/(b_star^4*(a - b_star + b_star*c)^2)];
f003 =[ 0; 0];  

f120 = [-(2*a^2*c^3*(3*c - 2))/(b_star^3*(a - b_star + b_star*c)^2);
         (2*a^2*c^3*(3*c - 2))/(b_star^3*(a - b_star + b_star*c)^2)];
f210 = [-(2*a^2*c^2*(3*c^2 - 4*c + 1))/(b_star^2*(a - b_star + b_star*c)^2);
         (2*a^2*c^2*(3*c^2 - 4*c + 1))/(b_star^2*(a - b_star + b_star*c)^2)];
f012 = [ 0; 0];  
f021 = [-(2*a*c^3*(2*a - 3*b_star + 3*b_star*c))/(b_star^3*(a - b_star + b_star*c)^2);
         (2*a*c^3*(2*a - 3*b_star + 3*b_star*c))/(b_star^3*(a - b_star + b_star*c)^2)];
f102 = [ 0; 0];  
f201 = [-(2*a*(- a^2 - 2*a*b_star*c + 2*a*b_star + b_star^2*c^4 - 3*b_star^2*c^3 + 2*b_star^2*c^2 + b_star^2*c - b_star^2))/(b_star^2*(a - b_star + b_star*c)^2);
         (2*a*c*(c - 1)^3)/(a - b_star + b_star*c)^2];

c_k210 = qk.'*(f300*pk1*abs(pk1)^2 + f030*pk2*abs(pk2)^2 + f210*(pk1^2*conj(pk2)+2*pk2*abs(pk1)^2) + ...
    f120*(pk2^2*conj(pk1)+2*pk1*abs(pk2)^2));
sum_k4 = sum(eigenvectors(:,k).^4);
C_k210 = c_k210 * 1/2 * sum_k4;

F_k200 = f200*pk1^2 + 2*f110*pk1*pk2 + f020*pk2^2;
F_k020 = conj(F_k200);
F_k110 = f200*abs(pk1)^2 + f110*(pk1*conj(pk2)+conj(pk1)*pk2) + f020*abs(pk2)^2;

d_k210 = (-(qk.'*F_k200)*(qk.'*F_k110) + 2*abs(qk.'*F_k110)^2 + 1/3*abs(qk.'*F_k020)^2);
sum_k3 = sum(eigenvectors(:,k).^3);
if abs(sum_k3) < 10e-7
    sum_k3 = 0;
end
D_k210 = d_k210 * 1/3 / i/omega_k * sum_k3^2;

h_kj200 = zeros(2,N);
h_kj110 = zeros(2,N);

for j=1:N
    Mj = [d1*eigenvalues(j)+b_star*(1-c^2)-a, -c^2;
            b_star*(1-c)^2,          d2*eigenvalues(j)-c*(1-c)];
    if j ~= k
        sum_k2_j = sum(eigenvectors(:,k).^2.*eigenvectors(:,j));
        if abs(sum_k2_j) < 10e-7
            sum_k2_j = 0;
        end
        h_kj200(:,j) = inv(2*i*omega_k*[1,0;0,1] - Mj)*F_k200*sum_k2_j;
        h_kj110(:,j) = - 2 * inv(Mj)*F_k110*sum_k2_j;
    else
        sum_k3 = sum(eigenvectors(:,k).^3);
        if abs(sum_k3) < 10e-7
            sum_k3 = 0;
        end
        %%%%%%%%%%%%%%
        h_kj200(:,j) = inv(2*i*omega_k*[1,0;0,1] - Mj)*( F_k200- (qk.'*F_k200*pk + conj(qk).'*F_k200*conj(pk)) )*sum_k3;
        %%%%%%%%%%%%%%
        h_kj110(:,j) = - 2 * inv(Mj)*( F_k110- (qk.'*F_k110*pk + conj(qk).'*F_k110*conj(pk)) )*sum_k3;
    end
end

E_k210_vector1 = zeros(1,N);
E_k210_vector2 = zeros(1,N);
for j=1:N
    sum_k2_j = sum(eigenvectors(:,k).^2.*eigenvectors(:,j));
    if abs(sum_k2_j) < 10e-7
        sum_k2_j = 0;
    end
    
    E_k210_vector1(j) = 1/3 * qk.'*((f200*pk1+f110*pk2)*h_kj110(1,j)+(f020*pk2+f110*pk1)*h_kj110(2,j))*sum_k2_j;
    E_k210_vector2(j) = 1/3 * qk.'*((f200*conj(pk1)+f110*conj(pk2))*h_kj200(1,j)+(f020*conj(pk2)+f110*conj(pk1))*h_kj200(2,j))*sum_k2_j;

end

E_k210 = sum(E_k210_vector1 + E_k210_vector2);

B_k210 = C_k210 + 3/2 * (D_k210 + E_k210);
nv_k2 = real(B_k210);

%%

b_star = b_star + weimu;
N = size(laplacian_matrix,1);
T_ks = zeros(N,1);
J_ks = zeros(N,1);

small_matrices_eigenvalues = zeros(N,2);
deltas = zeros(N,1);
for j=1:N
%     k
    T_ks(j) = -(d1+d2)*eigenvalues(j) - b_star* (1-c^2) + a + c*(1-c);
    J_ks(j) = d1*d2*eigenvalues(j)^2 + (d1*c*(c-1) + d2*(b_star*(1-c^2) - a))*eigenvalues(j) + c*(1-c)*(a + b_star*(c-1));
    one_delta = T_ks(j)^2 - 4*J_ks(j);
    deltas(j) = one_delta;
    small_matrices_eigenvalues(j,1) = (-T_ks(j) + sqrt(one_delta)) / 2;
    small_matrices_eigenvalues(j,2) = (-T_ks(j) - sqrt(one_delta)) / 2;
end

disp( ['T_k: ' num2str( [length(find(T_ks<0))] ) ' ---> ' num2str( [find(T_ks<0)'] ) ] )
disp( ['J_k: ' num2str( [length(find(J_ks<0))] ) ' ---> ' num2str( [find(J_ks<0)'] ) ] )
% small_matrices_eigenvalues(1,:)
% small_matrices_eigenvalues(2,:)
% small_matrices_eigenvalues(3,:)
% small_matrices_eigenvalues(4,:)

%%
% nv_k1*mu*rho + nv_k2*rho^3

mu = linspace(-weimu,weimu,1000);

gwLineWidth = 1.0;

figure
set(gcf,"Position",[300 300 500 320])
% axes1 = axes('Position',[0.16625 0.21375 0.75 0.75]);
% set(gcf,'Position',get(0,'ScreenSize'))
hold on
box on

mu0_left = linspace(-weimu,0,1000);
x0_left = zeros(size(mu0_left));
plot(mu0_left, x0_left,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle','-','Color','k');

mu0_right = linspace(0,weimu,1000);
x0_right = zeros(size(mu0_right));
plot(mu0_right, x0_right,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle','--','Color','k');

mu_up_right = linspace(0,weimu,1000);
rho_up_right = real(sqrt(-(nv_k1*mu_up_right)/nv_k2));
plot(mu_up_right, rho_up_right,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle','-','Color','k');

ax = gca;
% xticks = xticks(ax1);
xticks = [-weimu,0,weimu];
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
xticklabels(2) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-0.1,0,0.5];
yticklabels = arrayfun(@(x) sprintf('%.1f', x), yticks, 'UniformOutput', false);
yticklabels(2) = {'0'};
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

set(gca,'XColor','k');
set(gca,'YColor','k');

ylim([-0.15,0.65])
xlim([-weimu weimu])

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'linewidth',1.0,'layer','top');

figure_name = [results_path 'Figure Step02 HBNF diagram wsnetworks_n=10k=2p=0.05number=796.eps'];
saveas(gcf, figure_name, 'epsc');