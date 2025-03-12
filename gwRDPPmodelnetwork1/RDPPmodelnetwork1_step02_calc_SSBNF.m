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
%%

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

disp(['eigenvalues(k) = ' num2str(eigenvalues(k),'%.20f')])

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

phi = eigenvectors(:,k);
psi = eigenvectors(:,k);

qn_vector = kron(phi,MK_q);
qn_star_vector = kron(psi,MK_q_star);

%%
T_k = -(d1+d2)*eigenvalues(k) - b_star * (1-c^2) + a + c*(1-c);
J_k = d1*d2*eigenvalues(k)^2 + (d1*c*(c-1) + d2*(b_star*(1-c^2) - a))*eigenvalues(k) + c*(1-c)*(a + b_star*(c-1));
pk1 = 1; 
pk2 = (d1*eigenvalues(k) + b_star*(1-c^2) - a)/c^2;
pk = [pk1; pk2];

phi = eigenvectors(:,k);
psi = eigenvectors(:,k);

qn_vector = kron(phi,pk);

qk1 = (-d2*eigenvalues(k) + c*(1-c))/T_k;
qk2 = - c^2/T_k;
qk = [qk1; qk2];
qn_star_vector = kron(psi,qk);

f101 = [1 - c^2; (c - 1)^2];

f110 =[ (2*a*c^2*(c - 1))/(b_star*(a - b_star + b_star*c));
       -(2*a*c^2*(c - 1))/(b_star*(a - b_star + b_star*c))];
f011 = [0;0];

f200 =[-(2*a*(- b_star*c^3 + 2*b_star*c^2 + a - b_star))/(b_star*(a - b_star + b_star*c)); 
       -(2*a*c*(c - 1)^2)/(a - b_star + b_star*c)];
f020 =[ (2*a*c^3)/(b_star^2*(a - b_star + b_star*c));
       -(2*a*c^3)/(b_star^2*(a - b_star + b_star*c))];

F_k20 = f200*pk1^2 + 2*f110*pk1*pk2 + f020*pk2^2;

S_k11 = qk.'*(pk1*f101+pk2*f011);

sum_k3 = sum(eigenvectors(:,k).^3);
if abs(sum_k3) < 10e-7
    sum_k3 = 0;
end
S_k20 = 1/2 * qk.' * F_k20 * sum_k3;


f300 = [-(6*a^2*c^2*(c - 1)^2)/(b_star*(a - b_star + b_star*c)^2);
          (6*a^2*c^2*(c - 1)^2)/(b_star*(a - b_star + b_star*c)^2)];
f030 = [-(6*a^2*c^4)/(b_star^4*(a - b_star + b_star*c)^2);
         (6*a^2*c^4)/(b_star^4*(a - b_star + b_star*c)^2)];
f120 = [-(2*a^2*c^3*(3*c - 2))/(b_star^3*(a - b_star + b_star*c)^2);
         (2*a^2*c^3*(3*c - 2))/(b_star^3*(a - b_star + b_star*c)^2)];
f210 = [-(2*a^2*c^2*(3*c^2 - 4*c + 1))/(b_star^2*(a - b_star + b_star*c)^2);
         (2*a^2*c^2*(3*c^2 - 4*c + 1))/(b_star^2*(a - b_star + b_star*c)^2)];
f201 = [-(2*a*(- a^2 - 2*a*b_star*c + 2*a*b_star + b_star^2*c^4 - 3*b_star^2*c^3 + 2*b_star^2*c^2 + b_star^2*c - b_star^2))/(b_star^2*(a - b_star + b_star*c)^2);
         (2*a*c*(c - 1)^3)/(a - b_star + b_star*c)^2];
f021 = [-(2*a*c^3*(2*a - 3*b_star + 3*b_star*c))/(b_star^3*(a - b_star + b_star*c)^2);
          (2*a*c^3*(2*a - 3*b_star + 3*b_star*c))/(b_star^3*(a - b_star + b_star*c)^2)];
f111 = [-(2*a*c^2*(c - 1)*(a - 2*b_star + 2*b_star*c))/(b_star^2*(a - b_star + b_star*c)^2);
          (2*a*c^2*(c - 1)*(a - 2*b_star + 2*b_star*c))/(b_star^2*(a - b_star + b_star*c)^2)];
     
c_k30 = qk.' * (f300*pk1^3 + f030*pk2^3 + 3*f210*pk1^2*pk2 + 3*f120*pk1*pk2^2);
sum_k4 = sum(eigenvectors(:,k).^4);
C_k30 = c_k30 * sum_k4 * 1/6;

c_k21 = qk.' * ( f201*pk1^2 + 2*f111*pk1*pk2 + f021*pk2^2 );
sum_k3 = sum(eigenvectors(:,k).^3);
if abs(sum_k3) < 10e-7
    sum_k3 = 0;
end
C_k21 = c_k21 * 1/2 * sum_k3;

h_kj20 = zeros(2,N);
h_kj11 = zeros(2,N);

for j=1:N
    Mj = [d1*eigenvalues(j)+b_star*(1-c^2)-a, -c^2;
            b_star*(1-c)^2,          d2*eigenvalues(j)-c*(1-c)];
    if j ~= k
        sum_k2_j = sum(eigenvectors(:,k).^2.*eigenvectors(:,j));
        if abs(sum_k2_j) < 10e-7
            sum_k2_j = 0;
        end
        h_kj20(:,j) = - inv(Mj)*F_k20*sum_k2_j;
        h_kj11(:,j) = - inv(Mj)*[0;0];
    else
        M = [Mj,pk];
        M = [M;[qk;0]'];
        sum_k3 = sum(eigenvectors(:,k).^3);
        if abs(sum_k3) < 10e-7
            sum_k3 = 0;
        end
        %%%%%%%%%%%%%%
        b_kk20 = ( F_k20 - qk.' * F_k20 * pk ) * sum_k3;
        bar_h_kk20 = - inv(M) * [b_kk20;0];
        h_kj20(:,j) = bar_h_kk20([1,2]);
        %%%%%%%%%%%%%%
        bkk11 = 2 * ( (f101*pk1+f011*pk2) - qk.'*(f101*pk1+f011*pk2)*pk );
        bar_h_kk11 = - inv(M) * [bkk11;0];
        h_kj11(:,j) = bar_h_kk11([1,2]);
    end
end

E_k30_vector = zeros(1,N);
E_k21_vector1 = zeros(1,N);
E_k21_vector2 = zeros(1,N);
for j=1:N
    
    sum_k2_j = sum(eigenvectors(:,k).^2.*eigenvectors(:,j));
    
    if abs(sum_k2_j) < 10e-7
        sum_k2_j = 0;
    end
    E_k30_vector(j) = 1/3 * qk.'*((f200*pk1+f110*pk2)*h_kj20(1,j)+(f020*pk2+f110*pk1)*h_kj20(2,j))*sum_k2_j;
    sum_k_j = sum(eigenvectors(:,k).*eigenvectors(:,j));
    if abs(sum_k_j) < 10e-7
        sum_k_j = 0;
    end
    E_k21_vector1(j) = 1/3 * qk.'*((f200*pk1+f110*pk2)*h_kj11(1,j)+(f020*pk2+f110*pk1)*h_kj11(2,j))*sum_k2_j;
    E_k21_vector2(j) = 1/3 * qk.'*(f101*h_kj20(1,j)+f011*h_kj20(2,j))*sum_k_j;
end

E_k30 = sum(E_k30_vector);
E_k21 = sum(E_k21_vector1 + E_k21_vector2);

S_k21 = C_k21 + 3/2 * E_k21;

S_k30 = C_k30 + 3/2 * E_k30;

format long

% [S_k11,S_k20,S_k21,S_k30]



%%
% mu_star = - S_k20 / S_k21

% [S_k11*mu,S_k20+S_k21*mu]
% x = - S_k11*mu / (S_k20+S_k21*mu)


%%
[S_k30 S_k11; 
S_k20 S_k21;
S_k30*S_k11 S_k20*S_k21];

delta_mu = 16*S_k11*S_k30*(S_k30*S_k11-S_k20*S_k21);

% [fprintf('S_k30=%20.20f',S_k30),fprintf(',   S_k11=%20.20f\n',S_k11),...
%  fprintf('S_k20=%20.20f',S_k20),fprintf(',   S_k21=%20.20f\n',S_k21)];

disp(['S_k11= ', num2str(S_k11,'%.6f')])
disp(['S_k20= ', num2str(S_k20,'%.6f')])
disp(['S_k21= ', num2str(S_k21,'%.6f')])
disp(['S_k30= ', num2str(S_k30,'%.6f')])

% [fprintf('delta_mu=%20.20f',delta_mu),fprintf(',   S_k30*S_k11=%20.20f',S_k30*S_k11),fprintf(',   S_k20*S_k21=%20.20f\n',S_k20*S_k21)];
%%
% [S_k11*mu,S_k20+S_k21*mu,S_k30]
% delta = (S_k20+S_k21*mu)^2 - 4*S_k11*S_k30*mu;
% x1 = (-(S_k20+S_k21*mu) - sqrt(delta)) / (2*S_k30);
% x2 = (-(S_k20+S_k21*mu) + sqrt(delta)) / (2*S_k30);
% [x1, x2]

syms x
eq3 = (S_k20+S_k21*x)^2 - 4*S_k11*S_k30*x;

pde_sol = solve(eq3, x);
vpa_sol = [vpa(pde_sol,12)];



%
N = size(laplacian_matrix,1);
b_star;
b_star1=b_star;
b_star = roundn(b_star,-qu_n);
b_star2=b_star;
b_star = b_star + perturbating_mu;

% b_star = 0.302387 + perturbating_mu; % para2
b_star3=b_star;


T_ks = zeros(N,1);
J_ks = zeros(N,1);
for j=1:N
%     j
    T_ks(j) = -(d1+d2)*eigenvalues(j) - b_star* (1-c^2) + a + c*(1-c);
    J_ks(j) = d1*d2*eigenvalues(j)^2 + (d1*c*(c-1) + d2*(b_star*(1-c^2) - a))*eigenvalues(j) + c*(1-c)*(a + b_star*(c-1));
end

disp( ['T_k: ' num2str( [length(find(T_ks<0))] ) ' ---> ' num2str( [find(T_ks<0)'] ) ] )
disp( ['J_k: ' num2str( [length(find(J_ks<0))] ) ' ---> ' num2str( [find(J_ks<0)'] ) ] )
format short

%%
Bk30=S_k30;
Bk11=S_k11;
Bk20=S_k20;
Bk21=S_k21;

mu_star = - Bk20 / Bk21;
delta1 = 16*Bk11*Bk30 * (Bk11*Bk30 - Bk20*Bk21);
mu1 = (-(2*Bk20*Bk21-4*Bk11*Bk30) -sqrt(delta1) ) / 2 / Bk21^2;
mu2 = (-(2*Bk20*Bk21-4*Bk11*Bk30) +sqrt(delta1) ) / 2 / Bk21^2;

mu = linspace(-maxmu,maxmu,1000);

gwLineWidth = 1.0;

figure
set(gcf,"Position",[300 300 500 320])
% axes1 = axes('Position',[0.16625 0.21375 0.75 0.75]);
% set(gcf,'Position',get(0,'ScreenSize'))
hold on
box on

mu0_left = linspace(-maxmu,0,1000);
x0_left = zeros(size(mu0_left));
plot(mu0_left, x0_left,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle','-','Color','k');

mu0_right = linspace(0,maxmu,1000);
x0_right = zeros(size(mu0_right));
plot(mu0_right, x0_right,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle',':','Color','k');

mu_up_right = linspace(0,maxmu,1000);
delta = (Bk20 + Bk21*mu_up_right).^2 - 4*Bk11*Bk30*mu_up_right;
x_up_1 = (-(Bk20 + Bk21*mu_up_right) - sqrt(delta)) / 2/ Bk30;
plot(mu_up_right, x_up_1,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle','-','Color','k');

mu_down_right = linspace(0,maxmu,1000);
delta = (Bk20 + Bk21*mu_down_right).^2 - 4*Bk11*Bk30*mu_down_right;
x_down_1 = (-(Bk20 + Bk21*mu_down_right) + sqrt(delta)) / 2/ Bk30;
plot(mu_down_right, x_down_1,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle','-','Color','k');

ax = gca;
% xticks = xticks(ax1);
xticks = [-maxmu,0,maxmu];
xticklabels = arrayfun(@(x) sprintf('%.4f', x), xticks, 'UniformOutput', false);
xticklabels(2) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-max_y,0,max_y];
yticklabels = arrayfun(@(x) sprintf('%.2f', x), yticks, 'UniformOutput', false);
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

ylim([-max_y,max_y])


% set(gca,'XTick',linspace(-maxmu,maxmu,8+1));
% set(gca,'XTicklabel',linspace(-maxmu,maxmu,8+1), 'FontSize', 24,'FontWeight','bold');
% xtickangle(-30)
% xlabel('\mu','FontWeight','bold','FontName','Times New Roman', 'FontSize', 36)
% ylabel('z^{*}_{1}','FontWeight','bold','FontName','Times New Roman', 'FontSize', 36)
% set(gca,'LineWidth',gwLineWidth);

% xlabel('$\mu_{1}$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% ylabel('$z^{*}_{1}$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');


box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'linewidth',1.0,'layer','top');

figure_name = [results_path 'Figure Step02 SSBNF diagram npbcwsnetworks_n=10k=2p=0number=1.eps'];
saveas(gcf, figure_name, 'epsc');
