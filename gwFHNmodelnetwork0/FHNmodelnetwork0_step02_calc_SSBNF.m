clear all
close all
clc

tic
digits = 50;
format long

%%

perturbating_mu =  0.002;
qu_n = 6;
maxmu = 0.002;
max_y = 0.3;

%%

network_name = 'brain';
datanetwork_path = ['./network0 ' network_name '/'];
results_path = strrep(datanetwork_path,'network0','FHN results network0');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

mat_name = [datanetwork_path 'adj_matrix_brainnetwork'];
load(mat_name)
adjacent_matrix = full(adj_matrix);


cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); % 
disp(degree')

laplacian_matrix = adjacent_matrix - diag(degree);

N = size(laplacian_matrix,1);
[eigenvectors, eigenvalues] = eig(laplacian_matrix);
eigenvalues = diag(eigenvalues);
eigenvectors = fliplr(eigenvectors);
eigenvalues = flipud(eigenvalues);
eigenvalues(1) = 0;
eigenvectors(:,1) = - eigenvectors(:,1);

%% set diffusion rate
theta = -1;
du = 0.3; 
fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.8f\n', -theta*du/(-theta+du*eigenvalues(2)));

for i=1:N-1
    panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
end
disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])
% dv = panbie_THB(6);
% dv = (panbie_THB(5)+panbie_THB(6))/2;
dv = 0.8;
c_set = 1.8;
c = [0:0.001:3];  
ymax = 2;
% %% set diffusion rate
% theta = -2;
% du = 0.04; 
% fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.8f\n', -theta*du/(-theta+du*eigenvalues(2)));
% 
% for i=1:N-1
%     panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
% end
% disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])
% % dv = panbie_THB(4);
% dv = (panbie_THB(3)+panbie_THB(4))/2;
% % dv = 0.06;
% c_set = 6;
% c = [0:0.001:10];  
% ymax = 0.8;

% %% set diffusion rate
% theta = -1/2;
% du = 0.01; 
% fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.8f\n', -theta*du/(-theta+du*eigenvalues(2)));
% 
% for i=1:N-1
%     panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
% end
% disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])
% dv = panbie_THB(8);
% % dv = 0.06;
% c_set = 1;
% c = [0:0.001:2];  
% ymax = 3;
%%
panbie_H = -theta + eigenvalues * (du + dv);
p_H = find(panbie_H>0);
Ha = zeros(length(p_H),length(c));
for i=1:length(p_H)
    Ha(i,:) = (-theta + eigenvalues(p_H(i)) * (du + dv))./c; 
end
% figure
% hold on
% grid minor
% for i=1:length(p_H)
%     plot(c,Ha(i,:),'b-');
% end

%%
panbie_T = -theta + eigenvalues * du;
p_T = find(panbie_T>0);
Ta = zeros(length(p_T),length(c));
for i=1:length(p_T)
    Ta(i,:) = 1/(du*eigenvalues(p_T(i)) -theta) + 1./c * dv * eigenvalues(p_T(i));
end
% figure
% hold on
% grid minor
% 
% for i=1:length(p_H)
%     plot(c,Ha(i,:),'b-');
% end
% for i=1:length(p_T)
%     plot(c,Ta(i,:),'r-');
% end

%% determine intersecting points of T_{k} and T_{k+1}, its number = length(p_T)-1?

intersecting_points_TT = zeros(2,length(p_T)-1);
panbie_intersect_upper_x = zeros(1,length(p_T)-1);
panbie_intersect_upper_H1 = zeros(1,length(p_T)-1);
for i=1:length(p_T)-1
    intersecting_c_ii1 = dv/du * ( -theta + du*eigenvalues(p_T(i)) ) * ( -theta + du*eigenvalues(p_T(i+1)) );
    intersecting_a_ii1 = 1/(du*eigenvalues(p_T(i)) -theta) + 1/intersecting_c_ii1 * dv * eigenvalues(p_T(i));
    intersecting_points_TT(:,i) = [intersecting_c_ii1,intersecting_a_ii1]';
    panbie_intersect_upper_x(i) = -theta + du*eigenvalues(p_T(i)) + du*eigenvalues(p_T(i+1));
    panbie_intersect_upper_H1(i) = du * dv * (  eigenvalues(p_T(i)) + eigenvalues(p_T(i+1))  ) + theta*(du - dv);
end

p_intersect_upper_x = find(panbie_intersect_upper_x>=0);
intersecting_points_TT_upper_x = intersecting_points_TT(:,p_intersect_upper_x);

p_intersect_upper_H1 = find(panbie_intersect_upper_H1>=0);

% p_intersect_upper_H1 = [p_intersect_upper_H1,p_intersect_upper_H1(end)+1];

intersecting_points_TT_upper_H1 = intersecting_points_TT(:,p_intersect_upper_H1);

intersect_c_H1Tend = - du*dv*eigenvalues(p_intersect_upper_H1(end)+1)^2 - theta* (du-dv) * eigenvalues(p_intersect_upper_H1(end)+1) + theta^2 ;
intersect_a_H1Tend = 1/(du*eigenvalues(p_intersect_upper_H1(end)+1) -theta) + 1/intersect_c_H1Tend * dv * eigenvalues(p_intersect_upper_H1(end)+1);


%% plot
figure
hold on
grid minor

plot(c,Ha(1,:),'b-');

for i=1:length(p_H)
    plot(c,Ha(i,:),'b-');
end
% for i=1:length(p_T)
% for i=1:length(p_intersect_upper_x)
for i=1:size(intersecting_points_TT_upper_H1,2)+1
    plot(c,Ta(i,:),'r-');
end
% plot(intersecting_points_TT(1,:), intersecting_points_TT(2,:), 'g*')
% plot(intersecting_points_TT_upper_x(1,:), intersecting_points_TT_upper_x(2,:), 'co')
plot(intersecting_points_TT_upper_H1(1,:), intersecting_points_TT_upper_H1(2,:), 'k*')
plot(intersect_c_H1Tend,intersect_a_H1Tend,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'Marker','pentagram','LineStyle','none','Color',[0 0 0]);
fprintf('critical THB point a=%.8f, c=%.8f, on Turing curve k=%d, Hopf curve k=%d.\n', intersect_a_H1Tend, intersect_c_H1Tend, size(intersecting_points_TT_upper_H1,2)+1,1);

plot(c,ones(size(c))*1/(-theta),'k-')
plot(c,zeros(size(c)),'k-')

for i=1:length(p_intersect_upper_H1)-1
%     [i p_intersect_upper_H1(i+1)]
    intersecting_c_ii1 = intersecting_points_TT_upper_H1(1,i);
    intersecting_c_ii2 = intersecting_points_TT_upper_H1(1,i+1);
    Turing_c_i = linspace(intersecting_c_ii2, intersecting_c_ii1, 100);
    Turing_a_i = 1/(du*eigenvalues(p_intersect_upper_H1(i+1)) -theta) + 1./Turing_c_i * dv * eigenvalues(p_intersect_upper_H1(i+1));
    plot(Turing_c_i,Turing_a_i,'k-');
end
i=length(p_intersect_upper_H1);
intersecting_c_ii1 = intersecting_points_TT_upper_H1(1,i);
intersecting_c_H1Tend = - du*dv*eigenvalues(p_intersect_upper_H1(end)+1)^2 - theta* (du-dv) * eigenvalues(p_intersect_upper_H1(end)+1) + theta^2 ;
Turing_c_i = linspace(intersecting_c_H1Tend, intersecting_c_ii1, 100);
Turing_a_i = 1/(du*eigenvalues(i+1) -theta) + 1./Turing_c_i * dv * eigenvalues(i+1);
plot(Turing_c_i,Turing_a_i,'k-');

% axis equal

% xlim([0 40])
ylim([0 ymax])

% figure
% hold on
% grid minor
% 
% plot(c,Ha(1,:),'b-');
% 
% plot(intersecting_points_TT_upper_H1(1,:), intersecting_points_TT_upper_H1(2,:), 'k*')
% plot(intersect_c_H1Tend,intersect_a_H1Tend,'k*')
% plot(c,ones(size(c)),'k-')
% plot(c,zeros(size(c)),'k-')
% 
% for i=1:length(p_intersect_upper_H1)-1
% %     [i p_intersect_upper_H1(i+1)]
%     intersecting_c_ii1 = intersecting_points_TT_upper_H1(1,i);
%     intersecting_c_ii2 = intersecting_points_TT_upper_H1(1,i+1);
%     Turing_c_i = linspace(intersecting_c_ii2, intersecting_c_ii1, 100);
%     Turing_a_i = 1/(du*eigenvalues(p_intersect_upper_H1(i+1)) -theta) + 1./Turing_c_i * dv * eigenvalues(p_intersect_upper_H1(i+1));
%     plot(Turing_c_i,Turing_a_i,'k-');
% end
% i=length(p_intersect_upper_H1);
% intersecting_c_ii1 = intersecting_points_TT_upper_H1(1,i);
% intersecting_c_H1Tend = - du*dv*eigenvalues(p_intersect_upper_H1(end)+1)^2 - theta* (du-dv) * eigenvalues(p_intersect_upper_H1(end)+1) + theta^2 ;
% Turing_c_i = linspace(intersecting_c_H1Tend, intersecting_c_ii1, 100);
% Turing_a_i = 1/(du*eigenvalues(i+1) -theta) + 1./Turing_c_i * dv * eigenvalues(i+1);
% plot(Turing_c_i,Turing_a_i,'k-');
% 
% 
% 
% % axis equal
% 
% % xlim([0 4])
% ylim([0 2.01])

%% 

c = c_set;
intersecting_points_c_TT_upper_H1 = intersecting_points_TT_upper_H1(1,:);
disp(['intersecting SSB points for c:   ' num2str(intersecting_points_c_TT_upper_H1)])
% [v p_critical_SSB] = min(abs(intersecting_points_c_TT_upper_H1-c))
p_critical_SSB = find(intersecting_points_c_TT_upper_H1-c > 0);
p_critical_SSB = p_critical_SSB(end);
a_critical_SSB = 1/(du*eigenvalues(p_critical_SSB + 1) -theta) + 1./c * dv * eigenvalues(p_critical_SSB + 1);



plot(c,a_critical_SSB,'g*')
fprintf('critical SSB point a=%.8f, with setting c=%.8f, on Turing curve k=%d.\n', a_critical_SSB, c, p_critical_SSB+1);


a_critical_HB = (-theta + eigenvalues(1) * (du + dv))./c; 
fprintf('critical HB point a=%.8f, with setting c=%.8f, on Hopf curve k=%d.\n', a_critical_HB, c, 1);
plot(c,a_critical_HB,'b*')

%% 计算 SSB 的规范型

a_star = a_critical_SSB;
format rat
a_star
format long

k = p_critical_SSB+1;
fprintf('确定的Turing线序号k= %d, 对应的 SSB 分支点 a* = %.10f, at c=%.10f.\n', k, a_star, c);
a_star_prepare = a_star;
% d1 = du; d2 =dv;
%
a11 = -theta;
a12 = -1;
a21 = c;
a22 = - a_star*c;

J = [a11,a12; a21,a22]; 
D = [du,0;0,dv];

T_k = (du+dv)*eigenvalues(k) + (a11 + a22);
J_k = du*dv*eigenvalues(k)^2 + (dv*a11 + du*a22)*eigenvalues(k) + a11*a22 - a12*a21;

pk1 = 1; 
pk2 = - (a11+du*eigenvalues(k))/a12;
pk = [pk1; pk2];

qk1 = (dv*eigenvalues(k) + a22)/T_k;
qk2 = - a12/T_k;
qk = [qk1; qk2];

pk'*qk % check!

f101 = [0; 0];
f110 = [0; 0];
f011 = [0; -c];

f200 = [2*theta + 2; 0];
f020 = [0; 0];

F_k20 = f200*pk1^2 + 2*f110*pk1*pk2 + f020*pk2^2;

S_k11 = qk.'*(pk1*f101+pk2*f011);

sum_k3 = sum(eigenvectors(:,k).^3);
if abs(sum_k3) < 10e-7
    sum_k3 = 0;
end
S_k20 = 1/2 * qk.' * F_k20 * sum_k3;


f300 = [-6; 0];
f030 = [0; 0];
f120 = [0; 0];
f210 = [0; 0];
f201 = [0; 0];
f021 = [0; 0];
f111 = [0; 0];
     
c_k30 = qk.' * (f300*pk1^3 + f030*pk2^3 + 3*f210*pk1^2*pk2 + 3*f120*pk1*pk2^2);
sum_k4 = sum(eigenvectors(:,k).^4);
C_k30 = c_k30 * 1/6 * sum_k4;

c_k21 = qk.' * ( f201*pk1^2 + 2*f111*pk1*pk2 + f021*pk2^2 );
sum_k3 = sum(eigenvectors(:,k).^3);
if abs(sum_k3) < 10e-7
    sum_k3 = 0;
end
C_k21 = c_k21 * 1/2 * sum_k3;

h_kj20 = zeros(2,N);
h_kj11 = zeros(2,N);

for j=1:N
    Mj = J + eigenvalues(j) * D;
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
% [S_k30 S_k11; 
% S_k20 S_k21;
% S_k30*S_k11 S_k20*S_k21];

delta_mu = 16*S_k11*S_k30*(S_k30*S_k11-S_k20*S_k21);

[fprintf('S_k30=%20.20f',S_k30),fprintf(',   S_k11=%20.20f\n',S_k11),...
 fprintf('S_k20=%20.20f',S_k20),fprintf(',   S_k21=%20.20f\n',S_k21)];
[fprintf('delta_mu=%20.20f',delta_mu),fprintf(',   S_k30*S_k11=%20.20f',S_k30*S_k11),fprintf(',   S_k20*S_k21=%20.20f\n',S_k20*S_k21)];
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

save([results_path 'SSBNormalForm'],'a_star','c','du','dv','S_k30','S_k11','S_k20','S_k21','k')

%%
a_star;
a_star1=a_star;
a_star = roundn(a_star,-qu_n);
a_star2=a_star;
a_star = a_star + perturbating_mu;

a_star3=a_star;
plot(c,a_star3,'r*')

T_ks = zeros(N,1);
J_ks = zeros(N,1);

a11 = -theta;
a12 = -1;
a21 = c;
a22 = - a_star*c;

for j=1:N
	Lambda = eigenvalues(j);
    T_ks(j) = -(du+dv)*Lambda - (a11 + a22);
    J_ks(j) = du*dv*Lambda^2 + (dv*a11 + du*a22)*Lambda + a11*a22 - a12*a21;
end

disp( ['T_k: ' num2str( [length(find(T_ks<0))] ) ' ---> ' num2str( [find(T_ks<0)'] ) ] )
disp( ['J_k: ' num2str( [length(find(J_ks<0))] ) ' ---> ' num2str( [find(J_ks<0)'] ) ] )
format short

close()

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
plot(mu0_right, x0_right,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle','--','Color','k');

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
xticklabels = arrayfun(@(x) sprintf('%.3f', x), xticks, 'UniformOutput', false);
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
ax.XAxis.FontSize = ax_FontSize;  % 
ax.XAxis.FontName = 'Times New Roman';  % 
ax.YAxis.FontSize = ax_FontSize;  % 
ax.YAxis.FontName = 'Times New Roman';  % 
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

figure_name = [results_path 'Figure Step02 SSBNF diagram brain.eps'];
saveas(gcf, figure_name, 'epsc');
