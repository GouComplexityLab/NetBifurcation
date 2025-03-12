clear all
close all
clc

tic
digits = 50;
format long

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

%% set diffusion rate
theta = -1;
du = 0.3; 
fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.6f\n', -theta*du/(-theta+du*eigenvalues(2)));

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
weimu = -0.05;
% %% set diffusion rate
% theta = -2;
% du = 0.04; 
% fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.6f\n', -theta*du/(-theta+du*eigenvalues(2)));
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
% weimu = -0.08;

% %% set diffusion rate
% theta = -1/2;
% du = 0.01; 
% fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.6f\n', -theta*du/(-theta+du*eigenvalues(2)));
% 
% for i=1:N-1
%     panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
% end
% disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])
% % dv = panbie_THB(8);
% dv = (panbie_THB(7)+panbie_THB(8))/2;
% % dv = 0.06;
% c_set = 1;
% c = [0:0.001:2];  
% ymax = 3;
% weimu = -0.08;
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
fprintf('critical THB point a=%.6f, c=%.6f, on Turing curve k=%d, Hopf curve k=%d.\n', intersect_a_H1Tend, intersect_c_H1Tend, size(intersecting_points_TT_upper_H1,2)+1,1);

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
fprintf('critical SSB point a=%.6f, with setting c=%.6f, on Turing curve k=%d.\n', a_critical_SSB, c, p_critical_SSB+1);


a_critical_HB = (-theta + eigenvalues(1) * (du + dv))./c; 
fprintf('critical HB point a=%.6f, with setting c=%.6f, on Hopf curve k=%d.\n', a_critical_HB, c, 1);
plot(c,a_critical_HB,'b*')

a_star = a_critical_HB;
format rat
a_star
format long

k = 1;
format long

a11 = -theta;
a12 = -1;
a21 = c;
a22 = - a_star*c;

J = [a11,a12; a21,a22]; 
D = [du,0;0,dv];

Mk = J + eigenvalues(k) * D;
% eig(Mk)
T_k = -(du+dv)*eigenvalues(k) - (a11 + a22);
J_k = du*dv*eigenvalues(k)^2 + (dv*a11 + du*a22)*eigenvalues(k) + a11*a22 - a12*a21;

omega_k = sqrt(J_k);

pk1 = 1; 
pk2 = (1i*omega_k - a11 - du*eigenvalues(k))/a12;
pk = [pk1; pk2];

qk1 = (1i*omega_k + a11 + du*eigenvalues(k))/(2*1i*omega_k);
qk2 = a12/(2*1i*omega_k);
qk = [qk1; qk2];

[pk, conj(pk)].'*[qk, conj(qk)] % check!

f200 = [2*theta + 2; 0];
f020 = [0; 0];
f002 =[ 0; 0];  

f110 = [0; 0];
f101 = [0; 0];
f011 = [0; -c];

%  B_k1 = ((1-c^2)*qk1 + (1-c)^2*qk2) * pk1
% B_k1 = pk1*(qk.'*f101) + pk2*(qk.'*f011)
B_k1 = qk.'*(f101*pk1 + f011*pk2);
nv_k1 = real(B_k1);
 
f300 = [-6; 0];
f030 = [0; 0];
f003 = [0; 0];

f120 = [0; 0];
f210 = [0; 0];
f012 = [0; 0];
f021 = [0; 0];
f102 = [0; 0];
f201 = [0; 0];


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
    Mj = J + eigenvalues(j) * D;
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
-nv_k1/nv_k2
[fprintf('nv_k1=%20.20f',nv_k1),fprintf(',    nv_k2=%20.20f\n',nv_k2)];
fprintf('-nv_k1/nv_k2=%.6f.\n', -nv_k1/nv_k2);

save([results_path 'Hom_HBNormalForm'],'a_star','c','du','dv','nv_k1','nv_k2','k')

%%
% weimu = -0.0002;
a_star = a_star + weimu;
plot(c,a_star,'r*')
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

% % small_matrices_eigenvalues(1,:)
% % small_matrices_eigenvalues(2,:)
% % small_matrices_eigenvalues(3,:)
% % small_matrices_eigenvalues(4,:)


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
plot(mu0_right, x0_right,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle',':','Color','k');

mu_up_right = linspace(0,weimu,1000);
rho_up_right = real(sqrt(-(nv_k1*mu_up_right)/nv_k2));
plot(mu_up_right, rho_up_right,'MarkerSize',8,'LineWidth',gwLineWidth,'LineStyle','-','Color','k');

ax = gca;
% xticks = xticks(ax1);
xticks = [weimu,0,-weimu];
xticklabels = arrayfun(@(x) sprintf('%.4f', x), xticks, 'UniformOutput', false);
xticklabels(2) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-0.1,0,0.6];
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

ylim([-0.1,0.6])
xlim([weimu -weimu])

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'linewidth',1.0,'layer','top');

figure_name = [results_path 'Figure Step02 HBNF diagram npbcwsnetworks_n=10k=2p=0number=1.eps'];
saveas(gcf, figure_name, 'epsc');

