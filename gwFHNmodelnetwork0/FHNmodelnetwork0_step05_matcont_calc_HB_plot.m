clear all
close all
clc

addpath('./system/');

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
%%
matcont_path = [results_path 'matcont NetFHNsystem ' network_name ' hb'];
if ~exist(matcont_path, 'dir')
    mkdir(matcont_path);
end
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

%% set diffusion rate
theta = -1;
du = 0.3; 
fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.8f\n', -theta*du/(-theta+du*eigenvalues(2)));

for i=1:N-1
    panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
end
disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])

dv = 0.8;
c_set = 1.8;
c = [0:0.001:3];  
ymax = 2;

%%
panbie_H = -theta + eigenvalues * (du + dv);
p_H = find(panbie_H>0);
Ha = zeros(length(p_H),length(c));
for i=1:length(p_H)
    Ha(i,:) = (-theta + eigenvalues(p_H(i)) * (du + dv))./c; 
end

%%
panbie_T = -theta + eigenvalues * du;
p_T = find(panbie_T>0);
Ta = zeros(length(p_T),length(c));
for i=1:length(p_T)
    Ta(i,:) = 1/(du*eigenvalues(p_T(i)) -theta) + 1./c * dv * eigenvalues(p_T(i));
end

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

intersecting_points_TT_upper_H1 = intersecting_points_TT(:,p_intersect_upper_H1);

intersect_c_H1Tend = - du*dv*eigenvalues(p_intersect_upper_H1(end)+1)^2 - theta* (du-dv) * eigenvalues(p_intersect_upper_H1(end)+1) + theta^2 ;
intersect_a_H1Tend = 1/(du*eigenvalues(p_intersect_upper_H1(end)+1) -theta) + 1/intersect_c_H1Tend * dv * eigenvalues(p_intersect_upper_H1(end)+1);

%% 

c = c_set;
intersecting_points_c_TT_upper_H1 = intersecting_points_TT_upper_H1(1,:);
disp(['intersecting SSB points for c:   ' num2str(intersecting_points_c_TT_upper_H1)])
% [v p_critical_SSB] = min(abs(intersecting_points_c_TT_upper_H1-c))
p_critical_SSB = find(intersecting_points_c_TT_upper_H1-c > 0);
p_critical_SSB = p_critical_SSB(end);
a_critical_SSB = 1/(du*eigenvalues(p_critical_SSB + 1) -theta) + 1./c * dv * eigenvalues(p_critical_SSB + 1);

fprintf('critical SSB point a=%.8f, with setting c=%.8f, on Turing curve k=%d.\n', a_critical_SSB, c, p_critical_SSB+1);


a_critical_HB = (-theta + eigenvalues(1) * (du + dv))./c; 
fprintf('critical HB point a=%.8f, with setting c=%.8f, on Hopf curve k=%d.\n', a_critical_HB, c, 1);


format rat
a_star = a_critical_HB
c
k = 2;
format long


%%
n = k;
phin = real(eigenvectors(:,n));

%%
mu0 = 0;
a0 = a_star + mu0;
u_star0 = 0; v_star0 = 0;

mu = -0.03; % mu = -0.05;
a_first = a_star + mu;
u_star = 0; v_star = 0;

% x0=zeros(1,2*N);
% x0(1:2:end-1) = u_star + phin* 0.05;
% x0(2:2:end  ) = v_star + phin* 0.05;
% 
% % x0(1:2:end-1) = u_star - phin* 0.05;
% % x0(2:2:end  ) = v_star - phin* 0.05;
% 
% hls=feval(@FHNSystemOnBrain);
% options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(2*N,1));
% % options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
% [t,y]=ode45(hls{2},[0:0.1:1e4],x0,options, a_first, c, du, dv);
% 
% simulated_u1 = y(:, 1:2:end);
% simulated_v1 = y(:, 2:2:end);


%%

% initial_x0=zeros(1,2*N)';
% initial_x0(1:2:end-1) = u_star + phin* 0;
% initial_x0(2:2:end  ) = v_star + phin* 0;

% initial_x0 = y(end,:)';


%%
%%


load([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
load([matcont_path '/data_matcontHB.mat'], 'xlc', 'vlc', 'slc', 'hlc', 'flc')


x=xlc;v=vlc;s=slc;

gwlinewidth = 1.0;
gwcolors = hsv(41);
gwcolors = [gwcolors;hsv(41)];
% gwcolors(3,:) = [];


figure
set(gcf,"Position",[300 300 880 380])
set(gca,'Position',[0.103807033601355 0.0849389137570983 0.85187801187373 0.884632927996743]);
hold on

%%
Node_pos=1;

for Node_pos=1:82


e=[size(xlc,1) 2*(Node_pos-1)+1 2*(Node_pos-1)+2];


gwnphase = N*2;
gwtps = (size(xlc,1)-2)/2/N;

x_p = zeros(1,size(x,2));
ymin_p = zeros(1,size(x,2));
ymax_p = zeros(1,size(x,2));

zmin_p = zeros(1,size(x,2));
zmax_p = zeros(1,size(x,2));

for j = 1:size(x,2)

    xx = x(e(1),j)*ones(gwtps,1);
    yy = x((0:gwtps-1)*gwnphase+e(2),j);
    zz = x((0:gwtps-1)*gwnphase+e(3),j);

    x_p(j) = xx(1);
    ymin_p(j) = min(yy);
    ymax_p(j) = max(yy);
    
    zmin_p(j) = min(zz);
    zmax_p(j) = max(zz);

end

right_para = linspace(a_star-0.1,a_star,100+1);
right_unode = zeros(size(right_para));
left_para = linspace(a_star,a_star+0.1,100+1);
left_unode = zeros(size(left_para));

% plot(right_para,right_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
% plot(left_para,left_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
% 
% plot(x_p,ymin_p,'r')
% plot(x_p,ymax_p,'r')

plot3(right_para,Node_pos*ones(size(right_para)),right_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
plot3(left_para,Node_pos*ones(size(left_para)),left_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color','k');

plot3(x_p,Node_pos*ones(size(x_p)),ymin_p,'linewidth',gwlinewidth,'Color',gwcolors(Node_pos,:))
plot3(x_p,Node_pos*ones(size(x_p)),ymax_p,'linewidth',gwlinewidth,'Color',gwcolors(Node_pos,:))


% xlim([a_star-0.05 a_star+0.05])
% zlim([-0.5 0.5])

xlim([a_star-0.05 a_star+0.05])
ylim([0 83])
zlim([-0.6 0.6])
view([75-180 25])

end


ax = gca;
xticks = xticks(ax);
xticks = [0.54,0.58];
xticklabels = arrayfun(@(x) sprintf('%.2f' , x), xticks, 'UniformOutput', false);
% xticklabels(1) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [1,10,20,30,50,60,70,82];
yticklabels = arrayfun(@(x) sprintf('%.0f', x), yticks, 'UniformOutput', false);
% yticklabels(1) = {'0'};
% yticklabels(end) = {'$a$'};
set(gca, 'ytick', yticks);
set(gca, 'yticklabel', yticklabels);

zticks = zticks(ax);
zticks = [-0.5,0,0.5];
zticklabels = arrayfun(@(x) sprintf('%.1f', x), zticks, 'UniformOutput', false);
zticklabels(find(zticks==0)) = {'0'};
% zticklabels(end) = {'$a$'};
set(gca, 'ztick', zticks);
set(gca, 'zticklabel', zticklabels);

ax_FontSize = 24;
ax.XAxis.FontSize = ax_FontSize;  % 
ax.XAxis.FontName = 'Times New Roman';  % 
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickLabelInterpreter = 'latex';


ax.YAxis.FontSize = ax_FontSize;  % 
ax.YAxis.FontName = 'Times New Roman';  % 
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickLabelInterpreter = 'latex';


ax.ZAxis.FontSize = ax_FontSize;  % 
ax.ZAxis.FontName = 'Times New Roman';  % 
ax.ZAxis.TickDirection = 'in';
ax.ZAxis.TickLabelInterpreter = 'latex';

set(gca, 'YDir', 'reverse');

set(gca,'Box','on');
% set(gca,'Box','on','BoxStyle','full');
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'linewidth',1.0,'layer','top');

figure_name = [results_path '/Figure step03 HB 3d plot.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = [results_path '/Figure step03 HB 3d plot.jpg'];
saveas(gcf, figure_name);





