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
matcont_path = [results_path 'matcont NetFHNsystem ' network_name ' ssb'];
if ~exist(matcont_path, 'dir')
    mkdir(matcont_path);
end
%%
cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); % 行和
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
a_star = a_critical_SSB
c
k = p_critical_SSB+1
format long


%%
n = k;
phin = real(eigenvectors(:,n));

%%
mu0 = 0;
a0 = a_star + mu0;
u_star0 = 0; v_star0 = 0;

mu = 0.001;
a_first = a_star + mu;
u_star = 0; v_star = 0;

% % % x0=zeros(1,2*N);
% % % x0(1:2:end-1) = u_star + phin* 0.05;
% % % x0(2:2:end  ) = v_star + phin* 0.05;
% % % 
% % % 
% % % hls=feval(@FHNSystemOnBrain);
% % % options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(2*N,1));
% % % % options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
% % % [t,y]=ode45(hls{2},[0:0.1:1e4],x0,options, a_first, c, du, dv);
% % % 
% % % simulated_u1 = y(:, 1:2:end);
% % % simulated_v1 = y(:, 2:2:end);
% % % 
% % % x0(1:2:end-1) = u_star - phin* 0.05;
% % % x0(2:2:end  ) = v_star - phin* 0.05;
% % % 
% % % hls=feval(@FHNSystemOnBrain);
% % % options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(2*N,1));
% % % % options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
% % % [t,y]=ode45(hls{2},[0:0.1:1e4],x0,options, a_first, c, du, dv);
% % % 
% % % simulated_u2 = y(:, 1:2:end);
% % % simulated_v2 = y(:, 2:2:end);
%%


%%

load([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
load([matcont_path '/data_matcont2.mat'], 'x2', 'v2', 's2', 'h2', 'f2')
load([matcont_path '/data_matcont3.mat'], 'x3', 'v3', 's3', 'h3', 'f3')
load([matcont_path '/data_matcont4.mat'], 'x4', 'v4', 's4', 'h4', 'f4')

% gwcplanim2(x1,v1,s1,'b',[2*N+1,1]);
% gwcplanim(x2,v2,s2,'c',[2*N+1, 1]);
% gwcplanim(x3,v3,s3,'r',[2*N+1,1]);
% gwcplanim(x4,v4,s4,'m',[2*N+1,1]);


gwlinewidth = 1.0;
gwcolors = hsv(41);
gwcolors = [gwcolors;hsv(41)];
% gwcolors(3,:) = [];


figure
set(gcf,"Position",[300 300 880 380])
set(gca,'Position',[0.103807033601355 0.0849389137570983 0.85187801187373 0.884632927996743]);
hold on


% plot(a0, u_star0, 'k.')

Node_pos = 10;
for Node_pos=[1:82]
% for Node_pos=[1:41]+41
% for Node_pos=[1:41]
% for Node_pos=[3, 8]


% gwcplanim2(x1,v1,s1,'b',[2*N+1,1]);
% gwcplanim(x2,v2,s2,'c',[2*N+1, 1]);
% gwcplanim(x3,v3,s3,'r',[2*N+1,1]);
% gwcplanim(x4,v4,s4,'m',[2*N+1,1]);

%% right-left
% gwcplanim2(x1,v1,s1,'b',[2*N+1,2*(Node_pos-1)+1]);
unode = x1(2*(Node_pos-1)+1,:); 
para = x1(2*N+1,:);
pos1_in = [s1(1).index:s1(2).index];
cont1_para = para(pos1_in); cont1_unode = unode(pos1_in);
% plot(cont1_para,cont1_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
% gwcplanim2(x2,v2,s2,'c',[2*N+1, 2*(Node_pos-1)+1]);
unode = x2(2*(Node_pos-1)+1,:); 
para = x2(2*N+1,:);
pos1_in = [s2(1).index:s2(3).index];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
% plot(cont2_para,cont2_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');

% cont_para = [cont1_para, cont2_para];
% cont_unode = [cont1_unode, cont2_unode];
% cont_para = [cont1_para, fliplr(cont2_para)];
% cont_unode = [cont1_unode, fliplr(cont2_unode)];
cont_para = [fliplr(cont1_para), cont2_para];
cont_unode = [fliplr(cont1_unode), cont2_unode];

% plot(cont_para,cont_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
pos_right_u_star = find(cont_para>=a0);
right_unode = [cont_unode(pos_right_u_star),u_star0,];
right_para = [cont_para(pos_right_u_star),a0,];
% % % plot(right_para,right_unode,'Marker','none',...
% % %     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

pos_left_u_star = find(cont_para<=a0);
left_unode = [cont_unode(pos_left_u_star),u_star0,];
left_para = [cont_para(pos_left_u_star),a0,];
% % % plot(left_para,left_unode,'Marker','none',...
% % %     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');


%% up-down
% gwcplanim2(x3,v3,s3,'r',[2*N+1,1]);
unode = x3(2*(Node_pos-1)+1,:); 
para = x3(2*N+1,:);
pos1_in = [s3(1).index+1:s3(2).index];
cont1_para = para(pos1_in); cont1_unode = unode(pos1_in);
% plot(cont1_para,cont1_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
% gwcplanim(x4,v4,s4,'m',[2*N+1,1]);
unode = x4(2*(Node_pos-1)+1,:); 
para = x4(2*N+1,:);
pos1_in = [s4(1).index+1:s4(2).index];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
% plot(cont2_para,cont2_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

% cont_para = [cont1_para, cont2_para];
% cont_unode = [cont1_unode, cont2_unode];
% cont_para = [cont1_para, fliplr(cont2_para)];
% cont_unode = [cont1_unode, fliplr(cont2_unode)];
cont_para = [fliplr(cont1_para), cont2_para];
cont_unode = [fliplr(cont1_unode), cont2_unode];
% % % plot(cont_para,cont_unode,'Marker','none',...
% % %     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));




% plot(right_para,right_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
% plot(left_para,left_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
% plot(cont_para,cont_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));
% 
% xlim([a_star-0.003 a_star+0.003])
% ylim([-0.1 0.1])

% plot3(Node_pos*ones(size(right_para)),right_para,right_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
% plot3(Node_pos*ones(size(left_para)),left_para,left_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
% plot3(Node_pos*ones(size(cont_para)),cont_para,cont_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));
% ylim([a_star-0.003 a_star+0.003])
% zlim([-0.1 0.1])

plot3(right_para,Node_pos*ones(size(right_para)),right_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
plot3(left_para,Node_pos*ones(size(left_para)),left_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
plot3(cont_para,Node_pos*ones(size(cont_para)),cont_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));

xlim([a_star-0.005 a_star+0.005])
% ylim([0 42])
ylim([41 83])
ylim([0 83])
zlim([-0.08 0.08])
view([75 30])

end



% xlabel('a','FontSize',18,'Interpreter','latex')
% ylabel('$u_{1}$','FontSize',18,'Interpreter','latex')

ax = gca;
xticks = xticks(ax);
xticks = [0.952,0.956];
xticklabels = arrayfun(@(x) sprintf('%.3f' , x), xticks, 'UniformOutput', false);
% xticklabels(1) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

yticks = yticks(ax);
% yticks = [1,10,20,30,41]+41;
yticks = [1,10,20,30,50,60,70,82];
yticklabels = arrayfun(@(x) sprintf('%.0f', x), yticks, 'UniformOutput', false);
% yticklabels(1) = {'0'};
% yticklabels(end) = {'$a$'};
set(gca, 'ytick', yticks);
set(gca, 'yticklabel', yticklabels);

zticks = zticks(ax);
% zticks = [0.04:0.02:0.15];
zticklabels = arrayfun(@(x) sprintf('%.2f', x), zticks, 'UniformOutput', false);
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

set(gca,'Box','on');
% set(gca,'Box','on','BoxStyle','full');
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'linewidth',1.0,'layer','top');


figure_name = [results_path '/Figure step03 SSB 3d plot.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = [results_path '/Figure step03 SSB 3d plot.jpg'];
saveas(gcf, figure_name);

