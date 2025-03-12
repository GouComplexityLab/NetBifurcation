clear all
% close all
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
% c_set = 1.5;
c = [0:0.001:3];  
ymax = 2;


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
% % grid minor
% 
% for i=1:length(p_H)
%     plot(c,Ha(i,:),'b-');
% end
% for i=1:5
%     plot(c,Ta(i,:),'r-');
% end
% 
% plot(linspace(0,3,100),ones(size(linspace(0,3,100))),'k')
% 
% plot(ones(size(linspace(0,ymax,100))),linspace(0,ymax,100),'k')
% ylim([0 ymax])
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

c_sk = zeros(length(p_T)-1,1);
mhs_k_kplus1 = zeros(length(p_T)-1,1);
css_k_kplus1 = zeros(length(p_T)-1,1);
ass_k_kplus1 = zeros(length(p_T)-1,1);
for i=1:length(p_T)-1
    c_sk(i) = -dv * eigenvalues(p_T(i)) * (1 + du*eigenvalues(p_T(i)));
    css_k_kplus1(i) = dv/du * (1 + du*eigenvalues(p_T(i))) * (1 + du*eigenvalues(p_T(i+1)));
    ass_k_kplus1(i) = (1 + du*( eigenvalues(p_T(i)) + eigenvalues(p_T(i+1)) ) ) / (1 + du*eigenvalues(p_T(i))) / (1 + du*eigenvalues(p_T(i+1)));
    mhs_k_kplus1(i) = (1 + du*( eigenvalues(p_T(i)) + eigenvalues(p_T(i+1)) ) );
end
%% plot
gwLineWidth = 0.5;
figure
% set(gcf,"Position",[300 300 400 400])
% axes1 = axes('Position',[0.16625 0.21375 0.75 0.75]);
set(gcf,"Position",[300 200 700 530])
axes1 = axes('Position',[0.118571428571429 0.16207527759708 0.845 0.781792646931222]);

hold on
% grid minor

% plot(c,Ha(1,:),'b-');

% for i=1:length(p_H)
%     plot(c,Ha(i,:),'b-');
% end
% % for i=1:length(p_T)
% for i=1:length(p_intersect_upper_x)
% % for i=1:size(intersecting_points_TT_upper_H1,2)+1
%     plot(c,Ta(i,:),'r-');
% end
% plot(intersecting_points_TT(1,:), intersecting_points_TT(2,:), 'g*')
% plot(intersecting_points_TT_upper_x(1,:), intersecting_points_TT_upper_x(2,:), 'co')
% plot(intersecting_points_TT_upper_H1(1,:), intersecting_points_TT_upper_H1(2,:), 'k*')

fprintf('critical THB point a=%.6f, c=%.6f, on Turing curve k=%d, Hopf curve k=%d.\n', intersect_a_H1Tend, intersect_c_H1Tend, size(intersecting_points_TT_upper_H1,2)+1,1);

plot(c,ones(size(c))*1/(-theta),'LineWidth',gwLineWidth,'Color','k')
% plot(c,zeros(size(c)),'k-')

%% SSB
p_SSB = find(mhs_k_kplus1>0);

M = 9;
dq=lines(M); 

k_SSB = 1;
c_SSB = linspace(css_k_kplus1(k_SSB), 3,100);
a_SSB = 1/(du*eigenvalues(p_T(k_SSB)) -theta) + 1./c_SSB * dv * eigenvalues(p_T(k_SSB));
plot(c_SSB,a_SSB,'LineWidth',gwLineWidth,'color',dq(k_SSB,:))

for k_SSB = 2:7
c_SSB = linspace(css_k_kplus1(k_SSB), css_k_kplus1(k_SSB-1),100);
a_SSB = 1/(du*eigenvalues(p_T(k_SSB)) -theta) + 1./c_SSB * dv * eigenvalues(p_T(k_SSB));
plot(c_SSB,a_SSB,'LineWidth',gwLineWidth,'color',dq(k_SSB,:))
end

%% HB
i=1; k=3;
c_HS_i_k = (1+(du+dv)*eigenvalues(i)-dv*eigenvalues(k))*(1+du*eigenvalues(k));
st = 1; et = 3;
x=linspace(st,et,100);
y = (-theta + eigenvalues(i) * (du + dv))./x; 
plot(x,y,'LineWidth',gwLineWidth,'LineStyle','--','Color','b')

i=1; k=3;
c_HS_i_k = (1+(du+dv)*eigenvalues(i)-dv*eigenvalues(k))*(1+du*eigenvalues(k));
st = c_HS_i_k; et = 3;
x=linspace(st,et,100);
y = (-theta + eigenvalues(i) * (du + dv))./x; 
plot(x,y,'LineWidth',gwLineWidth,'Color','b')

i=2; k=7;
c_HS_i_k = (1+(du+dv)*eigenvalues(i)-dv*eigenvalues(k))*(1+du*eigenvalues(k));
st = c_HS_i_k; et = 3;
x=linspace(st,et,100);
y = (-theta + eigenvalues(i) * (du + dv))./x; 
% plot(x,y,'LineWidth',gwLineWidth,'Color','b')

%% Turing
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

a_slice = linspace(0,2,100);
c_slice = 1.8 * ones(size(a_slice));
plot(c_slice,a_slice,'LineWidth',gwLineWidth,'LineStyle','--','Color','k')

plot(intersect_c_H1Tend,intersect_a_H1Tend, ...
    'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',12,...
    'Marker','pentagram','LineStyle','none','Color','r');
plot(css_k_kplus1,ass_k_kplus1,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,...
    'Marker','.','LineStyle','none','Color','k');
% axis equal
box on


ax = gca;
% xticks = xticks(ax);
xticks = [0,0.8,1.8,2.8];
xticklabels = arrayfun(@(x) sprintf('%.1f', x), xticks, 'UniformOutput', false);
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
ax.XAxis.FontSize = ax_FontSize;  %
ax.XAxis.FontName = 'Times New Roman';  %
ax.YAxis.FontSize = ax_FontSize;  %
ax.YAxis.FontName = 'Times New Roman';  %
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'in';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

% xlabel('$c$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% ylabel('$a$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');

ylim([0 ymax])

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'linewidth',1.0,'layer','top');


figure_name = [results_path 'Figure Step01 BifurcationDiagram npbcwsnetworks_n=10k=2p=0number=1.eps'];
saveas(gcf, figure_name, 'epsc');

