clear all
close all
clc

tic
network_name = 'ws_n=10k=2p=0.05number=796';
datanetwork_path = ['./network2 ' network_name '/'];
results_path = strrep(datanetwork_path,'network2','RDPPmodel results network2');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

load([results_path 'THBNormalForm'],'alpha11','alpha12','alpha21','alpha22','kappa30','kappa12','kappa21','kappa03','kappa021','kappa022','kappa023','n','c','d1','d2','a_star','b_star','laplacian_matrix')

format rat
b_star
a_star
format long

many = 1;

mu1 = 0.05; mu2 =  0.04;
mu1 = -0.06; mu2 =  -0.06;

mu1s = linspace(-0.05,0.05,5);
mu2s = linspace(-0.05,0.05,5);

%%
count_num = 5e3;


mu1_left = -0.07; mu1_right = 0.07;
mu2_left = -0.07; mu2_right = 0.07;

% mu1_left = -0.07; mu1_right = 0;
% mu2_left = -0.07; mu2_right = 0;

mu1_values = linspace(mu1_left,mu1_right,count_num+1);
mu2_values = linspace(mu2_left,mu2_right,count_num+1);
[mu1_values,mu2_values]=meshgrid(mu1_values,mu2_values);
mu10 = mu1_values(:);
mu20 = mu2_values(:);


%%

g1 = (kappa021+kappa022.*mu1_values+kappa023.*mu2_values).^2 - 4*kappa03*(alpha21.*mu1_values+alpha22.*mu2_values);
v = [-50000,0,50000];
p1 = find(g1>=0);
%%

g2 = - 4*kappa30*(alpha11.*mu1_values+alpha12.*mu2_values);
v = [-50000,0,50000];
p2 = find(g2>=0);

%% 
g3 =  kappa30^2*(kappa021+kappa022*mu1_values+kappa023*mu2_values).^2 - ...
    4*(kappa12*kappa21-kappa30*kappa03)*( kappa21*(alpha11*mu1_values+alpha12*mu2_values) - kappa30*(alpha21*mu1_values+alpha22*mu2_values) );
p3_zheng = find(g3>=0);

g41 =  kappa30*kappa12*(kappa021+kappa022*mu1_values+kappa023*mu2_values).^2 + 2 * (kappa12*kappa21-kappa30*kappa03) * ( kappa12*(alpha21*mu1_values+alpha22*mu2_values) - kappa03*(alpha11*mu1_values+alpha12*mu2_values) );
g42 = kappa12*(kappa021+kappa022*mu1_values+kappa023*mu2_values).*sqrt( kappa30^2*(kappa021+kappa022*mu1_values+kappa023*mu2_values).^2 - ...
    4*(kappa12*kappa21-kappa30*kappa03)*( kappa21*(alpha11*mu1_values+alpha12*mu2_values) - kappa30*(alpha21*mu1_values+alpha22*mu2_values) ) );
g43 = kappa12*(kappa021+kappa022*mu1_values+kappa023*mu2_values);

g45 = ( kappa30*kappa12*(kappa021+kappa022*mu1_values+kappa023*mu2_values).^2 + 2 * (kappa12*kappa21-kappa30*kappa03) * ( kappa12*(alpha21*mu1_values+alpha22*mu2_values) - kappa03*(alpha11*mu1_values+alpha12*mu2_values) ) ).^2 ...
    - ( kappa12*(kappa021+kappa022*mu1_values+kappa023*mu2_values).*sqrt( kappa30^2*(kappa021+kappa022*mu1_values+kappa023*mu2_values).^2 - ...
    4*(kappa12*kappa21-kappa30*kappa03)*( kappa21*(alpha11*mu1_values+alpha12*mu2_values) - kappa30*(alpha21*mu1_values+alpha22*mu2_values) ) ) ).^2;

p41_zheng = find(g41>0);
p41_fu = find(g41<0);

p43_zheng = find(g43>0);
p43_fu = find(g43<0);

p45_zheng = find(g45>0);
p45_fu = find(g45<0);

% case 1
zheng_p3_jiao_p41_fu_jiao_p43_fu = intersect(intersect(p3_zheng, p41_fu), p43_fu);
% case 2
zheng_p3_jiao_p41_fu_jiao_p43_zheng_p45_zheng = ...
    intersect(intersect(p3_zheng, p41_fu), intersect(p43_zheng, p45_zheng));
% case 3
zheng_p3_jiao_p41_zheng_jiao_p43_fu_p45_fu = ...
    intersect(intersect(p3_zheng, p41_zheng), intersect(p43_fu, p45_fu));

% case 1
fu_p3_jiao_p41_fu_jiao_p43_zheng = intersect(intersect(p3_zheng, p41_fu), p43_zheng);
% case 2
fu_p3_jiao_p41_fu_jiao_p43_fu_p45_zheng = ...
    intersect(intersect(p3_zheng, p41_fu), intersect(p43_fu, p45_zheng));
% case 3
fu_p3_jiao_p41_zheng_jiao_p43_zheng_p45_fu = ...
    intersect(intersect(p3_zheng, p41_zheng), intersect(p43_zheng, p45_fu));


%%
p_zheng1 = zheng_p3_jiao_p41_fu_jiao_p43_fu;% case 1
p_zheng2 = zheng_p3_jiao_p41_fu_jiao_p43_zheng_p45_zheng;% case 2
p_zheng3 = zheng_p3_jiao_p41_zheng_jiao_p43_fu_p45_fu;% case 3
p_fu1 = fu_p3_jiao_p41_fu_jiao_p43_zheng;% case 1
p_fu2 = fu_p3_jiao_p41_fu_jiao_p43_fu_p45_zheng;% case 2
% p_fu3 = fu_p3_jiao_p41_zheng_jiao_p43_zheng_p45_fu;% case 3
% 
p_zheng = union(union(p_zheng1,p_zheng2),p_zheng3);
p_fu = union(p_fu1,p_fu2);



%%
p_two = intersect(p_zheng, p_fu);
p_one = union(setdiff(p_zheng, p_two),setdiff(p_fu, p_two));


h41 = g42+g41;
h42 = -g42+g41;

h41(find(g3<0))=1;
h42(find(g3<0))=1;

gwLineWidth = 1.0;

figure
set(gcf,"Position",[300 200 700 530])
axes1 = axes('Position',[0.118571428571429 0.16207527759708 0.845 0.781792646931222]);

% set(fig,'visible','off')
hold on

contour(mu1_values,mu2_values,g1,'LevelList',0,'color','g','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,g2,'LevelList',0,'color','b','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,h41,'LevelList',0,'color','m','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,h42,'LevelList',0,'color','r', 'LineStyle','-','ShowText','off','LineWidth',gwLineWidth);


mu1 = 0.05; mu2 = 0.05; % 1
p1 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','k','LineStyle','none'); % 1
mu1 = 0.05; mu2 = 0.0474; % 2
p2 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','r','LineStyle','none'); % 2
mu1 = 0.05; mu2 = 0.046; % 3
p3 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','g','LineStyle','none'); % 3
mu1 = -0.05; mu2 = -0.0481; % 4
p4 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','b','LineStyle','none'); % 4
mu1 = -0.05; mu2 = -0.04795; % 5
p5 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','k','LineWidth',gwLineWidth,'LineStyle','none'); % 5
mu1 = -0.05; mu2 = -0.0478; % 6
p6 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','r','LineWidth',gwLineWidth,'LineStyle','none'); % 6
mu1 = -0.05; mu2 =  -0.0475; % 7
p7 = plot(mu1,mu2,'MarkerSize',15,'Marker','*','Color','m','LineWidth',gwLineWidth,'LineStyle','none'); % 6

x = [-0.07, 0.07];
y = zeros(size(x));
plot(x,y,'MarkerSize',15,'Marker','none','Color','k','LineWidth',0.5,'LineStyle',':'); % 6

y = [-0.07, 0.07];
x = zeros(size(y));
plot(x,y,'MarkerSize',15,'Marker','none','Color','k','LineWidth',0.5,'LineStyle',':'); % 6

ax = gca;
% xticks = xticks(ax1);
xticks = [-0.06:0.02:0.06];
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
xticklabels(find(xticks==0)) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-0.06:0.02:0.06];
yticklabels = arrayfun(@(x) sprintf('%.2f', x), yticks, 'UniformOutput', false);
yticklabels(find(yticks==0)) = {'0'};
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


xlim([-0.07 0.07])
ylim([-0.07 0.07])

% xlabel('\mu_{1}')
% ylabel('\mu_{2}')

% axis off

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.01 0.01],'linewidth',1.0,'layer','top');
%% 

%% 2
axes('Position',[0.727857142857143 0.400943396226415 0.216428571428571 0.262264150943398]);

hold on
contour(mu1_values,mu2_values,g1,'LevelList',0,'color','g','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,g2,'LevelList',0,'color','b','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,h41,'LevelList',0,'color','m','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,h42,'LevelList',0,'color','r', 'LineStyle','-','ShowText','off','LineWidth',gwLineWidth);


mu1 = 0.05; mu2 = 0.05; % 1
p1 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','k','LineStyle','none'); % 1
mu1 = 0.05; mu2 = 0.0474; % 2
p2 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','r','LineStyle','none'); % 2
mu1 = 0.05; mu2 = 0.046; % 3
p3 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','g','LineStyle','none'); % 3
mu1 = -0.05; mu2 = -0.0481; % 4
p4 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','b','LineStyle','none'); % 4
mu1 = -0.05; mu2 = -0.04795; % 5
p5 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','k','LineWidth',gwLineWidth,'LineStyle','none'); % 5
mu1 = -0.05; mu2 = -0.0478; % 6
p6 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','r','LineWidth',gwLineWidth,'LineStyle','none'); % 6
mu1 = -0.05; mu2 =  -0.0475; % 7
p7 = plot(mu1,mu2,'MarkerSize',15,'Marker','*','Color','m','LineWidth',gwLineWidth,'LineStyle','none'); % 6


ax = gca;
% xticks = xticks(ax1);
xticks = [0.05];
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
xticklabels(find(xticks==0)) = {''};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [0.043:0.001:0.066];
yticklabels = arrayfun(@(x) sprintf('%.3f', x), yticks, 'UniformOutput', false);
yticklabels(find(yticks==0)) = {''};
% yticklabels(end) = {'$a$'};
set(gca, 'ytick', yticks);
set(gca, 'yticklabel', yticklabels);

ax_FontSize = 12;
ax.XAxis.FontSize = ax_FontSize;  
ax.XAxis.FontName = 'Times New Roman';  
ax.YAxis.FontSize = ax_FontSize;  
ax.YAxis.FontName = 'Times New Roman';  
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

set(gca,'XAxisLocation','bottom');
% set(gca,'YAxisLocation','origin');

xlim([0.0480255917928736 0.0518645574974252]);
ylim([0.0452849962835435 0.0506959980261086]);

xtickangle(15)
ytickangle(15)

grid on
box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'linewidth',0.5,'layer','top');

%% 3
axes('Position',[0.142142857142857 0.410377358490566 0.212857142857142 0.270754716981132]);

hold on
contour(mu1_values,mu2_values,g1,'LevelList',0,'color','g','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,g2,'LevelList',0,'color','b','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,h41,'LevelList',0,'color','m','ShowText','off','LineWidth',gwLineWidth);
contour(mu1_values,mu2_values,h42,'LevelList',0,'color','r', 'LineStyle','-','ShowText','off','LineWidth',gwLineWidth);


mu1 = 0.05; mu2 = 0.05; % 1
p1 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','k','LineStyle','none'); % 1
mu1 = 0.05; mu2 = 0.0474; % 2
p2 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','r','LineStyle','none'); % 2
mu1 = 0.05; mu2 = 0.046; % 3
p3 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','g','LineStyle','none'); % 3
mu1 = -0.05; mu2 = -0.0481; % 4
p4 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','b','LineStyle','none'); % 4
mu1 = -0.05; mu2 = -0.04795; % 5
p5 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','k','LineWidth',gwLineWidth,'LineStyle','none'); % 5
mu1 = -0.05; mu2 = -0.0478; % 6
p6 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','r','LineWidth',gwLineWidth,'LineStyle','none'); % 6
mu1 = -0.05; mu2 =  -0.0475; % 7
p7 = plot(mu1,mu2,'MarkerSize',15,'Marker','*','Color','m','LineWidth',gwLineWidth,'LineStyle','none'); % 6


ax = gca;
% xticks = xticks(ax1);
xticks = [-0.05];
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
xticklabels(find(xticks==0)) = {''};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-0.06:0.0002:-0.02];
yticklabels = arrayfun(@(x) sprintf('%.4f', x), yticks, 'UniformOutput', false);
yticklabels(find(yticks==0)) = {''};
% yticklabels(end) = {'$a$'};
set(gca, 'ytick', yticks);
set(gca, 'yticklabel', yticklabels);

ax_FontSize = 12;
ax.XAxis.FontSize = ax_FontSize;  
ax.XAxis.FontName = 'Times New Roman';  
ax.YAxis.FontSize = ax_FontSize;  
ax.YAxis.FontName = 'Times New Roman';  
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

set(gca,'XAxisLocation','top');
set(gca,'YAxisLocation','right');


xlim([-0.0505168544610138 -0.0495165254909654]);
ylim([-0.0482726808867528 -0.0472723519167044]);

xtickangle(15)
ytickangle(35)
grid on
box on

set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'linewidth',0.5,'layer','top');

figure_name = [results_path 'Figure Step02 THB region.eps'];
saveas(gcf, figure_name, 'epsc');
