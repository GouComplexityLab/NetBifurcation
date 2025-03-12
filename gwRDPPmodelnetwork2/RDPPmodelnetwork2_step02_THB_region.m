clear all
close all
clc

network_name = 'ws_n=10k=4p=0.05number=69';
datanetwork_path = ['./network3 ' network_name '/'];
results_path = strrep(datanetwork_path,'network3','RDPPmodel results network3');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

load([results_path 'THBNormalForm'],'alpha11','alpha12','alpha21','alpha22','kappa30','kappa12','kappa21','kappa03','kappa021','kappa022','kappa023','n','c','d1','d2','a_star','b_star','laplacian_matrix')


% mu1_left = -0.066; mu1_right = 0.066;
% mu2_left = -0.066; mu2_right = 0.066;

mu1_left = -0.07; mu1_right = 0.07;
mu2_left = -0.07; mu2_right = 0.07;

gwLineWidth = 1.0;

figure
set(gcf,"Position",[300 200 700 530])
axes1 = axes('Position',[0.118571428571429 0.16207527759708 0.845 0.781792646931222]);
% set(gcf,'Position',get(0,'ScreenSize'))
hold on
% plot([mu1_left mu1_right],[0 0],'k-.')
% plot([0 0],[mu2_left mu2_right],'k-.')
f1 = @(mu1,mu2)  -(alpha11*mu1+alpha12*mu2)/kappa30;
f2 = @(mu1,mu2)  -(alpha21*mu1+alpha22*mu2)/kappa03;
f3 = @(mu1,mu2) (kappa03*(alpha11*mu1+alpha12*mu2)-kappa12*(alpha21*mu1+alpha22*mu2))...
    /(-kappa30*kappa03+kappa21*kappa12); 
f4 = @(mu1,mu2) (kappa21*(alpha11*mu1+alpha12*mu2)-kappa30*(alpha21*mu1+alpha22*mu2))...
    /(kappa30*kappa03-kappa21*kappa12);


h1 = fimplicit(f1,[mu1_left mu1_right mu2_left mu2_right],'-g','LineWidth',gwLineWidth);
h2 = fimplicit(f2,[mu1_left mu1_right mu2_left mu2_right],'-b','LineWidth',gwLineWidth);
h3 = fimplicit(f3,[mu1_left 0 mu2_left mu2_right],'-m','LineWidth',gwLineWidth);
h4 = fimplicit(f4,[mu1_left 0 mu2_left mu2_right],'-r','LineWidth',gwLineWidth);

mu1 = 0.05; mu2 = 0.05; % 1
p1 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','k','LineStyle','none'); % 1
mu1 = 0.05; mu2 = 0.0475; % 2
p2 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','r','LineStyle','none'); % 2
mu1 = 0.05; mu2 = 0.0464; % 3
p3 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','g','LineStyle','none'); % 3
mu1 = -0.05; mu2 = -0.049; % 4
p4 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','b','LineStyle','none'); % 4
mu1 = -0.05; mu2 = -0.0483; % 5
p5 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','k','LineWidth',gwLineWidth,'LineStyle','none'); % 5
mu1 = -0.05; mu2 = -0.0475; % 6
p6 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','r','LineWidth',gwLineWidth,'LineStyle','none'); % 6

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


xlim([-0.07 0.07])
ylim([-0.07 0.07])

% xlabel('$\mu_{1}$','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% ylabel('$\mu_{2}$','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');


% set(gca,'XAxisLocation','origin');
% set(gca,'YAxisLocation','origin');

% xtickangle(-15)
% ytickangle(-15)

% axis equal

% grid on
% box on
axis tight
box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.01 0.01],'linewidth',1.0,'layer','top');

%% 2
axes('Position',[0.727857142857143 0.400943396226415 0.216428571428571 0.262264150943398]);

hold on
h1 = fimplicit(f1,[mu1_left mu1_right mu2_left mu2_right],'-g','LineWidth',gwLineWidth);
h2 = fimplicit(f2,[mu1_left mu1_right mu2_left mu2_right],'-b','LineWidth',gwLineWidth);
h3 = fimplicit(f3,[mu1_left 0 mu2_left mu2_right],'-m','LineWidth',gwLineWidth);
h4 = fimplicit(f4,[mu1_left 0 mu2_left mu2_right],'-r','LineWidth',gwLineWidth);

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
mu1 = 0.05; mu2 = 0.05; % 1
p1 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','k','LineStyle','none'); % 1
mu1 = 0.05; mu2 = 0.0475; % 2
p2 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','r','LineStyle','none'); % 2
mu1 = 0.05; mu2 = 0.0464; % 3
p3 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','g','LineStyle','none'); % 3
mu1 = -0.05; mu2 = -0.049; % 4
p4 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','b','LineStyle','none'); % 4
mu1 = -0.05; mu2 = -0.0483; % 5
p5 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','k','LineWidth',gwLineWidth,'LineStyle','none'); % 5
mu1 = -0.05; mu2 = -0.0475; % 6
p6 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','r','LineWidth',gwLineWidth,'LineStyle','none'); % 6

set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'linewidth',0.5,'layer','top');


%% 3
axes('Position',[0.142142857142857 0.410377358490566 0.212857142857142 0.270754716981132]);

hold on
h1 = fimplicit(f1,[mu1_left mu1_right mu2_left mu2_right],'-g','LineWidth',gwLineWidth);
h2 = fimplicit(f2,[mu1_left mu1_right mu2_left mu2_right],'-b','LineWidth',gwLineWidth);
h3 = fimplicit(f3,[mu1_left 0 mu2_left mu2_right],'-m','LineWidth',gwLineWidth);
h4 = fimplicit(f4,[mu1_left 0 mu2_left mu2_right],'-r','LineWidth',gwLineWidth);

ax = gca;
% xticks = xticks(ax1);
xticks = [-0.05];
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
xticklabels(find(xticks==0)) = {''};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-0.06:0.001:-0.02];
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

set(gca,'XAxisLocation','top');
set(gca,'YAxisLocation','right');

xlim([-0.0513155656085387 -0.0474363310329603]);
ylim([-0.0504910951707173 -0.0466118605951388]);

xtickangle(15)
ytickangle(35)
grid on
box on


mu1 = 0.05; mu2 = 0.05; % 1
p1 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','k','LineStyle','none'); % 1
mu1 = 0.05; mu2 = 0.0475; % 2
p2 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','r','LineStyle','none'); % 2
mu1 = 0.05; mu2 = 0.0464; % 3
p3 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','g','LineStyle','none'); % 3
mu1 = -0.05; mu2 = -0.049; % 4
p4 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','b','LineStyle','none'); % 4
mu1 = -0.05; mu2 = -0.0483; % 5
p5 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','k','LineWidth',gwLineWidth,'LineStyle','none'); % 5
mu1 = -0.05; mu2 = -0.0475; % 6
p6 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','r','LineWidth',gwLineWidth,'LineStyle','none'); % 6

set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'linewidth',0.5,'layer','top');



figure_name = [results_path 'Figure Step02 THB region.eps'];
saveas(gcf, figure_name, 'epsc');