clear all
close all
clc

network_name = 'npbcws_n=10k=2p=0number=1';
datanetwork_path = ['./network1 ' network_name '/'];
results_path = strrep(datanetwork_path,'network1','FHN results network1');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

load([results_path 'THBNormalForm'],'a_star','c_star','du','dv','alpha11','alpha12','alpha21','alpha22','kappa30','kappa12','kappa21','kappa03','kappa021','kappa022','kappa023','n')


%%

mu1_left = -0.01; mu1_right = 0.01;
mu2_left = -0.01; mu2_right = 0.01;

gwLineWidth = 1.0;

figure
set(gcf,"Position",[300 200 700 530])
axes1 = axes('Position',[0.118571428571429 0.16207527759708 0.845 0.781792646931222]);
hold on

% plot([mu1_left mu1_right],[0 0],'k-.')
% plot([0 0],[mu2_left mu2_right],'k-.')
f1 = @(mu1,mu2)  -(alpha11*mu1+alpha12*mu2)/kappa30;
f2 = @(mu1,mu2)  -(alpha21*mu1+alpha22*mu2)/kappa03;
f3 = @(mu1,mu2) (kappa03*(alpha11*mu1+alpha12*mu2)-kappa12*(alpha21*mu1+alpha22*mu2))...
    /(-kappa30*kappa03+kappa21*kappa12); 
f4 = @(mu1,mu2) (kappa21*(alpha11*mu1+alpha12*mu2)-kappa30*(alpha21*mu1+alpha22*mu2))...
    /(kappa30*kappa03-kappa21*kappa12);

h1 = fimplicit(f1,[mu1_left mu1_right mu2_left mu2_right],'-g','LineWidth',1.5);
h2 = fimplicit(f2,[mu1_left mu1_right mu2_left mu2_right],'-b','LineWidth',1.5);
h3 = fimplicit(f3,[0 mu1_right mu2_left mu2_right],'-m','LineWidth',1.5);
h4 = fimplicit(f4,[0 mu1_right mu2_left mu2_right],'-r','LineWidth',1.5);

x = [-0.07, 0.07];
y = zeros(size(x));
plot(x,y,'MarkerSize',15,'Marker','none','Color','k','LineWidth',0.5,'LineStyle',':'); % 6

y = [-0.07, 0.07];
x = zeros(size(y));
plot(x,y,'MarkerSize',15,'Marker','none','Color','k','LineWidth',0.5,'LineStyle',':'); % 6

mu1 = 0.002; mu2 =  0.008;
p1 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','k','LineStyle','none'); % 1
mu1 = -0.003; mu2 =  0.002;
p2 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','r','LineStyle','none'); % 2
mu1 = -0.002; mu2 =  -0.008;
p3 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','g','LineStyle','none'); % 3
mu1 = 0.0035; mu2 =  -0.008;
p4 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','b','LineStyle','none'); % 4
mu1 = 0.0028; mu2 =  -0.005;
p5 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','k','LineWidth',1.5,'LineStyle','none'); % 5
mu1 = 0.005; mu2 =  0.002;
p6 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','r','LineWidth',1.5,'LineStyle','none'); % 6

ax = gca;
% xticks = xticks(ax1);
xticks = [-0.005:0.005:0.005];
xticklabels = arrayfun(@(x) sprintf('%.3f', x), xticks, 'UniformOutput', false);
xticklabels(find(xticks==0)) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-0.01:0.005:0.01];
yticklabels = arrayfun(@(x) sprintf('%.3f', x), yticks, 'UniformOutput', false);
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

xlim([mu1_left mu1_right])
ylim([mu2_left mu2_right])

% xlabel('$\mu_{1}$','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% ylabel('$\mu_{2}$','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.01 0.01],'linewidth',1.0,'layer','top');



%% 2
axes('Position',[0.725714285714286 0.339622641509434 0.211428571428571 0.250943396226415]);

hold on
h1 = fimplicit(f1,[mu1_left mu1_right mu2_left mu2_right],'-g','LineWidth',1.5);
h2 = fimplicit(f2,[mu1_left mu1_right mu2_left mu2_right],'-b','LineWidth',1.5);
h3 = fimplicit(f3,[0 mu1_right mu2_left mu2_right],'-m','LineWidth',1.5);
h4 = fimplicit(f4,[0 mu1_right mu2_left mu2_right],'-r','LineWidth',1.5);

x = [-0.07, 0.07];
y = zeros(size(x));
plot(x,y,'MarkerSize',15,'Marker','none','Color','k','LineWidth',0.5,'LineStyle',':'); % 6

y = [-0.07, 0.07];
x = zeros(size(y));
plot(x,y,'MarkerSize',15,'Marker','none','Color','k','LineWidth',0.5,'LineStyle',':'); % 6

mu1 = 0.002; mu2 =  0.008;
p1 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','k','LineStyle','none'); % 1
mu1 = -0.003; mu2 =  0.002;
p2 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','r','LineStyle','none'); % 2
mu1 = -0.002; mu2 =  -0.008;
p3 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','g','LineStyle','none'); % 3
mu1 = 0.0035; mu2 =  -0.008;
p4 = plot(mu1,mu2,'MarkerSize',15,'Marker','.','Color','b','LineStyle','none'); % 4
mu1 = 0.0028; mu2 =  -0.005;
p5 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','k','LineWidth',1.5,'LineStyle','none'); % 5
mu1 = 0.005; mu2 =  0.002;
p6 = plot(mu1,mu2,'MarkerSize',15,'Marker','x','Color','r','LineWidth',1.5,'LineStyle','none'); % 6

ax = gca;
% xticks = xticks(ax1);
xticks = [0.0035];
xticklabels = arrayfun(@(x) sprintf('%.4f', x), xticks, 'UniformOutput', false);
xticklabels(find(xticks==0)) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [-0.01:0.005:0.01];
yticklabels = arrayfun(@(x) sprintf('%.3f', x), yticks, 'UniformOutput', false);
yticklabels(find(yticks==0)) = {'0'};
% yticklabels(end) = {'$a$'};
set(gca, 'ytick', yticks);
set(gca, 'yticklabel', yticklabels);

ax_FontSize = 18;
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

xlim([0.00249720188808178 0.00435084329937695]);
ylim([-0.00931191855188248 -0.00677131249127641]);

% xlabel('$\mu_{1}$','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% ylabel('$\mu_{2}$','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.05 0.01],'linewidth',1.0,'layer','top');

figure_name = [results_path 'Figure Step02 THB region npbcwsnetworks_n=10k=2p=0number=1.eps'];
saveas(gcf, figure_name, 'epsc');

