close all

network_name = 'ws_n=10k=2p=0.05number=796';
datanetwork_path = ['./network2 ' network_name '/'];
results_path = strrep(datanetwork_path,'network2','RDPPmodel results network2');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

figure_path = [results_path 'figure step03 ssb nodes wsnetworks_n=10k=2p=0.05number=796'];
if ~exist(figure_path, 'dir')
    mkdir(figure_path);
end



% gwcolors = lines(10);
gwcolors = hsv(11);
gwcolors(3,:) = [];

Node_pos = 10;

for Node_pos=1:10

figure
set(gcf,"Position",[300 300 500 320])
set(gca,'Position',[0.146 0.129047619047619 0.814 0.745]);
hold on


% Node_pos = 1;

load([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
load([matcont_path '/data_matcont2.mat'], 'x2', 'v2', 's2', 'h2', 'f2')
load([matcont_path '/data_matcont3.mat'], 'x3', 'v3', 's3', 'h3', 'f3')
load([matcont_path '/data_matcont4.mat'], 'x4', 'v4', 's4', 'h4', 'f4')



[2*(Node_pos-1)+1 phin(Node_pos)]

% gwcplanim2(x1,v1,s1,'b',[2*N+1, 2*(Node_pos-1)+1]);
% gwcplanim2(x2,v2,s2,'c',[2*N+1, 2*(Node_pos-1)+1]);
% gwcplanim2(x3,v3,s3,'r',[2*N+1, 2*(Node_pos-1)+1]);
% gwcplanim2(x4,v4,s4,'m',[2*N+1, 2*(Node_pos-1)+1]);

% plot(b_first, simulated_u1(end,Node_pos), 'r*')
% plot(b_first, simulated_u2(end,Node_pos), 'bd')

%% right-left
% gwcplanim2(x1,v1,s1,'b',[2*N+1,2*(Node_pos-1)+1]);
unode = x1(2*(Node_pos-1)+1,:); 
para = x1(2*N+1,:);
pos1_in = [s1(1).index:s1(3).index];
cont1_para = para(pos1_in); cont1_unode = unode(pos1_in);
% plot(cont1_para,cont1_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));
% gwcplanim2(x2,v2,s2,'c',[2*N+1, 2*(Node_pos-1)+1]);
unode = x2(2*(Node_pos-1)+1,:); 
para = x2(2*N+1,:);
pos1_in = [s2(1).index:s2(4).index];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
% plot(cont2_para,cont2_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','b');

% cont_para = [cont1_para, cont2_para];
% cont_unode = [cont1_unode, cont2_unode];
% cont_para = [cont1_para, fliplr(cont2_para)];
% cont_unode = [cont1_unode, fliplr(cont2_unode)];
cont_para = [fliplr(cont1_para), cont2_para];
cont_unode = [fliplr(cont1_unode), cont2_unode];

% plot(cont_para,cont_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color',gwcolors(Node_pos,:));
pos_right_u_star = find(cont_para>=b0);
right_unode = [cont_unode(pos_right_u_star),u_star0,];
right_para = [cont_para(pos_right_u_star),b0,];
plot(right_para(1:5:end),right_unode(1:5:end),'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

pos_left_u_star = find(cont_para<=b0);
left_unode = [cont_unode(pos_left_u_star),u_star0,];
left_para = [cont_para(pos_left_u_star),b0,];
plot(left_para,left_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color','k');


%% up-down

% gwcplanim(x4,v4,s4,'m',[2*N+1,1]);
unode = x4(2*(Node_pos-1)+1,:); 
para = x4(2*N+1,:);
pos1_in = [43:100];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
plot(cont2_para,cont2_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color',gwcolors(Node_pos,:));

pos1_in = [101:170];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
plot(cont2_para,cont2_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));

s1 = 1097;
pos1_in = [s1:s1+58];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
plot(cont2_para,cont2_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color',gwcolors(Node_pos,:));

s = 1096;
pos1_in = [s-60:s];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
plot(cont2_para,cont2_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));

plot(b0, u_star0, 'k.')

xlim([b_star-0.0015 b_star+0.0008])
ylim([0.055 0.13])

ax = gca;
% xticks = xticks(ax);
xticks = [1.136:0.001:1.1380];
xticklabels = arrayfun(@(x) sprintf('%.3f' , x), xticks, 'UniformOutput', false);
% xticklabels(1) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [0.06:0.03:0.12];
yticklabels = arrayfun(@(x) sprintf('%.2f', x), yticks, 'UniformOutput', false);
% yticklabels(1) = {'0'};
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

% if Node_pos ~= 10
%     set(gca, 'xticklabel', []);
%     set(gca, 'yticklabel', []);
% end

% xlabel('$b$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% ylabel(['$u_{' num2str(Node_pos) '}$'],'FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
title(['node $i=' num2str(Node_pos) '$'],'FontWeight','bold', ...
       'FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex','Color',gwcolors(Node_pos,:));

if phin(Node_pos) > 0.001
    set(gca,'XColor','r','YColor','r');
end 
if phin(Node_pos) < -0.001
    set(gca,'XColor','b','YColor','b');
end 
if abs(phin(Node_pos)) < 0.001
    set(gca,'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5]);
end 
box on
set(gca,'TickLength',[0.02 0.05],'FontSize',24,'linewidth',1.0,'layer','top');

figure_name = [figure_path '/node i=' num2str(Node_pos) '.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = [figure_path '/node i=' num2str(Node_pos) '.jpg'];
saveas(gcf, figure_name);

% close()

end