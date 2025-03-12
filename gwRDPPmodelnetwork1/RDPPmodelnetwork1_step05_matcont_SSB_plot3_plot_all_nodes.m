close all
network_name = 'npbcws_n=10k=2p=0number=1';
datanetwork_path = ['./network1 ' network_name '/'];
results_path = strrep(datanetwork_path,'network1','RDPPmodel results network1');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

figure_path = [results_path 'figure step03 ssb nodes npbcwsnetworks_n=10k=2p=0number=1'];
if ~exist(figure_path, 'dir')
    mkdir(figure_path);
end



% gwcolors = lines(10);
gwcolors = hsv(11);
gwcolors(3,:) = [];

% Node_pos = 1;

load([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
load([matcont_path '/data_matcont2.mat'], 'x2', 'v2', 's2', 'h2', 'f2')
load([matcont_path '/data_matcont3.mat'], 'x3', 'v3', 's3', 'h3', 'f3')
load([matcont_path '/data_matcont4.mat'], 'x4', 'v4', 's4', 'h4', 'f4')


figure
set(gcf,"Position",[300 300 500 320])
set(gca,'Position',[0.172872039794922 0.164985119047619 0.792127960205078 0.764702380952381]);
hold on



Node_pos = 10;

for Node_pos=1:10



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
pos1_in = [s2(1).index:s2(2).index];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
% plot(cont2_para,cont2_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));

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
% % % plot(right_para,right_unode,'Marker','none',...
% % %     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

pos_left_u_star = find(cont_para<=b0);
left_unode = [cont_unode(pos_left_u_star),u_star0,];
left_para = [cont_para(pos_left_u_star),b0,];
% % % plot(left_para,left_unode,'Marker','none',...
% % %     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');


%% up-down
% gwcplanim2(x3,v3,s3,'r',[2*N+1,1]);
unode = x3(2*(Node_pos-1)+1,:); 
para = x3(2*N+1,:);
pos1_in = [s3(1).index+1:s3(2).index];
cont1_para = para(pos1_in); cont1_unode = unode(pos1_in);
% plot(cont1_para,cont1_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));
% gwcplanim(x4,v4,s4,'m',[2*N+1,1]);
unode = x4(2*(Node_pos-1)+1,:); 
para = x4(2*N+1,:);
pos1_in = [s4(1).index+1:s4(2).index];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
% plot(cont2_para,cont2_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color',gwcolors(Node_pos,:));

% cont_para = [cont1_para, cont2_para];
% cont_unode = [cont1_unode, cont2_unode];
cont_para = [cont1_para, fliplr(cont2_para)];
cont_unode = [cont1_unode, fliplr(cont2_unode)];

% % % plot(cont_para,cont_unode,'Marker','none',...
% % %     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));

% % % pos_up_u_star = find(cont_unode>=u_star0);
% % % up_unode = [cont_unode(pos_up_u_star),u_star0,];
% % % up_para = [cont_para(pos_up_u_star),b0,];
% % % % plot(up_para,up_unode,'Marker','none',...
% % % %     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));
% % % 
% % % pos_down_u_star = find(cont_unode<=u_star0);
% % % pos_down_u_star1 = intersect(find(cont_para<=b0+0.01), find(cont_para>=b0-0.01));
% % % pos_down_u_star = intersect(pos_down_u_star, pos_down_u_star1);
% % % down_unode = [cont_unode(pos_down_u_star)];
% % % down_para = [cont_para(pos_down_u_star)];
% % % plot(down_para,down_unode,'Marker','none',...
% % %     'LineWidth',gwlinewidth,'LineStyle','--','Color',gwcolors(Node_pos,:));
% % % ccc=cont_unode(pos_down_u_star);

% plot(right_para,right_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
% plot(left_para,left_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
% plot(cont_para,cont_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));

plot3(right_para,Node_pos*ones(size(right_para)),right_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
plot3(left_para,Node_pos*ones(size(left_para)),left_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
plot3(cont_para,Node_pos*ones(size(cont_para)),cont_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(Node_pos,:));

% close()

end

% plot(b0, u_star0, 'k.')

xlim([b_star-0.01 b_star+0.01])
ylim([0 11])
zlim([0.04 0.16])
% view([75 30])
view([80 20])

ax = gca;
% xticks = xticks(ax);
xticks = [1.13,1.14];
xticklabels = arrayfun(@(x) sprintf('%.2f' , x), xticks, 'UniformOutput', false);
% xticklabels(1) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [1,3,5,8,10];
yticklabels = arrayfun(@(x) sprintf('%.0f', x), yticks, 'UniformOutput', false);
% yticklabels(1) = {'0'};
% yticklabels(end) = {'$a$'};
set(gca, 'ytick', yticks);
set(gca, 'yticklabel', yticklabels);


% zticks = zticks(ax);
zticks = [0.05:0.05:0.15];
zticklabels = arrayfun(@(x) sprintf('%.2f' , x), zticks, 'UniformOutput', false);
% zticklabels(1) = {'0'};
% zticklabels(end) = {'$b$'};
set(gca, 'ztick', zticks);
set(gca, 'zticklabel', zticklabels);


ax_FontSize = 24;
ax.XAxis.FontSize = ax_FontSize;  
ax.XAxis.FontName = 'Times New Roman';  
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickLabelInterpreter = 'latex';


ax.YAxis.FontSize = ax_FontSize;  
ax.YAxis.FontName = 'Times New Roman';  
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickLabelInterpreter = 'latex';


ax.ZAxis.FontSize = ax_FontSize;  
ax.ZAxis.FontName = 'Times New Roman';  
ax.ZAxis.TickDirection = 'in';
ax.ZAxis.TickLabelInterpreter = 'latex';

% if Node_pos ~= 10
%     set(gca, 'xticklabel', []);
%     set(gca, 'yticklabel', []);
% end

% % xlabel('$b$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% % ylabel(['$u_{' num2str(Node_pos) '}$'],'FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% title(['node $i=' num2str(Node_pos) '$'],'FontWeight','bold', ...
%        'FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex','Color',gwcolors(Node_pos,:));
% 
% if phin(Node_pos) > 0.001
%     set(gca,'XColor','r','YColor','r');
% end 
% if phin(Node_pos) < -0.001
%     set(gca,'XColor','b','YColor','b');
% end 
% if abs(phin(Node_pos)) < 0.001
%     set(gca,'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5]);
% end 

set(gca,'Box','on');
% set(gca,'Box','on','BoxStyle','full');
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'linewidth',1.0,'layer','top');

figure_name = [figure_path '/all nodes.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = [figure_path '/all nodes.jpg'];
saveas(gcf, figure_name);

