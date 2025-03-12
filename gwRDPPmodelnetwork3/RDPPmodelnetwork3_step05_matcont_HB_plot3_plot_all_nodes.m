close all
network_name = 'ws_n=10k=4p=0.05number=69';
datanetwork_path = ['./network3 ' network_name '/'];
results_path = strrep(datanetwork_path,'network3','RDPPmodel results network3');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

figure_path = [results_path 'figure step03 hb nodes wsnetworks_n=10k=2p=0.05number=69'];
if ~exist(figure_path, 'dir')
    mkdir(figure_path);
end


% gwcolors = lines(10);
gwcolors = hsv(11);
gwcolors(3,:) = [];

Node_pos = 10;

figure
set(gcf,"Position",[300 300 500 320])
set(gca,'Position',[0.172872039794922 0.164985119047619 0.792127960205078 0.764702380952381]);
hold on


for Node_pos=1:10


e=[size(xlc,1) 2*(Node_pos-1)+1 2*(Node_pos-1)+2];

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

% gwcplanim2(x1,v1,s1,'b',[2*N+1,2*(Node_pos-1)+1]);
% gwcplanim2(x2,v2,s2,'k',[2*N+1,2*(Node_pos-1)+1]);
% gwcplanim2(x1,v1,s1,'b',[2*N+1,2*(Node_pos-1)+1]);
unode = x1(2*(Node_pos-1)+1,:); 
para = x1(2*N+1,:);
pos1_in = [s1(1).index:s1(3).index];
cont1_para = para(pos1_in); cont1_unode = unode(pos1_in);
% plot(cont1_para,cont1_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');

% gwcplanim2(x2,v2,s2,'c',[2*N+1, 2*(Node_pos-1)+1]);
unode = x2(2*(Node_pos-1)+1,:); 
para = x2(2*N+1,:);
pos1_in = [s2(1).index:s2(2).index];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
% plot(cont2_para,cont2_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

% cont_para = [cont1_para, cont2_para];
% cont_unode = [cont1_unode, cont2_unode];
% cont_para = [cont1_para, fliplr(cont2_para)];
% cont_unode = [cont1_unode, fliplr(cont2_unode)];
cont_para = [fliplr(cont1_para), cont2_para];
cont_unode = [fliplr(cont1_unode), cont2_unode];
% cont_para = [cont1_para];
% cont_unode = [cont1_unode];

% plot(cont_para,cont_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
pos_right_u_star = find(cont_para>=b0);
right_unode = [cont_unode(pos_right_u_star),u_star0,];
right_para = [cont_para(pos_right_u_star),b0,];
% plot(right_para,right_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

plot3(right_para,Node_pos*ones(size(right_para)),right_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

pos_left_u_star = find(cont_para<=b0);
left_unode = [cont_unode(pos_left_u_star),u_star0,];
left_para = [cont_para(pos_left_u_star),b0,];
% plot(left_para,left_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','k');
plot3(left_para,Node_pos*ones(size(left_para)),left_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color','k');

% plot(x_p,ymin_p,'linewidth',gwlinewidth,'Color',gwcolors(Node_pos,:))
% plot(x_p,ymax_p,'linewidth',gwlinewidth,'Color',gwcolors(Node_pos,:))


plot3(x_p,Node_pos*ones(size(x_p)),ymin_p,'linewidth',gwlinewidth,'Color',gwcolors(Node_pos,:))
plot3(x_p,Node_pos*ones(size(x_p)),ymax_p,'linewidth',gwlinewidth,'Color',gwcolors(Node_pos,:))




end

% 
% plot(b0, u_star0, 'k.')

xlim([b_star-0.008 b_star+0.008])
ylim([0 11])
zlim([0 1])
% view([75 30])
view([80 20])

ax = gca;
% xticks = xticks(ax);
xticks = [3.285,3.295];
xticklabels = arrayfun(@(x) sprintf('%.3f' , x), xticks, 'UniformOutput', false);
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
zticks = [0:0.5:1];
zticklabels = arrayfun(@(x) sprintf('%.1f' , x), zticks, 'UniformOutput', false);
zticklabels(1) = {'0'};
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


% set(gca,'Box','on');
set(gca,'Box','on','BoxStyle','full');
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.02],'linewidth',1.0,'layer','top');

figure_name = [figure_path '/all node.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = [figure_path '/all node.jpg'];
saveas(gcf, figure_name);
