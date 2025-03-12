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

for Node_pos=1:10

figure
set(gcf,"Position",[300 300 500 320])
set(gca,'Position',[0.146 0.129047619047619 0.814 0.745]);
hold on


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
plot(right_para(1:5:end),right_unode(1:5:end),'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

pos_left_u_star = find(cont_para<=b0);
left_unode = [cont_unode(pos_left_u_star),u_star0,];
left_para = [cont_para(pos_left_u_star),b0,];
plot(left_para,left_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color','k');

plot(x_p,ymin_p,'linewidth',gwlinewidth,'Color',gwcolors(Node_pos,:))
plot(x_p,ymax_p,'linewidth',gwlinewidth,'Color',gwcolors(Node_pos,:))

plot(b0, u_star0, 'k.')

xlim([b_star-0.025 b_star+0.025])
ylim([0 1.8])

ax = gca;
% xticks = xticks(ax);
xticks = [3.27:0.02:3.5];
xticklabels = arrayfun(@(x) sprintf('%.2f' , x), xticks, 'UniformOutput', false);
% xticklabels(1) = {'0'};
% xticklabels(end) = {'$b$'};
set(gca, 'xtick', xticks);
set(gca, 'xticklabel', xticklabels);

% yticks = yticks(ax);
yticks = [0:0.9:1.8];
yticklabels = arrayfun(@(x) sprintf('%.1f', x), yticks, 'UniformOutput', false);
yticklabels(1) = {'0'};
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

if Node_pos ~= 10
    set(gca, 'xticklabel', []);
    set(gca, 'yticklabel', []);
end

% xlabel('$b$','FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
% ylabel(['$u_{' num2str(Node_pos) '}$'],'FontWeight','bold','FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex');
title(['node $i=' num2str(Node_pos) '$'],'FontWeight','bold', ...
       'FontName','Times New Roman', 'FontSize', 28,'Interpreter','latex','Color',gwcolors(Node_pos,:));

% if phin(Node_pos) > 0.001
%     set(gca,'XColor','r','YColor','r');
% end 
% if phin(Node_pos) < -0.001
%     set(gca,'XColor','b','YColor','b');
% end 
% if abs(phin(Node_pos)) < 0.001
%     set(gca,'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5]);
% end 
set(gca,'XColor','k','YColor','k');
box on
set(gca,'TickLength',[0.02 0.05],'FontSize',24,'linewidth',gwlinewidth,'layer','top');

figure_name = [figure_path '/node i=' num2str(Node_pos) '.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = [figure_path '/node i=' num2str(Node_pos) '.jpg'];
saveas(gcf, figure_name);

% close()

end