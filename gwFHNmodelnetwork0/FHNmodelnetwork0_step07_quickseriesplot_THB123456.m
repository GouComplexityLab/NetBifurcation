clear all
close all
clc

network_name = 'brain';
datanetwork_path = ['./network0 ' network_name '/'];
results_path = strrep(datanetwork_path,'network0','FHN results network0');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

gwlinewidth = 0.5;

% gwcolors = lines(10);
gwcolors = hsv(41);
gwcolors = [gwcolors;hsv(41)];

L = 82;

domain_n = 1;
pianU = 0.01;
pianV = 0.01;

domain_n = 2;
pianU = 0.02;
pianV = 0.02;

domain_n = 3;
pianU = 0.02;
pianV = 0.02;

domain_n = 4;
pianU = 0.02;
pianV = 0.02;

domain_n = 5;
pianU = 0.02;
pianV = 0.02;

domain_n = 6;
pianU = 0.02;
pianV = 0.02;

figure_path = [results_path 'figure step06 THB series brain network' '/region ' num2str(domain_n)];
if ~exist(figure_path, 'dir')
    mkdir(figure_path);
end

%%

for init_n = [1,2,3,4]

    mat_pathname = [results_path 'data_THBFHN_region_i' num2str(domain_n) '_Init' num2str(init_n) '.mat'];
    load(mat_pathname,'region_i','init_n','x','mu1','mu2','time','U_series','V_series','Initial_U','Initial_V','Final_U','Final_V')

    
    cpos = [1:1:length(time)];

    ctime = time(cpos)';
    cU_series = U_series(:,cpos);
    cV_series = V_series(:,cpos);
    
    [X, Y] = meshgrid(ctime,x);
    
    figure('Name',['U line init n = ' num2str(init_n)])
    set(gcf,"Position",[300 200+init_n*100 700 200])
    axes1 = axes('Position',[0.15 0.1225 0.821428571428571 0.76]);
    hold on
    for i=1:L
        plot(ctime,cU_series(i,:),'Marker','none',...
        'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(i,:));
    end

    ax1 = gca;
    xticks1 = xticks(ax1);
    % xticks1 = [1.595:0.005:1.61];
    xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks1, 'UniformOutput', false);
    xticklabels(find(xticks1==0)) = {'0'};
    set(gca, 'xtick', xticks1);
    % set(gca, 'xticklabel', xticklabels);
    set(gca, 'xticklabel', []);
    
    yticks1 = yticks(ax1);
    step = fix(length(yticks1)/3);
    yticks1 = yticks1(1:step:end);
    if domain_n==1
        if init_n==1 | init_n ==2 
            yticks1 = [-0.01:0.01:0.01];
        end
    end
    yticklabels = arrayfun(@(x) sprintf('%.2f', x), yticks1, 'UniformOutput', false);
    [v p] = min(abs(yticks1));
    yticklabels(p) = {'0'};
    set(gca, 'ytick', yticks1);
    set(gca, 'yticklabel', yticklabels);
    
    ax_FontSize = 24;
    ax1.XAxis.FontSize = ax_FontSize;
    ax1.XAxis.FontName = 'Times New Roman'; 
    ax1.YAxis.FontSize = ax_FontSize;  
    ax1.YAxis.FontName = 'Times New Roman';  
    ax1.XAxis.TickDirection = 'in';
    ax1.YAxis.TickDirection = 'in';
    ax1.XAxis.TickLabelInterpreter = 'latex';
    ax1.YAxis.TickLabelInterpreter = 'latex';

    ylim([min(cU_series(:))-pianU max(cU_series(:))+pianU])

    annotation('textbox',...
    [0.73 0.1725 0.174285714285714 0.1725],...
    'Color','k',...
    'String',{['$T=' num2str(ctime(end)) '$']},...
    'Interpreter','latex',...
    'FontSize',22,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');

    box on
    set(gca,'XColor','k','YColor','k','TickLength',...
        [0.01 0.005],'linewidth',0.5,'layer','top');


    figure_name = [figure_path '/U line init n = ' num2str(init_n) '.eps'];
    saveas(gcf, figure_name, 'epsc');
    
    figure_name = [figure_path '/U line init n = ' num2str(init_n) '.jpg'];
    saveas(gcf, figure_name);
    %%
    figure('Name',['V line init n = ' num2str(init_n)])
    set(gcf,"Position",[500 200+init_n*100 700 200])
    axes1 = axes('Position',[0.15 0.1225 0.821428571428571 0.76]);
    hold on
    for i=1:L
        plot(ctime,cV_series(i,:),'Marker','none',...
        'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(i,:));
    end

    ax2 = gca;
    xticks2 = xticks(ax2);
    % xticks2 = [1.595:0.005:1.61];
    xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks2, 'UniformOutput', false);
    xticklabels(find(xticks2==0)) = {'0'};
    set(gca, 'xtick', xticks2);
    % set(gca, 'xticklabel', xticklabels);
    set(gca, 'xticklabel', []);
    
    yticks2 = yticks(ax2);
    step = fix(length(yticks2)/3);
    yticks2 = yticks2(1:step:end);
    
    if domain_n==1
        if init_n==1 | init_n ==2 
            yticks2 = [-0.01:0.01:0.01];
        end
    end
    yticklabels = arrayfun(@(x) sprintf('%.2f', x), yticks2, 'UniformOutput', false);
    [v p] = min(abs(yticks2));
    yticklabels(p) = {'0'};
    set(gca, 'ytick', yticks2);
    set(gca, 'yticklabel', yticklabels);
    
    ax2_FontSize = 24;
    ax2.XAxis.FontSize = ax2_FontSize;
    ax2.XAxis.FontName = 'Times New Roman'; 
    ax2.YAxis.FontSize = ax2_FontSize;  
    ax2.YAxis.FontName = 'Times New Roman';  
    ax2.XAxis.TickDirection = 'in';
    ax2.YAxis.TickDirection = 'in';
    ax2.XAxis.TickLabelInterpreter = 'latex';
    ax2.YAxis.TickLabelInterpreter = 'latex';

    ylim([min(cV_series(:))-pianV max(cV_series(:))+pianV])

    annotation('textbox',...
    [0.73 0.1725 0.174285714285714 0.1725],...
    'Color','k',...
    'String',{['$T=' num2str(ctime(end)) '$']},...
    'Interpreter','latex',...
    'FontSize',22,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');


    box on
    set(gca,'XColor','k','YColor','k','TickLength',...
        [0.01 0.005],'linewidth',0.5,'layer','top');


    figure_name = [figure_path '/V line init n = ' num2str(init_n) '.eps'];
    saveas(gcf, figure_name, 'epsc');
    
    figure_name = [figure_path '/V line init n = ' num2str(init_n) '.jpg'];
    saveas(gcf, figure_name);
    %%
    cpos = [length(time)-1000:1:length(time)];
    ctime = time(cpos)';
    cU_series = U_series(:,cpos);
    cV_series = V_series(:,cpos);
    figure('Name',['phase init n = ' num2str(init_n)])
    set(gcf,"Position",[800 200+init_n*30 300 300])
    axes1 = axes('Position',[0.223333333333333 0.235 0.65 0.65]);
    hold on
    for i=1:L
        plot(cU_series(i,:),cV_series(i,:),'MarkerSize',18,'Marker','.',...
        'LineWidth',gwlinewidth,'LineStyle','-','Color',gwcolors(i,:),'linewidth',1.5);
    end

    xlim_left = min(cU_series(:))-pianU;
    xlim_right = max(cU_series(:))+pianU;
    ylim_left = min(cV_series(:))-pianV;
    ylim_right = max(cV_series(:))+pianV;
    xlim([xlim_left, xlim_right])
    ylim([ylim_left, ylim_right])

    ax3 = gca;
%     xticks3 = xticks(ax3);
    xticks3 = [xlim_left, 0, xlim_right];
    xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks3, 'UniformOutput', false);
    xticklabels(find(xticks3==0)) = {'0'};
    set(gca, 'xtick', xticks3);
    set(gca, 'xticklabel', xticklabels);
%     set(gca, 'xticklabel', []);
    
%     yticks3 = yticks(ax3);
%     step = fix(length(yticks3)/3);
%     yticks3 = yticks3(1:step:end);
    yticks3 = [ylim_left, 0, ylim_right];
    yticklabels = arrayfun(@(x) sprintf('%.2f', x), yticks3, 'UniformOutput', false);
    yticklabels(find(yticks3==0)) = {'0'};
    set(gca, 'ytick', yticks3);
    set(gca, 'yticklabel', yticklabels);
    
    ax3_FontSize = 18;
    ax3.XAxis.FontSize = ax3_FontSize;
    ax3.XAxis.FontName = 'Times New Roman'; 
    ax3.YAxis.FontSize = ax3_FontSize;  
    ax3.YAxis.FontName = 'Times New Roman';  
    ax3.XAxis.TickDirection = 'in';
    ax3.YAxis.TickDirection = 'in';
    ax3.XAxis.TickLabelInterpreter = 'latex';
    ax3.YAxis.TickLabelInterpreter = 'latex';
    


%     annotation('textbox',...
%     [0.767142857142857 0.1725 0.174285714285714 0.1725],...
%     'Color','k',...
%     'String',{['$T=' num2str(ctime(end)) '$']},...
%     'Interpreter','latex',...
%     'FontSize',22,...
%     'FontName','Times New Roman',...
%     'FitBoxToText','off',...
%     'EdgeColor','none',...
%     'BackgroundColor','none');

    box on
    set(gca,'XColor','k','YColor','k','TickLength',...
        [0.05 0.05],'linewidth',0.5,'layer','top');

    figure_name = [figure_path '/phase line init n = ' num2str(init_n) '.eps'];
    saveas(gcf, figure_name, 'epsc');
    
    figure_name = [figure_path '/phase line init n = ' num2str(init_n) '.jpg'];
    saveas(gcf, figure_name);

    clear time U_series V_series

end

