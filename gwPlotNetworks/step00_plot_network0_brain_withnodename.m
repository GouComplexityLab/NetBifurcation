clear all
close all
clc

network_name = 'brain';
datanetwork_path = ['./network0 ' network_name '/'];
mat_name = [datanetwork_path 'adj_matrix_brainnetwork'];
load(mat_name)
adjacent_matrix = full(adj_matrix);

degree = sum(adjacent_matrix, 1);

LaplaceMatrix = adjacent_matrix - diag(degree);

Data = adjacent_matrix;

Class=[[1:41],[1:41]];

for i=1:82
    nodeName{i}=[num2str(i),'-',num2str(Class(i))];
end

nodeName = readcell([datanetwork_path 'brainregionnames.xlsx']);

for i=1:41
    nodeName{i}=[nodeName{i},', ', num2str(i)];
end

for i=42:82
    nodeName{i}=[num2str(i),', ',nodeName{i}];
end

figure
set(gcf,"Position",[60 100 700+60 700])
axes("Position",[0.15 0.17 0.7 0.7])
BD=bubbleGraph_BrainGou(Data,Class,'NodeName',nodeName);
BD=BD.draw();

colorclass = hsv(41);

BD.setBubbleColor(colorclass);
bubblesize([8,9]);

annotation('textbox',...
    [0.0721052631578946 0.333571428571429 0.412105263157895 0.0812886674049467],...
    'VerticalAlignment','middle',...
    'String',{'Left hemisphere'},...
    'Rotation',90,...
    'HorizontalAlignment','center',...
    'FontSize',28,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation('textbox',...
    [0.940394736842104 0.757142857142857 0.368815789473685 0.057717238833518],...
    'VerticalAlignment','middle',...
    'String',{'Right hemisphere'},...
    'Rotation',180*1+90,...
    'HorizontalAlignment','center',...
    'FontSize',28,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

figure_name = ['./Figure step00 plot brain with nodename.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = ['./Figure step00 plot brain with nodename.jpg'];
saveas(gcf, figure_name);
