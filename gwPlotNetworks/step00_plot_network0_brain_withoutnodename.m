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
set(gcf,"Position",[60 100 700 700])
axes("Position",[0.02 0.02 0.96 0.96])
BD=bubbleGraph_BrainGou(Data,Class);
BD=BD.draw();

colorclass = hsv(41);

BD.setBubbleColor(colorclass);
bubblesize([8,9]);
axis equal

figure_name = ['./Figure step00 plot brain without nodename.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = ['./Figure step00 plot brain without nodename.jpg'];
saveas(gcf, figure_name);

