clear all
close all
clc

network_name = 'ws_n=10k=2p=0.05number=796';
datanetwork_path = ['./network2 ' network_name '/'];

net_number = 796;
mat_name = [datanetwork_path 'adj_matrix_wsnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);

N = size(adjacent_matrix,2);
degree = sum(adjacent_matrix,2);

laplacian_matrix = adjacent_matrix - diag(degree);

[eigenvectors, matrix_eigenvalues] = eig(laplacian_matrix);
eigenvalues = diag(matrix_eigenvalues);
eigenvectors = fliplr(eigenvectors);
eigenvalues = flipud(eigenvalues);
eigenvalues(1) = 0;
eigenvectors(:,1) = - eigenvectors(:,1);

%%
gwcolors = hsv(11);
gwcolors(3,:) = [];
gwcolors = flipud(gwcolors);

Data = adjacent_matrix;
Class = [1:N];

nodeName = cell(1,N);
for i=1:N
    nodeName{i}=[num2str(i)];
end

figure
set(gcf,"Position",[300 300 400 300])
axes1 = axes('Position',[0.0749999828388293 -0.231666666666667 0.865000017161171 1.15333335621489]);
BD=bubbleGraph_gou(Data,Class);

BD=BD.draw();
colorclass = gwcolors;
BD.setBubbleColor(colorclass);
bubblesize([5,16]);


figure_name = ['./Figure Step00 plot wsnetworksn10k2p0dot05number796.eps'];
saveas(gcf, figure_name, 'epsc');

figure_name = ['./Figure Step00 plot wsnetworksn10k2p0dot05number796.jpg'];
saveas(gcf, figure_name);