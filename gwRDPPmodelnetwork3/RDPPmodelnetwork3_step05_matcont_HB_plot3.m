clear all
close all
clc

addpath('./Systems/');

network_name = 'ws_n=10k=4p=0.05number=69';
datanetwork_path = ['./network3 ' network_name '/'];
results_path = strrep(datanetwork_path,'network3','RDPPmodel results network3');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

net_number = 69;
mat_name = [datanetwork_path 'adj_matrix_wsnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);

matcont_path = [results_path 'matcont hb wsnetworks_n=10k=4p=0.05number=69'];
if ~exist(matcont_path, 'dir')
    mkdir(matcont_path);
end

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
a = 3.0;
c = 0.2;
d1 = 0.001; d2 = 10*d1;
Times = 0.01;

N = length(eigenvalues);
d1 = d1 / Times;
d2 = d2 / Times;

%%
k = 1;

b_star = a/(1-c^2) + c/(1+c) - (d1+d2)*eigenvalues(k)/(1-c^2);
T1 = -(d1+d2)*eigenvalues(1) - b_star * (1-c^2) + a + c*(1-c);

b_star2 = a/(1-c^2) + c/(1+c) - (d1+d2)*eigenvalues(k+1)/(1-c^2);

format rat
a
b_star
k
format long


%%
n = 2;
phin = real(eigenvectors(:,n));

%%
mu0 = 0;
b0 = b_star + mu0;
u_star0 = b0*(a+(c-1)*b0)/a; v_star0 = b0*(1-c)/c * u_star0;

% b1 = b_star2 + mu0;
% u_star1 = b1*(a+(c-1)*b1)/a; v_star1 = b1*(1-c)/c * u_star1;
% 
% mu = 0.01;
% b_first = b_star + mu;
% u_star = b_first*(a+(c-1)*b_first)/a; v_star = b_first*(1-c)/c * u_star;
% 
% x0=zeros(1,2*N);
% x0(1:2:end-1) = u_star + phin* 0.05;
% x0(2:2:end  ) = v_star + phin* 0.05;
% 
% % x0(1:2:end-1) = u_star - phin* 0.05;
% % x0(2:2:end  ) = v_star - phin* 0.05;
% 
% hls=feval(@RDwsnetworksn10k4p0dot05number69);
% options=odeset('RelTol',1e-8);
% % options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
% [t,y]=ode45(hls{2},[0:1:25e4],x0,options, a, b_first, c, d1, d2);
% 
% simulated_u1 = y(:, 1:2:end);
% simulated_v1 = y(:, 2:2:end);
% 
% figure
% subplot(2,1,1)
% hold on
% for i=1:size(laplacian_matrix,1)
%     plot(t,simulated_u1(:,i));
% end
% xlabel('time t', 'Interpreter','latex')
% ylabel('$u_{i}$', 'Interpreter','latex')
% 
% subplot(2,1,2)
% hold on
% for i=1:size(laplacian_matrix,1)
%     plot(t,simulated_v1(:,i));
% end
% xlabel('time t', 'Interpreter','latex')
% ylabel('$v_{i}$', 'Interpreter','latex')
% 
% %%
% initial_x0=zeros(1,2*N)';
% initial_x0(1:2:end-1) = u_star + phin* 0;
% initial_x0(2:2:end  ) = v_star + phin* 0;
% 
% % initial_x0 = y(end,:)';
% 
% figure
% hold on
% 
% xlim([b_star-0.03 b_star+0.03])
% ylim([0 1])
% 
% xlabel('b','FontSize',18,'Interpreter','latex')
% ylabel('$u_{1}$','FontSize',18,'Interpreter','latex')
% 
% box on
% set(gca,'XColor','k','YColor','k','TickLength',...
%     [0.02 0.05],'FontSize',18,'linewidth',1.0,'layer','top');
% 
% plot(b0, u_star0, 'k*')
% plot(b0, u_star0, 'ko')
% 
% plot(b1, u_star1, 'r*')
% plot(b1, u_star1, 'ro')
% 
% plot(b_first, u_star, 'r*')
% plot(b_first, simulated_u1(end,1), 'bd')

%%
load([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
load([matcont_path '/data_matcont2.mat'], 'x2', 'v2', 's2', 'h2', 'f2')
load([matcont_path '/data_matcont3.mat'], 'x3', 'v3', 's3', 'h3', 'f3')
load([matcont_path '/data_matcont4.mat'], 'x4', 'v4', 's4', 'h4', 'f4')
load([matcont_path '/data_matcontHB.mat'], 'xlc', 'vlc', 'slc', 'hlc', 'flc')

% gwcplanim2(x1,v1,s1,'b',[2*N+1,1]);
% gwcplanim(x2,v2,s2,'c',[2*N+1, 1]);
% gwcplanim(x3,v3,s3,'r',[2*N+1,1]);
% gwcplanim(x4,v4,s4,'m',[2*N+1,1]);

% figure 
% hold on
% xlim([b_star-0.03 b_star+0.03])
% ylim([0 1])
% 
% plotcycle(xlc,vlc,slc,[size(xlc,1) 1 2]);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% view(3)

%%

x=xlc;v=vlc;s=slc;

x_p = zeros(1,size(x,2));
ymin_p = zeros(1,size(x,2));
ymax_p = zeros(1,size(x,2));

zmin_p = zeros(1,size(x,2));
zmax_p = zeros(1,size(x,2));

gwnphase = N*2;
gwtps = (size(xlc,1)-2)/2/N;


Node_pos = 4;
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

gwlinewidth = 1.0;

figure 
hold on

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

plot(x_p,ymin_p,'r')
plot(x_p,ymax_p,'r')

plot(b0, u_star0, 'k.')

xlim([b_star-0.025 b_star+0.025])
ylim([0 2.0])



