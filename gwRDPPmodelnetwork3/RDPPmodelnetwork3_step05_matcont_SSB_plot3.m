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

matcont_path = [results_path 'matcont ssb wsnetworks_n=10k=4p=0.05number=69'];
if ~exist(matcont_path, 'dir')
    mkdir(matcont_path);
end
%%
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
a = 0.8; 
a = 1.0; 
c = 0.2;
d1 = 0.001; d2 = 10*d1;
Times = 0.01;

N = length(eigenvalues);
d1 = d1 / Times;
d2 = d2 / Times;

%%
if d1 < c*(1-c)*d2/(c*(1-c) - eigenvalues(2)*d2)
    disp('Network case: Turing instability occurs.')
else
    disp('Network case: Turing instability dose not occur.')
end

%%
%%
network_bk = [];
network_ak = [];
N = size(laplacian_matrix,1);
for k=1:N-1
    one_bk = d1/d2 - (eigenvalues(k)+eigenvalues(k+1))/c/(1-c)*d1 + eigenvalues(k)*eigenvalues(k+1)/c/c/(1-c)/(1-c)*d1*d2;
    network_bk = [network_bk,one_bk];
    
    s_1k = (-d2*(1-c^2)*eigenvalues(k)+c*(1-c)^2) / (-d2*eigenvalues(k)+c*(1-c));
    s_2k = - (d1*d2*eigenvalues(k)^2-d1*c*(1-c)*eigenvalues(k)) / (-d2*eigenvalues(k)+c*(1-c));
    one_ak = s_1k*one_bk + s_2k;
    network_ak = [network_ak,one_ak];

end

postion = find(network_ak>a);
k = postion(1);

%%
% k=19;
b_star = (-d2*eigenvalues(k)+c*(1-c)) / (-d2*(1-c^2)*eigenvalues(k)+c*(1-c)^2) *a ...
    + (d1*d2*eigenvalues(k)^2-d1*c*(1-c)*eigenvalues(k))/(-d2*(1-c^2)*eigenvalues(k)+c*(1-c)^2);
b_star_s=sprintf('%20.10f',b_star);
b_star_s=b_star;

format rat
a
b_star
k
format long


%%
n = k;
phin = real(eigenvectors(:,n));

%%
mu0 = 0;
b0 = b_star + mu0;
u_star0 = b0*(a+(c-1)*b0)/a; v_star0 = b0*(1-c)/c * u_star0;

mu = 0.0001;
b_first = b_star + mu;
u_star = b_first*(a+(c-1)*b_first)/a; v_star = b_first*(1-c)/c * u_star;

x0=zeros(1,2*N);
x0(1:2:end-1) = u_star + phin* 0.05;
x0(2:2:end  ) = v_star + phin* 0.05;

hls=feval(@RDPPSystemOnNetwork3);
options=odeset('RelTol',1e-8);
% options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
[t,y]=ode45(hls{2},[0:1:25e4],x0,options, a, b_first, c, d1, d2);

simulated_u1 = y(:, 1:2:end);
simulated_v1 = y(:, 2:2:end);

x0(1:2:end-1) = u_star - phin* 0.05;
x0(2:2:end  ) = v_star - phin* 0.05;

hls=feval(@RDPPSystemOnNetwork3);
options=odeset('RelTol',1e-8);
% options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
[t,y]=ode45(hls{2},[0:1:25e4],x0,options, a, b_first, c, d1, d2);

simulated_u2 = y(:, 1:2:end);
simulated_v2 = y(:, 2:2:end);

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

%%

initial_x0=zeros(1,2*N)';
initial_x0(1:2:end-1) = u_star + phin* 0;
initial_x0(2:2:end  ) = v_star + phin* 0;

% initial_x0 = y(end,:)';

gwlinewidth = 1.0;

figure
hold on

Node_pos = 10;

xlim([b_star-0.0001 b_star+0.0001])
ylim([0.045 0.051])

xlabel('b','FontSize',18,'Interpreter','latex')
ylabel('$u_{1}$','FontSize',18,'Interpreter','latex')

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'FontSize',18,'linewidth',1.0,'layer','top');

% gwcplanim(x2,v2,s2,'b',[3*N+1, 1]);

plot(b0, u_star0, 'k*')



plot(b_first, simulated_u1(end,Node_pos), 'r*')
plot(b_first, simulated_u2(end,Node_pos), 'bd')

% plot(b0, u_star0, 'ko')
% 
% plot(b_first, u_star, 'r*')
% plot(b_first, simulated_u1(end,1), 'bd')

load([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
load([matcont_path '/data_matcont2.mat'], 'x2', 'v2', 's2', 'h2', 'f2')
load([matcont_path '/data_matcont3.mat'], 'x3', 'v3', 's3', 'h3', 'f3')
load([matcont_path '/data_matcont4.mat'], 'x4', 'v4', 's4', 'h4', 'f4')

[2*(Node_pos-1)+1 phin(Node_pos)]

% gwcplanim(x1,v1,s1,'b',[2*N+1,2*(Node_pos-1)+1]);
% gwcplanim(x2,v2,s2,'c',[2*N+1,2*(Node_pos-1)+1]);
% gwcplanim(x3,v3,s3,'r',[2*N+1,2*(Node_pos-1)+1]);
% gwcplanim2(x4,v4,s4,'m',[2*N+1,2*(Node_pos-1)+1]);

plot(b_first, simulated_u1(end,Node_pos), 'r*')
plot(b_first, simulated_u2(end,Node_pos), 'bd')

%% right-left
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
%     'LineWidth',gwlinewidth,'LineStyle','--','Color','k');
pos_right_u_star = find(cont_para>=b0);
right_unode = [u_star0,cont_unode(pos_right_u_star),];
right_para = [b0,cont_para(pos_right_u_star),];
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
s = 100;
pos1_in = [s:s+26];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
plot(cont2_para,cont2_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color','k');

s = 126;
pos1_in = [s:s+32];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
plot(cont2_para,cont2_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','--','Color','k');

s1 = 126+33;
pos1_in = [s1:s1+50];
cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
plot(cont2_para,cont2_unode,'Marker','none',...
    'LineWidth',gwlinewidth,'LineStyle','-','Color','k');

% s = 1096;
% pos1_in = [s-60:s];
% cont2_para = para(pos1_in); cont2_unode = unode(pos1_in);
% plot(cont2_para,cont2_unode,'Marker','none',...
%     'LineWidth',gwlinewidth,'LineStyle','-','Color','r');