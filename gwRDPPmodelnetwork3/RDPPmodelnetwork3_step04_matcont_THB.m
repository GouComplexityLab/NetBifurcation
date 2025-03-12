clear all
close all
clc

init()

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

matcont_path = [results_path 'matcont wsnetworks_n=10k=4p=0.05number=69'];
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

k = 2;
special_eigenvector = eigenvectors(:,k);

%%
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
[network_v network_p]=min(abs((eigenvalues - ((d1-d2)*c*(1-c))/(2*d1*d2))));
% [network_v network_p]

n = network_p;
network_b_star = - d1*d2/(c^2*(1-c)^2)*eigenvalues(n)^2 + (d1-d2)/c/(1-c)*eigenvalues(n) + 1;
network_a_star = (1-c^2)*network_b_star + c*(c-1);



format rat
network_b_star;
network_a_star;
network_J_n = d1*d2*eigenvalues(n)^2 + (d1*c*(c-1) + d2*(network_b_star*(1-c^2) - network_a_star))*eigenvalues(n) + c*(1-c)*(network_a_star + network_b_star*(c-1));
network_T_1 = -(d1+d2)*eigenvalues(1) - network_b_star * (1-c^2) + network_a_star + c*(1-c);
format short

if network_a_star > network_b_star*(1-c)
    disp('there exists the positive equilibrium in Network.')
else
    disp('there exists no positive equilibrium in Network. Error for Network!')
end

% connecting points
network_bk = 0;
for k=1:N-1
    one_bk = d1/d2 - (eigenvalues(k)+eigenvalues(k+1))/c/(1-c)*d1 + eigenvalues(k)*eigenvalues(k+1)/c/c/(1-c)/(1-c)*d1*d2;
    network_bk = [network_bk,one_bk];
end

%%
Set_M = 50;
intersect_n = [];
for k=1:N
    
    b_HT1 = [];
    for j=1:N
    %     one_b_HT = - d1*d2/c^2/(1-c)^2 * eigenvalues(j)^2 + (d1 - d2)/c/(1-c) * eigenvalues(j) + 1;
        one_b_HT1 = (-d1*d2*eigenvalues(j)^2 + (d1-d2)*c*(1-c)*eigenvalues(j) - (d1+d2)*c*(1-c)*eigenvalues(k) ...
                      + d2*(d1+d2)*eigenvalues(j)*eigenvalues(k) + c^2*(1-c)^2) / c^2/(1-c)^2;
        b_HT1 = [b_HT1,one_b_HT1];
    end
    [max_b_HT intersect_n] = max(b_HT1);
    intersect_n = [intersect_n,intersect_n];

end

%%
choosed_k0 = 1;
b_HT0 = [];
for j=1:N
%     one_b_HT0 = - d1*d2/c^2/(1-c)^2 * eigenvalues(j)^2 + (d1 - d2)/c/(1-c) * eigenvalues(j) + 1;
    one_b_HT0 = (-d1*d2*eigenvalues(j)^2 + (d1-d2)*c*(1-c)*eigenvalues(j) - (d1+d2)*c*(1-c)*eigenvalues(choosed_k0) ...
                  + d2*(d1+d2)*eigenvalues(j)*eigenvalues(choosed_k0) + c^2*(1-c)^2) / c^2/(1-c)^2;
    b_HT0 = [b_HT0,one_b_HT0];
end
[intersect_kn_b0 intersect_kn0] = max(b_HT0);
intersect_kn_a0 = (1 - c^2) * intersect_kn_b0 - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k0);

%%
choosed_k1 = 2;
b_HT1 = [];
for j=1:N
%     one_b_HT1 = - d1*d2/c^2/(1-c)^2 * eigenvalues(j)^2 + (d1 - d2)/c/(1-c) * eigenvalues(j) + 1;
    one_b_HT1 = (-d1*d2*eigenvalues(j)^2 + (d1-d2)*c*(1-c)*eigenvalues(j) - (d1+d2)*c*(1-c)*eigenvalues(choosed_k1) ...
                  + d2*(d1+d2)*eigenvalues(j)*eigenvalues(choosed_k1) + c^2*(1-c)^2) / c^2/(1-c)^2;
    b_HT1 = [b_HT1,one_b_HT1];
end
[intersect_kn_b1 intersect_kn1] = max(b_HT1);
intersect_kn_a1 = (1 - c^2) * intersect_kn_b1 - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k1);

%%
choosed_k2 = 3;
b_HT2 = [];
for j=1:N
%     one_b_HT2 = - d1*d2/c^2/(1-c)^2 * eigenvalues(j)^2 + (d1 - d2)/c/(1-c) * eigenvalues(j) + 1;
    one_b_HT2 = (-d1*d2*eigenvalues(j)^2 + (d1-d2)*c*(1-c)*eigenvalues(j) - (d1+d2)*c*(1-c)*eigenvalues(choosed_k2) ...
                  + d2*(d1+d2)*eigenvalues(j)*eigenvalues(choosed_k2) + c^2*(1-c)^2) / c^2/(1-c)^2;
    b_HT2 = [b_HT2,one_b_HT2];
end
[intersect_kn_b2 intersect_kn2] = max(b_HT2);
intersect_kn_a2 = (1 - c^2) * intersect_kn_b2 - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k2);




%%
% M=n;
intersect_kn1 = 8;
M=intersect_kn1+1;
% M=8;

if M == N
    max_b = network_bk(M)*2;
else
    max_b = network_bk(M+1);
end
b_value = linspace(0,max_b,100);
T_a_values = zeros(M, length(b_value));

H_a_values = zeros(M, length(b_value));

points_number = 100;
segment_b_values = zeros(M, points_number);
segment_a_values = zeros(M, points_number);

s_1k = [];
s_2k = [];
connecting_b = 0;
connecting_a = 0;

for k=1:M-1
%     i
    one_s_1k = (-d2*(1-c^2)*eigenvalues(k)+c*(1-c)^2) / (-d2*eigenvalues(k)+c*(1-c));
    one_s_2k = - (d1*d2*eigenvalues(k)^2-d1*c*(1-c)*eigenvalues(k)) / (-d2*eigenvalues(k)+c*(1-c));
    
    T_a_values(k,:) = one_s_1k * b_value + one_s_2k;
    H_a_values(k,:) = (1 - c^2) * b_value - c*(1-c) + (d1 + d2)*eigenvalues(k);
    
    connecting_b = [connecting_b,network_bk(k+1)];
    one_connecting_a = one_s_1k * network_bk(k+1) + one_s_2k;
    connecting_a = [connecting_a,one_connecting_a];
    
    segment_b_values(k,:) = linspace(network_bk(k),network_bk(k+1),points_number);
    segment_a_values(k,:) = one_s_1k * segment_b_values(k,:) + one_s_2k;
    
    s_1k = [s_1k,one_s_1k];
    s_2k = [s_2k,one_s_2k];
end

intersect_1n_b = - d1*d2/c^2/(1-c)^2 * eigenvalues(n)^2 + (d1 - d2)/c/(1-c) * eigenvalues(n) + 1;
intersect_1n_a = (1 - c^2) * intersect_1n_b - c*(1-c) + (d1 + d2)*eigenvalues(1);

% linspace(network_bk(k),network_bk(k+1),points_number);
segment_H1_a_values = (1 - c^2) * segment_b_values(n,:) - c*(1-c) + (d1 + d2)*eigenvalues(1);
segment_b_values_to_end_1n = linspace(intersect_1n_b,max_b,points_number);
segment_H1_a_values_to_end = (1 - c^2) * segment_b_values_to_end_1n - c*(1-c) + (d1 + d2)*eigenvalues(1);

%%
segment_Hk_a_values1 = (1 - c^2) * segment_b_values(intersect_kn1,:) - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k1);
segment_b_values_to_end_kn1 = linspace(intersect_kn_b1,max_b,points_number);
segment_Hk_a_values1_to_end = (1 - c^2) * segment_b_values_to_end_kn1 - c*(1-c) + (d1 + d2)*eigenvalues(choosed_k1);


%%

gwlinewidth = 0.1;
figure
hold on
dq=lines(M); 
for i=1:M
    plot(segment_b_values(i,:),segment_a_values(i,:),'LineStyle','-','LineWidth',gwlinewidth,'color',dq(i,:))
end
plot(b_value,T_a_values(1,:),'LineStyle','-','color','r','LineWidth',gwlinewidth)

plot(segment_b_values_to_end_1n,segment_H1_a_values_to_end,'LineStyle','-.','color','b','LineWidth',gwlinewidth)
plot(segment_b_values_to_end_kn1,segment_Hk_a_values1_to_end,'LineStyle','-.','color','r','LineWidth',gwlinewidth)

s=10;scatter(connecting_b,connecting_a,s,'MarkerFaceColor','g','MarkerEdgeColor','g','Marker','pentagram')
s=10;scatter(intersect_1n_b,intersect_1n_a,s,'MarkerFaceColor','r','MarkerEdgeColor','r','Marker','pentagram')
s=10;scatter(intersect_kn_b1,intersect_kn_a1,s,'MarkerFaceColor','r','MarkerEdgeColor','r','Marker','pentagram')


% plot(1.16,0.5,'MarkerSize',3,'Marker','.','Color','k')
% plot(4.3,2.5,'MarkerSize',3,'Marker','.','Color','k')

plot([0:0.5:25],21.0*ones(size([0:0.5:25])),'LineStyle','-.','Color','k','LineWidth',gwlinewidth)
plot([0:0.5:25],1.65*ones(size([0:0.5:25])),'LineStyle','-.','Color','k','LineWidth',gwlinewidth)


% set(gca,'XTick',[0:5:25]);
% set(gca,'XTicklabel',[0:5:25], 'FontSize', 24,'FontWeight','bold');
% set(gca,'YTick',[0:5:25]);
% set(gca,'YTicklabel',[0:5:25], 'FontSize', 24,'FontWeight','bold');
xlim([0 25])
ylim([0 25])

ax = gca;
% xticks = xticks(ax1);
% xticks = [1.595:0.001:1.61];
% xticklabels = arrayfun(@(x) sprintf('%.3f', x), xticks, 'UniformOutput', false);
% % xticklabels(1) = {'0'};
% set(gca, 'xtick', xticks);
% set(gca, 'xticklabel', xticklabels);
% 
% % yticks = yticks(ax);
% yticks = [0.7:0.05:1.1];
% yticklabels = arrayfun(@(x) sprintf('%.2f', x), yticks, 'UniformOutput', false);
% % yticklabels(1) = {'0'};
% set(gca, 'ytick', yticks);
% set(gca, 'yticklabel', yticklabels);

ax_FontSize = 24;
ax.XAxis.FontSize = ax_FontSize;  
ax.XAxis.FontName = 'Times New Roman';  
ax.YAxis.FontSize = ax_FontSize;  
ax.YAxis.FontName = 'Times New Roman';  
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'in';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

box on

set(gca,'LineWidth',gwlinewidth);
xlabel('b','FontWeight','bold','FontName','Times New Roman', 'FontSize', 36,'Interpreter','latex');
ylabel('a','FontWeight','bold','FontName','Times New Roman', 'FontSize', 36,'Interpreter','latex');


%%

format rat
network_b_star
network_a_star
format long

b_TH = network_b_star;
a_TH = network_a_star;

% mu1 = -0.05; 
% mu2 = -0.0483; % 5

mu1 = -0.05; 
mu2 =  -0.0475; % 7

b = b_TH + mu1;
a = a_TH + mu2; 

% a = a_TH + 0.004; 
% b = b_TH + 0.00004;

u_star = b*(a+(c-1)*b)/a; v_star = b*(1-c)/c * u_star;


% X0=rand(1,10);
X0=zeros(1,2*size(laplacian_matrix,1));
Initial_U = u_star*ones(1,10) - 0.01*ones(1,10) + 0.05 * eigenvectors(:,n)';
Initial_V = v_star*ones(1,10) - 0.01*ones(1,10) + 0.05 * eigenvectors(:,n)';
% load data_Initial_UV
X0(1:2:end-1) = Initial_U;
X0(2:2:end  ) = Initial_V;

% X0(1:2:end-1) = u_star + eigenvectors(:,n)* 0.005;
% X0(2:2:end  ) = v_star + eigenvectors(:,n)* 0.005;

% X0(1:2:end-1) = u_star + rand(1,5)* 0.5;
% X0(2:2:end  ) = v_star + rand(1,5)* 0.5;

P0 =[a;b;c;d1;d2]; % par_a,par_b,par_c,par_d1,par_d2
hls=feval(@RDPPSystemOnNetwork3);
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(2*N,1));
% options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
[t,y]=ode45(hls{2},[0:1:10e4],X0',options,P0(1),P0(2),P0(3),P0(4),P0(5));
% size(t)
% size(y)

U_series = y(:,[1:2:end])'; V_series = y(:,[2:2:end])';

simulated_u1 = U_series';
simulated_v1 = V_series';

figure
subplot(2,1,1)
hold on
for i=1:size(laplacian_matrix,1)
    plot(t,simulated_u1(:,i));
end
xlabel('time t', 'FontSize',18, 'Interpreter','latex')
ylabel('$u_{i}$', 'FontSize',18, 'Interpreter','latex')

subplot(2,1,2)
hold on
for i=1:size(laplacian_matrix,1)
    plot(t,simulated_v1(:,i));
end
xlabel('time t', 'FontSize',18, 'Interpreter','latex')
ylabel('$v_{i}$', 'FontSize',18, 'Interpreter','latex')

%%

initial_x0 = y(end,:)';


%%

figure
hold on
xlim([b_TH-0.06 b_TH+0.06])
ylim([0.04 0.14])

% xlim([b_TH-0.06 b_TH+0.06])
% ylim([0.05 0.6])

b_TH = network_b_star;
a_TH = network_a_star;
c_TH = 0.2;
u_TH = b_TH*(a_TH+(c_TH-1)*b_TH)/a_TH; v_TH = b_TH*(1-c_TH)/c_TH * u_star;
plot(b_TH,u_TH,'kd')

xlabel('b')
ylabel('u1')


p = P0;
ap1=[2];

plot(p(ap1),initial_x0(1),'k.')


[x0,v0]=init_EP_EP(@RDPPSystemOnNetwork3,initial_x0,p,ap1);
opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'MaxNumPoints',1000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'InitStepsize',0.0001);
% % opt=contset(opt,'TSearchOrder',0);
[x1,v1,s1,h1,f1]=cont(@equilibrium,x0,[],opt);
gwcplanim2(x1,v1,s1,'c',[2*N+1,1]);

opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'MaxNumPoints',1000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'InitStepsize',0.0001);
opt=contset(opt,'Backward',1);
[x2,v2,s2,h2,f2]=cont(@equilibrium,x0,[],opt);
gwcplanim2(x2,v2,s2,'c',[2*N+1, 1]);

%% BP-EP 不需要
% BP_pos = 9;
% initial_x1=x2(1:2*N,s2(BP_pos).index);
% p(ap1)=x2(2*N+1,s2(BP_pos).index);
% 
% plot(p(ap1),initial_x1(1),'k*')
% 
% [x0,v0]=init_BP_EP(@RDPPSystemOnNetwork3,initial_x1,p,s2(BP_pos),0.01);
BP_pos = 3;
initial_x1=x1(1:2*N,s1(BP_pos).index);
p(ap1)=x1(2*N+1,s1(BP_pos).index);
[x0,v0]=init_BP_EP(@RDPPSystemOnNetwork3,initial_x1,p,s1(BP_pos),0.01);
opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.001);
opt=contset(opt,'MaxNumPoints',1000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'InitStepsize',0.0001);
[x3,v3,s3,h3,f3]=cont(@equilibrium,x0,v0,opt);
gwcplanim2(x3,v3,s3,'m',[2*N+1, 1]);
opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.001);
opt=contset(opt,'MaxNumPoints',1000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'InitStepsize',0.0001);
opt=contset(opt,'Backward',1);
[x4,v4,s4,h4,f4]=cont(@equilibrium,x0,v0,opt);
gwcplanim2(x4,v4,s4,'r',[2*N+1, 1]);

% % BP-EP 无效
% BP_pos = 3;
% initial_x1=x4(1:2*N,s4(BP_pos).index);
% p(ap1)=x4(2*N+1,s4(BP_pos).index);
% plot(p(ap1),initial_x1(1),'k*')
% [x0,v0]=init_BP_EP(@RDPPSystemOnNetwork3,initial_x1,p,s4(BP_pos),0.01);
% opt=contset;
% opt=contset(opt,'VarTolerance',1e-3);
% opt=contset(opt,'FunTolerance',1e-3);
% opt=contset(opt,'MaxStepsize',0.01);
% opt=contset(opt,'MaxNumPoints',300);
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'TSearchOrder',1);
% opt=contset(opt,'InitStepsize',0.0001);
% [x3,v3,s3,h3,f3]=cont(@equilibrium,x0,v0,opt);
% gwcplanim2(x3,v3,s3,'k',[2*N+1, 1]);
% opt=contset(opt,'TSearchOrder',0);
% opt=contset(opt,'Backward',1);
% [x4,v4,s4,h4,f4]=cont(@equilibrium,x0,v0,opt);
% gwcplanim2(x4,v4,s4,'k',[2*N+1, 1]);

% % %% BP-LP
% % initial_x1=x2(1:2*N,s2(2).index);
% % p(ap1)=x2(2*N+1,s2(2).index);
% % ap1=[2 1]
% % [x0,v0]=init_BP_LP(@RDPPSystemOnNetwork3,initial_x1,p,ap1);
% % opt=contset(opt,'InitStepsize',0.0001);
% % [x3,v3,s3,h3,f3]=cont(@limitpoint,x0,v0,opt);
% % gwcplanim(x3,v3,s3,[2*N+1, 1]);

% BP-LP
% BP_pos = 9;
% initial_x1=x2(1:2*N,s2(BP_pos).index);
% p(ap1)=x2(2*N+1,s2(BP_pos).index);
% [x0,v0]=init_BP_LP(@RDPPSystemOnNetwork3,initial_x1,p,[2 1]);

BP_pos = 3;
initial_x1=x1(1:2*N,s1(BP_pos).index);
p(ap1)=x1(2*N+1,s1(BP_pos).index);
[x0,v0]=init_BP_LP(@RDPPSystemOnNetwork3,initial_x1,p,[2 1]);

opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'MaxNumPoints',1000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'InitStepsize',0.0001);
[x5,v5,s5,h5,f5]=cont(@limitpoint,x0,v0,opt);
gwcplanim2(x5,v5,s5,'k',[2*N+1, 1]);

opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'InitStepsize',0.001);
opt=contset(opt,'TSearchOrder',0);
opt=contset(opt,'Backward',1);
[x6,v6,s6,h6,f6]=cont(@limitpoint,x0,v0,opt);
gwcplanim2(x6,v6,s6,'k',[2*N+1, 1]);



%% hopf
% H_p = 2;
% initial_xH=x3(1:2*N,s3(H_p).index);
% p(ap1)=x3(2*N+1,s3(2).index);
% ap1=[2 1];
% [x0,v0]=init_H_H(@RDPPSystemOnNetwork3,initial_xH,p,ap1);
% opt=contset;
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'MaxStepsize',0.01);
% opt=contset(opt,'backward',0);
% [xh1,vh1,sh1,hh1,fh1]=cont(@hopf,x0,v0,opt);
% gwcplanim(xh1,vh1,sh1,'r',[2*N+1 1]);
% 
% % opt=contset;
% % opt=contset(opt,'Singularities',1);
% % opt=contset(opt,'MaxStepsize',0.01);
% % % opt=contset(opt,'TSearchOrder',1);
% % opt=contset(opt,'backward',1);
% % [xh2,vh2,sh2,hh2,fh2]=cont(@hopf,x0,v0,opt);
% % gwcplanim(xh2,vh2,sh2,'r',[2*N+1 1]);
% % 

%% ZH-H
ZH_pos = 2;
initial_xZH=x5(1:2*N,s5(ZH_pos).index);
% p(ap1)=x5(2*N+1,s5(ZH_pos).index);
% plot(p(ap1), initial_xZH(1), 'ko')
% plot(p(ap1), initial_xZH(1), 'r*')
% ap1=[2 1];

ap1=[2 1];
p(ap1)=x5([2*N+1,2*N+2],s5(ZH_pos).index);
plot(p(ap1(1)), initial_xZH(1), 'ko')
plot(p(ap1(1)), initial_xZH(1), 'r*')
[x0,v0]=init_ZH_H(@RDPPSystemOnNetwork3,initial_xZH,p,ap1);
opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.0001);
opt=contset(opt,'MaxNumPoints',2000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'InitStepsize',0.0001);
[x7,v7,s7,h7,f7]=cont(@hopf,x0,v0,opt);
% gwcplanim2(x7,v7,s7,'k',[2*N+1, 1]);
gwcplanimnolabel(x7,v7,s7,'b','-',[2*N+1, 1]);

opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.0001);
opt=contset(opt,'MaxNumPoints',2000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'InitStepsize',0.0001);
opt=contset(opt,'backward',1);
[x8,v8,s8,h8,f8]=cont(@hopf,x0,v0,opt);
% gwcplanim2(x8,v8,s8,'k',[2*N+1, 1]);
gwcplanimnolabel(x8,v8,s8,'b','-',[2*N+1, 1]);

%% ZH-LP

ZH_pos = 2;
initial_xZH=x5(1:2*N,s5(ZH_pos).index);

ap1=[2 1];
p(ap1)=x5([2*N+1,2*N+2],s5(ZH_pos).index);
plot(p(ap1(1)), initial_xZH(1), 'go')
plot(p(ap1(1)), initial_xZH(1), 'g*')

% [x0,v0]=init_ZH_LP(@RDPPSystemOnNetwork3,initial_xZH,p,ap1);
[x0,v0]=init_LP_LP(@RDPPSystemOnNetwork3,initial_xZH,p,ap1,[1 2]);
opt=contset;
opt=contset(opt,'VarTolerance',1e-1);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.001);
opt=contset(opt,'MaxNumPoints',1000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'InitStepsize',0.001);
[x9,v9,s9,h9,f9]=cont(@limitpoint,x0,v0,opt);
gwcplanim2(x9,v9,s9,'r',[2*N+1, 1]);

opt=contset;
opt=contset(opt,'VarTolerance',1e-1);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'MaxNumPoints',300);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'InitStepsize',0.001);
opt=contset(opt,'backward',1);
[x10,v10,s10,h10,f10]=cont(@limitpoint,x0,v0,opt);
gwcplanim2(x10,v10,s10,'m',[2*N+1, 1]);

% % % % % % % % % [x0,v0]=init_ZH_LP(@RDPPSystemOnNetwork3,initial_xZH,p,ap1);
% % % % % % % % % opt=contset;
% % % % % % % % % opt=contset(opt,'Singularities',1);
% % % % % % % % % opt=contset(opt,'MaxStepsize',0.01);
% % % % % % % % % opt=contset(opt,'backward',0);
% % % % % % % % % [xzh1,vzh1,szh1,hzh1,fzh1]=cont(@limitpoint,x0,v0,opt);
% % % % % % % % % gwcplanim(xzh1,vzh1,szh1,[2*N+1 1]);
% % % 

save([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
save([matcont_path '/data_matcont2.mat'], 'x2', 'v2', 's2', 'h2', 'f2')
save([matcont_path '/data_matcont3.mat'], 'x3', 'v3', 's3', 'h3', 'f3')
save([matcont_path '/data_matcont4.mat'], 'x4', 'v4', 's4', 'h4', 'f4')
save([matcont_path '/data_matcont5.mat'], 'x5', 'v5', 's5', 'h5', 'f5')

save([matcont_path '/data_matcont6.mat'], 'x6', 'v6', 's6', 'h6', 'f6')
save([matcont_path '/data_matcont7.mat'], 'x7', 'v7', 's7', 'h7', 'f7')
save([matcont_path '/data_matcont8.mat'], 'x8', 'v8', 's8', 'h8', 'f8')
save([matcont_path '/data_matcont9.mat'], 'x9', 'v9', 's9', 'h9', 'f9')
save([matcont_path '/data_matcont10.mat'], 'x10', 'v10', 's10', 'h10', 'f10')

