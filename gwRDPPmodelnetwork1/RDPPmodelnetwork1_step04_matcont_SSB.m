clear all
close all
clc

network_name = 'npbcws_n=10k=2p=0number=1';
datanetwork_path = ['./network1 ' network_name '/'];
results_path = strrep(datanetwork_path,'network1','RDPPmodel results network1');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

net_number = 1;
mat_name = [datanetwork_path 'adj_matrix_nwnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);
%%
matcont_path = [results_path 'matcont ssb npbcwsnetworks_n=10k=2p=0number=1'];
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
n = 4;
phin = real(eigenvectors(:,n));

%%
mu0 = 0;
b0 = b_star + mu0;
u_star0 = b0*(a+(c-1)*b0)/a; v_star0 = b0*(1-c)/c * u_star0;

mu = 0.001;
b_first = b_star + mu;
u_star = b_first*(a+(c-1)*b_first)/a; v_star = b_first*(1-c)/c * u_star;

% x0=zeros(1,2*N);
% x0(1:2:end-1) = u_star + phin* 0.05;
% x0(2:2:end  ) = v_star + phin* 0.05;

x0=zeros(1,2*N);
x0(1:2:end-1) = u_star + randn(size(phin))* 0.05;
x0(2:2:end  ) = v_star + randn(size(phin))* 0.05;

% x0(1:2:end-1) = u_star - phin* 0.05;
% x0(2:2:end  ) = v_star - phin* 0.05;

hls=feval(@RDPPSystemOnNetwork1);
options=odeset('RelTol',1e-8);
% options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
[t,y]=ode45(hls{2},[0:0.1:5e4],x0,options, a, b_first, c, d1, d2);

simulated_u1 = y(:, 1:2:end);
simulated_v1 = y(:, 2:2:end);

figure
subplot(2,1,1)
hold on
for i=1:size(laplacian_matrix,1)
    plot(t,simulated_u1(:,i));
end
xlabel('time t', 'Interpreter','latex')
ylabel('$u_{i}$', 'Interpreter','latex')

subplot(2,1,2)
hold on
for i=1:size(laplacian_matrix,1)
    plot(t,simulated_v1(:,i));
end
xlabel('time t', 'Interpreter','latex')
ylabel('$v_{i}$', 'Interpreter','latex')

%%

initial_x0=zeros(1,2*N)';
initial_x0(1:2:end-1) = u_star + phin* 0;
initial_x0(2:2:end  ) = v_star + phin* 0;

% initial_x0 = y(end,:)';

figure
hold on

xlim([b_star-0.01 b_star+0.01])
ylim([0.02 0.12])

xlabel('b','FontSize',18,'Interpreter','latex')
ylabel('$u_{1}$','FontSize',18,'Interpreter','latex')

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'FontSize',18,'linewidth',1.0,'layer','top');

% gwcplanim(x2,v2,s2,'b',[3*N+1, 1]);

plot(b0, u_star0, 'k*')
plot(b0, u_star0, 'ko')

plot(b_first, u_star, 'r*')
plot(b_first, simulated_u1(end,1), 'bd')


MaxStepsize_ = 1e-3;

global cds
p=[a, b_first, c, d1, d2]; ap=[2];
[x0,v0]=init_EP_EP(@RDPPSystemOnNetwork1,initial_x0,p,ap);
opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',MaxStepsize_);
opt=contset(opt,'MaxNumPoints',3000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
[x1,v1,s1,h1,f1]=cont(@equilibrium,x0,[],opt);
gwcplanim2(x1,v1,s1,'b',[2*N+1,1]);


opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',MaxStepsize_);
opt=contset(opt,'MaxNumPoints',500);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'Backward',1);
[x2,v2,s2,h2,f2]=cont(@equilibrium,x0,[],opt);
gwcplanim(x2,v2,s2,'c',[2*N+1, 1]);

%%
%% BP-EP

% initial_x1=x2(1:2*N,s2(BP_pos).index);
% p(ap)=x2(2*N+1,s2(BP_pos).index);
% % plot(p(ap),initial_x1(1),'MarkerSize',8,'Marker','*','LineStyle','none','Color','m');
% [x0,v0]=init_BP_EP(@RDPPSystemOnNetwork1,initial_x1,p,s2(BP_pos),0.01);
BP_pos = 2;
initial_x1=x1(1:2*N,s1(BP_pos).index);
p(ap)=x1(2*N+1,s1(BP_pos).index);
[x0,v0]=init_BP_EP(@RDPPSystemOnNetwork1,initial_x1,p,s1(BP_pos),0.01);
opt=contset;
opt=contset(opt,'VarTolerance',1e-6);
opt=contset(opt,'FunTolerance',1e-6);
opt=contset(opt,'MaxStepsize',MaxStepsize_);
opt=contset(opt,'MaxNumPoints',3000);
opt=contset(opt,'InitStepsize',0.000001);
opt=contset(opt,'TSearchOrder',1);
% opt=contset(opt,'TSearchOrder',0);
[x3,v3,s3,h3,f3]=cont(@equilibrium,x0,v0,opt);
gwcplanim(x3,v3,s3,'r',[2*N+1,1]);
opt=contset;
opt=contset(opt,'VarTolerance',1e-6);
opt=contset(opt,'FunTolerance',1e-6);
opt=contset(opt,'MaxStepsize',MaxStepsize_);
opt=contset(opt,'MaxNumPoints',3000);
opt=contset(opt,'InitStepsize',0.000001);
opt=contset(opt,'TSearchOrder',1);
opt=contset(opt,'Backward',1);
[x4,v4,s4,h4,f4]=cont(@equilibrium,x0,v0,opt);
gwcplanim(x4,v4,s4,'m',[2*N+1,1]);


save([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
save([matcont_path '/data_matcont2.mat'], 'x2', 'v2', 's2', 'h2', 'f2')
save([matcont_path '/data_matcont3.mat'], 'x3', 'v3', 's3', 'h3', 'f3')
save([matcont_path '/data_matcont4.mat'], 'x4', 'v4', 's4', 'h4', 'f4')

