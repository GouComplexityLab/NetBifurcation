clear all
close all
clc

init()

network_name = 'ws_n=10k=2p=0.05number=796';
datanetwork_path = ['./network2 ' network_name '/'];
results_path = strrep(datanetwork_path,'network2','RDPPmodel results network2');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

net_number = 796;
mat_name = [datanetwork_path 'adj_matrix_wsnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);

matcont_path = [results_path 'matcont hb wsnetworks_n=10k=2p=0.05number=796'];
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

b1 = b_star2 + mu0;
u_star1 = b1*(a+(c-1)*b1)/a; v_star1 = b1*(1-c)/c * u_star1;

mu = 0.01;
b_first = b_star + mu;
u_star = b_first*(a+(c-1)*b_first)/a; v_star = b_first*(1-c)/c * u_star;

x0=zeros(1,2*N);
x0(1:2:end-1) = u_star + phin* 0.05;
x0(2:2:end  ) = v_star + phin* 0.05;

% x0(1:2:end-1) = u_star - phin* 0.05;
% x0(2:2:end  ) = v_star - phin* 0.05;

hls=feval(@RDPPSystemOnNetwork2);
options=odeset('RelTol',1e-8);
% options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
[t,y]=ode45(hls{2},[0:1:25e4],x0,options, a, b_first, c, d1, d2);

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

xlim([b_star-0.03 b_star+0.03])
ylim([0 1])

xlabel('b','FontSize',18,'Interpreter','latex')
ylabel('$u_{1}$','FontSize',18,'Interpreter','latex')

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'FontSize',18,'linewidth',1.0,'layer','top');

plot(b0, u_star0, 'k*')
plot(b0, u_star0, 'ko')

plot(b1, u_star1, 'r*')
plot(b1, u_star1, 'ro')

plot(b_first, u_star, 'r*')
plot(b_first, simulated_u1(end,1), 'bd')


MaxStepsize_ = 1e-3;

global cds
p=[a, b_first, c, d1, d2]; ap=[2];
plot(p(ap), initial_x0(1), 'm*')
plot(p(ap), initial_x0(1), 'mo')
[x0,v0]=init_EP_EP(@RDPPSystemOnNetwork2,initial_x0,p,ap);
opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxStepsize',1e-4);
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


%% BP-EP
BP_pos = 2;
initial_x1=x2(1:2*N,s2(BP_pos).index);
p(ap)=x2(2*N+1,s2(BP_pos).index);
% plot(p(ap),initial_x1(1),'MarkerSize',8,'Marker','*','LineStyle','none','Color','m');
[x0,v0]=init_BP_EP(@RDPPSystemOnNetwork2,initial_x1,p,s2(BP_pos),0.01);
% BP_pos = 2;
% initial_x1=x1(1:2*N,s1(BP_pos).index);
% p(ap)=x1(2*N+1,s1(BP_pos).index);
% [x0,v0]=init_BP_EP(@RDPPSystemOnNetwork2,initial_x1,p,s1(BP_pos),0.01);
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


%% H LC
% HP_pos = 2;
% initial_x1=x2(1:2*N,s2(HP_pos).index);
% p(ap)=x2(2*N+1,s2(HP_pos).index);
% [x0,v0]=init_H_LC(@RDPPSystemOnNetwork2,initial_x1,p,[2],1e-6,20,4);
HP_pos = 2;
initial_x1=x1(1:2*N,s1(HP_pos).index);
p(ap)=x1(2*N+1,s1(HP_pos).index);
[x0,v0]=init_H_LC(@RDPPSystemOnNetwork2,initial_x1,p,[2],1e-6,20,4);

opt = contset; 
% opt = contset(opt,'Singularities',1);
opt = contset(opt,'MaxNumPoints',1000);
opt = contset(opt,'Multipliers',1);
opt = contset(opt,'Adapt',1);

% opt = contset(opt,'Multipliers',0);
% opt = contset(opt,'Adapt',1);
% opt = contset(opt,'MaxStepsize',5);
% opt = contset(opt,'FunTolerance',1e-6);
% opt = contset(opt,'VarTolerance',1e-6);
% opt = contset(opt,'PRC',1);
% opt = contset(opt,'dPRC',1);
% opt = contset(opt,'Input',1);
% opt = contset(opt,'MaxNumPoints',100);

% opt=contset;
% opt=contset(opt,'IgnoreSingularity',1);
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'MaxNumPoints',50);
% opt=contset(opt,'FunTolerance',0.0000001);
% opt=contset(opt,'VarTolerance',0.0000001);

[xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);


figure 
hold on
xlim([b_star-0.03 b_star+0.03])
ylim([0 1])

plotcycle(xlc,vlc,slc,[size(xlc,1) 1 2]);
xlabel('x')
ylabel('y')
zlabel('z')

view(3)

%%


global lds

x=xlc;v=vlc;s=slc;e=[size(xlc,1) 1 2];

figure 
hold on

gwcplanim2(x1,v1,s1,'k',[2*N+1,1,2]);
gwcplanim2(x2,v2,s2,'k',[2*N+1,1,2]);
% gwcplanim2(x3,v3,s3,'k',[2*N+1,1,2]);
% gwcplanim2(x4,v4,s4,'k',[2*N+1,1,2]);

xlim([b_star-0.03 b_star+0.03])
ylim([0 1])

x_p = zeros(1,size(x,2));
ymin_p = zeros(1,size(x,2));
ymax_p = zeros(1,size(x,2));

zmin_p = zeros(1,size(x,2));
zmax_p = zeros(1,size(x,2));

for j = 1:size(x,2)

    xx = x(e(1),j)*ones(lds.tps,1);
    yy = x((0:lds.tps-1)*lds.nphase+e(2),j);
    zz = x((0:lds.tps-1)*lds.nphase+e(3),j);
    plot3(xx,yy,zz,'b-');

    x_p(j) = xx(1);
    ymin_p(j) = min(yy);
    ymax_p(j) = max(yy);
    
    zmin_p(j) = min(zz);
    zmax_p(j) = max(zz);

end

xlabel('x')
ylabel('y')
zlabel('z')
view(3)

figure 
hold on

gwcplanim2(x1,v1,s1,'k',[2*N+1,1]);
gwcplanim2(x2,v2,s2,'k',[2*N+1, 1]);
gwcplanim2(x3,v3,s3,'k',[2*N+1,1]);
gwcplanim2(x4,v4,s4,'k',[2*N+1,1]);

plot(x_p,ymin_p,'r')
plot(x_p,ymax_p,'r')

xlim([b_star-0.03 b_star+0.03])
ylim([0 1])




save([matcont_path '/data_matcont1.mat'], 'x1', 'v1', 's1', 'h1', 'f1')
save([matcont_path '/data_matcont2.mat'], 'x2', 'v2', 's2', 'h2', 'f2')
save([matcont_path '/data_matcont3.mat'], 'x3', 'v3', 's3', 'h3', 'f3')
save([matcont_path '/data_matcont4.mat'], 'x4', 'v4', 's4', 'h4', 'f4')
save([matcont_path '/data_matcontHB.mat'], 'xlc', 'vlc', 'slc', 'hlc', 'flc')
