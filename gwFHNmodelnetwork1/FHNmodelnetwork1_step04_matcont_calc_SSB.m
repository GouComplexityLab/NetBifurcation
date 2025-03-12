clear all
close all
clc

%%
network_name = 'npbcws_n=10k=2p=0number=1';
datanetwork_path = ['./network1 ' network_name '/'];
results_path = strrep(datanetwork_path,'network1','FHN results network1');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

net_number = 1;
mat_name = [datanetwork_path 'adj_matrix_nwnetwork_' num2str(net_number,'%02d')];
load(mat_name)
adjacent_matrix = full(adj_matrix);
%%
matcont_path = [results_path 'matcont NetFHNsystem ' network_name ' ssb'];
if ~exist(matcont_path, 'dir')
    mkdir(matcont_path);
end
%%
cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); 
disp(degree')

laplacian_matrix = adjacent_matrix - diag(degree);

N = size(laplacian_matrix,1);
[eigenvectors, eigenvalues] = eig(laplacian_matrix);
eigenvalues = diag(eigenvalues);
eigenvectors = fliplr(eigenvectors);
eigenvalues = flipud(eigenvalues);
eigenvalues(1) = 0;
eigenvectors(:,1) = - eigenvectors(:,1);

%% set diffusion rate
theta = -1;
du = 0.3; 
fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.8f\n', -theta*du/(-theta+du*eigenvalues(2)));

for i=1:N-1
    panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
end
disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])

dv = 0.8;
c_set = 1.8;
c = [0:0.001:3];  
ymax = 2;

%%
panbie_H = -theta + eigenvalues * (du + dv);
p_H = find(panbie_H>0);
Ha = zeros(length(p_H),length(c));
for i=1:length(p_H)
    Ha(i,:) = (-theta + eigenvalues(p_H(i)) * (du + dv))./c; 
end

%%
panbie_T = -theta + eigenvalues * du;
p_T = find(panbie_T>0);
Ta = zeros(length(p_T),length(c));
for i=1:length(p_T)
    Ta(i,:) = 1/(du*eigenvalues(p_T(i)) -theta) + 1./c * dv * eigenvalues(p_T(i));
end

%% determine intersecting points of T_{k} and T_{k+1}, its number = length(p_T)-1?

intersecting_points_TT = zeros(2,length(p_T)-1);
panbie_intersect_upper_x = zeros(1,length(p_T)-1);
panbie_intersect_upper_H1 = zeros(1,length(p_T)-1);
for i=1:length(p_T)-1
    intersecting_c_ii1 = dv/du * ( -theta + du*eigenvalues(p_T(i)) ) * ( -theta + du*eigenvalues(p_T(i+1)) );
    intersecting_a_ii1 = 1/(du*eigenvalues(p_T(i)) -theta) + 1/intersecting_c_ii1 * dv * eigenvalues(p_T(i));
    intersecting_points_TT(:,i) = [intersecting_c_ii1,intersecting_a_ii1]';
    panbie_intersect_upper_x(i) = -theta + du*eigenvalues(p_T(i)) + du*eigenvalues(p_T(i+1));
    panbie_intersect_upper_H1(i) = du * dv * (  eigenvalues(p_T(i)) + eigenvalues(p_T(i+1))  ) + theta*(du - dv);
end

p_intersect_upper_x = find(panbie_intersect_upper_x>=0);
intersecting_points_TT_upper_x = intersecting_points_TT(:,p_intersect_upper_x);

p_intersect_upper_H1 = find(panbie_intersect_upper_H1>=0);

intersecting_points_TT_upper_H1 = intersecting_points_TT(:,p_intersect_upper_H1);

intersect_c_H1Tend = - du*dv*eigenvalues(p_intersect_upper_H1(end)+1)^2 - theta* (du-dv) * eigenvalues(p_intersect_upper_H1(end)+1) + theta^2 ;
intersect_a_H1Tend = 1/(du*eigenvalues(p_intersect_upper_H1(end)+1) -theta) + 1/intersect_c_H1Tend * dv * eigenvalues(p_intersect_upper_H1(end)+1);

%% 

c = c_set;
intersecting_points_c_TT_upper_H1 = intersecting_points_TT_upper_H1(1,:);
disp(['intersecting SSB points for c:   ' num2str(intersecting_points_c_TT_upper_H1)])
% [v p_critical_SSB] = min(abs(intersecting_points_c_TT_upper_H1-c))
p_critical_SSB = find(intersecting_points_c_TT_upper_H1-c > 0);
p_critical_SSB = p_critical_SSB(end);
a_critical_SSB = 1/(du*eigenvalues(p_critical_SSB + 1) -theta) + 1./c * dv * eigenvalues(p_critical_SSB + 1);

fprintf('critical SSB point a=%.8f, with setting c=%.8f, on Turing curve k=%d.\n', a_critical_SSB, c, p_critical_SSB+1);


a_critical_HB = (-theta + eigenvalues(1) * (du + dv))./c; 
fprintf('critical HB point a=%.8f, with setting c=%.8f, on Hopf curve k=%d.\n', a_critical_HB, c, 1);


format rat
a_star = a_critical_SSB
c
k = p_critical_SSB+1
format long


%%
n = k;
phin = real(eigenvectors(:,n));

%%
mu0 = 0;
a0 = a_star + mu0;
u_star0 = 0; v_star0 = 0;

mu = 0.001;
a_first = a_star + mu;
u_star = 0; v_star = 0;

x0=zeros(1,2*N);
x0(1:2:end-1) = u_star + phin* 0.05;
x0(2:2:end  ) = v_star + phin* 0.05;

% x0(1:2:end-1) = u_star - phin* 0.05;
% x0(2:2:end  ) = v_star - phin* 0.05;

hls=feval(@FHNSystemOnNetwork1);
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(2*N,1));
% options = odeset('Jacobian',hls(3),'JacobianP',hls(4),'Hessians',hls(5),'HessiansP',hls(6));
[t,y]=ode45(hls{2},[0:0.1:1e4],x0,options, a_first, c, du, dv);

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

xlim([a_star-0.003 a_star+0.003])
ylim([-0.1 0.1])

xlabel('a','FontSize',18,'Interpreter','latex')
ylabel('$u_{1}$','FontSize',18,'Interpreter','latex')

box on
set(gca,'XColor','k','YColor','k','TickLength',...
    [0.02 0.05],'FontSize',18,'linewidth',1.0,'layer','top');

plot(a0, u_star0, 'k*')
plot(a0, u_star0, 'ko')

plot(a_first, u_star, 'r*')
plot(a_first, simulated_u1(end,1), 'bd')


MaxStepsize_ = 1e-3;

global cds
p=[a_first, c, du, dv]; ap=[1];
[x0,v0]=init_EP_EP(@FHNSystemOnNetwork1,initial_x0,p,ap);
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

%% BP-EP
BP_pos = 2;
initial_x1=x2(1:2*N,s2(BP_pos).index);
p(ap)=x2(2*N+1,s2(BP_pos).index);
% plot(p(ap),initial_x1(1),'MarkerSize',8,'Marker','*','LineStyle','none','Color','m');
[x0,v0]=init_BP_EP(@FHNSystemOnNetwork1,initial_x1,p,s2(BP_pos),0.01);
% BP_pos = 2;
% initial_x1=x1(1:2*N,s1(BP_pos).index);
% p(ap)=x1(2*N+1,s1(BP_pos).index);
% [x0,v0]=init_BP_EP(@FHNSystemOnNetwork1,initial_x1,p,s1(BP_pos),0.01);
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

