clear all

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

cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); 
disp(degree')

laplacian_matrix = adjacent_matrix - diag(degree);

%%
N = size(laplacian_matrix,1);

figure;
x = cos(linspace(pi/2, pi/2+2*pi, N+1)); 
y = sin(linspace(pi/2, pi/2+2*pi, N+1)); 
x = x(1:end-1);
y = y(1:end-1);

nodePositions = [x;y]'; 
str_NodeLabel = [repmat('',N,1) num2str([1:N]')];
G = graph(adjacent_matrix);
h = plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), 'NodeLabel', cellstr(str_NodeLabel));

h.NodeCData = ones(1,N); 
h.MarkerSize = 8; 

h.LineWidth = 1.5; 
h.EdgeColor = 'r'; 

xlim([-1.5 1.5])
ylim([-1.5 1.5])
set(gca, "XTick",[])
set(gca, "YTick",[])
axis equal;

%%
N = size(laplacian_matrix,1);
syms_uarray = sym('u', [1 N]);
syms_varray = sym('v', [1 N]);

%%
fn='unppp=un - un^3 - vn + d1*(Sumun)';
gn='vnppp=c*(un - a*vn) + d2*(Sumvn)';            

filename = [results_path 'FHNSystemOnNetwork1.txt'];
fid = fopen(filename, 'w');

fprintf(fid, '%s\n', ['FHNSystemOnNetwork1']);
fprintf(fid, '%s\n','');

%%
fprintf(fid, '%s\n','For matcont');
fprintf(fid, '%s\n','');
% P = [rhou, rhos, a, b];

variables_name_base = ['un,vn'];
variables_name = ['u1,v1'];
Xvariables_name = ['x1,x2'];
for i =2:N
    variables_name = [variables_name ',' strrep(variables_name_base,'n',num2str(i))];
    Xvariables_name = [Xvariables_name, ',' 'x' num2str(2*i-1) ',' 'x' num2str(2*i)];
end
disp(variables_name)

fprintf(fid, '%s\n', variables_name);
% fprintf(fid, '%s\n','');

parameters_name = ['a, c, d1, d2'];
fprintf(fid, '%s\n', parameters_name);
fprintf(fid, '%s\n','');

for i = 1:N

    syms_Sumun = sym('S');
    for j=1:N
        syms_Sumun = syms_Sumun(1) + syms_uarray(j) * laplacian_matrix(i,j);
    end
    str_Sumun = char(syms_Sumun - sym('S'));
    fni = strrep(fn,'n',num2str(i));
    fni = strrep(fni,['Sumu' num2str(i)],str_Sumun);
    fni = strrep(fni,'ppp',"'");
    fprintf(fid, '%s\n', fni);

    syms_Sumvn = sym('S');
    for j=1:N
        syms_Sumvn = syms_Sumvn(1) + syms_varray(j) * laplacian_matrix(i,j);
    end
    str_Sumvn = char(syms_Sumvn - sym('S'));
    gni = strrep(gn,'n',num2str(i));
    gni = strrep(gni,['Sumv' num2str(i)],str_Sumvn);
    gni = strrep(gni,'ppp',"'");
    fprintf(fid, '%s\n', gni);

%     fprintf(fid, '%s\n','');
end
fprintf(fid, '%s\n','');

%% 
% fprintf(fid, '%s\n','For ode45 without noise');
% fprintf(fid, '%s\n', 'dx = zeros(2*N,1);');

for i = 1:N

    syms_Sumun = sym('S');
    for j=1:N
        syms_Sumun = syms_Sumun(1) + syms_uarray(j) * laplacian_matrix(i,j);
    end
    str_Sumun = char(syms_Sumun - sym('S'));
    fni = strrep(fn,'n',num2str(i));
    fni = strrep(fni,['Sumu' num2str(i)],str_Sumun);

    for k=N:-1:1
        fni = strrep(fni,['u' num2str(k)],['x' num2str(2*k-1)]);
        fni = strrep(fni,['v' num2str(k)],['x' num2str(2*k  )]);
    end

    fni = strrep(fni,'ppp',"'");
    fprintf(fid, '%s\n', fni);

    syms_Sumvn = sym('S');
    for j=1:N
        syms_Sumvn = syms_Sumvn(1) + syms_varray(j) * laplacian_matrix(i,j);
    end
    str_Sumvn = char(syms_Sumvn - sym('S'));
    gni = strrep(gn,'n',num2str(i));
    gni = strrep(gni,['Sumv' num2str(i)],str_Sumvn);

    for k=N:-1:1
        gni = strrep(gni,['u' num2str(k)],['x' num2str(2*k-1)]);
        gni = strrep(gni,['v' num2str(k)],['x' num2str(2*k  )]);
    end

    gni = strrep(gni,'ppp',"'");
    fprintf(fid, '%s\n', gni);

%     fprintf(fid, '%s\n','');
end
fprintf(fid, '%s\n','');

%% ode45 without noise
fprintf(fid, '%s\n','%% For ode45 without noise');

fprintf(fid, '%s\n',['function dx = FHNSystemOnNetwork1(x, N, a, b, c, d1, d2)']);
fprintf(fid, '%s\n','');

fprintf(fid, '%s\n', 'dx = zeros(2*N,1);');


% fn='dunppp=a*un*(1-un/b) - b*un*vn/(b*un+vn) + d1*(Sumun);';
% gn='dvnppp=b*un*vn/(b*un+vn) - c*vn + d2*(Sumvn);'; 

fn='dunppp=un - un^3 - vn + d1*(Sumun);';
gn='dvnppp=c*(un - a*vn) + d2*(Sumvn);';    

for i = 1:N

    syms_Sumun = sym('S');
    for j=1:N
        syms_Sumun = syms_Sumun(1) + syms_uarray(j) * laplacian_matrix(i,j);
    end
    str_Sumun = char(syms_Sumun - sym('S'));
    fni = strrep(fn,'n',num2str(i));
    fni = strrep(fni,['Sumu' num2str(i)],str_Sumun);

    for k=N:-1:1
        fni = strrep(fni,['u' num2str(k)],['x(' num2str(2*k-1) ')']);
        fni = strrep(fni,['v' num2str(k)],['x(' num2str(2*k  ) ')']);
    end

    fni = strrep(fni,'ppp',"");
    fprintf(fid, '%s\n', fni);

    syms_Sumvn = sym('S');
    for j=1:N
        syms_Sumvn = syms_Sumvn(1) + syms_varray(j) * laplacian_matrix(i,j);
    end
    str_Sumvn = char(syms_Sumvn - sym('S'));
    gni = strrep(gn,'n',num2str(i));
    gni = strrep(gni,['Sumv' num2str(i)],str_Sumvn);

    for k=N:-1:1
        gni = strrep(gni,['u' num2str(k)],['x(' num2str(2*k-1) ')']);
        gni = strrep(gni,['v' num2str(k)],['x(' num2str(2*k  ) ')']);
    end

    gni = strrep(gni,'ppp',"");
    fprintf(fid, '%s\n', gni);

%     fprintf(fid, '%s\n','');
end
fprintf(fid, '%s\n','');
fprintf(fid, '%s\n','end');
fprintf(fid, '%s\n','');

%%

fprintf(fid, '%s\n',['function F = FHNSystemOnNetwork1(x, a, b, c, d1, d2)']);
fprintf(fid, '%s\n','');
% fn='a*un*(1-un/b) - b*un*vn/(b*un+vn) + d1*(Sumun);';
% gn='b*un*vn/(b*un+vn) - c*vn + d2*(Sumvn);'; 
fn='un - un^3 - vn + d1*(Sumun);';
gn='c*(un - a*vn) + d2*(Sumvn);';  
for i = 1:N

    syms_Sumun = sym('S');
    for j=1:N
        syms_Sumun = syms_Sumun(1) + syms_uarray(j) * laplacian_matrix(i,j);
    end
    str_Sumun = char(syms_Sumun - sym('S'));
    fni = strrep(fn,'n',num2str(i));
    fni = strrep(fni,['Sumu' num2str(i)],str_Sumun);

    for k=N:-1:1
        fni = strrep(fni,['u' num2str(k)],['x(' num2str(2*k-1) ')']);
        fni = strrep(fni,['v' num2str(k)],['x(' num2str(2*k  ) ')']);
    end

    fni = strrep(fni,'ppp',"'");
    fprintf(fid, '%s\n', fni);

    syms_Sumvn = sym('S');
    for j=1:N
        syms_Sumvn = syms_Sumvn(1) + syms_varray(j) * laplacian_matrix(i,j);
    end
    str_Sumvn = char(syms_Sumvn - sym('S'));
    gni = strrep(gn,'n',num2str(i));
    gni = strrep(gni,['Sumv' num2str(i)],str_Sumvn);

    for k=N:-1:1
        gni = strrep(gni,['u' num2str(k)],['x(' num2str(2*k-1) ')']);
        gni = strrep(gni,['v' num2str(k)],['x(' num2str(2*k  ) ')']);
    end

    gni = strrep(gni,'ppp',"'");
    fprintf(fid, '%s\n', gni);

%     fprintf(fid, '%s\n','');
end
fprintf(fid, '%s\n','');
fprintf(fid, '%s\n',['end']);

fprintf(fid, '%s\n', variables_name);
fprintf(fid, '%s\n', Xvariables_name);

fclose(fid);

disp(['Finish writing ' filename]);

