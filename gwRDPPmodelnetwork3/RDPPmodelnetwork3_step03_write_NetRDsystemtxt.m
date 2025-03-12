clear all
close all
clc

%%
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

cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2);
disp(degree')

laplacian_matrix = adjacent_matrix - diag(degree);


%%
N = size(laplacian_matrix,1);
syms_uarray = sym('u', [1 N]);
syms_varray = sym('v', [1 N]);

%%
%
fn='unppp=a*un*(1-un/b) - b*un*vn/(b*un+vn) + d1*(Sumun)';
gn='vnppp=b*un*vn/(b*un+vn) - c*vn + d2*(Sumvn)';            

filename = [results_path  'RDPPSystemOnNetwork3.txt'];
fid = fopen(filename, 'w');

fprintf(fid, '%s\n', ['RDPPSystemOnNetwork3']);
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

parameters_name = ['a, b, c, d1, d2'];
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

fprintf(fid, '%s\n',['function dx = RDPPSystemOnNetwork3(x, N, a, b, c, d1, d2)']);
fprintf(fid, '%s\n','');

fprintf(fid, '%s\n', 'dx = zeros(2*N,1);');



% fn='dunppp=a*un*(1-un/b) - b*un*vn/(b*un+vn) + d1*(Sumun)';
% gn='dvnppp=b*un*vn/(b*un+vn) - c*vn + d2*(Sumvn)';     

fn='dunppp=a*un*(1-un/b) - b*un*vn/(b*un+vn) + d1*(Sumun);';
gn='dvnppp=b*un*vn/(b*un+vn) - c*vn + d2*(Sumvn);'; 

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

fprintf(fid, '%s\n',['function F = RDPPSystemOnNetwork3(x, a, b, c, d1, d2)']);
fprintf(fid, '%s\n','');
fn='a*un*(1-un/b) - b*un*vn/(b*un+vn) + d1*(Sumun);';
gn='b*un*vn/(b*un+vn) - c*vn + d2*(Sumvn);'; 

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

