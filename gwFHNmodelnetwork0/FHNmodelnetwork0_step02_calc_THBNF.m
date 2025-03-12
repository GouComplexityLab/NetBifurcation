clear all
close all
clc

tic
digits = 50;
format long


network_name = 'brain';
datanetwork_path = ['./network0 ' network_name '/'];
results_path = strrep(datanetwork_path,'network0','FHN results network0');
if ~exist(results_path, 'dir')
    mkdir(results_path);
end

mat_name = [datanetwork_path 'adj_matrix_brainnetwork'];
load(mat_name)
adjacent_matrix = full(adj_matrix);


cha = isequal(adjacent_matrix, adjacent_matrix');
disp(['is symmetrical:' num2str(cha)])
degree = sum(adjacent_matrix,2); % 行和
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
fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.6f\n', -theta*du/(-theta+du*eigenvalues(2)));

for i=1:N-1
    panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
end
disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])
% dv = panbie_THB(6);
% dv = (panbie_THB(5)+panbie_THB(6))/2;
dv = 0.8;
% c_set = 1.8;
c = [0:0.001:3];  
ymax = 2;
% %% set diffusion rate
% theta = -2;
% du = 0.04; 
% fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.6f\n', -theta*du/(-theta+du*eigenvalues(2)));
% 
% for i=1:N-1
%     panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
% end
% disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])
% % dv = panbie_THB(4);
% dv = (panbie_THB(3)+panbie_THB(4))/2;
% % dv = 0.06;
% c_set = 6;
% c = [0:0.001:10];  
% ymax = 0.8;

% %% set diffusion rate
% theta = -1/2;
% du = 0.01; 
% fprintf('must > -theta*du/(-theta+du*eigenvalues(2)) = %.6f\n', -theta*du/(-theta+du*eigenvalues(2)));
% 
% for i=1:N-1
%     panbie_THB(i) = -theta*du / ( -theta + du*(eigenvalues(i)+eigenvalues(i+1)) );
% end
% disp([num2str(length(find(panbie_THB>0))) ' special dv for THB points:   ' num2str(panbie_THB(find(panbie_THB>0)))])
% % dv = panbie_THB(8);
% dv = (panbie_THB(7)+panbie_THB(8))/2;
% % dv = 0.06;
% c_set = 1;
% c = [0:0.001:2];  
% ymax = 3;
%%
panbie_H = -theta + eigenvalues * (du + dv);
p_H = find(panbie_H>0);
Ha = zeros(length(p_H),length(c));
for i=1:length(p_H)
    Ha(i,:) = (-theta + eigenvalues(p_H(i)) * (du + dv))./c; 
end
% figure
% hold on
% grid minor
% for i=1:length(p_H)
%     plot(c,Ha(i,:),'b-');
% end

%%
panbie_T = -theta + eigenvalues * du;
p_T = find(panbie_T>0);
Ta = zeros(length(p_T),length(c));
for i=1:length(p_T)
    Ta(i,:) = 1/(du*eigenvalues(p_T(i)) -theta) + 1./c * dv * eigenvalues(p_T(i));
end
% figure
% hold on
% grid minor
% 
% for i=1:length(p_H)
%     plot(c,Ha(i,:),'b-');
% end
% for i=1:length(p_T)
%     plot(c,Ta(i,:),'r-');
% end

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

% p_intersect_upper_H1 = [p_intersect_upper_H1,p_intersect_upper_H1(end)+1];

intersecting_points_TT_upper_H1 = intersecting_points_TT(:,p_intersect_upper_H1);

intersect_c_H1Tend = - du*dv*eigenvalues(p_intersect_upper_H1(end)+1)^2 - theta* (du-dv) * eigenvalues(p_intersect_upper_H1(end)+1) + theta^2 ;
intersect_a_H1Tend = 1/(du*eigenvalues(p_intersect_upper_H1(end)+1) -theta) + 1/intersect_c_H1Tend * dv * eigenvalues(p_intersect_upper_H1(end)+1);


%% plot
figure
hold on
grid minor

plot(c,Ha(1,:),'b-');

for i=1:length(p_H)
    plot(c,Ha(i,:),'b-');
end
% for i=1:length(p_T)
% for i=1:length(p_intersect_upper_x)
for i=1:size(intersecting_points_TT_upper_H1,2)+1
    plot(c,Ta(i,:),'r-');
end
% plot(intersecting_points_TT(1,:), intersecting_points_TT(2,:), 'g*')
% plot(intersecting_points_TT_upper_x(1,:), intersecting_points_TT_upper_x(2,:), 'co')
plot(intersecting_points_TT_upper_H1(1,:), intersecting_points_TT_upper_H1(2,:), 'k*')
plot(intersect_c_H1Tend,intersect_a_H1Tend,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'Marker','pentagram','LineStyle','none','Color',[0 0 0]);
fprintf('critical THB point a=%.6f, c=%.6f, on Turing curve k=%d, Hopf curve k=%d.\n', intersect_a_H1Tend, intersect_c_H1Tend, size(intersecting_points_TT_upper_H1,2)+1,1);

plot(c,ones(size(c))*1/(-theta),'k-')
plot(c,zeros(size(c)),'k-')

for i=1:length(p_intersect_upper_H1)-1
%     [i p_intersect_upper_H1(i+1)]
    intersecting_c_ii1 = intersecting_points_TT_upper_H1(1,i);
    intersecting_c_ii2 = intersecting_points_TT_upper_H1(1,i+1);
    Turing_c_i = linspace(intersecting_c_ii2, intersecting_c_ii1, 100);
    Turing_a_i = 1/(du*eigenvalues(p_intersect_upper_H1(i+1)) -theta) + 1./Turing_c_i * dv * eigenvalues(p_intersect_upper_H1(i+1));
    plot(Turing_c_i,Turing_a_i,'k-');
end
i=length(p_intersect_upper_H1);
intersecting_c_ii1 = intersecting_points_TT_upper_H1(1,i);
intersecting_c_H1Tend = - du*dv*eigenvalues(p_intersect_upper_H1(end)+1)^2 - theta* (du-dv) * eigenvalues(p_intersect_upper_H1(end)+1) + theta^2 ;
Turing_c_i = linspace(intersecting_c_H1Tend, intersecting_c_ii1, 100);
Turing_a_i = 1/(du*eigenvalues(i+1) -theta) + 1./Turing_c_i * dv * eigenvalues(i+1);
plot(Turing_c_i,Turing_a_i,'k-');

% xlim([0 40])
ylim([0 ymax])


%%
n = size(intersecting_points_TT_upper_H1,2)+1;
k = 1;

a_star = intersect_a_H1Tend;
c_star = intersect_c_H1Tend;
format rat
a_star
c_star
format long


fprintf('critical THB point a_star=%.6f, c_star=%.6f, on Turing curve n=%d, Hopf curve k=%d.\n', intersect_a_H1Tend, intersect_c_H1Tend, n,k);

a11 = -theta;
a12 = -1;
a21 = c_star;
a22 = - a_star*c_star;

J = [a11,a12; a21,a22]; 
D = [du,0;0,dv];

% Hopf, k=1

J_1 = a11*a22 - a12*a21;

omega_c = sqrt(J_1);

p11 = 1; 
p12 = (1i*omega_c - a11)/a12;
p1 = [p11; p12];

q11 = (1i*omega_c + a11)/(2*1i*omega_c);
q12 = a12/(2*1i*omega_c);
q1 = [q11; q12];

% [p1, conj(p1)].'*[q1, conj(q1)] % check!

% SSB, n
T_n = (du+dv)*eigenvalues(n) + (a11 + a22);
J_n = du*dv*eigenvalues(n)^2 + (dv*a11 + du*a22)*eigenvalues(n) + a11*a22 - a12*a21;

pn1 = 1; 
pn2 = - (a11+du*eigenvalues(n))/a12;
pn = [pn1; pn2];

qn1 = (dv*eigenvalues(n) + a22)/T_n;
qn2 = - a12/T_n;
qn = [qn1; qn2];

% pn'*qn % check!


f2000 = [2*theta + 2; 0];
f0200 = [0; 0];
f0020 = [0; 0];
f0002 = [0; 0];

f1100 = [0; 0];
f1010 = [0; 0];
f1001 = [0; 1];
f0110 = [0; -c_star];
f0101 = [0; -a_star];
f0011 = [0; 0];

f3000 = [-6; 0];
f0300 = [0; 0];
f2100 = [0; 0];
f1200 = [0; 0];
      
f2010 = [0; 0];
f2001 = [0; 0];
f0210 = [0; 0];
f0201 = [0; 0];
f1110 = [0; 0];
f1101 = [0; 0];


F200 = p11^2*f2000 + 2*p11*p12*f1100 + p12^2*f0200;
F020 = conj(F200);
F002 = pn1^2*f2000 + 2*pn1*pn2*f1100 + pn2^2*f0200;
F110 = abs(p11)^2*f2000 + 2*real(p11*conj(p12))*f1100 + abs(p12)^2*f0200;
F101 = p11*pn1*f2000 + (p11*pn2+p12*pn1)*f1100 + p12*pn2*f0200;
F011 = conj(F101);

B10010 = q1.'*(f1010*p11 + f0110*p12);
B10001 = q1.'*(f1001*p11 + f0101*p12);
B00110 = qn.'*(f1010*pn1 + f0110*pn2);
B00101 = qn.'*(f1001*pn1 + f0101*pn2);

% B10100 = q1.'*F101*sum(eigenvectors(:,1).^2.*eigenvectors(:,n));
% B11000 = qn.'*F110*sum((1/sqrt(N)*ones(N,1)).^2.*eigenvectors(:,n));
B10100 = 0;
B11000 = 0;
B00200 = 1/2 * qn.'*F002*sum(eigenvectors(:,n).^3);

if abs(B00200)<1e-8
    B00200 = 0;
end
F21000 = 1/2*( f3000*abs(p11)^2*p11 + f0300*abs(p12)^2*p12 ... 
    + f2100*(p11^2*conj(p12)+2*abs(p11)^2*p12) + f1200*(p12^2*conj(p11)+2*abs(p12)^2*p11) );
F10200 = 1/2*(f3000*p11*pn1^2 + f0300*p12*pn2^2 ... 
    + f2100*(p12*pn1^2+2*p11*pn1*pn2) + f1200*(p11*pn2^2+2*p12*pn1*pn2));
F11100 = ( f3000*abs(p11)^2*pn1 + f0300*abs(p12)^2*pn2 ... 
    + f2100*(abs(p11)^2*pn2+2*pn1*real(p11*conj(p12))) + f1200*(abs(p12)^2*pn1+2*pn2*real(p12*conj(p11))) );
F00300 = 1/6*(f3000*pn1^3 + f0300*pn2^3 + 3*f2100*pn1^2*pn2 + 3*f1200*pn1*pn2^2);
F00210 = 1/2*(f2010*pn1^2 + f0210*pn2^2 + 2*f1110*pn1*pn2);
F00201 = 1/2*(f2001*pn1^2 + f0201*pn2^2 + 2*f1101*pn1*pn2);

C21000 = q1.'*F21000/N;
C10200 = q1.'*F10200/N;
C11100 = qn.'*F11100/N;
C00300 = qn.'*F00300*sum(eigenvectors(:,n).^4); 
C00210 = qn.'*F00210*sum(eigenvectors(:,n).^3); 
C00201 = qn.'*F00201*sum(eigenvectors(:,n).^3); 
C10110 = 0;
C10101 = 0;
C11010 = 0;
C11001 = 0;


D21000 = 1/(3*omega_c*i)*( - (q1.'*F200)*(q1.'*F110)*1/N + 1/3*abs(q1.'*F020)^2*1/N ...
    + 2*abs(q1.'*F110)^2*1/N );
D10200 = 1/(3*omega_c*i)*( - (q1.'*F200)*(q1.'*F002)*1/N + (q1.'*F110)*(conj(q1).'*F002)*1/N ...
    + 2*(q1.'*F002)*(qn.'*F101)*1/N );
D11100 = - 4/(3*omega_c)*( imag((q1.'*F110)*(qn.'*F101))*1/N );
D00300 = - 2/(3*omega_c)*imag((q1.'*F002)*(qn.'*F101))*1/N;
D00210 = 0; 
D00201 = 0; 
D10110 = 0;
D10101 = 0;
D11010 = 0;
D11001 = 0;

%%
h_k20000 = zeros(2,N);
for k=1:N
    Mk = J + eigenvalues(k) * D;
    if k==1
        h_k20000(:,k) = inv(2*i*omega_c*[1,0;0,1] - Mk)*( F200- (q1.'*F200*p1 + conj(q1).'*F200*conj(p1)) )*1/sqrt(N);
    else
        h_k20000(:,k) = [0;0];
    end
end

k=1;
Mk = J + eigenvalues(k) * D;
h_120000 = inv(2*i*omega_c*[1,0;0,1] - Mk)*( F200- (q1.'*F200*p1 + conj(q1).'*F200*conj(p1)) )*1/sqrt(N);
clear k
% % % h_k02000 = zeros(2,N);
% % % for k=1:N
% % %     Mk = [d1*eigenvalues(k)+b_star*(1-c^2)-a_star, -c^2;
% % %         b_star*(1-c)^2,          d2*eigenvalues(k)-c*(1-c)];
% % %     if k==1
% % %         h_k02000(:,k) = - inv(2*i*omega_c*[1,0;0,1] + Mk)*( F020- (q1.'*F020*p1 + conj(q1).'*F020*conj(p1)) )*1/sqrt(N);
% % %     else
% % %         h_k02000(:,k) = [0;0];
% % %     end
% % % end

h_k00200 = zeros(2,N);
for k=1:N
    Mk = J + eigenvalues(k) * D;
    if k==1
        h_k00200(:,k) = - inv(Mk)*( F002- (q1.'*F002*p1 + conj(q1).'*F002*conj(p1)) )*1/sqrt(N);
    else if k==n
            sum_n3 = sum(eigenvectors(:,n).^3);
            if abs(sum_n3) < 10e-7
                sum_n3 = 0;
            end
            M = [Mk,pn];
            M = [M;[qn;0]'];
            bn00200 = ( F002- qn.'*F002*pn ) * sum_n3;
            bar_h_n00200 = - inv(M) * [bn00200;0];
            h_k00200(:,k) = bar_h_n00200([1,2]);
        else
            sum_n2_k = sum(eigenvectors(:,n).^2.*eigenvectors(:,k));
            if abs(sum_n2_k) < 10e-7
                sum_n2_k = 0;
            end
            h_k00200(:,k) = - inv(Mk)*F002*sum_n2_k;
        end
    end
end

h_k11000 = zeros(2,N);
for k=1:N
    Mk = J + eigenvalues(k) * D;
    if k==1
        h_k11000(:,k) = - 2*inv(Mk)*( F110- (q1.'*F110*p1 + conj(q1).'*F110*conj(p1)) )*1/sqrt(N);
    else
        h_k11000(:,k) = [0;0];
    end
end

h_k00110 = zeros(2,N);
for k=1:N
    Mk = J + eigenvalues(k) * D;
    if k==n
        M = [Mk,pn];
        M = [M;[qn;0]'];
        bn00110 = 2*( (f1010*pn1+f0110*pn2)- qn.'*(f1010*pn1+f0110*pn2)*pn );
        bar_h_k00110 = - inv(M) * [bn00110;0];
        h_k00110(:,k) = bar_h_k00110([1,2]);
    else
        h_k00110(:,k) = [0;0];
    end
end

h_k00101 = zeros(2,N);
for k=1:N
    Mk = J + eigenvalues(k) * D;
    if k==n
        M = [Mk,pn];
        M = [M;[qn;0]'];
        bn00101 = 2*( (f1001*pn1+f0101*pn2)- qn.'*(f1001*pn1+f0101*pn2)*pn );
        bar_h_k00101 = - inv(M) * [bn00101;0];
        h_k00101(:,k) = bar_h_k00101([1,2]);
    else
        h_k00101(:,k) = [0;0];
    end
end

h_k10100 = zeros(2,N);
for k=1:N
    Mk = J + eigenvalues(k) * D;
    if k==n
        h_k10100(:,k) = 2 * inv(i*omega_c*[1,0;0,1] - Mk)*( F101- qn.'*F101*pn )*1/sqrt(N);
    else
        h_k10100(:,k) = [0;0];
    end
end

k=n;
Mk = J + eigenvalues(k) * D;
h_n10100 = 2 * inv(i*omega_c*[1,0;0,1] - Mk)*( F101- qn.'*F101*pn )*1/sqrt(N);
clear k
%%%%%%%%%%%%%%%% 
% % % h_k10010 = zeros(2,N);
% % % for k=1:N
% % %     Mk = [d1*eigenvalues(k)+b_star*(1-c^2)-a_star, -c^2;
% % %         b_star*(1-c)^2,          d2*eigenvalues(k)-c*(1-c)];
% % %     if k==1
% % %         M = i*omega_c*[1,0;0,1] - Mk;
% % %         M = [M,p1];
% % %         M = [M;[q1;0]'];
% % %         b110010 = 2*( (f1010*p11+f0110*p12)- (q1.'*(f1010*p11+f0110*p12)*p1 + conj(q1).'*(f1010*p11+f0110*p12)*conj(p1)) );
% % %         bar_h_k10010 =  inv(M) * [b110010;0];
% % %         h_k10010(:,k) = bar_h_k10010([1,2]);
% % %     else
% % %         h_k10010(:,k) = [0;0];
% % %     end
% % % end
% % % h_k10001 = zeros(2,N);
% % % for k=1:N
% % %     Mk = [d1*eigenvalues(k)+b_star*(1-c^2)-a_star, -c^2;
% % %         b_star*(1-c)^2,          d2*eigenvalues(k)-c*(1-c)];
% % %     if k==1
% % %         M = i*omega_c*[1,0;0,1] - Mk;
% % %         M = [M,p1];
% % %         M = [M;[q1;0]'];
% % %         b110001 = 2*( (f1001*p11+f0101*p12)- (q1.'*(f1001*p11+f0101*p12)*p1 + conj(q1).'*(f1001*p11+f0101*p12)*conj(p1)) );
% % %         bar_h_k10001 =  inv(M) * [b110001;0];
% % %         h_k10001(:,k) = bar_h_k10001([1,2]);
% % %     else
% % %         h_k10001(:,k) = [0;0];
% % %     end
% % % end

%%%%%%%%%%%%%%%% 
h_k01100 = zeros(2,N);
for k=1:N
    Mk = J + eigenvalues(k) * D;
    if k==n
        h_k01100(:,k) = -2 * inv(i*omega_c*[1,0;0,1] + Mk)*( F011- qn.'*F011*pn )*1/sqrt(N);
    else
        h_k01100(:,k) = [0;0];
    end
end
k=n;
Mk = J + eigenvalues(k) * D;
h_n01100 = -2 * inv(i*omega_c*[1,0;0,1] + Mk)*( F011- qn.'*F011*pn )*1/sqrt(N);
clear k
% % h_k01010 = zeros(2,N);
% % for k=1:N
% %     Mk = [d1*eigenvalues(k)+b_star*(1-c^2)-a_star, -c^2;
% %         b_star*(1-c)^2,          d2*eigenvalues(k)-c*(1-c)];
% %     if k==1
% %         M = i*omega_c*[1,0;0,1] + Mk;
% %         M = [M,conj(p1)];
% %         M = [M;[conj(q1);0]'];
% %         b101010 = 2*( (f1010*conj(p11)+f0110*conj(p12))- (q1.'*(f1010*conj(p11)+f0110*conj(p12))*p1 + conj(q1).'*(f1010*conj(p11)+f0110*conj(p12))*conj(p1)) );
% %         bar_h_k01010 = - inv(M) * [b101010;0];
% %         h_k01010(:,k) = bar_h_k01010([1,2]);
% %     else
% %         h_k01010(:,k) = [0;0];
% %     end
% % end
% % h_k01001 = zeros(2,N);
% % for k=1:N
% %     Mk = [d1*eigenvalues(k)+b_star*(1-c^2)-a_star, -c^2;
% %         b_star*(1-c)^2,          d2*eigenvalues(k)-c*(1-c)];
% %     if k==1
% %         M = i*omega_c*[1,0;0,1] + Mk;
% %         M = [M,conj(p1)];
% %         M = [M;[conj(q1);0]'];
% %         b101001 = 2*( (f1001*conj(p11)+f0101*conj(p12))- (q1.'*(f1001*conj(p11)+f0101*conj(p12))*p1 + conj(q1).'*(f1001*conj(p11)+f0101*conj(p12))*conj(p1)) );
% %         bar_h_k01001 = - inv(M) * [b101001;0];
% %         h_k01001(:,k) = bar_h_k01001([1,2]);
% %     else
% %         h_k01001(:,k) = [0;0];
% %     end
% % end

%%
E21000_vector1 = zeros(1,N);
E21000_vector2 = zeros(1,N);
for k=1:N
    sum_12_k = sum(eigenvectors(:,1).^2.*eigenvectors(:,k));
    if abs(sum_12_k) < 10e-7
        sum_12_k = 0;
    end
    E21000_vector1(k) = 1/3 * q1.'*((f2000*p11+f1100*p12)*h_k11000(1,k)+(f0200*p12+f1100*p11)*h_k11000(2,k))*sum_12_k;
    E21000_vector2(k) = 1/3 * q1.'*((f2000*conj(p11)+f1100*conj(p12))*h_k20000(1,k)+(f0200*conj(p12)+f1100*conj(p11))*h_k20000(2,k))*sum_12_k;
end
% E21000 = sum(E21000_vector1 + E21000_vector2);

% E21000 = 1/3/sqrt(N) * q1.'*((f2000*p11+f1100*p12)*h_k11000(1,1)+(f0200*p12+f1100*p11)*h_k11000(2,1)) ...
%        + 1/3 * q1.'*((f2000*conj(p11)+f1100*conj(p12))*h_k20000(1,1)+(f0200*conj(p12)+f1100*conj(p11))*h_k20000(2,1));

E21000 = 1/3/sqrt(N) * q1.'*((f2000*p11+f1100*p12)*h_k11000(1,1)+(f0200*p12+f1100*p11)*h_k11000(2,1)) ...
   + 1/3 * q1.'*((f2000*conj(p11)+f1100*conj(p12))*h_120000(1)+(f0200*conj(p12)+f1100*conj(p11))*h_120000(2));
   
E10200_vector1 = zeros(1,N);
E10200_vector2 = zeros(1,N);
for k=1:N
    sum_12_k = sum(eigenvectors(:,1).^2.*eigenvectors(:,k));
    if abs(sum_12_k) < 10e-7
        sum_12_k = 0;
    end
    sum_1n_k = sum(eigenvectors(:,1).*eigenvectors(:,n).*eigenvectors(:,k));
    if abs(sum_1n_k) < 10e-7
        sum_1n_k = 0;
    end
    E10200_vector1(k) = 1/3 * q1.'*((f2000*p11+f1100*p12)*h_k00200(1,k)+(f0200*p12+f1100*p11)*h_k00200(2,k))*sum_12_k;
    E10200_vector2(k) = 1/3 * q1.'*((f2000*pn1+f1100*pn2)*h_k10100(1,k)+(f0200*pn2+f1100*pn1)*h_k10100(2,k))*sum_1n_k;
end
% E10200 = sum(E10200_vector1 + E10200_vector2);
% E10200 = 1/3/sqrt(N) * q1.'*((f2000*p11+f1100*p12)*h_k00200(1,1)+(f0200*p12+f1100*p11)*h_k00200(2,1)) ...
%        + 1/3/sqrt(N) * q1.'*((f2000*pn1+f1100*pn2)*h_k10100(1,n)+(f0200*pn2+f1100*pn1)*h_k10100(2,n));
E10200 = 1/3/sqrt(N) * q1.'*((f2000*p11+f1100*p12)*h_k00200(1,1)+(f0200*p12+f1100*p11)*h_k00200(2,1)) ...
       + 1/3/sqrt(N) * q1.'*((f2000*pn1+f1100*pn2)*h_n10100(1)+(f0200*pn2+f1100*pn1)*h_n10100(2));
   
E11100_vector1 = zeros(1,N);
E11100_vector2 = zeros(1,N);
E11100_vector3 = zeros(1,N);
for k=1:N
    sum_1n_k = sum(eigenvectors(:,1).*eigenvectors(:,n).*eigenvectors(:,k));
    if abs(sum_1n_k) < 10e-7
        sum_1n_k = 0;
    end
    sum_n2_k = sum(eigenvectors(:,n).^2.*eigenvectors(:,k));
    if abs(sum_n2_k) < 10e-7
        sum_n2_k = 0;
    end
    E11100_vector1(k) = 1/3 * qn.'*((f2000*p11+f1100*p12)*h_k01100(1,k)+(f0200*p12+f1100*p11)*h_k01100(2,k))*sum_1n_k;
    E11100_vector2(k) = 1/3 * qn.'*((f2000*conj(p11)+f1100*conj(p12))*h_k10100(1,k)+(f0200*conj(p12)+f1100*conj(p11))*h_k10100(2,k))*sum_1n_k;
    E11100_vector3(k) = 1/3 * qn.'*((f2000*pn1+f1100*pn2)*h_k11000(1,k)+(f0200*pn2+f1100*pn1)*h_k11000(2,k))*sum_n2_k;
end
% E11100 = sum(E11100_vector1 + E11100_vector2 + E11100_vector3);
% E11100 = 1/3/sqrt(N) * qn.'*((f2000*p11+f1100*p12)*h_k01100(1,n)+(f0200*p12+f1100*p11)*h_k01100(2,n)) ...
%        + 1/3/sqrt(N) * qn.'*((f2000*conj(p11)+f1100*conj(p12))*h_k10100(1,n)+(f0200*conj(p12)+f1100*conj(p11))*h_k10100(2,n)) ...
%        + sum(E11100_vector3);

E11100 = 1/3/sqrt(N) * qn.'*((f2000*p11+f1100*p12)*h_n01100(1)+(f0200*p12+f1100*p11)*h_n01100(2)) ...
   + 1/3/sqrt(N) * qn.'*((f2000*conj(p11)+f1100*conj(p12))*h_n10100(1)+(f0200*conj(p12)+f1100*conj(p11))*h_n10100(2)) ...
   + sum(E11100_vector3);

E00300_vector = zeros(1,N);
for k=1:N
    sum_n2_k = sum(eigenvectors(:,n).^2.*eigenvectors(:,k));
    if abs(sum_n2_k) < 10e-7
        sum_n2_k = 0;
    end
    E00300_vector(k) = 1/3 * qn.'*((f2000*pn1+f1100*pn2)*h_k00200(1,k)+(f0200*pn2+f1100*pn1)*h_k00200(2,k))*sum_n2_k;
end
E00300 = sum(E00300_vector);

% % % E10110_vectr1 = zeros(1,N);
% % % E10110_vectr2 = zeros(1,N);
% % % E10110_vectr3 = zeros(1,N);
% % % for k=1:N
% % %     sum_12_k = sum(eigenvectors(:,1).^2.*eigenvectors(:,k));
% % %     if abs(sum_12_k) < 10e-7
% % %         sum_12_k = 0;
% % %     end
% % %     sum_1n_k = sum(eigenvectors(:,1).*eigenvectors(:,n).*eigenvectors(:,k));
% % %     if abs(sum_1n_k) < 10e-7
% % %         sum_1n_k = 0;
% % %     end
% % %     sum_1_k = sum(eigenvectors(:,1).*eigenvectors(:,k));
% % %     if abs(sum_1_k) < 10e-7
% % %         sum_1_k = 0;
% % %     end
% % %     E10110_vectr1(k) = 1/3 * q1.'*((f2000*p11+f1100*p12)*h_k00110(1,k)+(f0200*p12+f1100*p11)*h_k00110(2,k))*sum_12_k;
% % %     E10110_vectr2(k) = 1/3 * q1.'*((f2000*pn1+f1100*pn2)*h_k10010(1,k)+(f0200*pn2+f1100*pn1)*h_k10010(2,k))*sum_1n_k;
% % %     E10110_vectr3(k) = 1/3 * q1.'*(f1010*h_k10100(1,k)+f0110*h_k10100(2,k))*sum_1_k;
% % % end
% % % E10110 = sum(E10110_vectr1 + E10110_vectr2 + E10110_vectr3);
E10110 = 0;

% % % E10101_vectr1 = zeros(1,N);
% % % E10101_vectr2 = zeros(1,N);
% % % E10101_vectr3 = zeros(1,N);
% % % for k=1:N
% % %     sum_12_k = sum(eigenvectors(:,1).^2.*eigenvectors(:,k));
% % %     if abs(sum_12_k) < 10e-7
% % %         sum_12_k = 0;
% % %     end
% % %     sum_1n_k = sum(eigenvectors(:,1).*eigenvectors(:,n).*eigenvectors(:,k));
% % %     if abs(sum_1n_k) < 10e-7
% % %         sum_1n_k = 0;
% % %     end
% % %     sum_1_k = sum(eigenvectors(:,1).*eigenvectors(:,k));
% % %     if abs(sum_1_k) < 10e-7
% % %         sum_1_k = 0;
% % %     end
% % %     E10101_vectr1(k) = 1/3 * q1.'*((f2000*p11+f1100*p12)*h_k00101(1,k)+(f0200*p12+f1100*p11)*h_k00101(2,k))*sum_12_k;
% % %     E10101_vectr2(k) = 1/3 * q1.'*((f2000*pn1+f1100*pn2)*h_k10001(1,k)+(f0200*pn2+f1100*pn1)*h_k10001(2,k))*sum_1n_k;
% % %     E10101_vectr3(k) = 1/3 * q1.'*(f1001*h_k10100(1,k)+f0101*h_k10100(2,k))*sum_1_k;
% % % end
% % % E10101 = sum(E10101_vectr1 + E10101_vectr2 + E10101_vectr3);
E10101 = 0;

% % % E11010_vectr1 = zeros(1,N);
% % % E11010_vectr2 = zeros(1,N);
% % % E11010_vectr3 = zeros(1,N);
% % % for k=1:N
% % %     sum_1n_k = sum(eigenvectors(:,1).*eigenvectors(:,n).*eigenvectors(:,k));
% % %     if abs(sum_1n_k) < 10e-7
% % %         sum_1n_k = 0;
% % %     end
% % %     sum_n_k = sum(eigenvectors(:,n).*eigenvectors(:,k));
% % %     if abs(sum_n_k) < 10e-7
% % %         sum_n_k = 0;
% % %     end
% % %     E11010_vectr1(k) = 1/3 * qn.'*((f2000*p11+f1100*p12)*h_k01010(1,k)+(f0200*p12+f1100*p11)*h_k01010(2,k))*sum_1n_k;
% % %     E11010_vectr2(k) = 1/3 * qn.'*((f2000*conj(p11)+f1100*conj(p12))*h_k10010(1,k)+(f0200*conj(p12)+f1100*conj(p11))*h_k10010(2,k))*sum_1n_k;
% % %     E11010_vectr3(k) = 1/3 * qn.'*(f1010*h_k11000(1,k)+f0110*h_k11000(2,k))*sum_n_k;
% % % end
% % % E11010 = sum(E11010_vectr1 + E11010_vectr2 + E11010_vectr3);
E11010 = 0;

% % % E11001_vectr1 = zeros(1,N);
% % % E11001_vectr2 = zeros(1,N);
% % % E11001_vectr3 = zeros(1,N);
% % % for k=1:N
% % %     sum_1n_k = sum(eigenvectors(:,1).*eigenvectors(:,n).*eigenvectors(:,k));
% % %     if abs(sum_1n_k) < 10e-7
% % %         sum_1n_k = 0;
% % %     end
% % %     sum_n_k = sum(eigenvectors(:,n).*eigenvectors(:,k));
% % %     if abs(sum_n_k) < 10e-7
% % %         sum_n_k = 0;
% % %     end
% % %     E11001_vectr1(k) = 1/3 * qn.'*((f2000*p11+f1100*p12)*h_k01001(1,k)+(f0200*p12+f1100*p11)*h_k01001(2,k))*sum_1n_k;
% % %     E11001_vectr2(k) = 1/3 * qn.'*((f2000*conj(p11)+f1100*conj(p12))*h_k10001(1,k)+(f0200*conj(p12)+f1100*conj(p11))*h_k10001(2,k))*sum_1n_k;
% % %     E11001_vectr3(k) = 1/3 * qn.'*(f1001*h_k11000(1,k)+f0101*h_k11000(2,k))*sum_n_k;
% % % end
% % % E11001 = sum(E11001_vectr1 + E11001_vectr2 + E11001_vectr3);
E11001 = 0;

E00210_vectr1 = zeros(1,N);
E00210_vectr2 = zeros(1,N);
for k=1:N
    sum_n2_k = sum(eigenvectors(:,n).^2.*eigenvectors(:,k));
    if abs(sum_n2_k) < 10e-7
        sum_n2_k = 0;
    end
    sum_n_k = sum(eigenvectors(:,n).*eigenvectors(:,k));
    if abs(sum_n_k) < 10e-7
        sum_n_k = 0;
    end
    E00210_vectr1(k) = 1/3 * qn.'*((f2000*pn1+f1100*pn2)*h_k00110(1,k)+(f0200*pn2+f1100*pn1)*h_k00110(2,k))*sum_n2_k;
    E00210_vectr2(k) = 1/3 * qn.'*(f1010*h_k00200(1,k)+f0110*h_k00200(2,k))*sum_n_k;
end
% E00210 = sum(E00210_vectr1 + E00210_vectr2);
E00210 = sum(E00210_vectr1) + 1/3 * qn.'*(f1010*h_k00200(1,n)+f0110*h_k00200(2,n));

E00201_vectr1 = zeros(1,N);
E00201_vectr2 = zeros(1,N);
for k=1:N
    sum_n2_k = sum(eigenvectors(:,n).^2.*eigenvectors(:,k));
    if abs(sum_n2_k) < 10e-7
        sum_n2_k = 0;
    end
    sum_n_k = sum(eigenvectors(:,n).*eigenvectors(:,k));
    if abs(sum_n_k) < 10e-7
        sum_n_k = 0;
    end
    E00201_vectr1(k) = 1/3 * qn.'*((f2000*pn1+f1100*pn2)*h_k00101(1,k)+(f0200*pn2+f1100*pn1)*h_k00101(2,k))*sum_n2_k;
    E00201_vectr2(k) = 1/3 * qn.'*(f1001*h_k00200(1,k)+f0101*h_k00200(2,k))*sum_n_k;
end
% E00201 = sum(E00201_vectr1 + E00201_vectr2);
E00201 = sum(E00201_vectr1) + 1/3 * qn.'*(f1001*h_k00200(1,n)+f0101*h_k00200(2,n));

B21000 = C21000 + 3/2*(D21000 + E21000);
B10200 = C10200 + 3/2*(D10200 + E10200);
B11100 = C11100 + 3/2*(D11100 + E11100);
B00300 = C00300 + 3/2*(D00300 + E00300);

B00210 = C00210 + 3/2*(D00210 + E00210);
B00201 = C00201 + 3/2*(D00201 + E00201);
B10110 = C10110 + 3/2*(D10110 + E10110);
B10101 = C10101 + 3/2*(D10101 + E10101);
B11010 = C11010 + 3/2*(D11010 + E11010);
B11001 = C11001 + 3/2*(D11001 + E11001);

%% 非特殊网络
format long
% B00200
if B00200~=0
    % format rat
    alpha11 = real(B10010);
    alpha12 = real(B10001) ;
    alpha21 = B00110;
    alpha22 = B00101;

    kappa30 = real(B21000);
    kappa12 = real(B10200);
    kappa21 = B11100;
    kappa03 = B00300;
    kappa021 = B00200;
    kappa022 = B00210;
    kappa023 = B00201;
    
%     format rat
    
    [alpha11 alpha12 alpha21 alpha22;
        kappa30 kappa12 kappa21 kappa03;
           kappa021 kappa022 kappa023 0]
    kappa30*kappa03
else
    disp('****************************** special case ******************************')

    % format rat
    alpha11 = real(B10010);
    alpha12 = real(B10001) ;
    alpha21 = B00110;
    alpha22 = B00101;

    kappa30 = real(B21000);
    kappa12 = real(B10200);
    kappa21 = B11100;
    kappa03 = B00300;

    kappa021 = B00200;
    kappa022 = B00210;
    kappa023 = B00201;

    % format rat
    % [alpha11 alpha12 alpha21 alpha22;
    %     kappa30 kappa12 kappa21 kappa03]

    %%%

    % H
    res1 = -alpha11/alpha12;
    % T
    res2 = -alpha21/alpha22;
    % T1
    % (alpha11-kappa12/kappa22*alpha21)/(kappa12/kappa22*kappa21 - kappa11)
    % (alpha12-kappa12/kappa22*alpha22)/(kappa12/kappa22*kappa21 - kappa11)

    res3 = - (alpha11-kappa12/kappa03*alpha21)/(kappa12/kappa03*kappa21 - kappa30) / ((alpha12-kappa12/kappa03*alpha22)/(kappa12/kappa03*kappa21 - kappa30));
    % T2
    % -(alpha21/kappa22 + kappa21/(kappa12*kappa21-kappa22*kappa11)*(alpha11-kappa12/kappa22*alpha21))
    % -(alpha22/kappa22 + kappa21/(kappa12*kappa21-kappa22*kappa11)*(alpha12-kappa12/kappa22*alpha22))

    res4 = - (alpha21/kappa03 + kappa21/(kappa12*kappa21-kappa03*kappa30)*(alpha11-kappa12/kappa03*alpha21)) / ((alpha22/kappa03 + kappa21/(kappa12*kappa21-kappa03*kappa30)*(alpha12-kappa12/kappa03*alpha22)));

    % [res1 res2 res3 res4]

    [alpha11 alpha12 alpha21 alpha22;
        kappa30 kappa12 kappa21 kappa03;
           res1 res2 res3 res4]
    roundn([alpha11 alpha12 alpha21 alpha22;
        kappa30 kappa12 kappa21 kappa03;
           res1 res2 res3 res4],-6)
    kappa30*kappa03

end
format long
            
save([results_path 'THBNormalForm'],'a_star','c_star','du','dv','alpha11','alpha12','alpha21','alpha22','kappa30','kappa12','kappa21','kappa03','kappa021','kappa022','kappa023','n')

format short



