%VI-based method
clc; clear all; close all; 
%% load data from system identification file
A = [0 1;
     -2*rand(1) 5*rand(1)];
B = [0;
     1];
D = [0;
     rand(1)];

N=100;
x = zeros(2,N);
x(:,1) = [10;0];

K0 = place(A,B,[0.9+0.1i,0.9-0.1i]);
Q = eye(2);
R = 1;
[Kstar,Pstar,~] = dlqr(A,B,Q,R);
for i =1 : N
    u(:,i) =  -K0*x(:,i) + rand(1);
    x(:,i+1) = A*x(:,i) + B*u(:,i) + D; 
end
figure(1);
plot(1:length(x), x(:,:));


rho = ones(1,N+1);
Nsample_total = length(x);
Number_sample_need = 20;
Rand_samples = randsample(Nsample_total-1,Number_sample_need);

%% learning phase 1: compute K P
K=K0;P_1=zeros(2);index=1;errorK=[];errorP=[];INDEX=[];
while(1) 
Psi=[];Phi=[];YY=[];test=[];TEST=[];

%random select j
%rng(0,'twister');
for i=1:length(Rand_samples)
    j= Rand_samples(i);

Gamma_Kx=(K*x(:,j))^2;
Gamma_u=u(:,j)^2;
Gamma_xu=transpose(kron(x(:,j),u(:,j)));
Gamma_xx=transpose(kron(x(:,j),x(:,j)));
temp2=kron(x(:,j)',x(:,j)')-kron(x(:,j+1)',x(:,j+1)'); %sign is opposite to the definition 

%dim(delta_xx) = 2*(2+1)/2 = 3
delta_xx=temp2(1,[1:2,4]);

Gamma_rhox=transpose(kron(rho(:,j),x(:,j)));
Gamma_rhou=transpose(kron(rho(:,j),u(:,j)));
Gamma_rho=rho(:,j)^2;

Psi=[Psi;[Gamma_u-Gamma_Kx,2*(Gamma_xu+Gamma_xx*(kron(eye(2),K'))),delta_xx,2*Gamma_rhox,2*Gamma_rhou,Gamma_rho]];
%Psi=[Psi;[Gamma_u-Gamma_Kx,2*(Gamma_xu+Gamma_xx*(kron(eye(2),K'))),delta_xx]];
Phi=[Phi;x(:,j)'*Q*x(:,j)+x(:,j)'*K'*R*K*x(:,j)];
end

%rank(X) should be 15
rank(Psi);
pp=(Psi\Phi);
% pp=pinv(X)*Y;%dim(pp)

BPB = pp(1);
BPA = pp (2:3)';
K=inv((R+BPB))*BPA;




%temp5 = half Pbar
l=4;
for m = 1:2
    for n=m:2
        temp5(m,n)=pp(l)/2;
        l=l+1;
    end
end
P=temp5+temp5';

APD = pp(7:8);
% APD_star = A'*Pstar*D;
BPD = pp(9);
% BPD_star = B'*Pstar*D;
DPD = pp(10);
% DPD_star = D'*Pstar*D;

INDEX=[INDEX;[index]];
index;
errorK=[errorK;[norm(K-Kstar)]];
errorP=[errorP;[norm(P-Pstar)]];
if ((norm(P-P_1)<10e-5)||(index>50))    
    %disp('meet the converngent requirement and the norm of the P_k-P_{k-1} is')
    vpa(K,5);
    vpa(P,5) ; 
    indexf = index;
    break
end   
P_1=P;index=index+1;  
end
% K learning phase ends

%% output regulation
C = [1 0];
Null_C = null(C);
%% i = 1
Y1 = Null_C(:,1);
x1 = x - Y1*1*ones(1,length(x));
Psi=[];Phi=[];YY=[];
for i=1:length(Rand_samples)
    j= Rand_samples(i);

Gamma_Kx1=(K*x1(:,j))^2;
Gamma_u=u(:,j)^2;
Gamma_x1u=transpose(kron(x1(:,j),u(:,j)));
Gamma_x1x1=transpose(kron(x1(:,j),x1(:,j)));
temp2=kron(x1(:,j)',x1(:,j)')-kron(x1(:,j+1)',x1(:,j+1)'); %sign is opposite to the definition 

%dim(delta_x1x1) = 2*(2+1)/2 = 3
delta_x1x1=temp2(1,[1:2,4]);

Gamma_rhox1=transpose(kron(rho(:,j),x1(:,j)));
Gamma_rhou=transpose(kron(rho(:,j),u(:,j)));
Gamma_rho=rho(:,j)^2;

Psi=[Psi;[Gamma_u-Gamma_Kx1,2*(Gamma_x1u+Gamma_x1x1*(kron(eye(2),K'))),delta_x1x1,2*Gamma_rhox1,2*Gamma_rhou,Gamma_rho]];
%Psi=[Psi;[Gamma_u-Gamma_Kx1,2*(Gamma_x1u+Gamma_x1x1*(kron(eye(2),K'))),delta_x1x1]];
Phi=[Phi;x1(:,j)'*Q*x1(:,j)+x1(:,j)'*K'*R*K*x1(:,j)];
end
rank(Psi);
pp=(Psi\Phi);
%% check APE1
APE1 = pp(7:8);
APS1 = APD - APE1;
% E1 = D-Y1+A*Y1;
% APE1_star = A'*Pstar*E1;

K
Kstar
%%
APB = BPA';
AA = [-APS1 APB];
BB = -APD;
alphaU = inv(AA)*BB;
alpha1 = alphaU(1);
U = alphaU(2);
X = alpha1*Y1;
L = U+K*X

AA = [A-eye(2) B; C 0];
%AA = [A-eye(3) B; C 0];
BB = [-D ; 0];
XU = inv(AA)*BB;
Xstar = XU(1:2);
Ustar = XU(3);
Lstar = Ustar + Kstar*Xstar