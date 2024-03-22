%****************Section 6***********************
%Formation Control for 3 agents using Receding Horizon method
%  without Trajectory tracking
clear all
close all
clc

n = 3; %Dimensional Plane
N = 3; %number of agents
M = 2; %number of edges
ns = 12; % number of states
ni =4;  %number of inputs
tf = 20; % t final
dt = 0.1; %steps
eta =4;
Np = 0;

step = 100;

A11i =zeros(2*n);
A11 = kron(eye(N),A11i);
A12i =eye(2*n);
A12 = kron(eye(N),A12i);
A21i = [zeros(n),[0, -9.81, 0; 9.81, 0, 0; 0, 0, 0];zeros(n,2*n)];
A21 =  kron(eye(N),A21i);
A22i = zeros(2*n);
A22 = kron(eye(N),A22i);

%A  Matrix for whole system
A = [A11, zeros(2*n*N,1) , A12;
    zeros(1,2*n*N), 0, zeros(1,2*n*N);
    A21, zeros(2*n*N,1), A22];
F = eye(ns*N+1)+(dt*A) + (dt^2/2)*A^2 + (dt^3/6)*A^3;

m = 0.64;

B1i = [zeros(ni);
    0 0 1 0;
    0 0 0 1/m];
B2i = zeros(1,ni);
B3i = [zeros(n,ni);eye(n),zeros(n,1)];

% B matrix components for N agents
B1 = kron(eye(N), B1i);
B2 = zeros(1, N*ni);
B3 = kron(eye(N), B3i);

B = [B1;B2;B3]; %B for N agents

% individual B for
Bb(:,:,1) = B(:, 1:ni);
Bb(:,:,2) = B(:, ni+1:2*ni);
Bb(:,:,3) = B(:, (2*ni)+1:3*ni);

for i = 1:N
    Gb(:,:,i) = dt*Bb(:,:,i) + (dt^2)/2*A*Bb(:,:,i) + (dt^3)/6*A^2*(Bb(:,:,i)).^2 + (dt^4)/24*A^3*(Bb(:,:,i)).^3;
end

mui = [1 ,1;1,0;0,1];

Di = [-1,-1;1,0;0,1];
D= kron(Di,eye(2*n)); %incidence matrix for 3 agents

d12 = [-1.2;-0.4;0];
d13 = [0;-0.8;0];

d =[0.015; 0.055;0.05;0;0;0;0.055;0.055;0.055;0;0;0]; %formation parameter for 3 agents

for i = 1:N
    %Q matrix for all time
    Wtemp(:,:,i) = diag(mui(i,:));
    W(:,:,i) = kron(Wtemp(:,:,i),eye(2*n));

    L(:,:,i) = D*W(:,:,i)*D';
    Q(:,:,i) = [L(:,:,i), -D*W(:,:,i)*d, zeros(2*N*n);
        (-D*W(:,:,i)*d)', d'*W(:,:,i)*d, zeros(1,2*N*n);
        zeros(2*N*n), zeros(2*N*n,1), L(:,:,i)];
    %Terminal weights
    WTtemp(:,:,i) = eta*diag(mui(i,:));
    WT(:,:,i) = kron(WTtemp(:,:,i),eye(2*n));

    LT(:,:,i) = D*WT(:,:,i)*D';
    QT(:,:,i) = [LT(:,:,i), -D*WT(:,:,i)*d, zeros(2*N*n);
        (-D*WT(:,:,i)*d)', d'*WT(:,:,i)*d, zeros(1,2*N*n);
        zeros(2*N*n), zeros(2*N*n,1), LT(:,:,i)];
    Pb(:,:,i,(tf/dt)) = QT(:,:,i);

end
% weight on input
Rb = dt*eye(4);

aN = 1;
aE = 1;

p1 = [0;2;1;0;0;0];
p2 = [1;3;1;0;0;0];
p3 = [1;1;1;0;0;0];
po = [p1;p2;p3];
% Velocity vectors (linear and Angular)
v1 =[0;0;0.1;0;0;0];
v2 =[0;0;0.1;0;0;0];
v3 =[0;0;0.1;0;0;0];

x(:,1) = [p1; p2; p3; v1;v2;v3];
xb(:,1) = [p1; p2; p3;1; v1;v2;v3];

Ap1 = [A11i,A12i;
    A21i,A22i];
A11t = kron(eye(M),A11i);
A12t = kron(eye(M),A12i);
A21t = kron(eye(M),A21i);
A22t = kron(eye(M), A22i);

At = [A11t, A12t;
    A21t, A22t];
Ft = eye(4*n*M) + dt*At;
B1t = kron(eye(M),B1i);
B3t = kron(eye(M),B3i);
Bt = [B1t;B3t];

Fpinv = [F(1:18, 1:18), F(1:18,20:37);F(20:37, 1:18),F(20:37, 20:37)]
Gpinv =[Gb(1:18,1:12);Gb(20:37,1:12)];

Qt = dt*eye(4*M*n);
Qft   = eta*Qt;
Rt =eye(2*ni);
Gt    = dt*Bt + (dt^2)/2*At*Bt;
Phi   = -Di';
Phia  = kron(Phi,eye(4));
Phiai = pinv(Phia);

Iter    = 200;
In      = eye(2*M*n);
u0      = zeros(2*M*n,1);
eps     = 1e-3;
udata   = [];
alpha = 0.1;


for idx=1:Np+1
    if idx==1
        z(:,1) = [p1(1,1)-p2(1,1)-d(1,1);p1(2,1)-p2(2,1)-d(2,1);p1(3,1)-p2(3,1)-d(3,1);p1(4,1)-p2(4,1);p1(5,1)-p2(5,1);p1(6,1)-p2(6,1);p1(1,1)-p3(1,1)-d(7,1);p1(2,1)-p3(2,1)-d(8,1);p1(3,1)-p3(3,1)-d(9,1);p1(4,1)-p3(4,1);p1(5,1)-p3(5,1);p1(6,1)-p3(6,1);v1(1,1)-v2(1,1);v1(2,1)-v2(2,1);v1(3,1)-v2(3,1);v1(4,1)-v2(4,1);v1(5,1)-v2(5,1);v1(6,1)-v2(6,1);v1(1,1)-v3(1,1);v1(2,1)-v3(2,1);v1(3,1)-v3(3,1);v1(4,1)-v3(4,1);v1(5,1)-v3(5,1);v1(6,1)-v3(6,1)];
        %have to initialize the position and velocity vector
        xpinv(:,idx) = x(:,1);
        xhat(:,idx)  = x(:,1);
    else
        idxn = ((idx-1) * 10) + 1;
        z(:,idxn) = [xf{idx-1}(1,idxn)- xf{idx-1}(7,idxn); xf{idx-1}(2,idxn)- xf{idx-1}(8,idxn);xf{idx-1}(3,idxn)- xf{idx-1}(9,idxn);xf{idx-1}(4,idxn)- xf{idx-1}(10,idxn);xf{idx-1}(5,idxn)- xf{idx-1}(11,idxn); xf{idx-1}(6,idxn)- xf{idx-1}(12,idxn); xf{idx-1}(1,idxn)- xf{idx-1}(13,idxn); xf{idx-1}(2,idxn)- xf{idx-1}(14,idxn);xf{idx-1}(3,idxn)- xf{idx-1}(15,idxn);xf{idx-1}(4,idxn)- xf{idx-1}(16,idxn),xf{idx-1}(5,idxn)- xf{idx-1}(17,idxn),xf{idx-1}(6,idxn)- xf{idx-1}(18,idxn),xf{idx-1}(19,idxn)- xf{idx-1}(25,idxn),xf{idx-1}(20,idxn)- xf{idx-1}(26,idxn),xf{idx-1}(21,idxn)- xf{idx-1}(27,idxn),xf{idx-1}(22,idxn)- xf{idx-1}(28,idxn),xf{idx-1}(23,idxn)- xf{idx-1}(29,idxn);xf{idx-1}(24,idxn)- xf{idx-1}(30,idxn);xf{idx-1}(19,idxn)- xf{idx-1}(31,idxn);xf{idx-1}(20,idxn)- xf{idx-1}(32,idxn);xf{idx-1}(21,idxn)- xf{idx-1}(33,idxn);xf{idx-1}(22,idxn)- xf{idx-1}(34,idxn);xf{idx-1}(23,idxn)- xf{idx-1}(35,idxn);xf{idx-1}(24,idxn)- xf{idx-1}(36,idxn)];
        %if reced the horizon then x(l) be the initial state
        %x(:,idxn)    = xd{idx-1}(:,idxn);
        xpinv(:,idxn)= xe{idx-1}(:,idxn);
        xhat(:,idxn) = xf{idx-1}(:,idxn);
    end


    T = (idx-1)*10 + step;
    Pt(:,:,T+1) = Qft; %initial Riccati

    for k = 1:T
        for j = T-1:-1:1
            Pt(:,:,j) = Qt + Ft'*Pt(:,:,j+1)*Ft - Ft'*Pt(:,:,j+1)*Gt*inv(Rt+Gt'*Pt(:,:,j+1)*Gt)*Gt'*Pt(:,:,j+1)*Ft;
        end
        Kt(:,:,k) = -inv(Rt+Gt'*Pt(:,:,k+1)*Gt)*Gt'*Pt(:,:,k+1)*Ft;
        a(:,k)    = Kt(:,:,k) * z(:,k);
        z(:,k+1)  = Ft * z(:,k) + Gt * a(:,k);

        for v = 1:N
            upinv(:,k) = Phiai * a(:,k);
        end

        xpinv(:,k+1) =  Fpinv * xpinv(:,k)+Gpinv*upinv(:,k);
        %iterative method to get back uhat from stored a
        %tic;
        t = 1;
        if k>1
            u0 = usol;
        end
        while (t <= Iter)
            u0_temp = (In - 2*alpha*Phia'*Phia)*u0(:,t) + 2*alpha*Phia'*a(:,k);
            u0 = [u0 u0_temp];
            t = t+1;
        end
        udata = [udata,u0];
        %toc;
        usol      = u0(:,end);
        uhat(:,k) = usol;
        xhat(:,k+1) = Fpinv * xhat(:,k) + Gpinv*uhat(:,k);
        %calculate the receding stability equation
        term1 = (Ft + Gt*Kt(:,:,end))' * Qft * (Ft + Gt*Kt(:,:,end)) - Qft;
        term2 = -Qt - Kt(:,:,end)'*Rt*Kt(:,:,end);
    end
    xd{idx} = x;
    xe{idx} = xpinv;
    xf{idx} = xhat;
    zf{idx} = z;
end
%

figure('Name', 'formation', 'NumberTitle', 'off')
plot3(xhat(1,:),xhat(2,:),xhat(3,:),'-b',xhat(7,:),xhat(8,:),xhat(9,:),'-r',...
    xhat(13,:),xhat(14,:),xhat(15,:),'-m','linewidth',3.5)
hold on

plot3(xpinv(1,:),xpinv(2,:),xpinv(3,:),':b',xpinv(7,:),xpinv(8,:),xpinv(9,:),':r',...
    xpinv(13,:),xpinv(14,:),xpinv(15,:),':m','linewidth',3.5)
hold on

%Plot to show the formation for 3 agent system.

legend('$p^{1N}$','$p^{2N}$','$p^{3N}$','$p^{4N}$','$p^{1D}$','$p^{2D}$','$p^{3D}$','$p^{4D}$','fontweight','bold','fontsize',14,'interpreter','latex')
xlabel('x-axis','fontweight','bold','fontsize',12)
ylabel('y-axis','fontweight','bold','fontsize',12)
zlabel('z-axis','fontweight','bold','fontsize',12)
title('Motion trajectories of each agent','fontweight','bold','fontsize',12)
grid on
set(gca,'color',[0.9,0.9,0.9]);
