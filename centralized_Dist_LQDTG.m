
%****************Section 1***********************
%Formation Control for 3 agents without Trajectory tracking
%This section has two problems 1. Corresponding to Nash
%Equilibrium and 2. The distributed problem Inorder to solve
%Problem 2 First part of section must be run as problem 2
%takes values from problem 1.
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

%**************************start***********************
%************************Problem1**********************

% individual elements for A matrix
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
d =[0.5; -0.5;2;0;0;0;0;-0.5;-1;-1;0;0];
for i = 1:N
    %Q matrix for all time
    Wtemp(:,:,i) = diag(mui(i,:));
    W(:,:,i) = kron(Wtemp(:,:,i),eye(2*n));

    L(:,:,i) = D*W(:,:,i)*D';
    Q(:,:,i) = [L(:,:,i), -D*W(:,:,i)*d, zeros(2*N*n);
        (-D*W(:,:,i)*d)', d'*W(:,:,i)*d, zeros(1,2*N*n);
        zeros(2*N*n), zeros(2*N*n,1), L(:,:,i)];
    Qn(:,:,i) = [Q(1:18,1:18,i), Q(20:37,20:37,i);Q(20:37,1:18,i), Q(20:37,20:37. ...
        ,i)]
    %Terminal weights
    WTtemp(:,:,i) = eta*diag(mui(i,:));
    WT(:,:,i) = kron(WTtemp(:,:,i),eye(2*n));

    LT(:,:,i) = D*WT(:,:,i)*D';
    QT(:,:,i) = [LT(:,:,i), -D*WT(:,:,i)*d, zeros(2*N*n);
        (-D*WT(:,:,i)*d)', d'*WT(:,:,i)*d, zeros(1,2*N*n);
        zeros(2*N*n), zeros(2*N*n,1), LT(:,:,i)];
    Pb(:,:,i,(tf/dt)) = QT(:,:,i);
    QbT(:,:,i) = [QT(1:18,1:18,i), QT(20:37,20:37,i);QT(20:37,1:18,i), QT(20:37,20:37,i)]

end
% weight on input
Rb = dt*eye(4);




% linear and angular Position Vector
p1 =[0;0;0;0;0;0];
p2=[-4;5;0;0;0;0];
p3=[4;5;0;0;0;0];
aN = 1;
aE = 1;
% Velocity vectors (linear and Angular)
v1 =[0;0;0.1;0;0;0];
v2 =[0;0;0.1;0;0;0];
v3 =[0;0;0.1;0;0;0];

x(:,1) = [p1; p2; p3; v1;v2;v3];
xb(:,1) = [p1; p2; p3;1; v1;v2;v3];

%%
%backward computation to solve single ricatti equation
for k=1:100

    for j = (tf/dt)-1:-1:1
        for v = 1:N
            S(:,:,v) = Gb(:,:,v)*inv(Rb)*Gb(:,:,v)';
        end
        Lam(:,:,k) = eye((ns*N)+1);
        for v = 1:N
            Lam(:,:,k) = Lam(:,:,k) + S(:,:,v)*Pb(:,:,v,j+1);
        end
        for v = 1:N
            Pb(:,:,v,j) = Q(:,:,v) + (F'*Pb(:,:,v,j+1) * inv(Lam(:,:,k))*F);
        end
    end
    for v = 1:N
        ub(:,:,v,k) = -inv(Rb)*Gb(:,:,v)' * Pb(:,:,v,k+1) * inv(Lam(:,:,k)) * F * xb(:,k);
    end
    uball(:,k) = [ub(:,:,1,k);ub(:,:,2,k);ub(:,:,3,k)];
    xb(:,k+1)  = F * xb(:,k) +  Gb(:,:,1)*ub(:,:,1,k) + Gb(:,:,2)*ub(:,:,2,k) + Gb(:,:,3)*ub(:,:,3,k);
    % x(:,k+1) = inv(Lam(:,:,k)) * F * x(:,k);
    for l =1:N
        Jn(:,l,k) =  xb(:,k)'*Q(:,:,l)*xb(:,k) + uball(1:4,k)'*Rb*uball(1:4,k)+ uball(5:8,k)'*Rb*uball(5:8,k) + uball(9:12,k)'*Rb*uball(9:12,k);
    end

end

Jn1 = xb(:,100)'*QT(:,:,1)*xb(:,100) + sum(Jn(:,1,:));
Jn2 = xb(:,100)'*QT(:,:,2)*xb(:,100) + sum(Jn(:,2,:));
Jn3 = xb(:,100)'*QT(:,:,3)*xb(:,100) + sum(Jn(:,3,:));
Jnall = Jn1+Jn2+Jn3;




%
%Plot to show the formation control using Nash Equilibrium Coupled problem
figure('Name', 'formation', 'NumberTitle', 'off')
plot3(xb(1,:),xb(2,:),xb(3,:),'-b',xb(7,:),xb(8,:),xb(9,:),'-r',...
    xb(13,:),xb(14,:),xb(15,:),'-m','linewidth',3.5)
hold on
legend('$p^{1N}$','$p^{2N}$','$p^{3N}$','$p^{4N}$','$p^{1D}$','$p^{2D}$','$p^{3D}$','$p^{4D}$','fontweight','bold','fontsize',14,'interpreter','latex')
xlabel('x-axis','fontweight','bold','fontsize',12)
ylabel('y-axis','fontweight','bold','fontsize',12)
zlabel('z-axis','fontweight','bold','fontsize',12)
title('Motion trajectories of each agent','fontweight','bold','fontsize',12)
grid on
set(gca,'color',[0.9,0.9,0.9]);
%%
%*******************start************************
%*****************problem 2 Distributed*********

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
% Bt    = [zeros(2*n*M);eye(2*n*M)];

Qt = dt*eye(4*M*n);
Qft   = eta*Qt;
Rt =eye(2*ni);
Gt    = dt*Bt + (dt^2)/2*At*Bt;
% Dt= kron(Di,eye(n));
Phi   = -Di';
Phia  = kron(Phi,eye(4));
Phiai = pinv(Phia);

Iter    = 200;
In      = eye(2*M*n);
u0      = zeros(2*M*n,1);
eps     = 1e-3;
udata   = [];
alpha = 0.1;

z(:,1) = [p1(1,1)-p2(1,1)-d(1,1);p1(2,1)-p2(2,1)-d(2,1);p1(3,1)-p2(3,1)-d(3,1);p1(4,1)-p2(4,1);p1(5,1)-p2(5,1);p1(6,1)-p2(6,1);p1(1,1)-p3(1,1)-d(7,1);p1(2,1)-p3(2,1)-d(8,1);p1(3,1)-p3(3,1)-d(9,1);p1(4,1)-p3(4,1);p1(5,1)-p3(5,1);p1(6,1)-p3(6,1);v1(1,1)-v2(1,1);v1(2,1)-v2(2,1);v1(3,1)-v2(3,1);v1(4,1)-v2(4,1);v1(5,1)-v2(5,1);v1(6,1)-v2(6,1);v1(1,1)-v3(1,1);v1(2,1)-v3(2,1);v1(3,1)-v3(3,1);v1(4,1)-v3(4,1);v1(5,1)-v3(5,1);v1(6,1)-v3(6,1)];


for idx=1:Np+1
    if idx==1

        %have to initialize the position and velocity vector
        %x(:,idx)     = x(:,1);
        xpinv(:,idx) = x(:,1);
        xhat(:,idx)  = x(:,1);
    else
        idxn = ((idx-1) * 10) + 1;
        z(:,idxn) = [xf{idx-1}(1:2,idxn)- xf{idx-1}(3:4,idxn); xf{idx-1}(3:4,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(5:6,idxn)- xf{idx-1}(7:8,idxn);xf{idx-1}(1:2,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(9:10,idxn)- xf{idx-1}(11:12,idxn); xf{idx-1}(11:12,idxn)- xf{idx-1}(13:14,idxn); xf{idx-1}(13:14,idxn)- xf{idx-1}(15:16,idxn); xf{idx-1}(9:10,idxn)- xf{idx-1}(13:14,idxn)];

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
        %xpinv(:,k+1) = F * xpinv(:,k) +G(:,1:2)*upinv(1:2,k) + G(:,3:4)*upinv(3:4,k) + G(:,5:6)*upinv(5:6,k) + G(:,7:8)*upinv(7:8,k) ;

        xpinv(:,k+1) =  Fpinv * xpinv(:,k)+Gpinv*upinv(:,k);
        %iterative method to get back uhat from stored a
        %tic;
        t = 1;
        if k>1
            u0 = usol;
        end
        %             while (norm(a(:,k)-Phia*u0(:,t)) > eps) && (t <= Iter)%for receding horizon
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
        for l = 1:N
            Jt(:,l,k) =  xpinv(:,k)'*Qn(:,:,l)*xpinv(:,k) + upinv(1:4,k)'*Rb*upinv(1:4,k)+ upinv(5:8,k)'*Rb*upinv(5:8,k)+ upinv(9:12,k)'*Rb*upinv(9:12,k);
        end

    end
end
Jt1 = xpinv(:,100)'*QbT(:,:,1)*xpinv(:,100) + sum(Jt(:,1,:));
Jt2 = xpinv(:,100)'*QbT(:,:,2)*xpinv(:,100) + sum(Jt(:,2,:));
Jt3 = xpinv(:,100)'*QbT(:,:,3)*xpinv(:,100) + sum(Jt(:,3,:));
Jtall = Jt1+Jt2+Jt3;


%Plot to show the Formation control for distributed system
figure('Name', 'formation', 'NumberTitle', 'off')
plot3(xhat(1,:),xhat(2,:),xhat(3,:),'-b',xhat(7,:),xhat(8,:),xhat(9,:),'-r',...
    xhat(13,:),xhat(14,:),xhat(15,:),'-m','linewidth',3.5)
hold on

plot3(xpinv(1,:),xpinv(2,:),xpinv(3,:),':b',xpinv(7,:),xpinv(8,:),xpinv(9,:),':r',...
    xpinv(13,:),xpinv(14,:),xpinv(15,:),':m','linewidth',3.5)
hold on


plot3(xhat(1,1:1),xhat(2,1:1),xhat(3,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#2f54eb')
plot3(xhat(1,100:100),xhat(2,100:100),xhat(3,100:100),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#2f54eb')

plot3(xhat(7,1:1),xhat(8,1:1),xhat(9,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#f51146')
plot3(xhat(7,100:100),xhat(8,100:100),xhat(9,100:100),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#f51146')

plot3(xhat(13,1:1),xhat(14,1:1),xhat(15,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#384d3f')
plot3(xhat(13,100:100),xhat(14,100:100),xhat(15,100:100),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#384d3f')

legend('$p^{1N}$','$p^{2N}$','$p^{3N}$','$p^{4N}$','$p^{1D}$','$p^{2D}$','$p^{3D}$','$p^{4D}$','fontweight','bold','fontsize',14,'interpreter','latex')
xlabel('x-axis','fontweight','bold','fontsize',12)
ylabel('y-axis','fontweight','bold','fontsize',12)
zlabel('z-axis','fontweight','bold','fontsize',12)
title('Motion trajectories of each agent','fontweight','bold','fontsize',12)
grid on
set(gca,'color',[0.9,0.9,0.9]);
