                             %Project_Arbeit_Rohith
% Linear Quadratic Distributed-Time Game for Quadrocopter Formation Control

% There are totally 6 sections in the code divided according to the problem
%statement and based on number of agents. Each section can be run
%seperately to get the desired results

%Section 1: Helical Trajectory tracking with formation for 3 agents.
%Section 2: Helical Trajectory tracking with formation for 4 agents.
%Section 3: Infinity Trajectory tracking with formation for 4 agents.
%Section 4: Infinity Trajectory tracking with formation for 4 agents.
%Section 5: Formation Control for 3 agents without Trajectory tracking
%Section 6: Formation Control for 3 agents using Receding Horizon method
            % without Trajectory tracking

%%
                 %****************Section 1***********************
             %Helical Trajectory tracking with formation for 3 agents 
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

d =[1; 0.5;-1;0;0;0;1;-1;1;0;0;0]; %formation parameter for 3 agents


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

Qn = [Q(1:18,1:18), Q(20:37,20:37);Q(20:37,1:18), Q(20:37,20:37)]
% weight on input
Rb = dt*eye(4);

% linear and angular Position Vector
aN = 1;
aE = 1;
p1 = [1;4;1;0;0;0];
p2 = [2;4;1;0;0;0];
p3 = [0;4;1;0;0;0];

% Velocity vectors (linear and Angular)
v1 =[0;0;0.1;0;0;0];
v2 =[0;0;0.1;0;0;0];
v3 =[0;0;0.1;0;0;0];

x(:,1) = [p1; p2; p3; v1;v2;v3];
xb(:,1) = [p1; p2; p3;1; v1;v2;v3]; %State Vector

Ap1 = [A11i,A12i;
    A21i,A22i];
A11t = kron(eye(M),A11i);
A12t = kron(eye(M),A12i);
A21t = kron(eye(M),A21i);
A22t = kron(eye(M), A22i);

%A matrix for $N$ agents
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
% Dt= kron(Di,eye(n));
Phi   = -Di';
Phia  = kron(Phi,eye(4));
Phiai = pinv(Phia);

Iter    = 50;
In      = eye(2*M*n);
u0      = zeros(2*M*n,1);
eps     = 1e-3;
udata   = [];
alpha = 0.1;

z(:,1) = [p1(1,1)-p2(1,1)-d(1,1);p1(2,1)-p2(2,1)-d(2,1);p1(3,1)-p2(3,1)-d(3,1);p1(4,1)-p2(4,1);p1(5,1)-p2(5,1);p1(6,1)-p2(6,1);p1(1,1)-p3(1,1)-d(7,1);p1(2,1)-p3(2,1)-d(8,1);p1(3,1)-p3(3,1)-d(9,1);p1(4,1)-p3(4,1);p1(5,1)-p3(5,1);p1(6,1)-p3(6,1);v1(1,1)-v2(1,1);v1(2,1)-v2(2,1);v1(3,1)-v2(3,1);v1(4,1)-v2(4,1);v1(5,1)-v2(5,1);v1(6,1)-v2(6,1);v1(1,1)-v3(1,1);v1(2,1)-v3(2,1);v1(3,1)-v3(3,1);v1(4,1)-v3(4,1);v1(5,1)-v3(5,1);v1(6,1)-v3(6,1)];


for idx=1:Np+1
    if idx==1
        %have to initialize the position and velocity vector
        xpinv(:,idx) = x(:,1);
        xhat(:,idx)  = x(:,1);
    else
        idxn = ((idx-1) * 10) + 1;
        z(:,idxn) = [xf{idx-1}(1:2,idxn)- xf{idx-1}(3:4,idxn); xf{idx-1}(3:4,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(5:6,idxn)- xf{idx-1}(7:8,idxn);xf{idx-1}(1:2,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(9:10,idxn)- xf{idx-1}(11:12,idxn); xf{idx-1}(11:12,idxn)- xf{idx-1}(13:14,idxn); xf{idx-1}(13:14,idxn)- xf{idx-1}(15:16,idxn); xf{idx-1}(9:10,idxn)- xf{idx-1}(13:14,idxn)];
        %if reced the horizon then x(l) be the initial state
        xpinv(:,idxn)= xe{idx-1}(:,idxn);
        xhat(:,idxn) = xf{idx-1}(:,idxn);
    end
    T = (idx-1)*10 + step;
    Pt(:,:,201) = Qft; %initial Riccati
    aN = 5;
    aE =5;
    %Discrete
    for k=1:200%horizon time
        p1d     = [aN*cos((k*pi)/20);aE*sin((k*pi)/20);1*dt*k;0;0;0];
        p2d     = [aN*cos((k*pi)/20);aE*sin((k*pi)/20);1*dt*k;0;0;0];
        p3d     = [aN*cos((k*pi)/20);aE*sin((k*pi)/20);1*dt*k;0;0;0];
        v1d   = [0;0;0.1;0;0;0];
        v2d = [0;0;0.1;0;0;0];
        v3d = [0;0;0.1;0;0;0];
        r(:,k)  = [p1d;p2d;p3d;v1d;v2d;v3d];
        e(:,1)  = abs(r(:,1) - x(:,1));

        for j = T-1:-1:1
            Pt(:,:,j) = Qt + Ft'*Pt(:,:,j+1)*Ft - Ft'*Pt(:,:,j+1)*Gt*inv(Rt+Gt'*Pt(:,:,j+1)*Gt)*Gt'*Pt(:,:,j+1)*Ft;
        end
        Kt(:,:,k) = -inv(Rt+Gt'*Pt(:,:,k+1)*Gt)*Gt'*Pt(:,:,k+1)*Ft;
        a(:,k)    = Kt(:,:,k) * z(:,k);
        z(:,k+1)  = Ft * z(:,k) + Gt * a(:,k);

        for v = 1:N
            upinv(:,k) = Phiai * a(:,k);
        end

        e(:,k+1) =  Fpinv * e(:,k)+Gpinv*upinv(:,k);
        xpinv(:,k+1) = e(:,k)+r(:,k);
        %iterative method to get back uhat from stored a
        %tic;
        t = 1;
        if k>1
            u0 = usol;
        end
        %
        while (t <= Iter)
            u0_temp = (In - 2*alpha*Phia'*Phia)*u0(:,t) + 2*alpha*Phia'*a(:,k);
            u0 = [u0 u0_temp];
            t = t+1;
        end
        udata = [udata,u0];
        %toc;
        usol      = u0(:,end);
        uhat(:,k) = usol;
        e(:,k+1) = Fpinv * e(:,k) + Gpinv*uhat(:,k);
        xhat(:,k+1) = e(:,k)+r(:,k);
        
    end

    end


%plot to get the formation and TRajectory tracking.
figure(1);
plot3(xhat(1,2:60),xhat(2,2:60),xhat(3,2:60),'-r',xhat(7,2:60),xhat(8,2:60),xhat(9,2:60),'-b',xhat(13,2:60),xhat(14,2:60),xhat(15,2:60),'-g','linewidth',1.5);
hold on
grid on
set(gca,'color',[0.9,0.9,0.9]);
title('Motion of Positions in Fomation','fontweight','bold')

plot3(xhat(1,2:2),xhat(2,2:2),xhat(3,2:2),'-o','Color','b','MarkerSize',9,...
    'MarkerFaceColor','#2f54eb')
plot3(xhat(1,60:60),xhat(2,60:60),xhat(3,60:60),'-o','Color','b','MarkerSize',9,...
    'MarkerFaceColor','#2f54eb')

plot3(xhat(7,2:2),xhat(8,2:2),xhat(9,2:2),'-o','Color','b','MarkerSize',9,...
    'MarkerFaceColor','#f51146')
plot3(xhat(7,60:60),xhat(8,60:60),xhat(9,60:60),'-o','Color','b','MarkerSize',9,...
    'MarkerFaceColor','#f51146')

plot3(xhat(13,2:2),xhat(14,2:2),xhat(15,2:2),'-o','Color','b','MarkerSize',9,...
    'MarkerFaceColor','#384d3f')
plot3(xhat(13,60:60),xhat(14,60:60),xhat(15,60:60),'-o','Color','b','MarkerSize',9,...
    'MarkerFaceColor','#384d3f')

legend('Agent 1','Agent 2','Agent 3','Location','Best')
xlabel('x-axis');
ylabel('y-axis');

% Plot to show the convergence of Distributed to the centralized solution
fig = figure(96);clf;
ax  = axes;
plot(ax,udata(9,:),'*b')
hold on
for i = 1:length(upinv(9,:))
    plot(ax,(Iter+1)*i,upinv(9,i),'xr','linewidth',5);
end
grid on
xlabel('Time/ iterations')
ylabel('Control input of Agent 3')
legend('Distributed approach','Centralized solution','fontsize',12,'location','northwest');
xlim([0 600])


%zoom window
% a2 = axes();
% a2.Position = [0.500 0.480 0.35 0.42]; % xlocation, ylocation, xsize, ysize
% plot(a2,udata(9,:),'*b');
% hold on
% for i = 1:100
%     plot(a2,(Iter+1)*i,upinv(1,i),'xr','linewidth',5);
%     hold on
% end


%%
                   %****************Section 3***********************
             %Helical Trajectory tracking with formation for 4 agents 

clear all
close all
clc

n = 3; %Dimensional Plane
N = 4; %number of agents
M = 4; %number of edges
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
Bb(:,:,4) = B(:, (3*ni)+1:4*ni);

for i = 1:N
    Gb(:,:,i) = dt*Bb(:,:,i) + (dt^2)/2*A*Bb(:,:,i) + (dt^3)/6*A^2*(Bb(:,:,i)).^2 + (dt^4)/24*A^3*(Bb(:,:,i)).^3;
end

mui = [-1 0 0 -1;
    1 -1 0 0;
    0 1 -1 1;
    0 0 1 0];

Di = [-1 0 0 -1;
    1 -1 0 0;
    0 1 -1 1;
    0 0 1 0];
D= kron(Di,eye(2*n)); %incidence matrix for 3 agents

d12 = [-1.2;-0.4;0];
d13 = [0;-0.8;0];
d =[1; 3;1;0;0;0;1;4;-1;0;0;0;-1; 1;-1;0;0;0;-1;1;1;0;0;0]; %formation parameter for 4 agents


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

% linear and angular Position Vector
p1 =[0;1;1;0;0;0];
p2=[0;0;1;0;0;0];
p3=[1;0;1;0;0;0];
p4 = [0;-1;1;0;0;0];
aN = 1;
aE = 1;

% Velocity vectors (linear and Angular)
v1 =[0;0;1;0;0;0];
v2 =[0;0;1;0;0;0];
v3 =[0;0;1;0;0;0];
v4 =[0;0;1;0;0;0];
x(:,1) = [p1; p2; p3;p4; v1;v2;v3;v4];
xb(:,1) = [p1; p2; p3;p4;1; v1;v2;v3;v4];

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

Fpinv = [F(1:24, 1:24), F(1:24,26:49);F(26:49, 1:24),F(26:49, 26:49)]
Gpinv =[Gb(1:24,1:16);Gb(26:49,1:16)];
% Bt    = [zeros(2*n*M);eye(2*n*M)];

Qt = dt*eye(4*M*n);
Qft   = eta*Qt;
Rt =eye(4*ni);
Gt    = dt*Bt + (dt^2)/2*At*Bt;


% Dt= kron(Di,eye(n));
Phi   = -Di';
Phia  = kron(Phi,eye(4));
Phiai = pinv(Phia);

Iter    = 200;
In      = eye(16);
u0      = zeros(4*M,1);
eps     = 1e-3;
udata   = [];
alpha = 0.1;

z(:,1) = [p1(1,1)-p2(1,1)-d(1,1);p1(2,1)-p2(2,1)-d(2,1);p1(3,1)-p2(3,1)-d(3,1);p1(4,1)-p2(4,1);p1(5,1)-p2(5,1);p1(6,1)-p2(6,1);p2(1,1)-p3(1,1)-d(7,1);p2(2,1)-p3(2,1)-d(8,1);p2(3,1)-p3(3,1)-d(9,1);p2(4,1)-p3(4,1);p2(5,1)-p3(5,1);p2(6,1)-p3(6,1);p1(1,1)-p3(1,1)-d(13,1);p1(2,1)-p3(2,1)-d(14,1);p1(3,1)-p3(3,1)-d(15,1);p1(4,1)-p3(4,1);p1(5,1)-p1(5,1);p1(6,1)-p3(6,1);p3(1,1)-p4(1,1)-d(13,1);p3(2,1)-p4(2,1)-d(14,1);p3(3,1)-p4(3,1)-d(15,1);p3(4,1)-p4(4,1);p3(5,1)-p4(5,1);p3(6,1)-p4(6,1);v1(1,1)-v2(1,1);v1(2,1)-v2(2,1);v1(3,1)-v2(3,1);v1(4,1)-v2(4,1);v1(5,1)-v2(5,1);v1(6,1)-v2(6,1);v2(1,1)-v3(1,1);v2(2,1)-v3(2,1);v2(3,1)-v3(3,1);v2(4,1)-v3(4,1);v2(5,1)-v3(5,1);v2(6,1)-v3(6,1);v1(1,1)-v3(1,1);v1(2,1)-v3(2,1);v1(3,1)-v3(3,1);v1(4,1)-v3(4,1);v1(5,1)-v3(5,1);v1(6,1)-v3(6,1);v3(1,1)-v4(1,1);v3(2,1)-v4(2,1);v3(3,1)-v4(3,1);v3(4,1)-v4(4,1);v3(5,1)-v4(5,1);v3(6,1)-v4(6,1)];

for idx=1:Np+1
    if idx==1

        %have to initialize the position and velocity vector
        %x(:,idx)     = x(:,1);
        xpinv(:,idx) = x(:,1);
        xhat(:,idx)  = x(:,1);
    else
        idxn = ((idx-1) * 10) + 1;
        z(:,idxn) = [xf{idx-1}(1:2,idxn)- xf{idx-1}(3:4,idxn); xf{idx-1}(3:4,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(5:6,idxn)- xf{idx-1}(7:8,idxn);xf{idx-1}(1:2,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(9:10,idxn)- xf{idx-1}(11:12,idxn); xf{idx-1}(11:12,idxn)- xf{idx-1}(13:14,idxn); xf{idx-1}(13:14,idxn)- xf{idx-1}(15:16,idxn); xf{idx-1}(9:10,idxn)- xf{idx-1}(13:14,idxn)];
        %if reced the horizon then x(l) be the initial state
        %x(:,idxn)    = xd{idx-1}(:,idxn);
        xpinv(:,idxn)= xe{idx-1}(:,idxn);
        xhat(:,idxn) = xf{idx-1}(:,idxn);
    end
    %

    T = (idx-1)*10 + step;
    Pt(:,:,101) = Qft; %initial Riccati
    aN = 1;
    aE =1;
    %Discrete
    for k=1:100 %horizon time
        p1d     = [aN*pi*cos((k*pi)/20);-aE*pi*sin((k*pi)/20);1*dt*k;0;0;0];
        p2d  = p1d;
        p3d  = p1d;
        p4d  = p1d;
        v1d   = [0;0;0.1;0;0;0];
        v2d = [0;0;0.1;0;0;0];
        v3d = [0;0;0.1;0;0;0];
        v4d = [0;0;0.1;0;0;0];
        r(:,k)  = [p1d;p2d;p3d;p4d;v1d;v2d;v3d;v4d];
        e(:,1)  = abs(r(:,1) - x(:,1));

        for j = T-1:-1:1
            Pt(:,:,j) = Qt + Ft'*Pt(:,:,j+1)*Ft - Ft'*Pt(:,:,j+1)*Gt*inv(Rt+Gt'*Pt(:,:,j+1)*Gt)*Gt'*Pt(:,:,j+1)*Ft;
        end
        Kt(:,:,k) = -inv(Rt+Gt'*Pt(:,:,k+1)*Gt)*Gt'*Pt(:,:,k+1)*Ft;
        a(:,k)    = Kt(:,:,k) * z(:,k);
        z(:,k+1)  = Ft * z(:,k) + Gt * a(:,k);

        for v = 1:N
            upinv(:,k) = Phiai * a(:,k);
        end

        e(:,k+1) =  Fpinv * e(:,k)+Gpinv*upinv(:,k);
        xpinv(:,k+1) = e(:,k)+r(:,k);
        %iterative method to get back uhat from stored a
        %tic;
        t = 1;
        if k>1
            u0 = usol;
        end
        %
        while (t <= Iter)
            u0_temp = (In - 2*alpha*Phia'*Phia)*u0(:,t) + 2*alpha*Phia'*a(:,k);
            u0 = [u0 u0_temp];
            t = t+1;
        end
        udata = [udata,u0];
        %toc;
        usol      = u0(:,end);
        uhat(:,k) = usol;
        e(:,k+1) = Fpinv * e(:,k) + Gpinv*uhat(:,k);
        xhat(:,k+1) = e(:,k)+r(:,k);
    end
end

%Plot to show the formation and reference tracking for 4 agents and helical
%Trajectory.
figure(2);
plot3(xhat(1,2:50),xhat(2,2:50),xhat(3,2:50),'-r',xhat(7,2:50),xhat(8,2:50),xhat(9,2:50),'-b',xhat(13,2:50),xhat(14,2:50),xhat(15,2:50),'-g',xhat(19,2:50),xhat(20,2:50),xhat(21,2:50),'linewidth',1.5);
hold on
grid on
set(gca,'color',[0.9,0.9,0.9]);
title('Motion of Positions in Fomation for 4 agents','fontweight','bold')
legend('Agent 1','Agent 2','Agent 3','Agent 4','Location','Best')
xlabel('x-axis');
ylabel('y-axis');

plot3(xhat(1,2:2),xhat(2,2:2),xhat(3,2:2),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#2f54eb')
plot3(xhat(1,50:50),xhat(2,50:50),xhat(3,50:50),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#2f54eb')

plot3(xhat(7,2:2),xhat(8,2:2),xhat(9,2:2),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#f51146')
plot3(xhat(7,50:50),xhat(8,50:50),xhat(9,50:50),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#f51146')

plot3(xhat(13,2:2),xhat(14,2:2),xhat(15,2:2),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#384d3f')
plot3(xhat(13,50:50),xhat(14,50:50),xhat(15,50:50),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#384d3f')

plot3(xhat(19,2:2),xhat(20,2:2),xhat(21,2:2),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#a6e3b6')
plot3(xhat(19,50:50),xhat(20,50:50),xhat(21,50:50),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#a6e3b6')

legend('Agent 1','Agent 2','Agent 3','Agent 4','Location','Best')
xlabel('x-axis');
ylabel('y-axis');

%%

                %****************Section 3***********************
             %Infinity Trajectory tracking with formation for 3 agents 
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
d =[-1; -0.5;-1;0;0;0;1;-0.5;1;0;0;0]; %formation parameter for 3 agents


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
% p1 = [3;4;3;0;0;0];
% p2 = [2;4;3;0;0;0];
% p3 = [1;5;3;0;0;0];

p1 = [0;0;1;0;0;0];
p2 = [0;1;1;0;0;0];
p3 = [0;-1;1;0;0;0];

% Velocity vectors (linear and Angular)
v1 =[0;0;1;0;0;0];
v2 =[0;0;1;0;0;0];
v3 =[0;0;1;0;0;0];

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

z(:,1) = [p1(1,1)-p2(1,1)-d(1,1);p1(2,1)-p2(2,1)-d(2,1);p1(3,1)-p2(3,1)-d(3,1);p1(4,1)-p2(4,1);p1(5,1)-p2(5,1);p1(6,1)-p2(6,1);p1(1,1)-p3(1,1)-d(7,1);p1(2,1)-p3(2,1)-d(8,1);p1(3,1)-p3(3,1)-d(9,1);p1(4,1)-p3(4,1);p1(5,1)-p3(5,1);p1(6,1)-p3(6,1);v1(1,1)-v2(1,1);v1(2,1)-v2(2,1);v1(3,1)-v2(3,1);v1(4,1)-v2(4,1);v1(5,1)-v2(5,1);v1(6,1)-v2(6,1);v1(1,1)-v3(1,1);v1(2,1)-v3(2,1);v1(3,1)-v3(3,1);v1(4,1)-v3(4,1);v1(5,1)-v3(5,1);v1(6,1)-v3(6,1)];

for idx=1:Np+1
    if idx==1

        %have to initialize the position and velocity vector

        xpinv(:,idx) = x(:,1);
        xhat(:,idx)  = x(:,1);
    else
        idxn = ((idx-1) * 10) + 1;
        z(:,idxn) = [xf{idx-1}(1:2,idxn)- xf{idx-1}(3:4,idxn); xf{idx-1}(3:4,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(5:6,idxn)- xf{idx-1}(7:8,idxn);xf{idx-1}(1:2,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(9:10,idxn)- xf{idx-1}(11:12,idxn); xf{idx-1}(11:12,idxn)- xf{idx-1}(13:14,idxn); xf{idx-1}(13:14,idxn)- xf{idx-1}(15:16,idxn); xf{idx-1}(9:10,idxn)- xf{idx-1}(13:14,idxn)];
        %if reced the horizon then x(l) be the initial state
        xpinv(:,idxn)= xe{idx-1}(:,idxn);
        xhat(:,idxn) = xf{idx-1}(:,idxn);
    end
    %
    T = (idx-1)*10 + step;
    Pt(:,:,101) = Qft; %initial Riccati
    aN = 5;
    aE =5;
    %Discrete
    for k=1:100%horizon time
        p1d     = [p1(1,1)+aN*sin((k^1.5)/30);p1(2,1)+aE*sin((k^1.5)/60);p1(3,1)+0.1*dt*k;0;0;0];
        p2d     = [p2(1,1)+aN*sin((k^1.5)/30);p2(2,1)+aE*sin((k^1.5)/60);p2(3,1)+0.1*dt*k;0;0;0];
        p3d     = [p3(1,1)+aN*sin((k^1.5)/30);p3(2,1)+aE*sin((k^1.5)/60);p3(3,1)+0.1*dt*k;0;0;0];
        v1d   = [0;0;0.1;0;0;0];
        v2d = [0;0;0.1;0;0;0];
        v3d = [0;0;0.1;0;0;0];
        r(:,k)  = [p1d;p2d;p3d;v1d;v2d;v3d];
        e(:,1)  = abs(r(:,1) - x(:,1));

        for j = T-1:-1:1
            Pt(:,:,j) = Qt + Ft'*Pt(:,:,j+1)*Ft - Ft'*Pt(:,:,j+1)*Gt*inv(Rt+Gt'*Pt(:,:,j+1)*Gt)*Gt'*Pt(:,:,j+1)*Ft;
        end
        Kt(:,:,k) = -inv(Rt+Gt'*Pt(:,:,k+1)*Gt)*Gt'*Pt(:,:,k+1)*Ft;
        a(:,k)    = Kt(:,:,k) * z(:,k);
        z(:,k+1)  = Ft * z(:,k) + Gt * a(:,k);

        for v = 1:N
            upinv(:,k) = Phiai * a(:,k);
        end

        e(:,k+1) =  Fpinv * e(:,k)+Gpinv*upinv(:,k);
        xpinv(:,k+1) = e(:,k)+r(:,k);
        %iterative method to get back uhat from stored a
        %tic;
        t = 1;
        if k>1
            u0 = usol;
        end
        %
        while (t <= Iter)
            u0_temp = (In - 2*alpha*Phia'*Phia)*u0(:,t) + 2*alpha*Phia'*a(:,k);
            u0 = [u0 u0_temp];
            t = t+1;
        end
        udata = [udata,u0];
        %toc;
        usol      = u0(:,end);
        uhat(:,k) = usol;
        e(:,k+1) = Fpinv * e(:,k) + Gpinv*uhat(:,k);
        xhat(:,k+1) = e(:,k)+r(:,k);
    end
end
figure(3);
plot3(xhat(1,1:80),xhat(2,1:80),xhat(3,1:80),'-r',xhat(7,1:80),xhat(8,1:80),xhat(9,1:80),'-b',xhat(13,1:80),xhat(14,1:80),xhat(15,1:80),'-g','linewidth',1.5);
hold on
grid on
set(gca,'color',[0.9,0.9,0.9]);
title('Motion of Positions in Fomation','fontweight','bold')
legend('Agent 1','Agent 2','Agent 3','Location','Best')
xlabel('x-axis');
ylabel('y-axis');

plot3(xhat(1,1:1),xhat(2,1:1),xhat(3,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#2f54eb')
plot3(xhat(1,80:80),xhat(2,80:80),xhat(3,80:80),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#2f54eb')

plot3(xhat(7,1:1),xhat(8,1:1),xhat(9,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#f51146')
plot3(xhat(7,80:80),xhat(8,80:80),xhat(9,80:80),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#f51146')

plot3(xhat(13,1:1),xhat(14,1:1),xhat(15,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#384d3f')
plot3(xhat(13,80:80),xhat(14,80:80),xhat(15,80:80),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#384d3f')

legend('Agent 1','Agent 2','Agent 3','Location','Best')
xlabel('x-axis');
ylabel('y-axis');

%%
                %****************Section 4***********************
             %Infinity Trajectory tracking with formation for 4 agents 

clear all
close all
clc

n = 3; %Dimensional Plane
N = 4; %number of agents
M = 4; %number of edges
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
Bb(:,:,4) = B(:, (3*ni)+1:4*ni);

for i = 1:N
    Gb(:,:,i) = dt*Bb(:,:,i) + (dt^2)/2*A*Bb(:,:,i) + (dt^3)/6*A^2*(Bb(:,:,i)).^2 + (dt^4)/24*A^3*(Bb(:,:,i)).^3;
end

mui = [-1 0 0 -1;
    1 -1 0 0;
    0 1 -1 1;
    0 0 1 0];
%mu =kron(mui,eye(2*n));

Di = [-1 0 0 -1;
    1 -1 0 0;
    0 1 -1 1;
    0 0 1 0];
D= kron(Di,eye(2*n)); %incidence matrix for 3 agents

d12 = [-1.2;-0.4;0];
d13 = [0;-0.8;0];
d =[15; 0;-4;0;0;0;-5;-16;4;0;0;0;4; -15;16;0;0;0;2;7;-5;0;0;0]; %formation parameter for 3 agents


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
p1 = [14;3;1;0;0;0];
p2 = [2;4;1;0;0;0];
p3 = [5;15;1;0;0;0];
p4 = [16;4;1;0;0;0]
% Velocity vectors (linear and Angular)
v1 =[0;0;1;0;0;0];
v2 =[0;0;1;0;0;0];
v3 =[0;0;1;0;0;0];
v4 =[0;0;1;0;0;0];
x(:,1) = [p1; p2; p3;p4; v1;v2;v3;v4];
xb(:,1) = [p1; p2; p3;p4;1; v1;v2;v3;v4];

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

Fpinv = [F(1:24, 1:24), F(1:24,26:49);F(26:49, 1:24),F(26:49, 26:49)]
Gpinv =[Gb(1:24,1:16);Gb(26:49,1:16)];
% Bt    = [zeros(2*n*M);eye(2*n*M)];

Qt = dt*eye(4*M*n);
Qft   = eta*Qt;
Rt =eye(4*ni);
Gt    = dt*Bt + (dt^2)/2*At*Bt;


% Dt= kron(Di,eye(n));
Phi   = -Di';
Phia  = kron(Phi,eye(4));
Phiai = pinv(Phia);

Iter    = 200;
In      = eye(16);
u0      = zeros(4*M,1);
eps     = 1e-3;
udata   = [];
alpha = 0.1;
%Z(:,1) = [xb(1,:)-xb(7,:);xb(2,:)-xb(8,:);xb(3,:)-xb(9,:);xb(4,:)-xb(10,:);xb(5,:)-xb(11,:);xb(6,:)-xb(12,:);xb(1,:)-xb(13,:);xb(2,:)-xb(14,:);xb(3,:)-xb(15,:);xb(4,:)-xb(16,:);xb(5,:)-xb(17,:);xb(6,:)-xb(18,:);
z(:,1) = [p1(1,1)-p2(1,1)-d(1,1);p1(2,1)-p2(2,1)-d(2,1);p1(3,1)-p2(3,1)-d(3,1);p1(4,1)-p2(4,1);p1(5,1)-p2(5,1);p1(6,1)-p2(6,1);p2(1,1)-p3(1,1)-d(7,1);p2(2,1)-p3(2,1)-d(8,1);p2(3,1)-p3(3,1)-d(9,1);p2(4,1)-p3(4,1);p2(5,1)-p3(5,1);p2(6,1)-p3(6,1);p1(1,1)-p3(1,1)-d(13,1);p1(2,1)-p3(2,1)-d(14,1);p1(3,1)-p3(3,1)-d(15,1);p1(4,1)-p3(4,1);p1(5,1)-p1(5,1);p1(6,1)-p3(6,1);p3(1,1)-p4(1,1)-d(13,1);p3(2,1)-p4(2,1)-d(14,1);p3(3,1)-p4(3,1)-d(15,1);p3(4,1)-p4(4,1);p3(5,1)-p4(5,1);p3(6,1)-p4(6,1);v1(1,1)-v2(1,1);v1(2,1)-v2(2,1);v1(3,1)-v2(3,1);v1(4,1)-v2(4,1);v1(5,1)-v2(5,1);v1(6,1)-v2(6,1);v2(1,1)-v3(1,1);v2(2,1)-v3(2,1);v2(3,1)-v3(3,1);v2(4,1)-v3(4,1);v2(5,1)-v3(5,1);v2(6,1)-v3(6,1);v1(1,1)-v3(1,1);v1(2,1)-v3(2,1);v1(3,1)-v3(3,1);v1(4,1)-v3(4,1);v1(5,1)-v3(5,1);v1(6,1)-v3(6,1);v3(1,1)-v4(1,1);v3(2,1)-v4(2,1);v3(3,1)-v4(3,1);v3(4,1)-v4(4,1);v3(5,1)-v4(5,1);v3(6,1)-v4(6,1)];


for idx=1:Np+1
    if idx==1

        %have to initialize the position and velocity vector
        %x(:,idx)     = x(:,1);
        xpinv(:,idx) = x(:,1);
        xhat(:,idx)  = x(:,1);
    else
        idxn = ((idx-1) * 10) + 1;
        z(:,idxn) = [xf{idx-1}(1:2,idxn)- xf{idx-1}(3:4,idxn); xf{idx-1}(3:4,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(5:6,idxn)- xf{idx-1}(7:8,idxn);xf{idx-1}(1:2,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(9:10,idxn)- xf{idx-1}(11:12,idxn); xf{idx-1}(11:12,idxn)- xf{idx-1}(13:14,idxn); xf{idx-1}(13:14,idxn)- xf{idx-1}(15:16,idxn); xf{idx-1}(9:10,idxn)- xf{idx-1}(13:14,idxn)];
        %if reced the horizon then x(l) be the initial state
        %x(:,idxn)    = xd{idx-1}(:,idxn);
        xpinv(:,idxn)= xe{idx-1}(:,idxn);
        xhat(:,idxn) = xf{idx-1}(:,idxn);
    end
    %

    T = (idx-1)*10 + step;
    Pt(:,:,101) = Qft; %initial Riccati
    aN = 8;
    aE =8;
    %Discrete
    for k=1:100%horizon time

        p1d     = [aN*pi*sin((k^1.5)/30);aE*pi*sin((k^1.5)/60);1*dt*k;0;0;0];
        p2d = p1d;
        p3d = p1d;
        p4d = p1d;
        v1d   = [0;0;0.1;0;0;0];
        v2d = [0;0;0.1;0;0;0];
        v3d = [0;0;0.1;0;0;0];
        v4d = [0;0;0.1;0;0;0];
        r(:,k)  = [p1d;p2d;p3d;p4d;v1d;v2d;v3d;v4d];
        e(:,1)  = abs(r(:,1) - x(:,1));

        for j = T-1:-1:1
            Pt(:,:,j) = Qt + Ft'*Pt(:,:,j+1)*Ft - Ft'*Pt(:,:,j+1)*Gt*inv(Rt+Gt'*Pt(:,:,j+1)*Gt)*Gt'*Pt(:,:,j+1)*Ft;
        end
        Kt(:,:,k) = -inv(Rt+Gt'*Pt(:,:,k+1)*Gt)*Gt'*Pt(:,:,k+1)*Ft;
        a(:,k)    = Kt(:,:,k) * z(:,k);
        z(:,k+1)  = Ft * z(:,k) + Gt * a(:,k);

        for v = 1:N
            upinv(:,k) = Phiai * a(:,k);
        end

        e(:,k+1) =  Fpinv * e(:,k)+Gpinv*upinv(:,k);
        xpinv(:,k+1) = e(:,k)+r(:,k);
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
        e(:,k+1) = Fpinv * e(:,k) + Gpinv*uhat(:,k);
        xhat(:,k+1) = e(:,k)+r(:,k);
    end
end

%Plot to get theTrajectory tracking and Formation for 4 agent infinity loop
figure(4);

plot3(xhat(1,1:70),xhat(2,1:70),xhat(3,1:70),'-r',xhat(7,1:70),xhat(8,1:70),xhat(9,1:70),'-b',xhat(13,1:70),xhat(14,1:70),xhat(15,1:70),'-g',xhat(19,1:70),xhat(20,1:70),xhat(21,1:70),'linewidth',1.5);
hold on
grid on

set(gca,'color',[0.9,0.9,0.9]);
title('Motion of Positions in Fomation','fontweight','bold')

plot3(xhat(1,1:1),xhat(2,1:1),xhat(3,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#2f54eb')
plot3(xhat(1,70:70),xhat(2,70:70),xhat(3,70:70),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#2f54eb')

plot3(xhat(7,1:1),xhat(8,1:1),xhat(9,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#f51146')
plot3(xhat(7,70:70),xhat(8,70:70),xhat(9,70:70),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#f51146')

plot3(xhat(13,1:1),xhat(14,1:1),xhat(15,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#384d3f')
plot3(xhat(13,70:70),xhat(14,70:70),xhat(15,70:70),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#384d3f')

plot3(xhat(19,1:1),xhat(20,1:1),xhat(21,1:1),'-s','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#b521cc')
plot3(xhat(19,70:70),xhat(20,70:70),xhat(21,70:70),'-o','Color','b','MarkerSize',8,...
    'MarkerFaceColor','#b521cc')

legend('Agent 1','Agent 2','Agent 3','Agent 4','Location','Best')
xlabel('x-axis');
ylabel('y-axis');


%%

                %****************Section 5***********************
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

%%
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

%%
