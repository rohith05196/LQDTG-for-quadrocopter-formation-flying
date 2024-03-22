function fourAgentsInfinity(n,N,M,ns,ni,tf,dt,eta,Np,step)

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











end