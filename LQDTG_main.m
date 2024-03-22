clear all
close all
clc
n = 3;      %Dimensional Plane (3D) in this case
N = 3;      %number of agents 
M = 2;      %number of edges
ns = 12;    % number of states
ni =4;      %number of inputs
tf = 20;    % t final
dt = 0.1;   %steps
eta =4;
Np = 0;

step = 100;


threeAgentsHelical(n,N,M,ns,ni,tf,dt,eta,Np,step);
threeAgentsInfinity(n,N,M,ns,ni,tf,dt,eta,Np,step);
fourAgentsHelical(n,N,M,ns,ni,tf,dt,eta,Np,step);
threeAgentsInfinity(n,N,M,ns,ni,tf,dt,eta,Np,step);