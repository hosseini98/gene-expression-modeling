%%Solving dilmensionless PDE: diff(p,tau)=-diff(p,X)+1/Pe*diff(p,X,2)
%%exposed to the initial and boundary conditions:

%%p = 0           at t=0
%%p - 1/Pe*diff(p,x)=k_init/a*P_PIC     at x=0
%%diff(p,x) = 0   at x=1

clear;
clc;

%Rate constants
k_ae    = 0.001;
k_elong = 0.001;
k_init  = 0.6;
k_PIC   = 0.0029;
k_bind  = 0.0016;


%define variavbles
Xmin    = 0;           %minimum value of X
Xmax    = 1;           %maxmimum value of X
N       = 100;          %no. nodes - 1
a       = k_ae*k_elong/((k_ae+k_elong));           %velocity
Pe      = 2*N;
dtau    = 0.001 %timestep
tau     = 0;            %dimensionless time
tau_max = 20;        %maximum value of time


M=50
%descritize the domain
dX     = (Xmax-Xmin)/M
X=  Xmin:dX:Xmax; %includes ghost nodes
%x=X




%plot(X,p0)



r=dtau/(dX^2*Pe)




nsteps=tau_max/dtau
tau2=zeros(1,nsteps);
t=zeros(1,nsteps);
cg=zeros(1,nsteps+1);
%set initial condition
sz     = size(X);
p0     = zeros(sz(2),nsteps);
p    = p0;

PPIC = zeros(nsteps,1);
%loop through time

for n=1:nsteps
    
    PPIC(n)=-k_bind*k_PIC  * ((k_PIC - k_init)*exp(-tau/a*(k_bind))...
         +(-k_bind + k_init)*exp(-tau/a*(k_PIC)) - (k_PIC + k_bind)*exp(-k_init*tau/a))...
         /((k_PIC - k_bind) * (k_PIC - k_init)*(k_bind - k_init));
    Ftau(n)=k_init/a*PPIC(n);
    
    %calculate boundary conditions
    %p(1,n)=(Pe*dX*Ftau+p(2,n))/(1+Pe*dX);
    %p(N+1,n)=p(N,n);
    p(1,n+1) = Ftau(n) +(p(2,n)-p(1,n))/(Pe*dX);
    
    
    %first order Upwind scheme
    
    % for i = 2:N+2
    %  pnp1(i,n+1)=p(i-1,n)*(1+Pe*dX)*r + p(i,n)*(1-r*Pe*dX-2*r) + p(i+1,n)*r;
    %end
    
    p(2:M,n+1)=p(1:M-1,n)*(1+Pe*dX)*r + p(2:M,n)*(1-r*Pe*dX-2*r) + p(3:M+1,n)*r;
    p(M+1,n+1) = p(M,n+1);
    
    %update t and p
    tau = tau+dtau;
    %p = pnp1;
    
    %plot solution
    
    %plot(X,p,'bo', 'markerfacecolor', 'b');
    %shg
    %pause(dtau)
    
    tau2(n+1)=tau2(n)+dtau;
    
    
    
 
    
    plot(X,p(:,n+1))
    set(gca,'FontSize',10,...
        'TickDir','out',...
        'XLim',[0,1])
    
    shg;
end


figure


subplot(1, 3, 1);
plot(tau2/a,p(10,:))
xlabel('t')
ylabel('probability')
set(gca,'FontSize',12,...
    'TickDir','out')

subplot(1, 3, 2);
plot(X,p(:,1000))
xlabel('X')
ylabel('Probability')
set(gca,'FontSize',12,...
    'TickDir','out',...
    'XLim',[0,1])
%
subplot(1, 3, 3);
plot(X,p(:,800))
xlabel('X')
ylabel('Probability')
set(gca,'FontSize',12,...
    'TickDir','out',...
    'XLim',[0,1])
shg;
% pp=mesh(tau2,X,p)
% plot(X,p(:,3830),'r','LineWidth',2)
% xlabel('X')
% ylabel('Probability')
% set(gca,'FontSize',12,...
%     'TickDir','out',...
%     'XLim',[0,1])
% hold on
%  load 'set2.txt'
% plot(set2(:,1),set2(:,2),'LineWidth',2)









