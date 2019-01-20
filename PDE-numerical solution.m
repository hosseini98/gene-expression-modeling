%%Solving dilmensionless PDE: diff(p,tau)=-diff(p,X)+1/Pe*diff(p,X,2)
%%exposed to the initial and boundary conditions:

%%p = 0           at t=0
%%p - 1/Pe*diff(p,x)=k_init/a*P_PIC     at x=0
%%diff(p,x) = 0   at x=1

clear;
clc;

%Rate constants
k_ae    = 144;
k_elong = 144;
k_init  = 0.6;
k_PIC   = 0.0029;
k_bind  = 0.0016;


%define variavbles
Xmin    = 0;           %minimum value of X
Xmax    = 1;           %maxmimum value of X
N       = 1000;          %no. nodes - 1
a       = k_ae*k_elong/((k_ae+k_elong)*N);           %velocity
Pe      = 2*N;
dtau    = 0.004 %timestep
tau     = 0;            %dimensionless time
tau_max = 200;        %maximum value of time


M=200
%descritize the domain
dX     = (Xmax-Xmin)/M
X=  Xmin:dX:Xmax; %includes ghost nodes
x=X;




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
    p(1,n+1) = Ftau(n) + (p(2,n)-p(1,n))/(Pe*dX);
    
    
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
    
    
    
 
    
%     plot(X,p(:,n+1))
%     set(gca,'FontSize',10,...
%         'TickDir','out',...
%         'XLim',[0,1])
%     
%     shg;
end


figure


% subplot(1, 3, 1);
% plot(tau2,p(10,:))
% xlabel('tau')
% ylabel('probability')
% set(gca,'FontSize',12,...
%     'TickDir','out')
% 
% subplot(1, 3, 2);
% plot(X,p(:,3830))
% xlabel('X')
% ylabel('Probability')
% set(gca,'FontSize',12,...
%     'TickDir','out',...
%     'XLim',[0,1])
% %
% subplot(1, 3, 3);
% plot(X,p(:,4500))
% xlabel('X')
% ylabel('Probability')
% set(gca,'FontSize',12,...
%     'TickDir','out',...
%     'XLim',[0,1])
shg;
% pp=mesh(tau2,X,p)
plot(X,p(:,8400),'r','LineWidth',2)
xlabel('X')
ylabel('Probability')
set(gca,'FontSize',12,...
    'TickDir','out',...
    'XLim',[0,1])
hold on
%  load 'set9.txt'
% plot(set9(:,1),set9(:,2),'LineWidth',2)

% figure(1)
% filename = 'testnew51.gif';
% for n = 1:20:2000
%       
%       plot(X,p(:,n),'r','LineWidth',2)
% xlabel('X')
% ylabel('Probability')
% set(gca,'FontSize',12,...
%     'TickDir','out',...
%     'XLim',[0,1],...
%     'YLim',[0,1.52])
%       drawnow
%       frame = getframe(1);approximate-kae1
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if n == 1;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%       end
% end



% x=0:0.01:1;
% t=0:0.01:1;
% 
% ppp= -N * k_init * k_PIC * (k_ae + k_elong) * (heaviside(-t * k_ae * k_elong / ...
%     (k_ae + k_elong) / N + x) - 1) * k_bind * (-(k_PIC - k_init) * ((N * ...
%     k_bind ^ 2 * x + 4 * k_ae ^ 2 + 2 * k_ae * k_bind) * k_elong ^ 2 ...
%     + 2 * k_ae * k_bind * (N * k_bind * x + k_ae) * k_elong + N * x *...
%     k_ae ^ 2 * k_bind ^ 2) * exp((k_bind * ((N * x - k_ae * t) * k_elong ...
%     + x * N * k_ae) / k_ae / k_elong)) + (k_PIC - k_bind) * ((N * k_init...
%     ^ 2 * x + (4 * k_ae ^ 2) + 0.2e1 * k_ae * k_init) * (k_elong ^ 2) ...
%     + 0.2e1 * k_init * k_ae * (N * k_init * x + k_ae) * k_elong + N *...
%     x * (k_ae ^ 2) * k_init ^ 2) * exp(k_init * ((N * x - k_ae * t) * ...
%     k_elong + x * N * k_ae) / k_ae / k_elong) + exp(k_PIC * ((N * x -...
%     k_ae * t) * k_elong + x * N * k_ae) / k_ae / k_elong) * (k_bind ...
%     - k_init) * ((N * k_PIC ^ 2 * x + 0.2e1 * k_PIC * k_ae + (4 * k_ae ^ 2))...
%     * (k_elong ^ 2) + 0.2e1 * k_ae * k_PIC * (N * k_PIC * x + k_ae) * k_elong ...
%     + N * x * (k_ae ^ 2) * k_PIC ^ 2)) / (k_PIC - k_bind) / (k_PIC - k_init) /...
%     (k_bind - k_init) / (k_ae ^ 3) / (k_elong ^ 3) / 0.4e1;



% ppp(n)= -N * k_init * k_PIC * (k_ae + k_elong) * (heaviside(-t(n) * k_ae * k_elong / ...
%     (k_ae + k_elong) / N + x) - 1) * k_bind * (-(k_PIC - k_init) * ((N * ...
%     k_bind ^ 2 * x + 4 * k_ae ^ 2 + 2 * k_ae * k_bind) * k_elong ^ 2 ...
%     + 2 * k_ae * k_bind * (N * k_bind * x + k_ae) * k_elong + N * x *...
%     k_ae ^ 2 * k_bind ^ 2) * exp((k_bind * ((N * x - k_ae * t(n)) * k_elong ...
%     + x * N * k_ae) / k_ae / k_elong)) + (k_PIC - k_bind) * ((N * k_init...
%     ^ 2 * x + (4 * k_ae ^ 2) + 0.2e1 * k_ae * k_init) * (k_elong ^ 2) ...
%     + 0.2e1 * k_init * k_ae * (N * k_init * x + k_ae) * k_elong + N *...
%     x * (k_ae ^ 2) * k_init ^ 2) * exp(k_init * ((N * x - k_ae * t(n)) * ...
%     k_elong + x * N * k_ae) / k_ae / k_elong) + exp(k_PIC * ((N * x -...
%     k_ae * t(n)) * k_elong + x * N * k_ae) / k_ae / k_elong) * (k_bind ...
%     - k_init) * ((N * k_PIC ^ 2 * x + 0.2e1 * k_PIC * k_ae + (4 * k_ae ^ 2))...
%     * (k_elong ^ 2) + 0.2e1 * k_ae * k_PIC * (N * k_PIC * x + k_ae) * k_elong ...
%     + N * x * (k_ae ^ 2) * k_PIC ^ 2)) / (k_PIC - k_bind) / (k_PIC - k_init) /...
%     (k_bind - k_init) / (k_ae ^ 3) / (k_elong ^ 3) / 0.4e1;
% 
% t(n+1)=t(n)+dtau;
% end


save pde8.mat '-v7.3'
