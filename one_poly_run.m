% Control program for a one-polymerase model

more off
tic
%parpool(4);

% Number of nucleotides in the strand
n_nuc = int32(1000);

% States
statenames = struct('promoter_empty_state',int8(0),'TBPpro_state',int8(-2),...
		'PICstate',int8(-1),'Ostate',int8(3),'Astate',int8(4),...
		'Pstate',int8(5),'TCstate',int32(n_nuc+1),...
		'terminated_state',int32(n_nuc+2));

% Rate constants
c(1) = 0.0016;  % TBP + pro -> TBP.pro
c(2) = 0.0029;  % TBP.pro + RNAP -> PIC
c(3) = 0.6;     % PIC + U(Delta) -> O1
c(4) = 0;	% Unused in this simplified model
c(5) = 0.0;    % Ai -> RNA(aborted) + Ui for i= 4 to 15
c(6) = 144;     % Ai + U(i+Delta) -> O(i+1) + Ui
c(7) = 144;     % Oi -> Ai
c(8) = 0;     % Oi -> Oi(paused) for i = 16 to 50
c(9) = 0;   % Oi(paused) -> Oi for i = 16 to 50
c(10) = 0.0032; % An -> TC
c(11) = 0.0032; % TC -> RNA + Un

% Simulation parameters
RNAsynth_target = int32(40000);      % Number of RNAs to synthesize
tmax = 6000;
noutput = 12000;
t = [0:tmax/noutput:tmax];	% Requested output times

% Array in which to return time evolution of polymerase position
% Each row of this array corresponds to one simulation.
position = zeros(RNAsynth_target,noutput+1,'int32');

% Get the trajectories
for i=1:RNAsynth_target
    [position(i,:)] = one_poly_simplified(n_nuc,statenames,c,t);
end
disp('Trajectories computed')
toc

tic
% Compute probability distribution
% The next three variables give the probability of being in the TBPpro, PIC,
% TC and terminated states, respectively, as a function of time.
pTBPpro = sum(position == statenames.TBPpro_state)/double(RNAsynth_target);
pPIC = sum(position == statenames.PICstate)/double(RNAsynth_target);
pAn = sum(position==n_nuc)/double(RNAsynth_target);
pTC = sum(position==statenames.TCstate)/double(RNAsynth_target);
pterminated = sum(position==statenames.terminated_state)...
	      /double(RNAsynth_target);
% The distribution in the elongation compartment at a given t is in the
% corresponding column of elong_distrib.
% You can plot the elongation distribution at a given time point by
% plot(elong_distrib(:,i)), where i is the time index (corresponding to
% the values in the array t).
elong_distrib = zeros(n_nuc-1,noutput+1);
for i=1:n_nuc-1
    elong_distrib(i,:) = sum(position==i)/double(RNAsynth_target);
end
toc
% w=9;
%  course_grained_prob=zeros((n_nuc-1)/w,noutput+1,'int32');
% 
% 
% for j=1:noutput+1
%     course_grained_prob(:,j) = mean(reshape(elong_distrib(:,j),w,[]));
% end
% Save data.
save one_poly_simplified.mat -v7.3

dy=diff(pterminated)./diff(t);
 %plot(t(2:end),dy);
w=400;
tt = mean(reshape(t(2:end),w,[]));
     course_grained_prob = mean(reshape(dy,w,[]));
plot(tt,course_grained_prob);



