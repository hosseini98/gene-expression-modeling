function pos = one_poly_simplified(n_nuc,statenames,c,toutput)
% A simplified one-polymerase model
% Input arguments:
%	n_nuc:		number of nucleotides
%	statenames:	struct containing the names of the states
%	c:		rate constants
%	toutput:	output times
% Because this model reduces to a simple sequence of steps, the program
% short-circuits the Gillespie algorithm.

% Housekeeping variables
t = 0;
position = statenames.promoter_empty_state;    % Position of the polymerase
noutput = length(toutput);
pos = zeros(noutput,1,'int32');
next_output = 1;
if t == toutput(next_output)
    pos(next_output) = position;
    next_output = next_output + 1;
end

% TBP + pro -> TBP.pro
a = c(1);
r = rand(4,1);	% Random numbers needed for initiation steps.
t = t + log(1/r(1))/a;
while t > toutput(next_output)
    pos(next_output) = position;
    next_output = next_output + 1;
    if next_output > noutput
        return
    end
end
position = statenames.TBPpro_state;

% TBP.pro + RNAP -> PIC
a = c(2);
t = t + log(1/r(2))/a;
while t > toutput(next_output)
    pos(next_output) = position;
    next_output = next_output + 1;
    if next_output > noutput
        return
    end
end
position = statenames.PICstate;

% PIC + U(Delta) -> O1
a = c(3);
t = t + log(1/r(3))/a;
while t > toutput(next_output)
    pos(next_output) = position;
    next_output = next_output + 1;
    if next_output > noutput
        return
    end
end
position = 1;

% O1 -> A1
a = c(7);
t = t + log(1/r(4))/a;
while t > toutput(next_output)
    pos(next_output) = position;
    next_output = next_output + 1;
    if next_output > noutput
        return
    end
end

% Generate all random numbers required in elongation region.
r = rand(n_nuc-1,2);
% Elongate
for i=2:n_nuc
    % Ai-1 -> Oi
    a = c(6);
    t = t + log(1/r(i-1,1))/a;
    while t > toutput(next_output)
        pos(next_output) = position;
        next_output = next_output + 1;
        if next_output > noutput
            return
        end
    end
    position = i;
    
    % Oi -> Ai
    a = c(7);
    t = t + log(1/r(i-1,2))/a;
    while t > toutput(next_output)
        pos(next_output) = position;
        next_output = next_output + 1;
        if next_output > noutput
            return
        end
    end
end

% Terminate in two steps
r = rand(2,1);
% An -> TC
a = c(10);
t = t + log(1/r(1))/a;
while t > toutput(next_output)
    pos(next_output) = position;
    next_output = next_output + 1;
    if next_output > noutput
        return
    end
end
position = statenames.TCstate;

% TC -> RNA + Un
a = c(11);
t = t + log(1/r(2))/a;
while t > toutput(next_output)
    pos(next_output) = position;
    next_output = next_output + 1;
    if next_output > noutput
        return
    end
end
pos(next_output:end) = statenames.terminated_state;

end
