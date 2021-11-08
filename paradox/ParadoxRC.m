function [Votes] = ParadoxRC(RC,ptype,etype)
%% Strong No-Show
if strcmp(ptype, 'no-show') == 1
% If adding extra vote where A is ranked above B changes
% unique winner from A to B
% Need to consider all possible votes
[n, m] = size(RC);
ballot = perms(1:m);

% If any row vector appended to RC causes no-show paradox
% return 1

% winnera is the winner in the election assuming the n+1 
% voter abstains
Abstain = ElectionSimRC(RC, etype);
[maxVotesa, winnera] = max(Abstain);
if sum(Abstain == maxVotesa) > 1
        winnera = 0;
        Votes = zeros(1,6);
     return 
end
RCp = RC;

% winnerp(i) is the ith winner in the election assuming the
% n+1 voter votes the ith permutation in ballot(i)
for i = 1:factorial(m)
    RC = RCp;
    RC(n+1,:) = ballot(i,:);
    RC(n+2,:) = ballot(i,:);
    Participate = ElectionSimRC(RC,etype);
    [maxVotes, winnerp(i)] = max(Participate);
    if sum(Participate == maxVotes) > 1
        winnerp(i) = 0;
    end
    
    %winneraidx is the indices of the winner in the ballot
    %so the ith iteration of losera is the abstain election losers
    %that are ranked below the winnera
    %thus if losera becomes winnerp, no-show has occured
    winneraidx = find(ballot(i,:)==winnera);
    if winneraidx == m
        losera = NaN;
    else
        losera = ballot(i,winneraidx+1:m);
    end
    
    ParadoxPossible(i) = sum(winnerp(i)==losera);
    
    
end
end

Votes = ParadoxPossible






%RC = [1 2 3; 3 1 2; 2 3 1; 1 3 2; 3 1 2]

