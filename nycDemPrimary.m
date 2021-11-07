N = readmatrix('nyc2021DemMayorPrimaryDataOct28.csv');
Cnames = readtable('nyc2021DemMayorPrimaryCandidates.csv');
C = readmatrix('nyc2021DemMayorPrimaryCandidates.csv');

%%
[numVoters, numRankings] = size(N)

candID = C(:,1);
numCands = length(candID);
Nnormal = N;
for i = 1:numCands
    Nnormal(N==candID(i))=i;
end

Nnormal(Nnormal == 111111) = 0;

RC = Nnormal;
%%
for i = 1:length(RC)
RCrow = RC(i,:);
[a,c] = unique(RCrow,'first');
RCrow = zeros(size(RCrow));
RCrow(c) = a;

RCrow = move_me(RCrow);

RC(i,:) = RCrow;
sprintf('Did %d', i)

end
%%
ValidRanking = zeros(1,length(RC));
for i = 1:length(RC)
    RCrow = RC(i,:);
    
    nonzeroVotes = RCrow(~~RCrow);
    firstZero = min(find(RCrow==0));
    
    candsRanked = RCrow(RCrow~=0);
    
    if isempty(candsRanked)
        continue
    end
    
    lastRankedCand = candsRanked(end);
    lastRankedIdx = find(lastRankedCand==RC(i,:));
       
    if length(nonzeroVotes) ~= length(unique(nonzeroVotes))
       
        
    elseif isempty(firstZero)
            ValidRanking(i) = 1;
            
    elseif firstZero > lastRankedIdx
        ValidRanking(i) = 1;
                     
    end
    
    %sprintf('Did %d',i);
end

%% Types of Invalid Votes
% 0 before a ranking:[1 2 3 0 4], [0 1 2 3 4], [1 0 2 3 0]
% 2 or more of same candidate: [1 1 2 3 4], [1 1 1 1 1], [1 2 2 3 3]

invalidVote = find(ValidRanking == 0)
RC(invalidVote,:) = [];
RC = sortrows(RC);

%%
Nfinal = RC;


%%
candScale = [-2 2; -4 -1; -4 -3; -3 4; -2 1; -1 0; -5 5; -5 -5; -3 5; -5 2; -3 0; -3 1; -4 -2];
mu = candScale;
muShape = [-5 -5;-5 2;-5 5; -3 5; -3 4;-2 2; -1 0; -2 1; -3 1;-3 0; -4 -1; -4 -2; -4 -3]; 

NnormalSort = sortrows(Nnormal);

