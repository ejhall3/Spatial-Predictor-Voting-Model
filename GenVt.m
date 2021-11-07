%% 0. Plot

tic
RCorig = RC;
RC = sortrows(RC);
[n1,k1] = size(RC);
m1 = length(mu);
x = linspace(-5, 5, 1000);

%RC is n1 x m1 matrix where each row represents voters ranking
%mu is m1 x 2 where each row contains the ideology values of the ith
%candidate


%Need to plot all CI1 candidates on coordinate plane
%Randomly generate voters within sectors according to RC matrix

%Create polygon connecting m candidates
pgon = polyshape(mu);
plot(pgon)

% Create Sectors by drawing perpendicular bisectors
pair = nchoosek(1:m1,2);
[pairNum, ~] = size(pair);
hold on

for i = 1:pairNum
    [y(i,:), mperp(i), b(i), c(i)] = perpbi(mu(pair(i,1),:), mu(pair(i,2),:));
    eqn = line(x, y(i,:), 'LineWidth', 1);
    xlim([-5 5])
    ylim([-5 5])

end

% yfunc is represntation of y in slope intercept form
%1st column is slope and 2nd column is y-intercept
%i-th row is i-th y function corresponding to line segment in pair
yfunc = [mperp; b; c]';
numLn = length(yfunc);

Bfunc = [0 5 NaN; 0 -5 NaN; Inf Inf 5; Inf Inf -5];

 %% 1. Create boundaries
 xlimit = [-5, 5];
 ylimit = [-5, 5];
 xbox = xlimit([1 1 2 2 1]);
 ybox = ylimit([1 2 2 1 1]);
 mapshow(xbox, ybox)
 
 %% 2. Find ALL Intersection Points

pair2 = nchoosek(1:length(pair),2);

xi = zeros(1,length(pair2));
yi = zeros(1,length(pair2));

%Lines Interect with themselves (m1C2 C 2)
for i = 1:length(pair2)
    
[xi(i), yi(i)] = linEqnx(yfunc(pair2(i,1),:), yfunc(pair2(i,2),:));
   
end

%Lines Intersect with boundaries (4*m1C2)
for i = 1:length(pair)
    for j = 1:length(Bfunc)
        [xi(end+1), yi(end+1)] = linEqnx(yfunc(i,:), Bfunc(j,:));
    end
end

%Boundaries Intersect with themselves (4)
xi(end+1:end+4) = [xbox(1), xbox(2), xbox(3), xbox(4)];
yi(end+1:end+4) = [ybox(1), ybox(2), ybox(3), ybox(4)];
         

xiyi = [xi;yi]';

%Remove points not in box
inBox = xiyi <= 5 & xiyi >= -5;
pointInBox = inBox(:,1) & inBox(:,2);
pointNotInBox = ~pointInBox;
xiyi(pointNotInBox,:) = NaN;

%NOW: Find intersections needed for each voter



%% 3. Analyze RC
candList = 1:m1;


%Find the last ranked candidate for each voter
for i = 1:n1
    RCrow = RC(i,:);
    candsRanked = RCrow(RCrow~=0);
    lastRankedCand(i) = candsRanked(end);
    sprintf('Completed %d',i);
    
    % Find lastRankedCand combination with every non-ranked candidate
    notRankedCands = candList(~ismember(candList, RC(i,:)));
    
    notRankedCandsCell{i} = notRankedCands;
    
    %Find Intersections Corresponding to assumption that lastRankedCand
    %is above all notRankedcands
    %lastAboveNot = [zeros(1,length(notRankedCands))+lastRankedCand(i); notRankedCands]';

    
end


%% 4. Determine which intersections are in the sector for each voter
%% 5. Find Vertices for each voter and generate random points

d = pdist2(mu, xiyi);

%4a Find points satisfying lastRanked above NotRanked

numFullVotes = sum(~~RC,2);
lastAboveNot = ones(n1, length(d));

for i = 1:n1
    
    if i > 1 & isequal(RC(i,:), RC(i-1,:))
        lastAboveNot(i,:) = lastAboveNot(i-1,:);
        continue
    end
    
    notRankedCandsFori = notRankedCandsCell{i};
    
    if k1 == m1 && numFullVotes(i) == k1
        lastAboveNot(i,:) = ones(1, length(d));
        
    else

    for j = 1:length(d)
    
    %criteria1 = zeros(1,length(notRankedCandsFori));
        for k = 1:length(notRankedCandsFori)
            %criteria1 tells if the jth point in xiyi
            %satisfies the condition that the ith
            %voter's last ranked cand is above each
            %of the k non-ranked cand
           criteria1 = d(lastRankedCand(i),j) < d(notRankedCandsFori(k),j) ...
                | isalmost(d(lastRankedCand(i),j), d(notRankedCandsFori(k),j), 0.00001);
            if criteria1 == 0
                lastAboveNot(i,j) = 0;
                break
            end
        end
        
            %lastAboveNot tells if the jth point in xiyi
            %satisfies all conditions that the lastRankedCand
            %is above all nonRankedCands for voter i           
        %if sum(criteria1) == length(notRankedCandsFori)
            
        %end
        
    end
    
    end
 sprintf('Completed %d',i);
end
    
%% Find Points Satisfying the k-rankings
rankedCandsInOrder = ones(n1, length(d));
for i = 1:n1
    
    if i > 1 && sum(RC(i,:) == RC(i-1,:)) == k1
        rankedCandsInOrder(i,:) = rankedCandsInOrder(i-1,:);
        continue
    end
    
    numFullVotes = sum(~~RC(i,:));
    if numFullVotes == 1
        rankedCandsInOrder(i,:) = lastAboveNot(i,:);
        continue
    end
    
    for j = 1:length(d)
        %criteria2 tells if the jth point in xiyi
        %satisfies the condition that all k ranked 
        %pairs of candidates are ordered accordingly for voter i
         
        for k = 1:numFullVotes-1
                       
            criteria2 = d(RC(i,k),j) < d(RC(i,k+1),j) ...
                | isalmost(d(RC(i,k),j), d(RC(i,k+1),j), 0.00001);
            
            if criteria2 == 0
                rankedCandsInOrder(i,j) = 0;
                break
            end
            
        end
        
    end
    
    sprintf('Completed %d',i);
end
            
%% 5. Find Vertices for each voter and generate random points

xiyiInSec = lastAboveNot & rankedCandsInOrder;
%xiyiInSec = rankedCandsInOrder;


%NOTE: CONSIDER X=C CASE LATER
yBfunc = [yfunc;Bfunc];
lastAboveNot2 = ones(n1, length(d2));
rankedCandsInOrder2 = ones(n1, length(yBfunc));
votersNoHome = zeros(1,n1);
for i = 1:n1
    sprintf('%d',i)
    numFullVotes = sum(~~RC(i,:));
    
    vertForVoter = find(xiyiInSec(i,:)==1);
    
    if isempty(vertForVoter)
        xArr(i) = NaN;
        yArr(i) = NaN;
        
        continue
        
    end
        
    vertForVoterCell{i} = vertForVoter;
    
    pointsForVoter = xiyi(vertForVoter,:);
    
    xmin = min(pointsForVoter(:,1));
    xmax = max(pointsForVoter(:,1));
    
    xrand = (xmax-xmin)*rand + xmin;
    
    for j = 1:length(yBfunc)
    yPossible(j) = linEqn(yBfunc(j,1), yBfunc(j,2), xrand);
    end
    
   possiblePoints(1:length(yBfunc),1) = xrand;
   possiblePoints(1:length(yBfunc),2) = yPossible;
   
   d2 = pdist2(mu, possiblePoints);
   
%Find lastAbovenot2
notRankedCandsFori = notRankedCandsCell{i};

if k1 == m1 && numFullVotes == k1
    lastAboveNot2(i,:) = ones(1,length(d2));
    
else
    
    for j = 1:length(d2)
    
   
        for k = 1:length(notRankedCandsFori)
            %criteria1 tells if the jth point in xiyi
            %satisfies the condition that the ith
            %voter's last ranked cand is above each
            %of the k non-ranked cand
           criteria1 = d2(lastRankedCand(i),j) < d2(notRankedCandsFori(k),j) ...
                | isalmost(d2(lastRankedCand(i),j), d2(notRankedCandsFori(k),j), 0.00001);
            
            if criteria1 == 0
                lastAboveNot2(i,j) = 0;
                break
            end
        end
        
            %lastAboveNot tells if the jth point in xiyi
            %satisfies all conditions that the lastRankedCand
            %is above all nonRankedCands for voter i           
        
        
    end

end
   
   
   %Find rankedCandsinOrder2
   for j = 1:length(d2)
        %criteria2 tells if the jth point in xiyi
        %satisfies the condition that all k ranked 
        %pairs of candidates are ordered accordingly for voter i
        for k = 1:numFullVotes-1
                       
            criteria2(k) = d2(RC(i,k),j) < d2(RC(i,k+1),j) ...
                | isalmost(d2(RC(i,k),j), d2(RC(i,k+1),j), 0.00001);
            
            if criteria2 == 0
                rankedCandsInOrder2(i,j) = 0;
                break
            end
            
        end
        

        
   end
   
   lastAboveNot2 = logical(lastAboveNot2);
   rankedCandsInOrder2 = logical(rankedCandsInOrder2);
   possiblePointsInSec = lastAboveNot2 & rankedCandsInOrder2;
   
   yInSec = yPossible(possiblePointsInSec(i,:));
   
   %Remove points not in box
    inBox2 = yInSec >= xlimit(1) & yInSec <= xlimit(2);
    yInSec = yInSec(inBox2);
    
    ymax = max(yInSec);
    ymin = min(yInSec);
    
    yrand = (ymax-ymin)*rand + ymin;
    
    xArr(i) = xrand;
    yArr(i) = yrand;
    

end
   
xyrandMat = [xArr; yArr]'

[PermsUsed, ~] = size(unique(RC, 'rows'));
rgb = distinguishable_colors(PermsUsed, 'w');
RCunique = unique(RC, 'rows');

%% Find Colors

for i = 1:n1
    for j = 1:PermsUsed
    criteria3 = sum(RC(i,:) == RCunique(j,:));
    if criteria3 == k1
        RCcolor(i) = j;
    end
    end
end
%%

hold on
scatter(xArr, yArr, 10, rgb(RCcolor,:), 'filled') 
toc

