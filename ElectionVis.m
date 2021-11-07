function [Votes] = ElectionVis(voters, candidates, type)
%% Important Variables
%pdist2 finds the distances between EISI and mu, thus
%d is the distance matrix where each row is how far 
%the ith voter is from each candidate
d = pdist2(voters, candidates);

%Min is the minimum value in each row and I is the index
%of the minimum value in each row
%Thus I can be intepreted as the vector containing the
%candidate that is closest to the voter's views
[Min, I] = min(d');

%dsort is the sorted matrix d such that each row is
%sorted from least to greatest
%RC is the ranked choice matrix, it contains the
%indices of the sorted distances, i.e. it orders a
%voters candidates from least to most distance
[dsort, RC] = sort(d,2);

%VotesPN is 1x4 array of number of votes for each cand
VotesPN = [];
for j = 1:length(candidates)
    VotesPN(j) = sum(I == j);
end

% Votes starts with counting the simple plurality results
% thus setting Votes = VotesPN
Votes = VotesPN;
VotesR = VotesPN;

% RCt is the transpose of RC
% rowc and colc are the size of the rows and columns of the
% candiadates
% m is row-1 since the maximum possible iterations of 
% runoff needed is one less than the total number of 
% candidates
[rowc, colc] = size(candidates);
m = rowc-1;

if colc == 2
votersx = voters(:,1);
votersy = voters(:,2);
candidatesx = candidates(:,1);
candidatesy = candidates(:,2);
end

if colc == 3
votersx = voters(:,1);
votersy = voters(:,2);
votersz = voters(:,3);
candidatesx = candidates(:,1);
candidatesy = candidates(:,2);
candidatesz = candidates(:,3);
end

colors = input('Please list colors to represent candidates:');

% rgb function sourced from
% Author: Ben Mitch
% Website: MATLAB Central File Exchange
% Link: https://www.mathworks.com/matlabcentral/fileexchange/1805-rgb-m

for i = 1:length(colors)
   colorInd(i,:) = rgb(colors(i));
end
          
%% Plurality Algorithm
tic
if strcmp(type, 'plurality') == 1 | strcmp(type, 'main') == 1 | strcmp(type, 'runoff') == 1 


    
    VotesPN

% Scatter Plurality
hold on
grid on
for j = 1:length(candidates)
    [~, voteridx] = find(I==j);
    if colc == 2
    scatter(votersx(voteridx), votersy(voteridx),25, colorInd(j,:), 'filled')
    scatter(candidates(j,1), candidates(j,2), 80, colorInd(j,:), 'filled', 'd')
    title('Plurality')
    elseif colc == 3
    view(3)
    scatter3(votersx(voteridx), votersy(voteridx), votersz(voteridx), 25, colorInd(j,:),'filled')
    scatter3(candidatesx(j), candidatesy(j), candidatesz(j), 80, colorInd(j,:), 'filled', 'd')
    title('Plurality')
    end
end
hold off 
grid off

end
toc
%% Instant Runoff Algorithm
tic
if strcmp(type, 'runoff') == 1 | strcmp(type, 'main') == 1
    
% Variable k determines how many iterations are needed
% to complete the runoff algorithm
for k = 1:m




% Winner is the maximum number of votes a candidate
% received from plurality
Winner = max(VotesR);

% Loser is the least number of votes a candidate received
% from plurality and after each iteration of runoff
% Uses Votessort to account for eliminating candidates
% by giving them 0 votes
Votessort = sort(VotesR);
Loser = find(Votessort(k)==VotesR);
%Loser = find(min(Votes(Votes>0)))
[rowl, coll] = size(Loser);
for i = 1:coll
    RCt = RC';
    RC = RCt(RCt~=Loser(i));
end

    

% RC recalculates the ranked choice matrix after each 
% iteration
% After the loser is eliminated, all votes for loser
% are removed from RC while maintaing the correct 
% order of votes
% m is used as the number of rows since each iteration
% will have one less candidate as losers are eliminated
RCt = RC';

RC = reshape(RC, [] , length(voters)).';


% If the Winner has a majority of votes, then instant 
% runoff algorithm is ended and the candidate with the
% majority wins the election; display votes
if Winner/length(voters) > 0.5
    
    VotesR;
    
break
    
% If the winner doe not have a majority, instant runoff
% occurs
else

% This for loop just counts all the first choice votes 
% for a candidate and recalculates Votes array
for j = 1:length(candidates)
     VotesR(j) = sum(RC(:,1) == j);
end
end
VotesR
m = m-1;
figure(k+1)
hold on
grid on
for p = 1:length(candidates)
    [~, voteridx] = find(RC(:,1).'==p);
    if colc ==2
    scatter(votersx(voteridx), votersy(voteridx),25, colorInd(p,:), 'filled')
    scatter(candidates(p,1), candidates(p,2), 80, colorInd(p,:), 'filled', 'd')
    title(sprintf('Instant Runoff Iteration %d', k))
    elseif colc == 3
    view(3)
    scatter3(votersx(voteridx), votersy(voteridx), votersz(voteridx), 25, colorInd(p,:), 'filled')
    scatter3(candidatesx(p), candidatesy(p),candidatesz(p), 80, colorInd(p,:), 'filled', 'd')
     title(sprintf('Instant Runoff Iteration %d', k))
    end
end
grid off
hold off
end
end
toc
%% Approval Algorithm
if strcmp(type, 'approval') == 1 | strcmp(type, 'main') == 1
tic
r = (max(voters,[], 'all') - min(voters, [], 'all'))/4;
VotesA = sum(d <= r)
LV = (d <= r);


% Scatter Approval
for i = 1:length(LV)
    inrange = find(LV(i,:)==1);
    if isempty(inrange) == 1
    newcolor(i,:) = [0, 0, 0];
    else
    newcolor(i,:) = mean(colorInd(inrange,:),1);
    end
end

figure
hold on

if colc == 2
scatter(votersx, votersy, 25, newcolor, 'filled')

for k = 1:length(candidates)
    viscircles(candidates(k,:),r,'color', colorInd(k,:));
    scatter(candidates(k,1), candidates(k,2), 80, colorInd(k,:), 'filled', 'd')
end
title('Approval')
elseif colc == 3
    view(3)
    scatter3(votersx, votersy, votersz, 25, newcolor, 'filled')
    [X, Y, Z] = sphere;
for k = 1:length(candidates)
    %surf((X*r)+candidatesx(k),(Y*r)+candidatesy(k), (Z*r)+candidatesz(k),'FaceColor', colorInd(k));
    scatter3(candidatesx(k), candidatesy(k), candidatesz(k), 80, colorInd(k,:), 'filled', 'd')
end
title('Approval')
end


%xlim([-50, 50]);
%ylim([-50, 50]);
hold off
grid on
end

toc

%% Score Voting Algorithm
tic
if strcmp(type, 'score') == 1 | strcmp(type, 'main') == 1

% define score range for voters to assign to candidates
score_range = 1:10;
% maximum radii from candidates
great_range = (max(voters,[], 'all') - min(voters, [], 'all'))/2;
% amount radius decreases each iteration from max radii
inc_r = great_range/length(score_range);

% calculate score matrix to determin score for all cands for all voters
score_matrix = (d<=great_range);
score_matrix = double(score_matrix);
score_matrix(score_matrix==0) = NaN;

score_rngidx = 1;
for i = great_range:-inc_r:inc_r
score_matrix(d<=i) = score_range(score_rngidx);
score_rngidx = score_rngidx+1;
end

score_matrix(isnan(score_matrix)) = min(score_range);
% standardize score matrix on scale from 0 to 1
st_score_matrix = score_matrix./length(score_range);

% calculate mean votes
VotesS = mean(score_matrix)

% Scatter Score Voting

% calculate new color for all voters based on distances to candidates
for i = 1: length(voters) 
colormix = colorInd.*st_score_matrix(i,:).';
colorsum = sum(colormix);
nonzeros = sum(colormix~=0);
newcolor(i,:) = colorsum./nonzeros;
end

% create circles
total_circles = length(score_range)*length(candidates);
r = inc_r;
inc_ocolor = 1/length(score_range);
ocolor = 0;
hold on
figure

if colc == 2
    scatter(votersx, votersy, 25, newcolor, 'filled' )
    title('Score')
    
    %for k = 1:total_circles
   %modvalue = mod(k-1,length(candidates))+1;
    %newocolor = colors(modvalue,:)-ocolor;
    %newocolor(newocolor<0)=0;
    %viscircles(candidates(modvalue,:),r,'color', newocolor);
    %if mod(k, length(candidates)) == 0
    %r = r+inc_r;
    %ocolor = ocolor+inc_ocolor;
    %end
    
   % end

elseif colc == 3
    view(3)
    scatter3(votersx, votersy, votersz, 25, newcolor, 'filled')
    title('Score')
end
%xlim([m, 50]);
%ylim([-50, 50]);
hold off
end
toc

%% Condorcet Algorithm

if strcmp(type, 'condorcet') == 1 | strcmp(type, 'main') == 1

% h2h_no is the total number of head to head matchups between candidates
% h2h is a matrix listing all possible head to head matchups
h2h_no = nchoosek(length(candidates), 2);
h2h = nchoosek(1:length(candidates), 2);

% candlist is a vetor listing each candidates number
candlist = 1:length(candidates);

% i for loop runs through each h2h_itertion, which is each individual row
% of h2h, which corresponds to each individual matchup
% RC matrix is re-defined at the beginning of each iteration, then all
% candidates except for the 2 being considered are removed from RC in the j
% loop
for i = 1:h2h_no
    
    h2h_iteration = h2h(i,:);
    [dsort, RC] = sort(d,2);
    RCt =RC.';
    exclude = setdiff(candlist, h2h_iteration);
    
    for j = 1:length(exclude)
        RCt = RC.';
        RC = RCt(RCt ~= exclude(j));
    end

RC = reshape(RC, [], length(voters)).';

% Vote_count finds the number of votes for cand A vs cand B in each
% individual matchup to determine the winner of the matchup
Vote_count = [];
    for k = 1:length(candidates)
        Vote_count(k) = (sum(RC(:,1) == k));
    end
    
Vote_count;

% Graphs each matchup between voters
figure
hold on
grid on
    for p = 1:length(candidates)
        [~, voteridx] = find(RC(:,1).'==p);
        if colc == 2
            scatter(votersx(voteridx), votersy(voteridx),25, colorInd(p,:), 'filled')
            %scatter(candidates(p,1), candidates(p,2), 80, colorInd(p,:), 'filled', 'd')
            title(sprintf('Condorcet Pairwise Case %d', i));
        elseif colc == 3
            view(3)
            scatter3(votersx(voteridx), votersy(voteridx), votersz(voteridx), 25, colorInd(p,:), 'filled')
            %scatter3(candidatesx(p), candidatesy(p), candidatesz(p), 80, colorInd(p,:), 'filled', 'd')
            title(sprintf('Condorcet Pairwise Case %d', i));
        end
    end
grid off
hold off

% winner finds the winner of each matchup
[~, winner(i)] = max(Vote_count);

% Votes tallys the total wins for each candidate from all matchups
    for m = 1:length(candidates)
        VotesC(m) = sum(winner==m);
    end
end
VotesC
end
   
   
   
   
  




    
