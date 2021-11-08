function [dist, RC2] = VotingDistWhile(RC,CI1, CI2,k1, k2)
tic
%% 0. Define Important Variables
[n1,~] = size(RC);
[m1, ~] = size(CI1);
[m2, ~] = size(CI2);
x = linspace(-5, 5, 1000);

%RC is n1 x m1 matrix where each row represents voters ranking
%CI1 is m1 x 3 where each row contains the ideology values of the ith
%candidate
%CI2 is m2 x 3

%% 1. Plot m1 polygon and create m1C2 bisectors
%Need to plot all CI1 candidates on coordinate plane
%Randomly generate voters within sectors according to RC matrix

%Create polygon connecting m1 candidates
f1 = figure;
pgon1 = polyshape(CI1);
plot(pgon1)

%Find all pairs of candidates and create profile (all possible votes)
pair1 = nchoosek(1:m1,2);
[pairNum1, ~] = size(pair1);
profile1 = sortrows(combinator(m1, k1, 'p'));
isPermUsed1 = ismember(profile1, RC, 'rows');
PermsUsed1 = sum(isPermUsed1);
maxColors1 = distinguishable_colors(9000, 'w');

%Creat perpendicular bisectors
for i = 1:pairNum1
    [y(i,:), mperp(i), b(i)] = perpbi(CI1(pair1(i,1),:), CI1(pair1(i,2),:));
    eqn = line(x, y(i,:), 'LineWidth', 2);
end

%% 2. Generate voters within sectors accoring to RC
for i = 1:n1
   
    % Generate a random point on the plane, and calculate the vote perm
    % it corresponds to
        randpoint = (10*rand(1,2))-5;
        d = pdist2(CI1, randpoint);
        [dsort, randVote] = sort(d);
        randVote = randVote';
        randVote = randVote(1:k1)
        numZeros = sum(RC(i,:) == zeros(1,k1));
    
        %Generate the random point until the randomVote is the same as
        %the actual vote in the ith row of RC
        %Account for non-full votes by subtracting off the numer of zeros
        while sum(RC(i,:) == randVote) ~= k1-numZeros
            randpoint = (10*rand(1,2))-5;
            d = pdist2(CI1, randpoint);
            [dsort, randVote] = sort(d);
            randVote = randVote';
            randVote = randVote(1:k1);
        end
        
        sprintf('Finished %d', i)
        
        randVoteMat(i,:) = randVote;
        randpointMat(i,:) = randpoint;
        rgb(i,:) = maxColors1(mod(find(ismember(profile1, randVote, 'rows')),9000),:);
    
end
hold on
scatter(randpointMat(:,1), randpointMat(:,2),50, rgb, 'filled')
xlim([-6, 6])
ylim([-6, 6])

%% 3. Create m2 polygon and m2C2 bisectors

%Create m2 polygon
f2 = figure;
pgon2 = polyshape(CI2);
plot(pgon2)


%Find all pairs of candidates and create profile (all possible votes)
pair2 = nchoosek(1:m2,2);
[pairNum2, ~] = size(pair2);
profile2 = sortrows(combinator(m2, k2, 'p'));

%Create new Sectors
for i = 1:pairNum2
    [y(i,:), mperp(i), b(i)] = perpbi(CI2(pair2(i,1),:), CI2(pair2(i,2),:));
    eqn = line(x, y(i,:), 'LineWidth', 0.5);
end

%% 4. Dorrrv

%Find RC2 for new election based on previous plot
d = pdist2(randpointMat, CI2);
[~, RC2] = sort(d,2);
RC2 = RC2(:,1:k2);

for i = 1:n1
    permForVoter(i,:) = find(ismember(profile2, RC2(i,:), 'rows'));
end

possiblePerms2 = 1:length(profile2);
isPermUsed2 = ismember(possiblePerms2, permForVoter);
PermsUsed2 = sum(isPermUsed2);
maxColors2 = distinguishable_colors(PermsUsed2, 'w');

sortUniquePerms = sort(unique(permForVoter));
permForVoterNormal = permForVoter;
for i = 1:length(sortUniquePerms)
    permForVoterNormal(permForVoter==sortUniquePerms(i))=i;
end

%maxColors2 = distinguishable_colors(500, 'w');
rgb2 = maxColors2(permForVoterNormal,:);

hold on       
scatter(randpointMat(:,1), randpointMat(:,2), 50, rgb2, 'filled');
xlim([-6, 6])
ylim([-6, 6])
for j = 1:length(profile2)
    dist(j) = sum(permForVoter == j)/n1;
end
   
toc   
end