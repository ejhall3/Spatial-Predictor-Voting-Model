function [Votes] = ElectionSim(n,m,k,type, dist)
if strcmp(dist, 'IAC')
    
    [~, RC] = sort(rand(n,m),2);
   


else

 %if strcmp(dist, 'dist')
    
    profile = sortrows(combinator(m, k, 'p'));
    numSectors = length(profile);
    
    voterPerms = randsrc(1,n, [1:numSectors; dist]);
    
    %RC = zeros(n, m);
    RC(:,1:k) = profile(voterPerms,:);
    
    for i = 1:n
    residualCand(i,:) = find(~ismember(1:m, RC(i,:)));
    end
    
    
    
end
%% Plurality Algorithm

%RC is nxm matrix, where n is #voters, m is #candidates
%thus, each row represents each voters ranking of cands

% rowc and colc are the size of the rows and columns of the
% candiadates
% m is colc-1 since the maximum possible iterations of 
% runoff needed is one less than the total number of 
% candidates

RCstatic = RC

[~, m] = size(RC);
runoff_iterations = m-2;

if strcmp(type, 'plurality') == 1 | strcmp(type, 'main') == 1 | strcmp(type, 'runoff') == 1 


    %Vp is the plurality votes vector where the ith entry is the
    %number of plurality votes for the ith candidate
    Vp = zeros(m, 1);
    for j = 1:m
        Vp(j) = sum(RC(:,1) == j);
    end

    Votes = Vp
end



%% Instant Runoff Algorithm

if strcmp(type, 'runoff') == 1 | strcmp(type, 'main') == 1
    
    Vir = Vp;
    % Variable k determines how many iterations are needed
    % to complete the runoff algorithm
    for j = 1:runoff_iterations

        % Winner is the maximum number of votes a candidate
        % received from plurality
        Winner = max(Vir);

        % Loser is the least number of votes a candidate received
        % from plurality and after each iteration of runoff
        % Uses Votessort to account for eliminating candidates
        % by giving them 0 votes
        Votessort = sort(Vir);
        Loser = find(Votessort(j)==Vir);
        
        [rowl, ~] = size(Loser);
        for i = 1:rowl
                RC(RC==Loser(i))=0
        end
        
        [rowz, colz] = find(RC==0);
        zeroidx = [rowz, colz];
        rowz = unique(rowz);
        
        for i = 1:length(rowz)
            RC(rowz(i),:) = move_me(RC(rowz(i),:))
        end

            % RC recalculates the ranked choice matrix after each 
            % iteration
            % After the loser is eliminated, all votes for loser
            % are removed from RC while maintaing the correct 
            % order of votes
            % m is used as the number of rows since each iteration
            % will have one less candidate as losers are eliminated
            %RCt = RC';

            %RC = reshape(RC, [] , n).';


            % If the Winner has a majority of votes, then instant 
            % runoff algorithm is ended and the candidate with the
            % majority wins the election; display votes
                if Winner/sum(Vir) > 0.5

                    Vir;

                break

                % If the winner doe not have a majority, instant runoff
                % occurs
                %else

                % This for loop just counts all the first choice votes 
                % for a candidate and recalculates Votes array
                    for j = 1:m
                         Vir(j) = sum(RC(:,1) == j);
                    end
                end
            Vir
            j = j + length(Loser);
                
    end
    Votes = Vir

end  


%% Condorcet Algorithm

if strcmp(type, 'condorcet') == 1 | strcmp(type, 'main') == 1

% h2h_no is the total number of head to head matchups between candidates
% h2h is a matrix listing all possible head to head matchups
h2h_no = nchoosek(m, 2);
h2h = nchoosek(1:m, 2);

% candlist is a vetor listing each candidates number
candlist = 1:m;

% i for loop runs through each h2h_itertion, which is each individual row
% of h2h, which corresponds to each individual matchup
% RC matrix is re-defined at the beginning of each iteration, then all
% candidates except for the 2 being considered are removed from RC in the j
% loop

RC = RCstatic;
RCc = RC;

for i = 1:h2h_no
    
    h2h_iteration = h2h(i,:);
    RC = RCc;
    RCt =RC.';
    exclude = setdiff(candlist, h2h_iteration);
    
    for j = 1:length(exclude)
        RCt = RC.';
        RC = RCt(RCt ~= exclude(j));
    end

RC = reshape(RC, [], n).';

% Vote_count finds the number of votes for cand A vs cand B in each
% individual matchup to determine the winner of the matchup
Vote_count = [];
    for j = 1:m
        Vote_count(j) = (sum(RC(:,1) == j));
    end
    
Vote_count;

% winner finds the winner of each matchup
[~, winner(i)] = max(Vote_count);

% Votes tallys the total wins for each candidate from all matchups
    for i= 1:m
        Vc(i) = sum(winner==i);
    end
end
Votes = Vc;
end
