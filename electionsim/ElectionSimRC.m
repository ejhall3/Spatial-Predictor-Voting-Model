function [Votes] = ElectionSimRC(RCex,type)
%% Plurality Algorithm

%RC is nxm matrix, where n is #voters, m is #candidates
%thus, each row represents each voters ranking of cands

% rowc and colc are the size of the rows and columns of the
% candiadates
% m is colc-1 since the maximum possible iterations of 
% runoff needed is one less than the total number of 
% candidates

[n, ~] = size(RCex);
m = max(max(RCex));
runoff_iterations = m-2;

if strcmp(type, 'plurality') == 1 | strcmp(type, 'main') == 1 | strcmp(type, 'runoff') == 1 


    %Vp is the plurality votes vector where the ith entry is the
    %number of plurality votes for the ith candidate
    Vp = zeros(m, 1);
    for j = 1:m
        Vp(j) = sum(RCex(:,1) == j);
    end

    Votes = Vp;
end



%% Instant Runoff Algorithm
j = 1;

if strcmp(type, 'runoff') == 1 | strcmp(type, 'main') == 1
    
    Vir = Vp;
    
    while j <= runoff_iterations

        % numWinner is the maximum number of votes a candidate
        % received from plurality
        numWinner = max(Vir);

        % Loser is the least number of votes a candidate received
        % from plurality and after each iteration of runoff
        % Uses Votessort to account for eliminating candidates
        % by giving them 0 votes
        Votessort = sort(Vir);
        Loser = find(Votessort(j)==Vir);
        
        %Assign all Losers a value of 0
        [rowl, ~] = size(Loser);
        for i = 1:rowl
                RCex(RCex==Loser(i))=0;
        end
        
        %Find rowz, the rows in RC with a 0
        [rowz, colz] = find(RCex==0);
        zeroidx = [rowz, colz];
        rowz = unique(rowz);
        
        %Move all zeros to the end of the row
        for i = 1:length(rowz)
            RCex(rowz(i),:) = move_me(RCex(rowz(i),:));
        end

            % If the Winner has a majority of votes, then instant 
            % runoff algorithm is ended and the candidate with the
            % majority wins the election; display votes
                if numWinner/sum(Vir) > 0.5

                    Vir;

                break

                % If the winner doe not have a majority, instant runoff
                % occurs
                else

                % This for loop just counts all the first choice votes 
                % for a candidate and recalculates Votes array
                    for p = 1:m
                         Vir(p) = sum(RCex(:,1) == p);
                    end
                end
            Vir
            
            %Increase j by the size of loser to find next candidate(s)
            %with least votes that hasnt been eliminated yet
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

RCc = RCex;

for i = 1:h2h_no
    
    h2h_iteration = h2h(i,:);
    RCex = RCc;
    RCt =RCex.';
    exclude = setdiff(candlist, h2h_iteration);
    
    for j = 1:length(exclude)
        RCt = RCex.';
        RCex = RCt(RCt ~= exclude(j));
    end

RCex = reshape(RCex, [], n).';

% Vote_count finds the number of votes for cand A vs cand B in each
% individual matchup to determine the winner of the matchup
Vote_count = [];
    for k = 1:m
        Vote_count(k) = (sum(RCex(:,1) == k));
    end
    
Vote_count;

% winner finds the winner of each matchup
[~, winner(i)] = max(Vote_count);

% Votes tallys the total wins for each candidate from all matchups
    for i= 1:m
        Vc(i) = sum(winner==i);
    end
end
Votes = Vc
end