clear all
clc
clf


nAgents = 1000;
gridSize = 100;

moveProb = 0.7;
reproductionNumber = 10:10:120;
infectRates = [0.9:-0.1: 0.1];
initInfected = 0.01;

nRepeats = 10;
maxSteps = 50000;

finalRecovered = zeros(length(reproductionNumber),length(infectRates));

agentStatus = zeros(nAgents,3); % Holds x,y and susceptible/infected/resistant status for each agent

gridStatus = cell(gridSize);    % Lookup table for which agents are at each grid point

saveStatus = zeros(maxSteps,3);   % Save # susceptible/infected/recoverd for each iteration

colorCode = ['b' 'r' 'g'];

for iRecoverRate = 1:length(infectRates)
    infectProb = infectRates(iRecoverRate);
    for iInfectRate = 1:length(reproductionNumber)
        recoverProb = infectProb / reproductionNumber(iInfectRate);
        for iRepeat = 1:nRepeats
            disp(['infect ' num2str(infectProb) ' recover ' num2str(recoverProb) ' repeat ' num2str(iRepeat) ])
            agentStatus(:,3) = 1;   % All are susceptible

            % Initialize agent positions
            for iAgent = 1:nAgents
                agentStatus(iAgent,1:2) = [ceil(rand()*gridSize), ceil(rand()*gridSize)];
            end

            agentDistance = sqrt(sum((agentStatus(:,1:2)-gridSize*0.5).^2,2));
            [~, distanceList] = sort(agentDistance);

            for iAgent = 1:round(nAgents*initInfected)  % Set some of the pop to infected
                agentStatus(distanceList(iAgent),3) = 2;
            end


            zeroInfected = false;
            iStep = 0;
            while ~zeroInfected
                iStep = iStep + 1;
                % Update position
                for iAgent = 1:nAgents
                    if rand() < moveProb;
                        change = [round(rand())+1 round(2*rand()-1)]; % First number selects x or y, second selects in- or decrease
                        agentStatus(iAgent,change(1)) = agentStatus(iAgent,change(1)) + change(2);
                        agentStatus(iAgent,1:2) = mod(agentStatus(iAgent,1:2)-1,gridSize)+1;    % Periodic boundaries
                    end
                end

                % Update grid contents
                gridStatus = cell(gridSize);
                for iAgent = 1:nAgents
                    x = agentStatus(iAgent,1);
                    y = agentStatus(iAgent,2);
                    gridStatus{x,y} = [ gridStatus{x,y} iAgent ];   % Add agents to the grid
                end

                infectedAgents = agentStatus(agentStatus(:,3)==2,:);
                for iInfected = 1:size(infectedAgents,1)
                    x = infectedAgents(iInfected,1);
                    y = infectedAgents(iInfected,2);
                    for iOccupant = 1:length(gridStatus{x,y})
                        if agentStatus(gridStatus{x,y}(iOccupant),3) == 1
                                if rand() < infectProb
                                    agentStatus(gridStatus{x,y}(iOccupant),3) = 2;
                                end
                        end
                    end
                end

                recover = rand(nAgents,1);
                agentStatus(agentStatus(:,3)==2 & recover < recoverProb ,3) = 3;

                saveStatus(iStep,1) = sum(sum(agentStatus(:,3)==1));
                saveStatus(iStep,2) = sum(sum(agentStatus(:,3)==2));
                saveStatus(iStep,3) = sum(sum(agentStatus(:,3)==3));

                if saveStatus(iStep,2) == 0 | iStep == maxSteps
                    zeroInfected = true;
                    disp(['Steps ' num2str(iStep)]);
                end
            end
            saveStatus = saveStatus / nAgents;

            finalRecovered(iInfectRate,iRecoverRate) = finalRecovered(iInfectRate,iRecoverRate) + saveStatus(iStep,3);
        end
        finalRecovered(iInfectRate,iRecoverRate) = finalRecovered(iInfectRate,iRecoverRate) / nRepeats;
    end    
end
% hold off

plot(finalRecovered,'-x')
%title(['$$d = ' num2str(moveProb) ', \beta = ' num2str(infectProb) ', \gamma = ' num2str(recoverProb) '$$'],'Interpreter','latex');
hold off
filename = ['t4a' num2str(tic) '.mat'];
save(filename,'finalRecovered')

disp('Done!')

%% load and display

%load t4a7403758257454.mat
load t4a_full.mat
reproductionNumber = 5:5:140;
surf(infectRates, reproductionNumber, finalRecovered)
