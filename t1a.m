clear all
clc
clf

nAgents = 1;
gridSize = 100;

moveProb = 0.9;

nSteps = 1000;

agentStatus = zeros(nAgents,3); % Holds x,y and susceptible/infected/resistant status

gridStatus = cell(gridSize);

iAgent = 1; % Unneccessary for now, but will be needed later so might as well...
agentStatus(iAgent,1:2) = [ceil(rand()*gridSize), ceil(rand()*gridSize)];

for iStep = 1:nSteps
    % Update position
    for iAgent = 1:nAgents
        if rand() < moveProb;
            change = [round(rand())+1 round(2*rand()-1)]; % First number selects x or y, second selects in- or decrease
            agentStatus(iAgent,change(1)) = agentStatus(iAgent,change(1)) + change(2);
            agentStatus(iAgent,1:2) = mod(agentStatus(iAgent,1:2),gridSize);    % Periodic boundaries
        end
    end
    % Plot
    color = [iStep/nSteps 0 (nSteps-iStep)/nSteps];
    hold on
    plot(agentStatus(iAgent,1),agentStatus(iAgent,2),'.','Color',color)
    axis([0 gridSize 0 gridSize])
    
end

hold off
