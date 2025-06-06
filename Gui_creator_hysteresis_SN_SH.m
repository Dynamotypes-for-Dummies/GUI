createGUI();


function createGUI()
% Create the main figure window
addpath('')
fontsize = 9;
mindip = 0;
maxdip = 45;
textdip = 9.5;
maxLength = 120;
%%change name here
addpath('GUI helper files');
fig = uifigure('Position', [20, 50, 1300, 750], 'Name', 'Class SN/SH', 'Color',[1,1,1]);
Height = 735;
% A slider
%%change caption here 
%%to the left, height, length of text box, height of text box
uilabel(fig, 'Position', [250, Height, 500,15], 'Text', 'Model parameters', 'FontWeight', 'bold', 'FontSize', fontsize+1);
Height = Height - 20;
uilabel(fig, 'Position', [50,Height, 400, 22], 'Text', 'Onset curve index: Saddle node (DC+)', 'FontSize', fontsize, 'FontWeight', 'bold');
% Define the segments of the description
descriptionLines = {
    '';
    'Crossing the Saddle Node curve (DC+) will result in onset of a signal with a DC shift. ';

};

% Calculate the number of lines for positioning
combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into sub strings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 %Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)

uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
%% 
ASlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1, 62], 'Value', 1);
Height = Height - maxdip;

%B Slider
uilabel(fig, 'Position', [50,Height, 400, 22], 'Text', 'Offset curve index: Saddle Homoclinic (DC+)', 'FontSize', fontsize, 'FontWeight', 'bold');
descriptionLines = {
'';
    'Crossing the Saddle Homoclinic (DC+) curve will result in a signal terminating with a DC shift and decreasing frequency. ';

};
% Calculate the number of lines for positioning
combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 %Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)

uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
BSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1, 14], 'Value', 4);
Height = Height - maxdip;
descriptionLines = {
    'Dstar'
    '';
    'The dstar parameter is an excitability parameter that controls the ratio between ';
    '';
    'duration of seizure and duration of rest. When dstar is smaller or equal to zero, ';
    '';
    'no seizure activity is possible, and the system will always stay at rest. For small';
    '';
    'positive values of dstar, seizure and rest occur. For sufficiently big values of ';
    '';
    'dstar, only seizure activity is possible. '
};
% Calculate the number of lines for positioning
combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 %Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)

uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
DstarSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [0, 0.5], 'Value', 0.5);
%Height = Height - maxdip;
% % Add labels for the sliders
% descriptionLines = {
%     'N'
%     '';
%     'The N parameter controls solution of resting state. Upper Branch (Case 1): Smoother transitions, reduced hysteresis. Lower Branches (Cases 2 and 3): Potential for hysteresis, with complex, path-dependent responses and multiple equilibria. The systems state may not revert immediately when external conditions are reversed, creating the characteristic hysteresis loop.'
% };
% % Calculate the number of lines for positioning
% combinedText = strjoin(descriptionLines(2:end), ' ');
% 
% % Split the combined text into substrings of maximum length 600
% 
% splitText = {};
% currentIdx = 1;
% 
% while currentIdx <= length(combinedText)
%     if length(combinedText) - currentIdx + 1 <= maxLength
%         % If the remaining text is shorter than maxLength, take the rest of the text
%         splitText{end + 1} = combinedText(currentIdx:end);
%         break;
%     else
%         % Find the position to split the text without breaking words
%         endIdx = currentIdx + maxLength - 1;
%         spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
%         if isempty(spaceIdx)
%             % If no space is found, split at maxLength
%             splitText{end + 1} = combinedText(currentIdx:endIdx);
%             currentIdx = endIdx + 1;
%         else
%             % Split at the last space within the limit
%             endIdx = currentIdx + spaceIdx - 1;
%             splitText{end + 1} = combinedText(currentIdx:endIdx);
%             currentIdx = endIdx + 1;
%         end
%     end
% end
% numChunks = numel(splitText);
% 
% % Calculate the number of lines for positioning
%  %Height = Height - maxdip;
% %numLines = numel(descriptionLines);
% lineHeight = 22; % Height of each line (adjust as needed)
% 
% uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
%     Height = Height -  textdip;
% % Position and display each line as a separate uilabel
% for i = 1:numChunks
% 
%     uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
%     Height = Height -  textdip;
% 
% end
% Height = Height - mindip;
% NSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1, 3], 'Value', 1);
%NSlider = 1;

% Define the segments of the description
descriptionLines = {
    'K';

    'The parameter k determines how many oscillations in the burst by setting the speed.';
    '';
    'The faster the movement the the less the number of oscillations per burst, and vise versa.  ';
};
combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)

uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
%%change k value here
kkSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [0, 0.1], 'Value', 0.003);
% Add sliders for 'input pink noise' and 'tmax'
Height = Height - maxdip;
descriptionLines = {
    'Dynamical Pink Noise'

    'Dynamical pink noise, or parametric noise, is added to the fast variable (x) of';
    '';
    'the governing model equations. This represents noise in the brain';
    '';
    '(i.e. random voltage fluctuations) that creates small perturbations, some of which';
    '';
    'may push the system into or out of the seizure state.'   
};

combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 %Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)

uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
pinkNoiseSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [0, 2000], 'Value', 00);

Height = Height - maxdip;


descriptionLines = {
    'Tmax'

    'How long simulation will run for.';
 
};

combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 %Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)

uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
TmaxSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [100, 10000], 'Value', 500);
% Add labels for the new sliders

Height = Height - maxdip;


descriptionLines = {
    'Signal Acquisition Pink Noise'

    'This filter adds “pink noise”, or 1/f spatial noise with a normal error distribution on top of the simulation to emulate signal acquisition noise. In human data, signal acquisition noise is picked up';
    '';
    'by the electrodes used to record brain activity, and the amplitude depends on the';
    '';
    'signal-to-noise ratio. Most human data has a 20 to 40 percent signal to noise ratio.'   
};

combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
% Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)
uilabel(fig, 'Position', [230, Height, 500,15], 'Text', 'Post-processing parameters', 'FontWeight', 'bold', 'FontSize', fontsize+1);
Height = Height - 20;
uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
acqSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [0, 100], 'Value', 30);

Height = Height - maxdip;

descriptionLines = {
    'Spiking Frequency adjustment factor (Hz)'

    'This filter adjusts the sampling rate of simulated seizures such that spike frequency is consistent with clinical guidelines of rhythmic bursting activity between 1-30 Hz by artificially changing the mean spike rate to a value between 1 and 30 Hz.';
    
 
};

combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 %Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)

uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
FrequencySlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1, 30], 'Value', 20);

Height = Height - maxdip;
descriptionLines = {
    'Electrode Drift Filter'
    'To emulate the residual baseline drift' ;
     '';
    'observed from platinum-iridium electrodes in human EEG, this filter applies a 2nd order digital high-pass filters with' ;
     '';
    'cutoff frequencies between 0.1 and 0.5 Hz. Use 0 - 0.1 Hertz for no correction/minimal correction.';
 
};

% Calculate the number of lines for positioning
combinedText = strjoin(descriptionLines(2:end), ' ');

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 %Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)

uilabel(fig, 'Position', [50,Height, 500, lineHeight], 'Text', descriptionLines{1}, 'FontWeight', 'bold', 'FontSize', fontsize);
    Height = Height -  textdip;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [50,Height, 600, lineHeight], 'Text', splitText{i}, 'FontSize', fontsize);
    Height = Height -  textdip;
    
end
Height = Height - mindip;
DriftSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [0, 0.5], 'Value', 0);

%uilabel(fig, 'Position', [50,240, 100, 22], 'Text', 'Dynamical noise:');
%%change title
titleLabel = uilabel(fig, 'Position', [800, 700, 500, 40], 'Text', 'Class SN/SH - Parameter Control and Visualization', 'FontSize', 16, 'FontWeight', 'bold');
   %change description
descriptionLines = {

    'Class SN/SH seizures have a saddle node (SN) onset bifurcation and a saddle homoclinic (SH) offset bifurcation. In the x time series, this appears as a DC shift that begins when the seizure starts and ends when the seizure stops. Spikes also slow, logarithmically, in frequency at seizure offset. Class SN/SH seizures are found in the Epileptor, which only revealed the most dominant dynamotype: SN/SH';
};



% Calculate the number of lines for positioning
combinedText = descriptionLines{1};

% Split the combined text into substrings of maximum length 600

splitText = {};
currentIdx = 1;

while currentIdx <= length(combinedText)
    if length(combinedText) - currentIdx + 1 <= maxLength
        % If the remaining text is shorter than maxLength, take the rest of the text
        splitText{end + 1} = combinedText(currentIdx:end);
        break;
    else
        % Find the position to split the text without breaking words
        endIdx = currentIdx + maxLength - 1;
        spaceIdx = find(combinedText(currentIdx:endIdx) == ' ', 1, 'last');
        if isempty(spaceIdx)
            % If no space is found, split at maxLength
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        else
            % Split at the last space within the limit
            endIdx = currentIdx + spaceIdx - 1;
            splitText{end + 1} = combinedText(currentIdx:endIdx);
            currentIdx = endIdx + 1;
        end
    end
end
numChunks = numel(splitText);

% Calculate the number of lines for positioning
 %Height = Height - maxdip;
%numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)
Height = 660;
maxLength = 300;
% Position and display each line as a separate uilabel
for i = 1:numChunks
 
    uilabel(fig, 'Position', [620,Height, 850, lineHeight], 'Text', splitText{i}, 'FontSize', 10);
    Height = Height -  20;
    
end
% Add a button to run the simulation
runButton = uibutton(fig, 'Position', [630, 550, 150, 22], 'Text', 'Run Simulation', ...
    'ButtonPushedFcn', @(btn, event) runSimulation(ASlider.Value, BSlider.Value, DstarSlider.Value,...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value, TmaxSlider.Value, FrequencySlider.Value, DriftSlider.Value, fig));
runButton = uibutton(fig, 'Position', [805, 550, 250, 22], 'Text', 'Run Simulation with postprocessing', ...
    'ButtonPushedFcn', @(btn, event) runppSimulation(ASlider.Value, BSlider.Value, DstarSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value, TmaxSlider.Value, FrequencySlider.Value, DriftSlider.Value, fig));
runButton = uibutton(fig, 'Position', [1080, 550, 150, 22], 'Text', 'Run Video Simulation', ...
    'ButtonPushedFcn', @(btn, event) runVideoSimulation(ASlider.Value, BSlider.Value, DstarSlider.Value,...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value, TmaxSlider.Value, FrequencySlider.Value, DriftSlider.Value, fig));
buttomCaption = uilabel(fig, 'Position', [990, 520, 300, 22], 'Text', 'Note: May take up to 20 seconds to run', 'FontSize', fontsize);
runButton = uibutton(fig, 'Position', [610, 50, 100, 22], 'Text', 'Save Simulation', ...
    'ButtonPushedFcn', @(btn, event) get_mat_file(ASlider.Value, BSlider.Value, DstarSlider.Value,...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value, TmaxSlider.Value, FrequencySlider.Value, DriftSlider.Value, fig));
buttomCaption = uilabel(fig, 'Position', [610, 20, 300, 22], 'Text', 'Save generated timeseries as a .mat', 'FontSize', fontsize);
% Add the first plot area
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';

% Add the second plot area
ax2 = uiaxes(fig, 'Position', [950, 100, 350,350]);
ax2.Tag = 'PlotArea2';

% Add a plot area to display the results

   savefig(fig, 'SN_SH_GUI.fig'); 
end

function runSimulation(AIndex, BIndex,  dstar, k, sigma, acq_noise, tmax, frequency, drift, fig)
    addpath('GUI helper files');
    % Get the plot area
    load('bifurcation_crossing.mat')
    load('sphere_mesh.mat')
    load('curves2.mat')
    tstep = 0.01;
    load('curves.mat')
 ax1 = findobj(fig, 'Tag', 'PlotArea1');
ax2 = findobj(fig, 'Tag', 'PlotArea2');
 if ishandle(ax1)
    delete(ax1);
end
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';

if ishandle(ax2)
    delete(ax2);
end
    ax2 = uiaxes(fig, 'Position', [950, 100, 350,350]);
ax2.Tag = 'PlotArea2';
    
N = 1;

    
    % Clear the existing plots
    cla(ax1);
    cla(ax2);
     
    
    % Clear the existing plots
   



A = SHl(:,floor(AIndex));
B = SNr_LCs(:,floor(BIndex));
A = A';
B = B';



     
    
    % Clear the existing plots
   
    
x0=[0;0;0]; % initial conditions (must be a column)

% Settings - Model Sphere
b = 1.0; % focus
R = 0.4; % radius

% Class Information (paths)
% cSN/SH 

tmax = floor(tmax);
% Class specific timespan

tspan = 0:tstep:tmax;

% Points on onset and offset defining path
 % A = AA(CL,:); % canonical path
 % B = BB(CL,:); % canonical path



% Class specific timescale separation


% Class specific threshold

N = floor(N);
% Equilibrium Branch for Resting State


   
% Create circular path based 3 defining points
[E, F] = Parametrization_2PointsArc(A,B,R);

N_t = length(tspan);
X = zeros(3,N_t);
xx = x0;
Rn =  [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 0);pinknoise([1,N_t],-1, 0)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));

    
for n = 1:N_t

    % Euler-Meruyama method
    
    [Fxx, mu2,mu1, nu] = HysteresisLoop_Model(tspan(n),xx,b,k,R,dstar,E,F,N);
    xx = xx + tstep*Fxx + sqrt(tstep)*Rn(:,n);
    X(:,n) = xx;

    mu2_big(n) = mu2;
    mu1_big(n) = mu1;
    nu_big(n) = nu;
    

end

x = X';
z = x(:,3);


% x = smoothdata(x, 'gaussian', onset_offset_range0);
% 
% 
% 
% x_HP = highpass(xtemp(:,1), 0.1, 100);
% Amp2 = movmean(abs(x_HP),5000);
% Amp2=normalize(Amp2);
% Amp2 = smoothdata(Amp2, 'gaussian', sigma*8+onset_offset_range);
% min_val = min(Amp2);
% max_val = max(Amp2);
% Amp2  = 2 * ((Amp2  - min_val) / (max_val - min_val)) - 1;
% Amp2_base = -0.95;
% Add_term = -sigma/50000- Amp2_base;
% Amp2 = Amp2 + Add_term;
% zeroCrossings = find(Amp2(1:end-1) .* Amp2(2:end) < 0);
% 
% 
% 
% 
% if length(zeroCrossings) < 3
%     onset = 1;
% 
% % Compute points along the arc
% offset = 1000;
% else
% if(Amp2(1) < 0.0)
% onset_time = zeroCrossings(1:2:end)*0.01;  % Even indices
% offset_time = zeroCrossings(2:2:end)*0.01; % Odd indices
% 
% else
% onset_time = zeroCrossings(2:2:end)*0.01;  % Even indices
% offset_time = zeroCrossings(3:2:end)*0.01; % Odd indices
% end
% onset = onset_time(1)*100 - 10000;
% if onset < 1
%     onset = 1;
% end
% offset = offset_time(1)*100 + 10000;
% end
onset = 1;
offset = tmax*100;
seizure = x(:,1);
z = x(:,3);

% noise_percentage= acq_noise/100;
%      L = length(seizure);
%     noise = pinknoise([1,L],-1,5.69e3*noise_percentage); 
%    seizure = seizure + noise';

hold(ax1,'on')
% plot(ax1,xtemp(:,1));
% plot(ax1,Amp2);
plot(ax1,tspan,seizure,'color','#696969','LineWidth',1)
plot(ax1,tspan, x(:,3), 'r');
xlabel(ax1, 'Time')
ylabel(ax1, "Simulated Seizure")
title(ax1, 'Timeseries');
 


linewidth = 2;

  % % openfig('Map.fig');
 
 
hold(ax2, 'on')
% load("map_regions.mat");


  

%[Left,right  up/down]
%%%
lineVector = [0,1,0];
load('region_mesh.mat')
vertices = BCSmesh.vertices;
faces = BCSmesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.973, 0.965, 0.722], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off' );

vertices = Active_restmesh.vertices;
faces = Active_restmesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.9216  ,  0.9216  ,  0.9216], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Active rest region');
vertices = Seizure_mesh.vertices;
faces = Seizure_mesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [0.894, 0.706, 0.831], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Seizure region');
faces = Bistable_Lcb_mesh.faces;
vertices = Bistable_Lcb_mesh.vertices;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.973, 0.965, 0.722], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceLighting', 'none', 'BackFaceLighting', 'reverselit', 'DisplayName', "Bistability region");
   
scale_array = [0.4];
scale_array = scale_array/0.4;
num_scales = length(scale_array);


% Scale factors to adjust radius from 0.4 to 0.39 and 0.38
for i = 1:length(scale_array)
% Scale the coordinates of the points for radius 0.39
Fold_of_cycles_scaled = scale_array(i) * Fold_of_cycles;
Homoclinic_to_saddle3_scaled = scale_array(i) * Homoclinic_to_saddle3;
Homoclinic_to_saddle2_scaled = scale_array(i) * Homoclinic_to_saddle2;
Homoclinic_to_saddle1_scaled = scale_array(i) * Homoclinic_to_saddle1;
Homoclinic_to_saddle_scaled = scale_array(i) * Homoclinic_to_saddle;
Fold_scaled = scale_array(i) * Fold;
Hopf_scaled = scale_array(i) * Hopf;
SNIC_scaled = scale_array(i) * SNIC;
plot3(ax2,Fold_of_cycles_scaled(1, :), Fold_of_cycles_scaled(2, :), Fold_of_cycles_scaled(3, :), 'Color', [0.9725,0.2667,0.5843], 'LineWidth', linewidth, 'HandleVisibility', 'off' );
plot3(ax2,Homoclinic_to_saddle3_scaled(1, :), Homoclinic_to_saddle3_scaled(2, :), Homoclinic_to_saddle3_scaled(3, :), 'Color', [0.404, 0.702, 0.851], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle2_scaled(1, :), Homoclinic_to_saddle2_scaled(2, :), Homoclinic_to_saddle2_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle1_scaled(1, :), Homoclinic_to_saddle1_scaled(2, :), Homoclinic_to_saddle1_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth,  'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'DisplayName', 'Offset curve');
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'DisplayName', 'Onset curve');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--', 'DisplayName', 'Onset and Offset curve', 'HandleVisibility', 'off');
end
%_half_sphere_lcs.mat');
%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.8, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

scale = 0.41/0.4;

scale = 0.42/0.4;
plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scatter3(ax2, A(1)*scale,A(2)*scale,A(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker', 'diamond') % point on the offset curve
scatter3(ax2, B(1)*scale,B(2)*scale,B(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker','square')
%scatter3(ax2, p3(1)*scale,p3(2)*scale,p3(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker','pentagram')
lgd = legend(ax2, 'show');
lgd.FontSize = 6;
%lgd.Position = [0.35, 0.75, 0.1, 0.1];

   xlabel(ax2, 'mu2');
    ylabel(ax2, '-mu1');
    zlabel(ax2, 'nu');
   
  az = atan2d(lineVector(2), lineVector(1)); % Azimuth angle
    el = atan2d(lineVector(3), norm(lineVector(1:2))); % Elevation angle
 view(ax2, az, el);
  title(ax2, 'Parameter space');


end

function runVideoSimulation(AIndex, BIndex, dstar,  k, sigma, acq_noise, tmax, frequency, drift, fig)
addpath('GUI helper files');
    % Get the plot area
    load('bifurcation_crossing.mat')
    load('sphere_mesh.mat')
    load('curves2.mat')
    
    load('curves.mat')
 ax1 = findobj(fig, 'Tag', 'PlotArea1');
ax2 = findobj(fig, 'Tag', 'PlotArea2');

tstep = 0.01;
N = 1;
 if ishandle(ax1)
    delete(ax1);
end


if ishandle(ax2)
    delete(ax2);
end




img = imread('hysterisis_legend.png'); % Replace 'your_image.jpg' with the path to your JPG file
[img_height, img_width, ~] = size(img);

  

    % Create an axes component that covers the entire figure window
    ax2 = uiaxes(fig, 'Position', [950, 120, img_width/1.5, img_height/1.5]);

    % Display the image in the axes component
    imshow(img, 'Parent', ax2);
img = imread('new_hysterisis_paths_2d.png'); % Replace 'your_image.jpg' with the path to your JPG file
[img_height, img_width, ~] = size(img);
ax1 = uiaxes(fig, 'Position', [600, 80, img_width/2, img_height/2]);


    % Display the image in the axes component
    imshow(img, 'Parent', ax1);
    %axis(ax2, 'off'); % Turn off the axes
    % Display the image in the axes component
 
    pause(2);

if ishandle(ax2)
    delete(ax2);
end


    ax2 = uiaxes(fig, 'Position', [950, 100, 350,350]);
ax2.Tag = 'PlotArea2';


A = SHl(:,floor(AIndex));
B = SNr_LCs(:,floor(BIndex));
A = A';
B = B';



     
    
    % Clear the existing plots
   
    
x0=[0;0;0]; % initial conditions (must be a column)

% Settings - Model Sphere
b = 1.0; % focus
R = 0.4; % radius

% Class Information (paths)
% cSN/SH 

tmax = floor(tmax);
% Class specific timespan

tspan = 0:tstep:tmax;

% Points on onset and offset defining path
 % A = AA(CL,:); % canonical path
 % B = BB(CL,:); % canonical path



% Class specific timescale separation


% Class specific threshold

N = floor(N);
% Equilibrium Branch for Resting State


   
% Create circular path based 3 defining points
[E, F] = Parametrization_2PointsArc(A,B,R);

N_t = length(tspan);
X = zeros(3,N_t);
xx = x0;
Rn =  [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 0);pinknoise([1,N_t],-1, 0)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));

    
for n = 1:N_t

    % Euler-Meruyama method
    
    [Fxx, mu2,mu1, nu] = HysteresisLoop_Model(tspan(n),xx,b,k,R,dstar,E,F,N);
    xx = xx + tstep*Fxx + sqrt(tstep)*Rn(:,n);
    X(:,n) = xx;

    mu2_big(n) = mu2;
    mu1_big(n) = mu1;
    nu_big(n) = nu;
    

end
x = X';

xtemp = x;


linewidth = 2;
 mu2_big = mu2_big';
    mu1_big = mu1_big';
    nu_big = nu_big';

  % % openfig('Map.fig');
 
hold(ax2, 'on')
% load("map_regions.mat");


  

%[Left,right  up/down]
%%%

initialLineVector = [0,1,0];
targetLineVector =  [0,1,0.6];

load('region_mesh.mat')
vertices = BCSmesh.vertices;
faces = BCSmesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.973, 0.965, 0.722], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off' );

vertices = Active_restmesh.vertices;
faces = Active_restmesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.9216  ,  0.9216  ,  0.9216], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Active rest region');
vertices = Seizure_mesh.vertices;
faces = Seizure_mesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [0.894, 0.706, 0.831], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Seizure region');
faces = Bistable_Lcb_mesh.faces;
vertices = Bistable_Lcb_mesh.vertices;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.973, 0.965, 0.722], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceLighting', 'none', 'BackFaceLighting', 'reverselit', 'DisplayName', "Bistability region");
   
scale_array = [0.4];
scale_array = scale_array/0.4;
num_scales = length(scale_array);


% Scale factors to adjust radius from 0.4 to 0.39 and 0.38
for i = 1:length(scale_array)
% Scale the coordinates of the points for radius 0.39
Fold_of_cycles_scaled = scale_array(i) * Fold_of_cycles;
Homoclinic_to_saddle3_scaled = scale_array(i) * Homoclinic_to_saddle3;
Homoclinic_to_saddle2_scaled = scale_array(i) * Homoclinic_to_saddle2;
Homoclinic_to_saddle1_scaled = scale_array(i) * Homoclinic_to_saddle1;
Homoclinic_to_saddle_scaled = scale_array(i) * Homoclinic_to_saddle;
Fold_scaled = scale_array(i) * Fold;
Hopf_scaled = scale_array(i) * Hopf;
SNIC_scaled = scale_array(i) * SNIC;
plot3(ax2,Fold_of_cycles_scaled(1, :), Fold_of_cycles_scaled(2, :), Fold_of_cycles_scaled(3, :), 'Color', [0.9725,0.2667,0.5843], 'LineWidth', linewidth, 'HandleVisibility', 'off' );
plot3(ax2,Homoclinic_to_saddle3_scaled(1, :), Homoclinic_to_saddle3_scaled(2, :), Homoclinic_to_saddle3_scaled(3, :), 'Color', [0.404, 0.702, 0.851], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle2_scaled(1, :), Homoclinic_to_saddle2_scaled(2, :), Homoclinic_to_saddle2_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle1_scaled(1, :), Homoclinic_to_saddle1_scaled(2, :), Homoclinic_to_saddle1_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth,  'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'DisplayName', 'Offset curve');
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'DisplayName', 'Onset curve');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--', 'DisplayName', 'Onset and Offset curve', 'HandleVisibility', 'off');
end
%_half_sphere_lcs.mat');
%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

scale = 0.41/0.4;

scale = 0.4/0.4;
%plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scatter3(ax2, A(1)*scale,A(2)*scale,A(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker', 'diamond') % point on the offset curve
scatter3(ax2, B(1)*scale,B(2)*scale,B(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker','square')
%scatter3(ax2, p3(1)*scale,p3(2)*scale,p3(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker','pentagram')
lgd = legend(ax2, 'show');
lgd.FontSize = 6;
%lgd.Position = [0.35, 0.75, 0.1, 0.1];



   xlabel(ax2, 'mu2');
    ylabel(ax2, '-mu1');
    zlabel(ax2, 'nu');
  
 %  az = atan2d(lineVector(2), lineVector(1)); % Azimuth angle
 %    el = atan2d(lineVector(3), norm(lineVector(1:2))); % Elevation angle
 % view(ax2, az, el);
  title(ax2, 'Parameter space');

numSteps = 20;
pauseDuration = 0.05; % Pause duration for effect
currentLineVector = initialLineVector  ;
    
    % Compute azimuth and elevation angles
    az = atan2d(currentLineVector(2), currentLineVector(1)); % Azimuth angle
     el = atan2d(currentLineVector(3), norm(currentLineVector(1:2))); % Elevation angle
    
    % Update the view
    view(ax2, az, el);

pause(2);
% Number of steps for rotation

% Linearly interpolate the vector from initial to target
for step = 1:numSteps
    % Interpolated vector
    currentLineVector = initialLineVector + step / numSteps * (targetLineVector - initialLineVector);
    
    % Compute azimuth and elevation angles
    az = atan2d(currentLineVector(2), currentLineVector(1)); % Azimuth angle
     el = atan2d(currentLineVector(3), norm(currentLineVector(1:2))); % Elevation angle
    
    % Update the view
    view(ax2, az, el);
    
 
    
    % Pause for effect
    pause(pauseDuration);
end

onset = 1;
offset = tmax*100;
seizure = xtemp(onset:offset,1);

 if ishandle(ax1)
    delete(ax1);
end
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';


n = length(seizure);
%hold(ax2,'off');
interval = 500;
hold(ax1,'on')
% Calculate Onset Times
[pks,times]=findpeaks(x(:,3),'MinPeakProminence',0.005);
onset_time = times*tstep;

% Calculate Offset Times
[pks2,times2]=findpeaks(-x(:,3),'MinPeakProminence',0.005);
offset_time = times2*tstep;







onset_plot = false;
offset_plot = false;
seizure_plot = false;
plot_size = 2;
increase_plot_size = floor(((2*pi/k)/tstep)/interval)*interval +1;

sigma = 20; % Adjust sigma based on the desired level of smoothing
z = x(:,3);
% Smooth the time series data
z_smooth = smoothdata(z, 'gaussian', 5000);

dt = diff(tspan'); % Difference between consecutive time points
dz = diff(z_smooth); % Difference between consecutive data points

% Calculate the slope at each point (excluding the last point since diff reduces the array size by 1)
slope = dz ./ dt; 

% If you want to keep the slope array the same size as the original time series:
slope = [slope; slope(end)];

heatmap = zeros(size(slope));

% Assign values based on the slope conditions
heatmap(slope > 0) = 1;          % Assign 1 if slope > 0.1
heatmap(slope <= 0.000001 & slope >= -0.000001) = 0;  % Assign 0 if -0.1 <= slope <= 0.1
heatmap(slope <= 0) = -1;        % Assign -1 if slope < -0.1


xlim(ax2, [-0.5 0.5]);
ylim(ax2, [-0.5 0.5]);
zlim(ax2, [-0.5 0.5]);
xticks(ax2, []);
yticks(ax2, []);
zticks(ax2, []);


plot_size = 20;
scale = 1;
% Calculate the slope of the smoothed data
dzdt = gradient(z_smooth) ./ gradient(tspan');
% Loop through the variable with the specified interval
for i = 1:interval:n
    % Handle the case for the last interval
    if i + interval - 1 <= n
        % if i == interval_onset_start
        %  onset_plot = true;
        % end
        % if i == interval_onset_stop
        %    onset_plot = false;
        %    interval_onset = floor((interval_onset + (2*pi/k)/tstep)/interval)*interval +1;
        %    interval_onset_start = interval_onset - onset_offset_range;
        %    interval_onset_stop = interval_onset + onset_offset_range;
        %    seizure_plot = true;
        % end
        % if i == interval_offset_start
        % 
        %     seizure_plot = false;
        %  offset_plot = true;
        % end
        % if i == interval_offset_stop
        %    offset_plot = false;
        %    interval_offset = floor((interval_offset + (2*pi/k)/tstep)/interval)*interval + 1;
        %     interval_offset_start = interval_offset - onset_offset_range;
        %     interval_offset_stop = interval_offset + onset_offset_range;
        % 
        % end
        % if i == increase_plot_size
        % increase_plot_size = floor((increase_plot_size+(2*pi/k)/tstep)/interval)*interval +1;
        % plot_size = plot_size+2;
        % end
        % if onset_plot == true
        % plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'Color',[0.957, 0.612, 0.204],  'LineStyle', '--', 'HandleVisibility', 'off');
        % plot3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,'k','LineWidth',plot_size,  'HandleVisibility', 'off');
        % elseif offset_plot == true
        % plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'Color',[0.957, 0.612, 0.204],  'LineStyle', '--', 'HandleVisibility', 'off');
        % plot3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,'k','LineWidth',plot_size,  'HandleVisibility', 'off');
        % elseif seizure_plot == true
        % plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'Color',[0.894, 0.706, 0.831],  'HandleVisibility', 'off');
        % plot3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,'k','LineWidth',plot_size,  'HandleVisibility', 'off');
        % else
            
        plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'k');
        %plot(ax1,tspan(i:i+interval-1), x(i:(i+interval-1),3), 'r');
        scatter(ax1,tspan(i:i+interval-1),z(i:i+interval-1), 50, heatmap(i:i+interval-1), 'filled');
        scatter3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,10, heatmap(i:i+interval-1) , 'HandleVisibility', 'off');
       scale = scale+0.001;

    else
        plot(ax1,tspan(i:n), seizure(i:n), 'k');
        scatter3(ax2, mu2_big(i:n,1)*scale,-mu1_big(i:n,1)*scale,nu_big(i:n,1)*scale,plot_size, heatmap(i:n) , 'HandleVisibility', 'off');
    end
    pause(0.0000001)
end
%%%%change here
%plot(ax1,tspan,dzdt);
hold(ax2,'off');

if((length(onset_time) >= 1 && length(offset_time) >= 1 && offset_time(1)>onset_time(1)) || length(onset_time) >= 1 && length(offset_time) >= 2) 

% Single seizure
if offset_time(1)>onset_time(1) % if system starts at rest
    start_index = times(1);
    stop_index = times2(1);

else % if system starts in a seizure
    start_index = times(1);
    stop_index = times2(2);

end

[pks,locs] = findpeaks(seizure(floor(start_index):floor(stop_index)), 'MinPeakProminence', 0.15);




sampling_rate = 10000; % Hz
time_in_seconds = locs / sampling_rate; % Convert peak indices to seconds

% Calculate spike rates
spike_rates = diff(time_in_seconds);
% disp('Interspike interval');
% disp(spike_rates);

% Adjust spike rates to achieve a mean average spiking rate of 5 Hz
target_avg_spike_rate = (1/frequency); % Hz

% Calculate the current average spiking rate
current_avg_spike_rate = mean(spike_rates);

% Calculate the adjustment factor
adjustment_factor = target_avg_spike_rate / current_avg_spike_rate;

% Adjust spike rates
adjusted_spike_rates = spike_rates * adjustment_factor;

% Calculate the mean average spiking rate after adjustment
mean_avg_spike_rate = mean(adjusted_spike_rates);

% Calculate the new sampling frequency
new_sampling_frequency = sampling_rate / adjustment_factor;


num_samples = length(seizure);

time = (0:num_samples-1) / new_sampling_frequency; % Time array in seconds
hold(ax1,'off');

if ishandle(ax1)
    delete(ax1);
end
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';
% Add text to the axis
x_position = 0.05; % X-coordinate of the text
y_position = 0.5; % Y-coordinate of the text
str = 'Adjusting sampling frequency to match human data'; % The text to display

% Use the text function to place the text on the plot
text(ax1, x_position, y_position, str, 'FontSize', 12, 'Color', 'k');
%title(ax1, 'Adjusting sampling frequency to match human data');
pause(2);
cla(ax1);
plot(ax1,time, seizure,'color','#696969','LineWidth',1)
xlabel(ax1, 'Adjusted time');
ylabel(ax1, "Simulated Seizure")
title(ax1, 'Timeseries');
pause(3);


hpf = drift;
if drift > 0
HPF1 = designfilt('highpassiir', ...       % Response type
    'FilterOrder',2, ...     % Frequency constraints
       'HalfPowerFrequency',hpf, ...     % Frequency constraints
       'DesignMethod','butter', ...     % Design method
       'SampleRate',new_sampling_frequency);               % Sample rate



seizure= filter(HPF1,seizure);
end
if ishandle(ax1)
    delete(ax1);
end
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';
x_position = 0.25; % X-coordinate of the text
y_position = 0.5; % Y-coordinate of the text
str = 'Applying electrode drift filter'; % The text to display

% Use the text function to place the text on the plot
text(ax1, x_position, y_position, str, 'FontSize', 12, 'Color', 'k');
pause(2);
plot(ax1,time, seizure,'color','#696969','LineWidth',1)
xlabel(ax1, 'Adjusted time');
ylabel(ax1, "Simulated Seizure")
title(ax1, 'Timeseries');
pause(3);

noise_percentage= acq_noise/100;
     L = length(seizure);
    noise = pinknoise([1,L],-1,5.69e3*noise_percentage); 
   seizure = seizure + noise';


if ishandle(ax1)
    delete(ax1);
end
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';
x_position = 0.23; % X-coordinate of the text
y_position = 0.5; % Y-coordinate of the text
str = 'Adding signal acquisition noise'; % The text to display

% Use the text function to place the text on the plot
text(ax1, x_position, y_position, str, 'FontSize', 12, 'Color', 'k');
pause(2);
plot(ax1,time, seizure,'color','#696969','LineWidth',1)
xlabel(ax1, 'Adjusted time');
ylabel(ax1, "Simulated Seizure")
title(ax1, 'Timeseries');
else
title(ax1, 'No post-processing done, simulation does not have at least one seizure');
end





end


function runppSimulation(AIndex, BIndex, dstar,  k, sigma, acq_noise, tmax, frequency, drift, fig)
addpath('GUI helper files');
    % Get the plot area
    load('bifurcation_crossing.mat')
    load('sphere_mesh.mat')
    load('curves2.mat')
    
    load('curves.mat')
 ax1 = findobj(fig, 'Tag', 'PlotArea1');
ax2 = findobj(fig, 'Tag', 'PlotArea2');

tstep = 0.01;

 if ishandle(ax1)
    delete(ax1);
end
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';

if ishandle(ax2)
    delete(ax2);
end
    ax2 = uiaxes(fig, 'Position', [950, 100, 350,350]);
ax2.Tag = 'PlotArea2';


N = 1;


A = SHl(:,floor(AIndex));
B = SNr_LCs(:,floor(BIndex));
A = A';
B = B';



     
    
    % Clear the existing plots
   
    
x0=[0;0;0]; % initial conditions (must be a column)

% Settings - Model Sphere
b = 1.0; % focus
R = 0.4; % radius

% Class Information (paths)
% cSN/SH 

tmax = floor(tmax);
% Class specific timespan

tspan = 0:tstep:tmax;

% Points on onset and offset defining path
 % A = AA(CL,:); % canonical path
 % B = BB(CL,:); % canonical path



% Class specific timescale separation


% Class specific threshold


% Equilibrium Branch for Resting State


   
% Create circular path based 3 defining points
[E, F] = Parametrization_2PointsArc(A,B,R);

N_t = length(tspan);
X = zeros(3,N_t);
xx = x0;
Rn =  [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 0);pinknoise([1,N_t],-1, 0)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));

    
for n = 1:N_t

    % Euler-Meruyama method
    
    [Fxx, mu2,mu1, nu] = HysteresisLoop_Model(tspan(n),xx,b,k,R,dstar,E,F,N);
    xx = xx + tstep*Fxx + sqrt(tstep)*Rn(:,n);
    X(:,n) = xx;

    mu2_big(n) = mu2;
    mu1_big(n) = mu1;
    nu_big(n) = nu;
    

end
x = X';

xtemp = x;


linewidth = 2;
 mu2_big = mu2_big';
    mu1_big = mu1_big';
    nu_big = nu_big';

  % % openfig('Map.fig');
 
hold(ax2, 'on')
% load("map_regions.mat");


  

%[Left,right  up/down]
%%%

initialLineVector = [0,1,0];
targetLineVector = [0,1,0];

load('region_mesh.mat')
vertices = BCSmesh.vertices;
faces = BCSmesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.973, 0.965, 0.722], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off' );

vertices = Active_restmesh.vertices;
faces = Active_restmesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.9216  ,  0.9216  ,  0.9216], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Active rest region');
vertices = Seizure_mesh.vertices;
faces = Seizure_mesh.faces;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [0.894, 0.706, 0.831], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Seizure region');
faces = Bistable_Lcb_mesh.faces;
vertices = Bistable_Lcb_mesh.vertices;
patch(ax2,'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [ 0.973, 0.965, 0.722], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceLighting', 'none', 'BackFaceLighting', 'reverselit', 'DisplayName', "Bistability region");
   
scale_array = [0.4];
scale_array = scale_array/0.4;
num_scales = length(scale_array);


% Scale factors to adjust radius from 0.4 to 0.39 and 0.38
for i = 1:length(scale_array)
% Scale the coordinates of the points for radius 0.39
Fold_of_cycles_scaled = scale_array(i) * Fold_of_cycles;
Homoclinic_to_saddle3_scaled = scale_array(i) * Homoclinic_to_saddle3;
Homoclinic_to_saddle2_scaled = scale_array(i) * Homoclinic_to_saddle2;
Homoclinic_to_saddle1_scaled = scale_array(i) * Homoclinic_to_saddle1;
Homoclinic_to_saddle_scaled = scale_array(i) * Homoclinic_to_saddle;
Fold_scaled = scale_array(i) * Fold;
Hopf_scaled = scale_array(i) * Hopf;
SNIC_scaled = scale_array(i) * SNIC;
plot3(ax2,Fold_of_cycles_scaled(1, :), Fold_of_cycles_scaled(2, :), Fold_of_cycles_scaled(3, :), 'Color', [0.9725,0.2667,0.5843], 'LineWidth', linewidth, 'HandleVisibility', 'off' );
plot3(ax2,Homoclinic_to_saddle3_scaled(1, :), Homoclinic_to_saddle3_scaled(2, :), Homoclinic_to_saddle3_scaled(3, :), 'Color', [0.404, 0.702, 0.851], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle2_scaled(1, :), Homoclinic_to_saddle2_scaled(2, :), Homoclinic_to_saddle2_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle1_scaled(1, :), Homoclinic_to_saddle1_scaled(2, :), Homoclinic_to_saddle1_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth,  'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'DisplayName', 'Offset curve');
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'DisplayName', 'Onset curve');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--', 'DisplayName', 'Onset and Offset curve', 'HandleVisibility', 'off');
end
%_half_sphere_lcs.mat');
%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

scale = 0.41/0.4;

scale = 0.4/0.4;
plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scatter3(ax2, A(1)*scale,A(2)*scale,A(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker', 'diamond') % point on the offset curve
scatter3(ax2, B(1)*scale,B(2)*scale,B(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker','square')
%scatter3(ax2, p3(1)*scale,p3(2)*scale,p3(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker','pentagram')
lgd = legend(ax2, 'show');
lgd.FontSize = 6;
%lgd.Position = [0.35, 0.75, 0.1, 0.1];



   xlabel(ax2, 'mu2');
    ylabel(ax2, '-mu1');
    zlabel(ax2, 'nu');
  
 %  az = atan2d(lineVector(2), lineVector(1)); % Azimuth angle
 %    el = atan2d(lineVector(3), norm(lineVector(1:2))); % Elevation angle
 % view(ax2, az, el);
  title(ax2, 'Parameter space');

numSteps = 20;
pauseDuration = 0.05; % Pause duration for effect
currentLineVector = [0,1,0] ;
    
    % Compute azimuth and elevation angles
    az = atan2d(currentLineVector(2), currentLineVector(1)); % Azimuth angle
     el = atan2d(currentLineVector(3), norm(currentLineVector(1:2))); % Elevation angle
    
    % Update the view
    view(ax2, az, el);



    
 
    

onset = 1;
offset = tmax*100;
seizure = xtemp(onset:offset,1);




n = length(seizure);
%hold(ax2,'off');
interval = 500;
hold(ax1,'on')
% Calculate Onset Times
[pks,times]=findpeaks(x(:,3),'MinPeakProminence',0.005);
onset_time = times*tstep;

% Calculate Offset Times
[pks2,times2]=findpeaks(-x(:,3),'MinPeakProminence',0.005);
offset_time = times2*tstep;







onset_plot = false;
offset_plot = false;
seizure_plot = false;
plot_size = 2;
increase_plot_size = floor(((2*pi/k)/tstep)/interval)*interval +1;

sigma = 20; % Adjust sigma based on the desired level of smoothing
z = x(:,3);
% Smooth the time series data
z_smooth = smoothdata(z, 'gaussian', 5000);

dt = diff(tspan'); % Difference between consecutive time points
dz = diff(z_smooth); % Difference between consecutive data points

% Calculate the slope at each point (excluding the last point since diff reduces the array size by 1)
slope = dz ./ dt; 

% If you want to keep the slope array the same size as the original time series:
slope = [slope; slope(end)];

heatmap = zeros(size(slope));

% Assign values based on the slope conditions
heatmap(slope > 0.000001) = 1;          % Assign 1 if slope > 0.1
heatmap(slope <= 0.000001 & slope >= -0.000001) = 0;  % Assign 0 if -0.1 <= slope <= 0.1
heatmap(slope < -0.000001) = -1;        % Assign -1 if slope < -0.1


xlim(ax2, [-0.5 0.5]);
ylim(ax2, [-0.5 0.5]);
zlim(ax2, [-0.5 0.5]);
xticks(ax2, []);
yticks(ax2, []);
zticks(ax2, []);


plot_size = 20;
scale = 1;
% Calculate the slope of the smoothed data
dzdt = gradient(z_smooth) ./ gradient(tspan');
% Loop through the variable with the specified interval


if((length(onset_time) >= 1 && length(offset_time) >= 1 && offset_time(1)>onset_time(1)) || length(onset_time) >= 1 && length(offset_time) >= 2) 

% Single seizure
if offset_time(1)>onset_time(1) % if system starts at rest
    start_index = times(1);
    stop_index = times2(1);

else % if system starts in a seizure
    start_index = times(1);
    stop_index = times2(2);

end

[pks,locs] = findpeaks(seizure(floor(start_index):floor(stop_index)), 'MinPeakProminence', 0.15);




sampling_rate = 10000; % Hz
time_in_seconds = locs / sampling_rate; % Convert peak indices to seconds

% Calculate spike rates
spike_rates = diff(time_in_seconds);
% disp('Interspike interval');
% disp(spike_rates);

% Adjust spike rates to achieve a mean average spiking rate of 5 Hz
target_avg_spike_rate = (1/frequency); % Hz

% Calculate the current average spiking rate
current_avg_spike_rate = mean(spike_rates);

% Calculate the adjustment factor
adjustment_factor = target_avg_spike_rate / current_avg_spike_rate;

% Adjust spike rates
adjusted_spike_rates = spike_rates * adjustment_factor;

% Calculate the mean average spiking rate after adjustment
mean_avg_spike_rate = mean(adjusted_spike_rates);

% Calculate the new sampling frequency
new_sampling_frequency = sampling_rate / adjustment_factor;


num_samples = length(seizure);

time = (0:num_samples-1) / new_sampling_frequency; % Time array in seconds



hpf = drift;
if drift > 0
HPF1 = designfilt('highpassiir', ...       % Response type
    'FilterOrder',2, ...     % Frequency constraints
       'HalfPowerFrequency',hpf, ...     % Frequency constraints
       'DesignMethod','butter', ...     % Design method
       'SampleRate',new_sampling_frequency);               % Sample rate



seizure= filter(HPF1,seizure);
end

noise_percentage= acq_noise/100;
     L = length(seizure);
    noise = pinknoise([1,L],-1,5.69e3*noise_percentage); 
   seizure = seizure + noise';




plot(ax1,time, seizure,'color','#696969','LineWidth',1)
xlabel(ax1, 'Adjusted time');
ylabel(ax1, "Simulated Seizure")
title(ax1, 'Timeseries');
else
title(ax1, 'No post-processing done, simulation does not have at least one seizure');
end





end
function get_mat_file(AIndex, BIndex, dstar, k, sigma, acq_noise, tmax, frequency, drift, fig)
    ax1 = findobj(fig, 'Tag', 'PlotArea1');
    lines = findobj(ax1, 'Type', 'Line');

% Initialize arrays to hold the data
Time = {};
Signal = {};
if length(lines) == 0
    xlabel(ax1, 'Simulation has not been run. Run simulation to save file');
else
% Loop through each Line object and store the data
for i = 1:length(lines)
    Time{i} = lines(i).XData;
    Signal{i} = lines(i).YData;
end

xlabel_text = ax1.XLabel.String;
if strcmp(xlabel_text, 'Adjusted time')
    filename = sprintf('SN_SH_timeseries_Offset_index_%d_Onset_index_%d_dstar_%.2f_k_%.2f_dynamical_noise_%.2f_tMax_%.2f_adjusted_frequency_%.2f_electrode_drift_filter_%.2f.mat', ...
        AIndex, BIndex, dstar, k, sigma, tmax, acq_noise,  frequency, drift);


else

%change here to title for class
filename = sprintf('SN_SH_timeseries_Offset_index_%d_Onset_index_%d_dstar_%.2f_k_%.2f_dynamical_noise_%.2f_tMax_%.2f.mat', ...
        AIndex, BIndex, dstar, k, sigma, tmax);
end


% Save the data to a .mat file
save(filename, 'Time', 'Signal');
title(ax1, 'Timeseries has saved succesfully!');

end

end



function [Xdot, mu2, mu1,nu] = HysteresisLoop_Model(~,x,~,k,R,dstar,E,F,N)
    
    % Parametrization of the path in the spherical parameter space in terms of great
    % circles      
    mu2=R*(E(1)*cos(x(3))+F(1)*sin(x(3)));
    mu1=-R*(E(2)*cos(x(3))+F(2)*sin(x(3)));
    nu=R*(E(3)*cos(x(3))+F(3)*sin(x(3)));

    % x coordinate of resting state (i.e. upper branch of eq)     
    x_rs=real(Resting_State(mu2,mu1,nu,N));
    
    % equations    
    xdot = - x(2);
    ydot = x(1)^3 - mu2*x(1) - mu1 - x(2)*( nu + x(1) + x(1)^2);
    
    zdot =  -k*(sqrt((x(1)-x_rs)^2+x(2)^2)-dstar);
    
    Xdot = [xdot;ydot;zdot];

end  

function [E,F] = Parametrization_2PointsArc(A,B,R)

    E = A./R;

    F=cross(cross(A,B),A);
    F=F./norm(F);

end

function x = pinknoise(DIM,BETA, MAG),
% This function generates 1/f spatial noise, with a normal error
% distribution 
%
% DIM is a two component vector that sets the size of the spatial pattern
%       (DIM=[10,5] is a 10x5 spatial grid)
%
% BETA defines the spectral distribution.
%      Spectral density S(f) = N f^BETA
%      (f is the frequency, N is normalisation coeff).
%           BETA = 0 is random white noise.  
%           BETA  -1 is pink noise
%           BETA = -2 is Brownian noise
%      The fractal dimension is related to BETA by, D = (6+BETA)/2
%
% MAG is the scaling variable for the noise amplitude
%
% The method is briefly descirbed in Lennon, J.L. "Red-shifts and red
% herrings in geographical ecology", Ecography, Vol. 23, p101-113 (2000)

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
u = repmat(u,1,DIM(2));
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]/DIM(2);
v = repmat(v,DIM(1),1);
S_f = (u.^2 + v.^2).^(BETA/2);
S_f(S_f==inf) = 0;
phi = rand(DIM);
y= S_f.^0.5 .* (cos(2*pi*phi)+i*sin(2*pi*phi));
y=y.*MAG/max(abs(y));  %set the mag to the level you want
x= ifft2(y);
x = real(x);

end

function x_rs = Resting_State(mu2,mu1,nu,N)

    switch N
        case 1 % resting state on upper branch
            x_rs=mu2/(3*(mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3)) + (mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3);
    
        case 2 % resting state on lower branch
            x_rs=- mu2/(6*(mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3)) - (mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3)/2 - (3^(1/2)*i*(mu2/(3*(mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3)) - (mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3)))/2;
     
        case 3
            x_rs= (3^(1/2)*i*(mu2/(3*(mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3)) - (mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3)))/2 - (mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3)/2 - mu2/(6*(mu1/2 + (mu1^2/4 - mu2^3/27)^(1/2))^(1/3));
    end

end
function [mu2,mu1,nu] = sphere(point1, point2, numPoints)
    % sphereArcPath - Generates an arc path between two points on a sphere
    %
    % Syntax: arcPath = sphereArcPath(point1, point2, numPoints)
    %
    % Inputs:
    %    point1 - [x1, y1, z1] Coordinates of the first point on the sphere
    %    point2 - [x2, y2, z2] Coordinates of the second point on the sphere
    %    numPoints - Number of points along the arc
    %
    % Outputs:
    %    arcPath - An Nx3 matrix containing the coordinates of points along the arc
    
    % Check the input points
    radius = 0.4;
    % if norm(point1) ~= radius || norm(point2) ~= radius
    %     error('The points must lie on the sphere of radius 0.4.');
    % end
    
    % Normalize the input points to make sure they are on the sphere
    point1 = point1 / norm(point1) * radius;
    point2 = point2 / norm(point2) * radius;
    
    % Compute the quaternion for rotation
    theta = acos(dot(point1, point2) / (radius^2));
    axis = cross(point1, point2);
    if norm(axis) == 0
        error('The points are the same or antipodal.');
    end
    axis = axis / norm(axis);

    % Compute points along the arc
    arcPath = zeros(numPoints, 3);
    for i = 0:numPoints-1
        t = i / (numPoints - 1);
        angle = t * theta;
        R = rotationMatrix(axis, angle);
        arcPath(i+1, :) = (R * point1')';
    end
    mu2 = arcPath(:,1)';
    mu1 = arcPath(:,2)';
    nu = arcPath(:,3)';
end

function R = rotationMatrix(axis, angle)
    % rotationMatrix - Generates a rotation matrix given an axis and an angle
    %
    % Syntax: R = rotationMatrix(axis, angle)
    %
    % Inputs:
    %    axis - A 3-element vector representing the axis of rotation
    %    angle - The angle of rotation in radians
    %
    % Outputs:
    %    R - A 3x3 rotation matrix

    ux = axis(1);
    uy = axis(2);
    uz = axis(3);

    c = cos(angle);
    s = sin(angle);
    t = 1 - c;

    R = [t*ux*ux + c,    t*ux*uy - s*uz, t*ux*uz + s*uy;
         t*ux*uy + s*uz, t*uy*uy + c,    t*uy*uz - s*ux;
         t*ux*uz - s*uy, t*uy*uz + s*ux, t*uz*uz + c];
end
function movedPoints = movePointsFromCenter(points, distance)
    % Function to move each point in Cartesian coordinates away from the center
    % by a specified distance
    %
    % Parameters:
    % points - A 3xN matrix where each column is a point in 3D Cartesian space
    % distance - The distance to move each point away from the center

    % Calculate the current distances from the origin
    currentDistances = sqrt(sum(points.^2, 1));
    
    % Calculate the scale factors
    scaleFactors = (currentDistances + distance) ./ currentDistances;
    
    % Scale the points
    movedPoints = points .* scaleFactors;
end

function [Xdot, mu2,mu1,nu] = SlowWave_Model(~,x,~,k,E,F,C,r)
    
    % Parametrization of the path in the spherical parameter space in terms
    % of a circle defined by 3 points
    mu2=C(1)+r*(E(1)*cos(x(3))+F(1)*sin(x(3)));
    mu1=-(C(2)+r*(E(2)*cos(x(3))+F(2)*sin(x(3))));
    nu=C(3)+r*(E(3)*cos(x(3))+F(3)*sin(x(3)));
 
    % System
    xdot = - x(2);
    ydot = x(1)^3 - mu2*x(1) - mu1 - x(2)*( nu + x(1) + x(1)^2);
    zdot = k;
   
    Xdot = [xdot;ydot;zdot];
end

function [E,F,C,r]=Parametrization_3PointsCircle(p1,p2,p3)
    
    syms x y z
    
    V12=(p1-p2)/norm(p1-p2);
    V13=(p1-p3)/norm(p1-p3);
    
    n=cross(V12,V13);
    n=n/norm(n);

    dalpha=sum(p1.*n);
    alpha=n(1)*x+n(2)*y+n(3)*z-dalpha;

    hp1p2=p1+(p2-p1)/2; % punto medio
    dbeta=sum(V12.*hp1p2);
    beta=V12(1)*x+V12(2)*y+V12(3)*z-dbeta;

    hp1p3=p1+(p3-p1)/2; % punto medio
    dgamma=sum(V13.*hp1p3);
    gamma=V13(1)*x+V13(2)*y+V13(3)*z-dgamma;
    
    sx=solve(alpha==0,x);
    s1=subs(beta,x,sx);
    sy=solve(s1==0,y);
    s2=subs(gamma,x,sx);
    s2=subs(s2,y,sy);
    Cz=solve(s2==0,z);
    Cy=subs(sy,z,Cz);
    Cx=subs(sx,y,Cy);
    Cx=subs(Cx,z,Cz);
    C=[Cx, Cy, Cz];
    C=eval(C);     
    
    E=(p1-C)/norm(p1-C);
    F=-cross(E,n);
    r=norm(p1-C);

end

