createGUI();


function createGUI()
% Create the main figure window
fontsize = 10;
mindip = 0;
maxdip = 55;
textdip = 11;
maxLength = 125;
%%change name here

fig = uifigure('Position', [20, 50, 1300, 750], 'Name', 'SNIC/SH','Color', [1 1 1]);

Height = 735;
% A slider
%%change caption here 
%%to the left, height, length of text box, height of text box
uilabel(fig, 'Position', [250, Height, 500,15], 'Text', 'Model parameters', 'FontWeight', 'bold', 'FontSize', fontsize+1);
Height = Height - 20;
uilabel(fig, 'Position', [50,Height, 400, 22], 'Text', 'Onset curve index: SNIC', 'FontSize', fontsize, 'FontWeight', 'bold');
% Define the segments of the description
descriptionLines = {
    '';
    'Crossing the SNIC curve will result in onset of a signal with increasing frequency and fixed amplitude. ';

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
%% 
ASlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1, 44], 'Value', 1);
Height = Height - maxdip;

%B Slider
uilabel(fig, 'Position', [50,Height, 400, 22], 'Text', 'Offset curve index: SH', 'FontSize', fontsize, 'FontWeight', 'bold');
descriptionLines = {
'';
    'Crossing the SH curve will result in a signal terminating with decreasing frequency and DC shift. ';

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
BSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1, 54], 'Value', 44);



% Add labels for the sliders

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
kkSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [0.00, 0.1], 'Value', 0.007);
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
pinkNoiseSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [0, 2000], 'Value', 100);

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
TmaxSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1000, 10000], 'Value', 1300);
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
titleLabel = uilabel(fig, 'Position', [830, 700, 500, 40], 'Text', 'SNIC/SH - Parameter Control and Visualization', 'FontSize', 16, 'FontWeight', 'bold');
   %change description
descriptionLines = {

    'SNIC/SH seizures have a saddle node invariant cycle (SNIC) onset and a saddle homoclinic (SH) offset. In the x time series, the SNIC onset creates a square root scaling of the frequency and the SH offset creates a logarithmic scaling of the frequency. At seizure onset, no DC shift occurs and the spiking rate increases in frequency. At seizure offset, there is a DC shift and the spiking rate decreases in frequency. '
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
 
    uilabel(fig, 'Position', [620,Height, 850, lineHeight], 'Text', splitText{i}, 'FontSize', 12);
    Height = Height -  20;
    
end

% Add a button to run the simulation
runButton = uibutton(fig, 'Position', [630, 550, 150, 22], 'Text', 'Run Simulation', ...
    'ButtonPushedFcn', @(btn, event) runSimulation(ASlider.Value, BSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value, TmaxSlider.Value, FrequencySlider.Value, DriftSlider.Value, fig));
runButton = uibutton(fig, 'Position', [805, 550, 250, 22], 'Text', 'Run Simulation with postprocessing', ...
    'ButtonPushedFcn', @(btn, event) runppSimulation(ASlider.Value, BSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value, TmaxSlider.Value, FrequencySlider.Value, DriftSlider.Value, fig));
runButton = uibutton(fig, 'Position', [1080, 550, 150, 22], 'Text', 'Run Video Simulation', ...
    'ButtonPushedFcn', @(btn, event) runVideoSimulation(ASlider.Value, BSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value, TmaxSlider.Value, FrequencySlider.Value, DriftSlider.Value, fig));
buttomCaption = uilabel(fig, 'Position', [940, 520, 300, 22], 'Text', 'Note: May take up to 20 seconds to run', 'FontSize', fontsize);
runButton = uibutton(fig, 'Position', [610, 50, 100, 22], 'Text', 'Save Simulation', ...
    'ButtonPushedFcn', @(btn, event) get_mat_file(ASlider.Value, BSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value,TmaxSlider.Value,  FrequencySlider.Value, DriftSlider.Value, fig));
buttomCaption = uilabel(fig, 'Position', [610, 20, 300, 22], 'Text', 'Save generated timeseries as a .mat', 'FontSize', fontsize);
% Add the first plot area
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';

% Add the second plot area
ax2 = uiaxes(fig, 'Position', [950,100, 350,350]);
ax2.Tag = 'PlotArea2';

% Add a plot area to display the results

   savefig(fig, 'SNIC_SH_GUI.fig'); 
end

function runSimulation(AIndex, BIndex,  k, sigma, acq_noise, tmax, frequency, drift, fig)
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
    ax2 = uiaxes(fig, 'Position', [950,100, 350,350]);
ax2.Tag = 'PlotArea2';
    
    %%change matrix for onset offset based on what is in bifurcation one
    %%path


p2 = SNIC(:,floor(AIndex))';
p3 = [0.33, 0.11, 0.18];
p1 = SHl(:,floor(BIndex+50))';

if(p1 == p2)
title(ax1, 'Cannot create path if Onset curve are at same point');
    return;
end


    
    % Clear the existing plots
    cla(ax1);
    cla(ax2);
     
    
    % Clear the existing plots
   
    
x0=[0;0;0]; % initial conditions (must be a column)

% Settings - Model Sphere
b = 1.0; % focus
R = 0.4; % radius

% Class Information (paths)
% c2s 
N= 1;
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
[E,F,C,r]=Parametrization_3PointsCircle(p1,p2,p3);

N_t = length(tspan);
X = zeros(3,N_t);
xx = x0;
Rn =  [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 0);pinknoise([1,N_t],-1, 0)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));
for n = 1:N_t

    % Euler-Meruyama method
    [Fxx, mu2, mu1, nu] = SlowWave_Model(0,xx,0,k,E,F,C,r);
    xx = xx + tstep*Fxx + sqrt(tstep)*Rn(:,n);
    X(:,n) = xx;
    mu2_big(n) = mu2;
    mu1_big(n) = mu1;
    nu_big(n) = nu;
    

end

x = X';

xtemp = x;





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
seizure = xtemp(onset:offset,1);
z = xtemp(onset:offset,3);

% noise_percentage= acq_noise/100;
%      L = length(seizure);
%     noise = pinknoise([1,L],-1,5.69e3*noise_percentage); 
%    seizure = seizure + noise';

hold(ax1,'on')
% plot(ax1,xtemp(:,1));
% plot(ax1,Amp2);
plot(ax1,seizure,'color','#696969','LineWidth',1)
xlabel(ax1, 'Time')
ylabel(ax1, "Simulated Seizure")
title(ax1, 'Timeseries');
 


linewidth = 2;

  % % openfig('Map.fig');
 
 
hold(ax2, 'on')
% load("map_regions.mat");


  

%[Left,right  up/down]
%%%
lineVector = [-0.19,0.6,-0.07];
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
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'DisplayName', 'Offset curve' );
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--', 'DisplayName', 'Onset curve');
end
%load('top_half_sphere_lcs.mat');
%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.8, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

scale = 0.41/0.4;

scale = 0.42/0.4;
plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')


scatter3(ax2, p1(1)*scale,p1(2)*scale,p1(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker','square') % point on the offset curve
scatter3(ax2, p2(1)*scale,p2(2)*scale,p2(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker','diamond')
scatter3(ax2, p3(1)*scale,p3(2)*scale,p3(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker', 'pentagram')
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

function runVideoSimulation(AIndex, BIndex,  k, sigma, acq_noise, tmax, frequency, drift, fig)
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


if ishandle(ax2)
    delete(ax2);
end
ax2 = uiaxes(fig, 'Position', [950,100, 350,350]);
ax2.Tag = 'PlotArea2';

ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';
    
  if(k == 0)
    title(ax1, 'K or tmax value is too small for video animation, path does not traverse a full circle');
    return;

end
if(floor(((2*pi)/k)) > tmax)
      title(ax1, 'K or tmax value is too small for video animation, path does not traverse a full circle');
return;
end  

    %%change matrix for onset offset based on what is in bifurcation one
    %%path
p2 = SNIC(:,floor(AIndex))';
p3 = [0.33, 0.11, 0.18];
p1 = SHl(:,floor(BIndex+50))';
if(p1 == p2)
title(ax1, 'Cannot create path if Onset curve are at same point');
    return;
end

 if ishandle(ax1)
    delete(ax1);
end


if ishandle(ax2)
    delete(ax2);
end

img = imread('slow_wave_legend.png'); % Replace 'your_image.jpg' with the path to your JPG file
[img_height, img_width, ~] = size(img);

  

    % Create an axes component that covers the entire figure window
    ax2 = uiaxes(fig, 'Position', [850, 120, img_width/1.5, img_height/1.5]);

    % Display the image in the axes component
    imshow(img, 'Parent', ax2);
img = imread('new_slow_wave_paths_2d.png'); % Replace 'your_image.jpg' with the path to your JPG file
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


    ax2 = uiaxes(fig, 'Position', [950,100, 350,350]);
ax2.Tag = 'PlotArea2';


    % Clear the existing plots

     
    
    % Clear the existing plots
   
    
x0=[0;0;0]; % initial conditions (must be a column)

% Settings - Model Sphere
b = 1.0; % focus
R = 0.4; % radius

% Class Information (paths)
% c2s 
N= 1;
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
[E,F,C,r]=Parametrization_3PointsCircle(p1,p2,p3);

N_t = length(tspan);
X = zeros(3,N_t);
xx = x0;
Rn =  [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 0);pinknoise([1,N_t],-1, 0)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));
for n = 1:N_t

    % Euler-Meruyama method
    [Fxx, mu2, mu1, nu] = SlowWave_Model(0,xx,0,k,E,F,C,r);
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
targetLineVector = [-0.19,0.6,-0.07];

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
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'DisplayName', 'Offset curve' );
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--', 'DisplayName', 'Onset curve');
end
%load('top_half_sphere_lcs.mat');
%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

scale = 0.41/0.4;

scale = 0.42/0.4;
%plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scatter3(ax2, p1(1)*scale,p1(2)*scale,p1(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker','square') % point on the offset curve
scatter3(ax2, p2(1)*scale,p2(2)*scale,p2(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker','diamond')
scatter3(ax2, p3(1)*scale,p3(2)*scale,p3(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker', 'pentagram')
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
class = 6;
if k > 0
            plot_onset_offset = 1;
        point1 = p1' - C;
        point2 = p2' - C;
        point3 = p3'-C;
        point1 = point1 / norm(point1) * r;
        point2 = point2 / norm(point2) * r;
        
        point3 = point3 / norm(point3) * r;
        % Compute the quaternion for rotation
        theta1 = acos(dot(point1, point2) / (r^2));
        %%%change here
        numPoints1 = floor((theta1/k)/tstep);
        point = [mu2_big(numPoints1), -mu1_big(numPoints1), nu_big(numPoints1)];
        if round(point,2) == round(p2,2)'
        onset_index = numPoints1;
        else
        numPoints1 = floor(((2*pi - theta1)/k)/tstep);
        onset_index = numPoints1;
        theta1 = theta1-2*pi;
        end
        theta2 = acos(dot(point1, point3) / (r^2));
        numPoints2 = floor(((theta2)/k)/tstep);
        point = [mu2_big(numPoints2), -mu1_big(numPoints2), nu_big(numPoints2)];
        if round(point,2) == round(p2,2)'
        %offset_index = numPoints2;
        else
        numPoints2 = floor(((2*pi  - theta2)/k)/tstep);
        theta2 = 2*pi - theta2;
        offset_index = numPoints2;
        end
        theta3 = 2*pi;
        numPoints3 = floor(((theta3)/k)/tstep);
        point = [mu2_big(numPoints3), -mu1_big(numPoints3), nu_big(numPoints3)];
        offset_index = numPoints3;
        end
        if class == 15 || class == 12
             onset_index_temp = onset_index;
        onset_index = offset_index;
        offset_index = onset_index_temp + floor((((2*pi)/k)/tstep));
        end

onset_offset_range = ((2*pi/k)/tstep)/20;
onset_offset_range = floor((onset_offset_range)/interval)*interval;
interval_onset = floor(onset_index/interval)*interval + 1;
interval_onset_start = interval_onset - onset_offset_range;
interval_onset_stop = interval_onset + onset_offset_range;

interval_offset = floor(offset_index/interval)*interval + 1;
interval_offset_start = interval_offset - onset_offset_range;
interval_offset_stop = interval_offset + onset_offset_range;
onset_plot = false;
offset_plot = false;
seizure_plot = false;

plot_size = 2;
increase_plot_size = floor(((2*pi/k)/tstep)/interval)*interval +1;
for i = 1:interval:n
    % Handle the case for the last interval
    if i + interval - 1 <= n
        if i == interval_onset_start
         onset_plot = true;
        end
        if i == interval_onset_stop
           onset_plot = false;
           interval_onset = floor((interval_onset + (2*pi/k)/tstep)/interval)*interval +1;
           interval_onset_start = interval_onset - onset_offset_range;
           interval_onset_stop = interval_onset + onset_offset_range;
           seizure_plot = true;
        end
        if i == interval_offset_start
            
            seizure_plot = false;
         offset_plot = true;
        end
        if i == interval_offset_stop
           offset_plot = false;
           interval_offset = floor((interval_offset + (2*pi/k)/tstep)/interval)*interval + 1;
            interval_offset_start = interval_offset - onset_offset_range;
            interval_offset_stop = interval_offset + onset_offset_range;
            
        end
        if i == increase_plot_size
        increase_plot_size = floor((increase_plot_size+(2*pi/k)/tstep)/interval)*interval +1;
        plot_size = plot_size+2;
        end
        if onset_plot == true
        plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'Color',[0.957, 0.612, 0.204],  'LineStyle', '--', 'HandleVisibility', 'off');
        plot3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,'k','LineWidth',plot_size,  'HandleVisibility', 'off');
        elseif offset_plot == true
        plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'Color',[0.404, 0.702, 0.851],   'HandleVisibility', 'off');
        plot3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,'k','LineWidth',plot_size,  'HandleVisibility', 'off');
        elseif seizure_plot == true
        plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'Color',[0.894, 0.706, 0.831],  'HandleVisibility', 'off');
        plot3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,'k','LineWidth',plot_size,  'HandleVisibility', 'off');
        else
        plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'k');
        plot3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,'k','LineWidth',plot_size,  'HandleVisibility', 'off');
        end

    else
        plot(ax1,tspan(i:n), seizure(i:n), 'k');
        plot3(ax2, mu2_big(i:n,1)*scale,-mu1_big(i:n,1)*scale,nu_big(i:n,1)*scale,'k','LineWidth',plot_size,  'DisplayName', 'Bursting path')
    end
    pause(0.0000001)
end
%%%%change here
hold(ax2,'off');




[pks,locs] = findpeaks(seizure(floor(onset_index):floor(offset_index)), 'MinPeakProminence', 0.15);




sampling_rate = 100; % Hz
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
time = linspace(0, (num_samples-1)/new_sampling_frequency, num_samples); % Assuming time starts from 0 and ends at 1
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
text(ax1, x_position, y_position, str, 'Fontsize', 12, 'Color', 'k');
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
text(ax1, x_position, y_position, str, 'Fontsize', 12, 'Color', 'k');
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
text(ax1, x_position, y_position, str, 'Fontsize', 12, 'Color', 'k');
pause(2);
plot(ax1,time, seizure,'color','#696969','LineWidth',1)
xlabel(ax1, 'Adjusted time');
ylabel(ax1, "Simulated Seizure")
title(ax1, 'Timeseries');
 





end


function runppSimulation(AIndex, BIndex,  k, sigma, acq_noise, tmax, frequency, drift, fig)
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
    ax2 = uiaxes(fig, 'Position', [950,100, 350,350]);
ax2.Tag = 'PlotArea2';
      if(k == 0)
    title(ax1, 'K or tmax value is too small for post-processing, path does not traverse a full circle');
    return;

end



    %%change matrix for onset offset based on what is in bifurcation one
    %%path
p2 = SNIC(:,floor(AIndex))';
p3 = [0.33, 0.11, 0.18];
p1 = SHl(:,floor(BIndex+50))';
if(p1 == p2)
title(ax1, 'Cannot create path if Onset curve are at same point');
    return;
end


    
    % Clear the existing plots

     
    
    % Clear the existing plots
   
    
x0=[0;0;0]; % initial conditions (must be a column)

% Settings - Model Sphere
b = 1.0; % focus
R = 0.4; % radius

% Class Information (paths)
% c2s 
N= 1;
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
[E,F,C,r]=Parametrization_3PointsCircle(p1,p2,p3);

N_t = length(tspan);
X = zeros(3,N_t);
xx = x0;
Rn =  [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 0);pinknoise([1,N_t],-1, 0)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));
for n = 1:N_t

    % Euler-Meruyama method
    [Fxx, mu2, mu1, nu] = SlowWave_Model(0,xx,0,k,E,F,C,r);
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
lineVector = [-0.19,0.6,-0.07];
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
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'DisplayName', 'Offset curve' );
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth,'HandleVisibility', 'off');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--', 'DisplayName', 'Onset curve');
end
%load('top_half_sphere_lcs.mat');
%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

scale = 0.41/0.4;

scale = 0.42/0.4;
plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scatter3(ax2, p1(1)*scale,p1(2)*scale,p1(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker','square') % point on the offset curve
scatter3(ax2, p2(1)*scale,p2(2)*scale,p2(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker','diamond')
scatter3(ax2, p3(1)*scale,p3(2)*scale,p3(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker', 'pentagram')
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



onset = 1;
offset = tmax*100;
seizure = xtemp(onset:offset,1);

n = length(seizure);
%hold(ax2,'off');

% Compute the quaternion for rotation
point1 = p1 - C;
point2 = p2 - C;
point3 = p3-C;
point1 = point1 / norm(point1) * r;
point2 = point2 / norm(point2) * r;

point3 = point3 / norm(point3) * r;
% Compute the quaternion for rotation
theta1 = acos(dot(point1, point2) / (r^2));
%%%change here
numPoints1 = floor((theta1/k)/tstep);
point = [mu2_big(numPoints1), -mu1_big(numPoints1), nu_big(numPoints1)];
if round(point,2) == round(p2,2)
onset_index = numPoints1;
else
numPoints1 = floor(((2*pi - theta1)/k)/tstep);
onset_index = numPoints1;
theta1 = theta1-2*pi;
end
theta2 = acos(dot(point1, point3) / (r^2));
numPoints2 = floor(((theta2)/k)/tstep);
point = [mu2_big(numPoints2), -mu1_big(numPoints2), nu_big(numPoints2)];
if round(point,2) == round(p2,2)
%onset_index = numPoints;
else
numPoints2 = floor(((2*pi  - theta2)/k)/tstep);
theta2 = 2*pi - theta2;
%onset_index = numPoints;
end
theta3 = 2*pi;
numPoints3 = floor(((theta3)/k)/tstep);
point = [mu2_big(numPoints3), -mu1_big(numPoints3), nu_big(numPoints3)];

onset_index = numPoints3;
offset_index = numPoints2 + floor(((2*pi)/k)/tstep);


%%%%change here
%hold(ax2,'off');

if(offset_index > tmax*100)
      title(ax1, 'K or tmax value is too small for post-processing, path does not traverse a full circle');
return;
end 


[pks,locs] = findpeaks(seizure(floor(onset_index):floor(offset_index)), 'MinPeakProminence', 0.35);




sampling_rate = 100; % Hz
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
time = linspace(0, (num_samples-1)/new_sampling_frequency, num_samples); % Assuming time starts from 0 and ends at 1
hold(ax1,'off');
hpf = drift;
if drift > 0
HPF1 = designfilt('highpassiir', ...       % Response type
    'FilterOrder',1, ...     % Frequency constraints
       'HalfPowerFrequency',hpf, ...     % Frequency constraints
       'DesignMethod','butter', ...     % Design method
       'SampleRate',new_sampling_frequency);               % Sample rate



seizure= filter(HPF1,seizure);
end

noise_percentage= acq_noise/100;
     L = length(seizure);
    noise = pinknoise([1,L],-1,5.69e3*noise_percentage); 
   seizure = seizure + noise';




title(ax1, 'Adding signal acquisition noise');
plot(ax1,time, seizure,'color','#696969','LineWidth',1);
xlabel(ax1, 'Adjusted time');
ylabel(ax1, "Simulated Seizure");
title(ax1, 'Timeseries');
 





end
function get_mat_file(AIndex, BIndex,  k, sigma, acq_noise, tmax, frequency, drift, fig)
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
%change here to title for class


xlabel_text = ax1.XLabel.String;
if strcmp(xlabel_text, 'Adjusted time')
    param_text = sprintf('_Offset_index_%d_Onset_index_%d_k_%.2f_dynamical_noise_%.2f_acquisition_noise_%.2f_tMax_%.2f_adjusted_frequency_%.2f_electrode_drift_filter_%.2f', ...
        AIndex, BIndex, k, sigma, acq_noise, tmax, frequency, drift);
    

else

%change here to title for class
param_text = sprintf('_Offset_index_%d_Onset_index_%d_k_%.2f_dynamical_noise_%.2f_acquisition_noise_%.2f_tMax_%.2f', ...
        AIndex, BIndex, k, sigma, acq_noise, tmax);
end


% Save the data to a .mat file
save('SNIC_SH_timeseries.mat', 'Time', 'Signal', 'param_text')


title(ax1, 'Timeseries has saved succesfully!');

end
end



function Xdot = HysteresisLoop_Model(~,x,~,k,R,dstar,E,F,N)
    
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

function [E, F, C, r] = Parametrization_3PointsCircle(p1, p2, p3)
    % Calculate unit direction vectors
    p1 = p1';
    p2 = p2';
    p3 = p3';
    V12 = (p1 - p2) / norm(p1 - p2);
    V13 = (p1 - p3) / norm(p1 - p3);
    
    % Compute the normal vector to the plane defined by the points
    n = cross(V12, V13);
    n = n / norm(n); % Normalize the normal vector

    % Calculate the coefficients for the plane equations
    dalpha = dot(p1, n);
    dbeta = dot(V12, p1 + (p2 - p1) / 2);
    dgamma = dot(V13, p1 + (p3 - p1) / 2);

    % Set up the linear equations to find the center C
    A = [n(1), n(2), n(3);
         V12(1), V12(2), V12(3);
         V13(1), V13(2), V13(3)];

    b = [dalpha; dbeta; dgamma];

    % Solve for C using least squares
    C = A\b;

    % Calculate E (unit vector from C to p1)
    E = (p1 - C) / norm(p1 - C);

    % Calculate F (perpendicular vector)
    F = -cross(E, n);
    
    % Calculate the radius r
    r = norm(p1 - C);
end



