createGUI();


function createGUI()
% Create the main figure window
fontsize = 10;
mindip = 0;
maxdip = 55;
textdip = 11;
maxLength = 125;
%%change name here
addpath(fullfile('GUI helper files'));
fig = uifigure('Position', [20, 50, 1300, 750], 'Name', 'SupH/SupH', 'Color',[1,1,1]);
Height = 735;
% A slider
%%change caption here 
%%to the left, height, length of text box, height of text box
uilabel(fig, 'Position', [250, Height, 500,15], 'Text', 'Model parameters', 'FontWeight', 'bold', 'FontSize', fontsize+1);
Height = Height - 20;
uilabel(fig, 'Position', [50,Height, 400, 22], 'Text', 'Onset curve index: SupH', 'FontSize', fontsize, 'FontWeight', 'bold');
% Define the segments of the description
descriptionLines = {
    '';
    'Crossing the SupH curve will result in onset of a signal with fixed frequency and increasing amplitude. ';

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
ASlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1, 100], 'Value', 34);
Height = Height - maxdip;

%B Slider
uilabel(fig, 'Position', [50,Height, 400, 22], 'Text', 'Offset curve index: SupH', 'FontSize', fontsize, 'FontWeight', 'bold');
descriptionLines = {
'';
    'Crossing the SupH curve will result in a signal terminating with fixed frequency and decreasing amplitude. ';

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
BSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [1, 100], 'Value', 44);



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
kkSlider = uislider(fig, 'Position', [50,Height, 500, 3], 'Limits', [0.00, 0.1], 'Value', 0.004);
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
titleLabel = uilabel(fig, 'Position', [780, 700, 400, 40], 'Text', 'SupH/SupH - Parameter Control and Visualization', 'FontSize', 16, 'FontWeight', 'bold');
   %change description
descriptionLines = {

    'SupH/SupH seizures have a supercritical Hopf (SupH) onset bifurcation and a supercritical Hopf (SupH)  offset bifurcation. The path';
    '';
    'traverses though the seizure region and the rest region. The x time series does exhibits increasing amplitude at onset and decreasing'
    '';
   ' amplitude at offset. ';
   
};

% Calculate the number of lines for positioning
numLines = numel(descriptionLines);
lineHeight = 22; % Height of each line (adjust as needed)


    Height = 670 -  10;
% Position and display each line as a separate uilabel
for i = 1:numLines
 
    uilabel(fig, 'Position', [620, Height, 850, lineHeight], 'Text', descriptionLines{i}, 'FontSize', 10);
    Height = Height -  10;
    
end

% Add a button to run the simulation
runButton = uibutton(fig, 'Position', [630, 550, 150, 22], 'Text', 'Run Simulation', ...
    'ButtonPushedFcn', @(btn, event) runSimulation(ASlider.Value, BSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value,  FrequencySlider.Value, DriftSlider.Value, fig));
runButton = uibutton(fig, 'Position', [805, 550, 250, 22], 'Text', 'Run Simulation with postprocessing', ...
    'ButtonPushedFcn', @(btn, event) runppSimulation(ASlider.Value, BSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value,  FrequencySlider.Value, DriftSlider.Value, fig));
runButton = uibutton(fig, 'Position', [1080, 550, 150, 22], 'Text', 'Run Video Simulation', ...
    'ButtonPushedFcn', @(btn, event) runVideoSimulation(ASlider.Value, BSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value,  FrequencySlider.Value, DriftSlider.Value, fig));
buttomCaption = uilabel(fig, 'Position', [940, 520, 300, 22], 'Text', 'Note: May take up to 20 seconds to run', 'FontSize', fontsize);
runButton = uibutton(fig, 'Position', [610, 50, 100, 22], 'Text', 'Save Simulation', ...
    'ButtonPushedFcn', @(btn, event) get_mat_file(ASlider.Value, BSlider.Value, ...
     kkSlider.Value, pinkNoiseSlider.Value, acqSlider.Value,  FrequencySlider.Value, DriftSlider.Value, fig));
buttomCaption = uilabel(fig, 'Position', [610, 20, 300, 22], 'Text', 'Save generated timeseries as a .mat', 'FontSize', fontsize);
% Add the first plot area
ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';

% Add the second plot area
ax2 = uiaxes(fig, 'Position', [950, 100, 350,350]);
ax2.Tag = 'PlotArea2';

% Add a plot area to display the results

   savefig(fig, 'SupH_SupH_GUI.fig'); 
end

function runSimulation(AIndex, BIndex,  k, sigma, acq_noise,  frequency, drift, fig)
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
    
    %%change matrix for onset offset based on what is in bifurcation one
    %%path




    
    % Clear the existing plots
    cla(ax1);
    cla(ax2);
     
    
tstep = 0.01;
load('bifurcation_crossing.mat');
load("curves2.mat");

x0=[0;0;0];
index_array = [3,7,9,10,11];

b = 0;
    bifurcation = index_array(5);
   [p0,p1,p1_5,p2,p3]=Random_bifurcation_path(bifurcation);
%p1_5 =  Homoclinic_to_saddle2(:,20)';
p1 = Hopf(:,floor(650+AIndex))';
p2 = Hopf(:,floor(650+BIndex))';
% [r1,t1] = calculateArcLength(p0,p1,0.4);
% [r2,d2] = calculateArcLength(p1,p1_5,0.4);
% [r3, t2] = calculateArcLength(p1_5,p2,0.4);
% [r4,t4] = calculateArcLength(p2,p3,0.4);
% total_len = r1+r2+r3+r4;
% total_rad = t1+d2+t2+t4;
% % disp(p0)
% % % disp(p1);

% % disp(randomNumber2);
% % disp(p3)
% Class specific timescale separation


stall_val = 30000;
[mu2_straight_path0,mu1_straight_path0,nu_straight_path0,rad1] = sphereArcPath(k,tstep,p0,p1);
[mu2_straight_path0_5,mu1_straight_path0_5,nu_straight_path0_5,rad2] = sphereArcPath(k,tstep,p1,p1_5);
points = repmat(p1_5, stall_val, 1)';
sigma_temp = 0;
Rn = [pinknoise([1,length(points)],-1, sigma_temp);pinknoise([1,length(points)],-1, sigma_temp);pinknoise([1,length(points)],-1, sigma_temp)];
points = points + Rn;
[mu2_straight_path,mu1_straight_path,nu_straight_path,rad3] = sphereArcPath(k,tstep,p1_5,p2);
[mu2_straight_path1,mu1_straight_path1,nu_straight_path1,rad4] = sphereArcPath(k,tstep,p2,p3);
mu2_all = [mu2_straight_path0, mu2_straight_path0_5, points(1, :), mu2_straight_path, mu2_straight_path1];
mu1_all = [mu1_straight_path0, mu1_straight_path0_5, points(2, :), mu1_straight_path, mu1_straight_path1];
mu1_all = -mu1_all;
nu_all = [nu_straight_path0, nu_straight_path0_5, points(3,:), nu_straight_path, nu_straight_path1];

onset_time = [];
offset_time = [];
mu2_straight_path = [];
mu1_straight_path = [];
nu_straight_path = [];
num_seizures = 1;
N_t = length(mu1_all);
tmax_small = N_t*tstep;
for i = 1:num_seizures
 mu2_straight_path = [mu2_straight_path, mu2_all];
 mu1_straight_path = [mu1_straight_path, mu1_all];
 nu_straight_path = [nu_straight_path, nu_all];
 % onset_time = [onset_time, onset_index + (N_t*(i-1))];
 % offset_time = [offset_time, offset_index + (N_t*(i-1))];
end
tspan = 0:tstep:((tmax_small*num_seizures-0.01));
N_t = length(tspan);
tmax = tmax_small*num_seizures-0.01;
X = zeros(3,N_t);
xx = x0;

Rn = [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 00);pinknoise([1,N_t],-1, 00)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));
%[E,F] = Parametrization_2PointsArc(p2,p3,0.4);
for n = 1:N_t
[Fxx,mu2,mu1,nu] = SlowWave_Model2(tspan(n),xx,b,k,mu2_straight_path(n), mu1_straight_path(n),nu_straight_path(n));
xx = xx + tstep*Fxx + sqrt(tstep)*Rn(:,n);
X(:,n) = xx;
mu2_big(n) = mu2;
mu1_big(n) = mu1;
nu_big(n) = nu;
end
x = X';
t = tspan;

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
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'HandleVisibility', 'off' );
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth ,'DisplayName', 'Onset and Offset curve');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--','HandleVisibility', 'off');
end

%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.8, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

scale = 0.41/0.4;

scale = 0.42/0.4;
plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scatter3(ax2, p1(1)*scale,p1(2)*scale,p1(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker', 'diamond') % point on the offset curve
scatter3(ax2, p2(1)*scale,p2(2)*scale,p2(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker','square')
scatter3(ax2, p1_5(1)*scale,p1_5(2)*scale,p1_5(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker','pentagram')
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

function runVideoSimulation(AIndex, BIndex,  k, sigma, acq_noise,  frequency, drift, fig)
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
ax2 = uiaxes(fig, 'Position', [950, 100, 350,350]);
ax2.Tag = 'PlotArea2';

ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';
    


    %%change matrix for onset offset based on what is in bifurcation one
    %%path


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
img = imread('pw_paths2.png'); % Replace 'your_image.jpg' with the path to your JPG file
[img_height, img_width, ~] = size(img);
ax1 = uiaxes(fig, 'Position', [600, 80, img_width/1.5, img_height/1.5]);


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


    % Clear the existing plots

     
    
    % Clear the existing plots
   
    
x0=[0;0;0]; % initial conditions (must be a column)

% Settings - Model Sphere
b = 1.0; % focus
R = 0.4; % radius

% Class Information (paths)
% c2s 
N= 1;

% Class specific timespan



tstep = 0.01;
load('bifurcation_crossing.mat');
load("curves2.mat");

x0=[0;0;0];
index_array = [3,7,9,10,11];

b = 0;
    bifurcation = index_array(5);
   [p0,p1,p1_5,p2,p3]=Random_bifurcation_path(bifurcation);
%p1_5 =  Homoclinic_to_saddle2(:,20)';
p1 = Hopf(:,floor(800))';
p2 = Hopf(:,floor(650+BIndex))';
% [r1,t1] = calculateArcLength(p0,p1,0.4);
% [r2,d2] = calculateArcLength(p1,p1_5,0.4);
% [r3, t2] = calculateArcLength(p1_5,p2,0.4);
% [r4,t4] = calculateArcLength(p2,p3,0.4);
% total_len = r1+r2+r3+r4;
% total_rad = t1+d2+t2+t4;
% % disp(p0)
% % % disp(p1);

% % disp(randomNumber2);
% % disp(p3)
% Class specific timescale separation


stall_val = 30000;
[mu2_straight_path0,mu1_straight_path0,nu_straight_path0,rad1] = sphereArcPath(k,tstep,p0,p1);
[mu2_straight_path0_5,mu1_straight_path0_5,nu_straight_path0_5,rad2] = sphereArcPath(k,tstep,p1,p1_5);
points = repmat(p1_5, stall_val, 1)';
sigma_temp = 0;
Rn = [pinknoise([1,length(points)],-1, sigma_temp);pinknoise([1,length(points)],-1, sigma_temp);pinknoise([1,length(points)],-1, sigma_temp)];
points = points + Rn;
[mu2_straight_path,mu1_straight_path,nu_straight_path,rad3] = sphereArcPath(k,tstep,p1_5,p2);
[mu2_straight_path1,mu1_straight_path1,nu_straight_path1,rad4] = sphereArcPath(k,tstep,p2,p3);
mu2_all = [mu2_straight_path0, mu2_straight_path0_5, points(1, :), mu2_straight_path, mu2_straight_path1];
mu1_all = [mu1_straight_path0, mu1_straight_path0_5, points(2, :), mu1_straight_path, mu1_straight_path1];
mu1_all = -mu1_all;
nu_all = [nu_straight_path0, nu_straight_path0_5, points(3,:), nu_straight_path, nu_straight_path1];

onset_time = [];
offset_time = [];
mu2_straight_path = [];
mu1_straight_path = [];
nu_straight_path = [];
num_seizures = 1;
N_t = length(mu1_all);
tmax_small = N_t*tstep;
for i = 1:num_seizures
 mu2_straight_path = [mu2_straight_path, mu2_all];
 mu1_straight_path = [mu1_straight_path, mu1_all];
 nu_straight_path = [nu_straight_path, nu_all];
 % onset_time = [onset_time, onset_index + (N_t*(i-1))];
 % offset_time = [offset_time, offset_index + (N_t*(i-1))];
end
tspan = 0:tstep:((tmax_small*num_seizures-0.01));
N_t = length(tspan);
tmax = tmax_small*num_seizures-0.01;
X = zeros(3,N_t);
xx = x0;

Rn = [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 00);pinknoise([1,N_t],-1, 00)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));
%[E,F] = Parametrization_2PointsArc(p2,p3,0.4);
for n = 1:N_t
[Fxx,mu2,mu1,nu] = SlowWave_Model2(tspan(n),xx,b,k,mu2_straight_path(n), mu1_straight_path(n),nu_straight_path(n));
xx = xx + tstep*Fxx + sqrt(tstep)*Rn(:,n);
X(:,n) = xx;
mu2_big(n) = mu2;
mu1_big(n) = mu1;
nu_big(n) = nu;
end
x = X';
t = tspan;

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
targetLineVector = [-0.89,0.6,-0.07];

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
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'HandleVisibility', 'off' );
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth,  'DisplayName', 'Onset and Offset curve');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--','HandleVisibility', 'off');
end

%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

scale = 0.41/0.4;

scale = 0.42/0.4;
%plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scatter3(ax2, p1(1)*scale,p1(2)*scale,p1(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker', 'diamond') % point on the offset curve
scatter3(ax2, p2(1)*scale,p2(2)*scale,p2(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker','square')
scatter3(ax2, p1_5(1)*scale,p1_5(2)*scale,p1_5(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker','pentagram')
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
if k > 0
% Compute the quaternion for rotation

numPoints1 = floor((rad1/k)/tstep);
point = [mu2_big(numPoints1), -mu1_big(numPoints1), nu_big(numPoints1)];

onset_index = numPoints1;

theta2 = rad1+rad2;
numPoints2 = floor(((theta2)/k)/tstep);

offset_index = numPoints2;

theta3 = rad1+rad2+rad3;
numPoints3 = floor(((theta3)/k)/tstep)+length(points);
point = [mu2_big(numPoints3), -mu1_big(numPoints3), nu_big(numPoints3)];

%%comment in or out here based on fixed point location
offset_index = numPoints3;

% Define the interval

onset_offset_range = ((2*pi/k)/tstep)/20;
onset_offset_range = floor((onset_offset_range)/interval)*interval;
interval_onset = floor(onset_index/interval)*interval + 1;
interval_onset_start = interval_onset - onset_offset_range;
interval_onset_stop = interval_onset + onset_offset_range;

interval_offset = floor(offset_index/interval)*interval + 1;
interval_offset_start = interval_offset - onset_offset_range;
interval_offset_stop = interval_offset + onset_offset_range;

else
interval_onset_start = 0;
interval_onset_stop = 0;
interval_offset_start = 0;

interval_offset_stop = 0;
end
onset_plot = false;
offset_plot = false;
seizure_plot = false;
plot_size = 2;
increase_plot_size = floor(((2*pi/k)/tstep)/interval)*interval +1;
% Loop through the variable with the specified interval
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
        plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'Color',[0.4549 ,0.7490  ,  0.2706],  'LineStyle', '--', 'HandleVisibility', 'off');
        plot3(ax2, mu2_big(i:i+interval-1,1)*scale,-mu1_big(i:i+interval-1,1)*scale,nu_big(i:i+interval-1,1)*scale,'k','LineWidth',plot_size,  'HandleVisibility', 'off');
        elseif offset_plot == true
        plot(ax1,tspan(i:i+interval-1), seizure(i:i+interval-1), 'Color',[0.4549 ,0.7490  ,  0.2706],  'LineStyle', '--', 'HandleVisibility', 'off');
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

data = seizure;
min_val = min(data(:));
max_val = max(data(:));
data = (data - min_val) / (max_val - min_val);
rms_signal = get_amp(data);
normalized_data = data;
seizure = add_pink_noise(normalized_data, rms_signal, acq_noise/100);


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
 





end


function runppSimulation(AIndex, BIndex,  k, sigma, acq_noise,  frequency, drift, fig)
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
ax2 = uiaxes(fig, 'Position', [950, 100, 350,350]);
ax2.Tag = 'PlotArea2';

ax1 = uiaxes(fig, 'Position', [600, 100, 350,350]);
ax1.Tag = 'PlotArea1';
    


    
x0=[0;0;0]; % initial conditions (must be a column)

% Settings - Model Sphere
b = 1.0; % focus
R = 0.4; % radius

% Class Information (paths)
% c2s 
N= 1;



tstep = 0.01;
load('bifurcation_crossing.mat');
load("curves2.mat");

x0=[0;0;0];
index_array = [3,7,9,10,11];

b = 0;
    bifurcation = index_array(5);
   [p0,p1,p1_5,p2,p3]=Random_bifurcation_path(bifurcation);
%p1_5 =  Homoclinic_to_saddle2(:,20)';
p1 = Hopf(:,floor(650+AIndex))';
p2 = Hopf(:,floor(650+BIndex))';
% [r1,t1] = calculateArcLength(p0,p1,0.4);
% [r2,d2] = calculateArcLength(p1,p1_5,0.4);
% [r3, t2] = calculateArcLength(p1_5,p2,0.4);
% [r4,t4] = calculateArcLength(p2,p3,0.4);
% total_len = r1+r2+r3+r4;
% total_rad = t1+d2+t2+t4;
% % disp(p0)
% % % disp(p1);

% % disp(randomNumber2);
% % disp(p3)
% Class specific timescale separation


stall_val = 30000;
[mu2_straight_path0,mu1_straight_path0,nu_straight_path0,rad1] = sphereArcPath(k,tstep,p0,p1);
[mu2_straight_path0_5,mu1_straight_path0_5,nu_straight_path0_5,rad2] = sphereArcPath(k,tstep,p1,p1_5);
points = repmat(p1_5, stall_val, 1)';
sigma_temp = 0;
Rn = [pinknoise([1,length(points)],-1, sigma_temp);pinknoise([1,length(points)],-1, sigma_temp);pinknoise([1,length(points)],-1, sigma_temp)];
points = points + Rn;
[mu2_straight_path,mu1_straight_path,nu_straight_path,rad3] = sphereArcPath(k,tstep,p1_5,p2);
[mu2_straight_path1,mu1_straight_path1,nu_straight_path1,rad4] = sphereArcPath(k,tstep,p2,p3);
mu2_all = [mu2_straight_path0, mu2_straight_path0_5, points(1, :), mu2_straight_path, mu2_straight_path1];
mu1_all = [mu1_straight_path0, mu1_straight_path0_5, points(2, :), mu1_straight_path, mu1_straight_path1];
mu1_all = -mu1_all;
nu_all = [nu_straight_path0, nu_straight_path0_5, points(3,:), nu_straight_path, nu_straight_path1];

numPoints1 = floor((rad1/k)/tstep);
onset_index = numPoints1;
theta3 = rad1+rad2+rad3;
numPoints3 = floor(((theta3)/k)/tstep);
offset_index = numPoints3;
onset_time = [];
offset_time = [];
mu2_straight_path = [];
mu1_straight_path = [];
nu_straight_path = [];
num_seizures = 1;
N_t = length(mu1_all);
tmax_small = N_t*tstep;
for i = 1:num_seizures
 mu2_straight_path = [mu2_straight_path, mu2_all];
 mu1_straight_path = [mu1_straight_path, mu1_all];
 nu_straight_path = [nu_straight_path, nu_all];
 % onset_time = [onset_time, onset_index + (N_t*(i-1))];
 % offset_time = [offset_time, offset_index + (N_t*(i-1))];
end
tspan = 0:tstep:((tmax_small*num_seizures-0.01));
N_t = length(tspan);
tmax = tmax_small*num_seizures-0.01;
X = zeros(3,N_t);
xx = x0;
Rn = [pinknoise([1,N_t],-1, sigma);pinknoise([1,N_t],-1, 00);pinknoise([1,N_t],-1, 00)];
mu2_big = zeros(1, length(tspan));
mu1_big = zeros(1, length(tspan));
nu_big = zeros(1, length(tspan));
%[E,F] = Parametrization_2PointsArc(p2,p3,0.4);
for n = 1:N_t
[Fxx,mu2,mu1,nu] = SlowWave_Model2(tspan(n),xx,b,k,mu2_straight_path(n), mu1_straight_path(n),nu_straight_path(n));
xx = xx + tstep*Fxx + sqrt(tstep)*Rn(:,n);
X(:,n) = xx;
mu2_big(n) = mu2;
mu1_big(n) = mu1;
nu_big(n) = nu;
end
x = X';
t = tspan;

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
targetLineVector = [-0.89,0.6,-0.07];

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
plot3(ax2,Homoclinic_to_saddle_scaled(1, :), Homoclinic_to_saddle_scaled(2, :), Homoclinic_to_saddle_scaled(3, :), 'Color', [0.404, 0.702, 0.851],  'LineWidth', linewidth, 'HandleVisibility', 'off' );
plot3(ax2,Fold_scaled(1, 140:564), Fold_scaled(2, 140:564), Fold_scaled(3, 140:564), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 575:end), Fold_scaled(2, 575:end), Fold_scaled(3, 575:end), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 1:80), Fold_scaled(2, 1:80), Fold_scaled(3, 1:80), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'HandleVisibility', 'off');
plot3(ax2,Fold_scaled(1, 565:575), Fold_scaled(2,  565:575), Fold_scaled(3,  565:575), 'Color', [0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', ':', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 1:400), Hopf_scaled(2, 1:400), Hopf_scaled(3, 1:400), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth, 'LineStyle', '--', 'HandleVisibility', 'off');
plot3(ax2,Hopf_scaled(1, 400:973), Hopf_scaled(2, 400:973), Hopf_scaled(3, 400:973), 'Color',[0.4549 ,0.7490  ,  0.2706], 'LineWidth', linewidth,  'DisplayName', 'Onset and Offset curve');
plot3(ax2,SNIC_scaled(1, :), SNIC_scaled(2, :), SNIC_scaled(3, :), 'Color',[0.957, 0.612, 0.204], 'LineWidth', linewidth, 'LineStyle', '--','HandleVisibility', 'off');
end

%surf(ax2,X, Y,Z,'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none',  'DisplayName', 'LCs region');
surf(ax2,(X_sphere), (Y_sphere), (Z_sphere), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scale = 0.41/0.4;

scale = 0.42/0.4;
plot3(ax2, mu2_big*scale,-mu1_big*scale,nu_big*scale,'k','LineWidth',5,  'DisplayName', 'Bursting path')
scatter3(ax2, p1(1)*scale,p1(2)*scale,p1(3)*scale,100, 'MarkerFaceColor', 'k', 'DisplayName', 'Onset point', 'Marker', 'diamond') % point on the offset curve
scatter3(ax2, p2(1)*scale,p2(2)*scale,p2(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Offset point', 'Marker','square')
scatter3(ax2, p1_5(1)*scale,p1_5(2)*scale,p1_5(3)*scale,100,  'MarkerFaceColor', 'k', 'DisplayName', 'Fixed point', 'Marker','pentagram')

lgd = legend(ax2, 'show');
lgd.FontSize = 6;
%lgd.Position = [0.35, 0.75, 0.1, 0.1];

currentLineVector = targetLineVector;

   xlabel(ax2, 'mu2');
    ylabel(ax2, '-mu1');
    zlabel(ax2, 'nu');
  

    
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









hpf = drift;
if drift > 0
HPF1 = designfilt('highpassiir', ...       % Response type
    'FilterOrder',2, ...     % Frequency constraints
       'HalfPowerFrequency',hpf, ...     % Frequency constraints
       'DesignMethod','butter', ...     % Design method
       'SampleRate',new_sampling_frequency);               % Sample rate
seizure= filter(HPF1,seizure);
end




data = seizure;
min_val = min(data(:));
max_val = max(data(:));
data = (data - min_val) / (max_val - min_val);
rms_signal = get_amp(data);
normalized_data = data;
seizure = add_pink_noise(normalized_data, rms_signal, acq_noise/100);




plot(ax1,time, seizure,'color','#696969','LineWidth',1)
xlabel(ax1, 'Adjusted time');
ylabel(ax1, "Simulated Seizure")
title(ax1, 'Timeseries');
 










end
function get_mat_file(AIndex, BIndex,  k, sigma, acq_noise,  frequency, drift, fig)
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
    param_text = sprintf('Offset_index_%d_Onset_index_%d_k_%.2f_dynamical_noise_%.2f_acquisition_noise_%.2f_tMax_%.2f_adjusted_frequency_%.2f_electrode_drift_filter_%.2f', ...
        BIndex, AIndex, k, sigma, acq_noise,  frequency, drift);
    

else

%change here to title for class
param_text = sprintf('Offset_index_%d_Onset_index_%d_k_%.2f_dynamical_noise_%.2f_acquisition_noise_%.2f', ...
        BIndex, AIndex, k, sigma, acq_noise);
end


% Save the data to a .mat file
save('Sup_Sup_timeseries.mat', 'Time', 'Signal', 'param_text')



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

function [Xdot, mu2,mu1,nu] = SlowWave_Model2(~,x,~,k,mu2,mu1,nu)
% Parametrization of the path in the spherical parameter space in terms
% of a circle defined by 3 points
% System
xdot = - x(2);
ydot = x(1)^3 - mu2*x(1) - mu1 - x(2)*( nu + x(1) + x(1)^2);
zdot = k;
Xdot = [xdot;ydot;zdot];
end
function [mu2,mu1,nu, theta] = sphereArcPath(k, tstep,point1, point2)
% sphereArcPath - Generates an arc path between two points on a sphere
%
% Syntax: arcPath = sphereArcPath(point1, point2, numPoints)
%
% Inputs:
% point1 - [x1, y1, z1] Coordinates of the first point on the sphere
% point2 - [x2, y2, z2] Coordinates of the second point on the sphere
% numPoints - Number of points along the arc
%
% Outputs:
% arcPath - An Nx3 matrix containing the coordinates of points along the arc
% Check the input points
radius = 0.4;
% if norm(point1) ~= radius || norm(point2) ~= radius
% error('The points must lie on the sphere of radius 0.4.');
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
numPoints = floor((theta/k)/tstep);
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

function point= get_random_point
radius = 0.4;
% Generate two random numbers
theta = 2 * pi * rand(); % Random angle between 0 and 2*pi
phi = acos(2 * rand() - 1); % Random angle between 0 and pi
% Convert spherical coordinates to Cartesian coordinates
x = radius * sin(phi) * cos(theta);
y = radius * sin(phi) * sin(theta);
z = radius * cos(phi);
point = [x,y,z];
% Display the point
end
function point= get_random_point_hopf
%load("map_regions.mat")
radius = 0.4;
% Loop until a valid point with y > 0 is found
while true
% Generate two random numbers
theta = 2 * pi * rand(); % Random angle between 0 and 2*pi
phi = acos(2 * rand() - 1); % Random angle between 0 and pi
% Convert spherical coordinates to Cartesian coordinates
x = radius * sin(phi) * cos(theta);
y = radius * sin(phi) * sin(theta);
z = radius * cos(phi);
% Check if y is positive
if y > 0
point = [x, y, z];
break;
end
end
% Display the point
end
% 
function point= get_random_point_c9
data = [
%0.157119269030135	0.161091875334199	0.330700684911216
-0.159101878627765	0.366197019600565	0.0242143563362679
0.240325089686399	0.304071689348835	-0.0989154133781482
0.380078172355370	0.116531602029670	-0.0442828254225105
0.376085815202690	0.0636356567306722	-0.120457306942271
0.166174992977878	0.327294142268742	0.158947840958507
0.172012914508465	0.357581320788856	-0.0504693596670386
-0.0758915517458272	0.306208966174370	0.245919786532182
-0.0486739184007403	0.108823597501720	-0.381822307224580
0.0796165762275659	0.325569951814250	-0.218324087689533
0.251981235011323	0.310614829711489	0.00488720425919465
-0.232117006686916	0.282224648390582	-0.162698933762698
%0.220275294465348	0.201237672469794	0.266424837106348
-0.332207536237643	0.183689008934820	0.126081326386816
%0.147579315557777	0.0266736449114354	0.370821874067101
%-0.243972195896027	0.0461366053360331	0.313606411410446
0.312005760944042	0.232343423323037	-0.0931071360114856
%-0.183233489251201	0.115164818276131	0.336396422466494
%0.161849841530673	0.129909400887954	0.341947622242689
-0.285057819447844	0.276172527226907	-0.0497068886246541
0.293516972528577	0.262050399776243	0.0719539770601309
0.100452183568069	0.133430682841316	-0.363463356740273
%0.0587657128391141	0.273279906407037	0.286120051287343
0.135209543713507	0.118213069098463	-0.357412995822910
%-0.399343120628290	0.0156914048904897	-0.0166988568297806
-0.208541275002113	0.323171146651893	0.109867859686961
-0.0449185896441651	0.258992664847436	-0.301504759265261
0.0449357963445848	0.160975079585810	-0.363411334384634
0.0462987815061658	0.0359526389944827	-0.395681476165333
-0.287432027169657	0.213474519904710	-0.178357671852362
%-0.0334309018751504	0.251955318843665	0.308870348376469
0.211875145592833	0.322495970331763	-0.105381553413254
0.350411490578410	0.149947999830475	0.121356436242404
-0.117148812770085	0.265217051374884	-0.275565003813205
%-0.190205660482966	0.250539203555510	0.247086855583233
0.232946913731378	0.315457300993538	0.0788823594531643
0.0427771439339962	0.351925488543204	-0.185252709752050
0.296164473879145	0.244426175265017	0.111993076827116
0.142393319279842	0.371620981937393	0.0402739172214049
0.111421877013698	0.0969703711017663	-0.371728277712538
%0.384674464390556	0.0201207759882954	0.107799400830745
0.0946871119121513	0.381474943125708	0.0742375821600734
%-0.359866162303914	0.0769286748857566	0.156774756287403
0.0536941087296937	0.387454928395370	0.0836398299248805
%0.0544849486836204	0.0436414639382243	0.393861413434055
%-0.128741399122706	0.0344709001636801	0.377143751365224
0.288413853864687	0.240285781103050	0.138131069273411
%-0.0175431145198549	0.0209970586931555	0.399063106111276
%-0.0133870427997347	0.262216240418017	0.301767179040596
-0.0940903089878144	0.357288847540074	0.153270000907104
];
idx = randi([1,length(data)]);
point = data(idx,:);
end
function point= get_random_point_c10
data = [
    0.139443244717430	0.180180368497788	-0.328771373922177;
-0.123686721647726	0.338825918756816	-0.172912092308889;
-0.100168545573643	0.0377628990918715	-0.385409166899074;
0.377957700206015	0.0522413891185729	-0.120078366569263;
0.0670217468691031	0.393741391170461	-0.0218128935669752;
0.165551316967523	0.0423021044581705	0.361667379519813;
-0.190629713140052	0.100956538347682	-0.336850248379001;
0.237742614103209	0.317715346392282	-0.0503528361287097;
0.133352560213168	0.205062052013036	0.316491152338887;
0.116930014179559	0.380594524578872	0.0384080674409766;
    0.247294091034076	0.313367115323397	0.0254299739190403;
0.319668729544415	0.149389595714773	0.188400244277008;
%-0.265630076158132	0.280893029748762	-0.102663374573342;
0.206794351734318	0.328104632230012	0.0978950785279422;
    0.0778563146462470	0.373124628962106	-0.121311192936165;
0.0369732031123261	0.391492997572173	0.0732544544964284;
       -0.1050    0.3835   -0.0436;
        0.2092    0.3125    0.1365;
        0.3456    0.1805    0.0896;
    0.3221    0.2370    0.0071;
%0.0963164157568826	0.388100733721629	0.0100483103296496;
%-0.0137890033618096	0.234736675106870	-0.323587015725384;
-0.198974692220617	0.0411752373414647	-0.344548504111663;
];
idx = randi([1,length(data)]);
point = data(idx,:);
end
function point= get_random_point_SupH_SupH
data = [
0.0373, 0.2497, -0.3102;
-0.0441, 0.2591, -0.3015;
-0.2104, 0.3180, -0.1209;
% 0.3475, 0.1347, 0.1452;
0.0806, 0.1560, -0.3594;
% 0.1743, 0.2562, 0.2529;
-0.0564, 0.3820, 0.1042;
0.3196, 0.2389, -0.0279;
-0.0564, 0.3820, 0.1042;
0.3196, 0.2389, -0.0279;
0.2394, 0.3202, -0.0114;
0.0612, 0.3715, 0.1350;
0.3186, 0.1795, -0.1620;
-0.1300, 0.2586, -0.2761;
0.1833, 0.2218, -0.2779;
-0.3059, 0.1997, -0.1629;
0.2071, 0.2684, -0.2124;
-0.0467, 0.3294, -0.2220;
-0.2104, 0.3180, -0.1209;
-0.2316, 0.3055, 0.1142;
0.1314, 0.3298, -0.1843
];
idx = randi([1,length(data)]);
point = [-0.2104, 0.3180, -0.1209];
end
% Display the point
function point= get_random_point_fixed
radius = 0.4;
% Loop until a valid point with y > 0 is found
while true
% Generate two random numbers
theta = 2 * pi * rand(); % Random angle between 0 and 2*pi
phi = acos(2 * rand() - 1); % Random angle between 0 and pi
% Convert spherical coordinates to Cartesian coordinates
x = radius * sin(phi) * cos(theta);
y = radius * sin(phi) * sin(theta);
z = radius * cos(phi);
% Check if y is positive
if y < 0 
point = [x, y, z];
break;
end
end
% Display the point
end
function bool= Bif_bool(mu2,mu1,nu)
load('curves2.mat');
load('bifurcation_crossing.mat');
% plot3(FLC(1, :), FLC(2, :), FLC(3, :), 'm-', 'LineWidth', 2);
% %
% plot3(Fold(1, 1:95), Fold(2, 1:95), Fold(3, 1:95), 'm-', 'LineWidth', 2);
% plot3(Fold(1, 145:end), Fold(2, 145:end), Fold(3, 145:end), 'm-', 'LineWidth', 2);
% plot3(SNIC(1, :), SNIC(2, :), SNIC(3, :), 'g-', 'LineWidth', 2);
% plot3(Hopf(1, :), Hopf(2, :), Hopf(3, :), 'm-', 'LineWidth', 2);
% %
%
% plot3(Homoclinic_to_saddle(1, :), Homoclinic_to_saddle(2, :), Homoclinic_to_saddle(3, :), 'm-', 'LineWidth', 2);
% plot3(Homoclinic_to_saddle1(1, :), Homoclinic_to_saddle1(2, :), Homoclinic_to_saddle1(3, :), 'm-', 'LineWidth', 2);
% plot3(Homoclinic_to_saddle2(1, :), Homoclinic_to_saddle2(2, :), Homoclinic_to_saddle2(3, :), 'm-', 'LineWidth', 2);
% plot3(Homoclinic_to_saddle3(1, 1:90), Homoclinic_to_saddle3(2, 1:90), Homoclinic_to_saddle3(3, 1:90), 'm-', 'LineWidth', 2);
bool = false;
%no_hopf
matrix_cell = [FLC, Fold(:,1:95), Fold(:,145:682), SNIC, Homoclinic_to_saddle, Homoclinic_to_saddle1, Homoclinic_to_saddle2, Homoclinic_to_saddle3(:,1:90) ]';
points = [mu2;-mu1;nu]';
%Homoclinic_to_saddle3, Homoclinic_to_saddle2, Homoclinic_to_saddle1, Homoclinic_to_saddle
for i = 1:size(points, 1)
point_to_check = points(i, :);
% Initialize a flag for match
match_found = false;
% Loop through matrix_cell to check for a match
for j = 1:size(matrix_cell, 1)
if all(round(matrix_cell(j, :), 2) == round(point_to_check, 2))
match_found = true;
break; % Exit the loop if a match is found
end
end
if match_found
bool = true;
break; % Exit the outer loop if a match is found
end
end
end
function point= get_nearest_hopf(mu2,mu1,nu)
load('curves2.mat');
load('bifurcation_crossing.mat');
% plot3(FLC(1, :), FLC(2, :), FLC(3, :), 'm-', 'LineWidth', 2);
% %
% plot3(Fold(1, 1:95), Fold(2, 1:95), Fold(3, 1:95), 'm-', 'LineWidth', 2);
% plot3(Fold(1, 145:end), Fold(2, 145:end), Fold(3, 145:end), 'm-', 'LineWidth', 2);
% plot3(SNIC(1, :), SNIC(2, :), SNIC(3, :), 'g-', 'LineWidth', 2);
% plot3(Hopf(1, :), Hopf(2, :), Hopf(3, :), 'm-', 'LineWidth', 2);
% %
%
% plot3(Homoclinic_to_saddle(1, :), Homoclinic_to_saddle(2, :), Homoclinic_to_saddle(3, :), 'm-', 'LineWidth', 2);
% plot3(Homoclinic_to_saddle1(1, :), Homoclinic_to_saddle1(2, :), Homoclinic_to_saddle1(3, :), 'm-', 'LineWidth', 2);
% plot3(Homoclinic_to_saddle2(1, :), Homoclinic_to_saddle2(2, :), Homoclinic_to_saddle2(3, :), 'm-', 'LineWidth', 2);
% plot3(Homoclinic_to_saddle3(1, 1:90), Homoclinic_to_saddle3(2, 1:90), Homoclinic_to_saddle3(3, 1:90), 'm-', 'LineWidth', 2);
Heatmap = 0;
matrix_cell = {Hopf(:,450:650)};
%Homoclinic_to_saddle3, Homoclinic_to_saddle2, Homoclinic_to_saddle1, Homoclinic_to_saddle
Proximity_vals = zeros(1, length(matrix_cell));
point = [mu2,-mu1,nu];
for i = 1:length(matrix_cell)
matrix = matrix_cell{i};
% Calculate the Euclidean distance from the point to each point in the matrix
%distances = sqrt(sum((matrix - point').^2, 1));
distances = sqrt((matrix(1,:)- point(1)).^2 + (matrix(2,:)- point(2)).^2 + (matrix(3,:)- point(3)).^2);
% Find the minimum distance
Proximity_vals(i) =min(distances);
% Output the minimum distance and the corresponding closest point
%closest_point = matrix(:, idx);
end
[min_dist, idx] = min(Proximity_vals);
point = Hopf(:,idx+450);
end
function [arcLength, theta] = calculateArcLength(P1, P2, radius)
 % calculateArcLength computes the arc length and central angle between two points on a sphere.
 % 
 % Input:
 % P1 - First point [x1, y1, z1]
 % P2 - Second point [x2, y2, z2]
 % radius - Radius of the sphere (default: 0.4 if not provided)
 %
 % Output:
 % arcLength - Arc length between the two points
 % theta - Central angle between the two points in radians
 
 if nargin < 3
 radius = 0.4;
 end
 
 % Compute the dot product of P1 and P2
 dotProduct = dot(P1, P2);
 
 % Compute the magnitudes of P1 and P2
 magnitudeP1 = norm(P1);
 magnitudeP2 = norm(P2);
 
 % Compute the cosine of the central angle
 cosTheta = dotProduct / (magnitudeP1 * magnitudeP2);
 
 % Compute the central angle in radians
 theta = acos(cosTheta);
 
 % Compute the arc length
 arcLength = radius * theta;
end

function [mu2,mu1,nu] = get_seizure_points(P1,P2,P3,k,tstep)
[E,F,C,r]=Parametrization_3PointsCircle(P1,P2,P3);
point1 = P1 - C;
point2 = P2 - C;
point3 = P3-C;
point1 = point1 / norm(point1) * r;
point2 = point2 / norm(point2) * r;

point3 = point3 / norm(point3) * r;
% Compute the quaternion for rotation
theta = acos(dot(point1, point3) / (r^2));
theta = 2*pi - theta;
numPoints = floor((theta/k)/tstep);
% theta = acos(dot(point2, point3) / (r^2));
% numPoints2 = floor((theta/k)/tstep);
N_t = floor(numPoints);
xx = 0;


X= zeros(1, N_t);
mu2_big = zeros(1, N_t);
mu1_big = zeros(1, N_t);
nu_big = zeros(1, N_t);
%[E,F] = Parametrization_2PointsArc(p2,p3,0.4);
for n = 1:N_t
[Fxx, mu2,mu1,nu] = SlowWave_Model(0,xx,0,k,E,F,r);
xx = xx + tstep*Fxx ;
X(:,n) = xx;

mu2_big(n) = mu2;
mu1_big(n) = mu1;
nu_big(n) = nu;
end
mu2_big = mu2_big';
mu1_big = mu1_big';
nu_big = nu_big';
z =X';
mu2=C(1)+r*(E(1)*cos(z)+F(1)*sin(z));
mu1=-(C(2)+r*(E(2)*cos(z)+F(2)*sin(z)));
nu=C(3)+r*(E(3)*cos(z)+F(3)*sin(z));
% mu1 = -mu1;
end

function   [p0,p1,p1_5,p2,p3]=Random_bifurcation_path(bifurcation)

load('curves.mat');
load('bifurcation_crossing.mat')
load("curves2.mat")

if bifurcation==3
    %fixed rest point
    p0 = Hopf(:,930)';
    %bifurcation curve
    randomNumber = randi([145,170]);
    p1 = Fold(:,randomNumber)';
    randomNumber2 = randi([600,750]);
    p1_5 = get_random_point_hopf();
    %bifurcation curve
    p2 = Hopf(:,randomNumber2)' ;
    %fixed rest
    p3 = Hopf(:,930)';

end

if bifurcation==7
    randomNumber = randi([600,750]);
 randomNumber2 = randi([1,44]);
%fixed rest point
p0 = Homoclinic_to_saddle2(:,30)';
%bifurcation curve
p1 =  SNIC(:,randomNumber2)' ;
%random point in limit cycle
p1_5 = get_random_point_c9();
%bifurcation curve
p2 = Hopf(:,randomNumber)';
%fixed rest
p3 = [ 0.1944 , 0.0893 , 0.3380];
end

if bifurcation==9
    randomNumber = randi([600,750]);
    %fixed rest point
    p0 = [ 0.1944 , 0.0893 , 0.3380];
    %bifurcation curve
    p1 = Hopf(:,randomNumber)';
    %change here
    randomNumber2 = randi([1,44]);
    %random point in limit cycle
    p1_5 = get_random_point_c9();
    %bifurcation curve
    p2 = SNIC(:,randomNumber2)' ;
    %fixed rest
    p3 = Fold(:,450)';
end

if bifurcation==10
    randomNumber = randi([600,750]);
    %fixed rest point
    p0 = [ 0.1944 , 0.0893 , 0.3380];
    %bifurcation curve
    p1 = Hopf(:,randomNumber)';%get_nearest_hopf(p0(1),p0(2),p0(3))';
    %change here
    randomNumber2 = randi([1,124]);
    %random point in limit cycle
    p1_5 = get_random_point_c10();
    %bifurcation curve
    p2 = Homoclinic_to_saddle(:,randomNumber2)' ;
    %fixed rest
    p3 = SNIC(:,30)';
end

if bifurcation==11
     randomNumber = randi([600,750]);
%fixed rest point
p0 = [ 0.1944 , 0.0893 , 0.3380];
%bifurcation curve
p1 = Hopf(:,randomNumber)';
%change here
randomNumber2 = randi([600,750]);
%random point in limit cycle
p1_5 = get_random_point_SupH_SupH();
%bifurcation curve
p2 = Hopf(:,randomNumber2)' ;
%fixed rest
p3 = [ 0.1944 , 0.0893 , 0.3380];
end



end






function noisy_signal = add_pink_noise(signal, rms_signal, noise_amplitude_ratio)
    % Inputs:
    % signal - input signal (1D array)
    % noise_amplitude_ratio - fraction of signal amplitude for noise (e.g., 0.4 for 40%)
    
    % Compute the RMS amplitude of the signal
    
    
    % Generate pink noise of the same length as the signal
    % Pink noise can be generated using dsp.ColoredNoise in MATLAB
    L = length(signal);
    pink_noise = pinknoise([1,L],-1,10000)';
    
    % % Compute the RMS amplitude of the pink noise
     rms_noise = get_amp(pink_noise);
    % 
    % % % Scale the noise so its amplitude is noise_amplitude_ratio of the signal's amplitude
    % scaling_factor = noise_amplitude_ratio * (1 / rms_noise);
    % scaled_noise = pink_noise * scaling_factor;
    min_val = min(pink_noise(:));
max_val = max(pink_noise(:));
scaled_noise = 2*noise_amplitude_ratio*(pink_noise - min_val) / (max_val - min_val);
    % Add the scaled noise to the original signal
    noisy_signal = signal + scaled_noise;
end

function amp = get_amp(signal)

[peaks,locs] = findpeaks(signal ,'MinPeakProminence', 0.55);
[troughs_neg,locs_troughs] = findpeaks(signal, 'MinPeakProminence', 0.55);
troughs = -1*troughs_neg;

newnew = [];
len = 0;
if length(troughs) > length(peaks)
    len = length(peaks);
else
    len = length(troughs);
end
for i = 1:len
    newnew = [newnew; -1*troughs(i) + peaks(i)];
end
amp = mean(newnew);
end