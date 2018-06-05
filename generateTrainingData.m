% Generate data for the CNN
%
% Required: isetbio
%
%%
ieInit;

%% Parameter initialization
sFreq         = 4; 
nPCs          = 2;
fov           = 1;
sContrast     = 1;

freqRange      = [0.1, 1.5441]; % log space
contrastRange  = [-3.5, -1]; % log space
nImages = 3000; % For each data point, we generate two images, one with the harmonic and one without.

% Save the generated data?
saveFlag = true;
saveName = 'sampleSet';
% accuracy = zeros(numel(scanFreq), numel(scanContrast));

%% Set up the stimulus parameters
clear hparams
hparams(2)           = harmonicP;
hparams(2).freq      = sFreq;  % Set the Frequency
hparams(2).contrast  = sContrast;
hparams(1)           = hparams(2);
hparams(1).contrast  = 0;

sparams.fov = 1.5;

nTimeSteps = 20;
stimWeights = ones(1, nTimeSteps);

%% Set up cone mosaic parameters

integrationTime = 0.005;
sampleTimes = ((1:nTimeSteps) - 1) * integrationTime;   % Five ms integration time
nTrials    = 100;

%% Randomly generate parameters:
% scanFreq
% scanContrast
% aberrations
rng(1); % Set the random seed

% Sample in log space

scanFreq = 10.^((freqRange(2)-freqRange(1)).*rand(nImages,1) + freqRange(1));
scanContrast = 10.^((contrastRange(2)-contrastRange(1)).*rand(nImages,1) + contrastRange(1));

% Aberrations
measPupilMM = 4.5; % 4.5 mm pupil size (Thibos data)
calcPupilMM = 3; % 3 mm pupil size (what we want to calculate with)
zCoeffs = VirtualEyes(nImages,measPupilMM);

%% Generate images

samplesTemp = zeros(249,249,2*nImages);
samplesNoNoise = zeros(249,249,2*nImages);
trainLabels = zeros(2*nImages,1);
trainContrasts = zeros(2*nImages,1);
trainFreqs = zeros(2*nImages,1);

k = 1;

for ii = 1 : nImages
    
    fprintf('Generating image: %i \n',ii)
    %% Change the frequency and contrast for the stimulus
    hparams(2).freq      = scanFreq(ii);  % Set the Frequency
    hparams(2).contrast  = scanContrast(ii);
    
    %% Create the oi with aberrations
    
    %z = zeros(65,1);
    z(1:13) = zCoeffs(ii,1:13);
    
    % Create the example subject
    sbjWvf = wvfCreate;                                     % Initialize
    sbjWvf = wvfSet(sbjWvf,'zcoeffs',z);                    % Zernike
    sbjWvf = wvfSet(sbjWvf,'measured pupil',measPupilMM);   % Data
    sbjWvf = wvfSet(sbjWvf,'calculated pupil',calcPupilMM); % What we calculate
    sbjWvf = wvfSet(sbjWvf,'measured wavelength',550);
    sbjWvf = wvfSet(sbjWvf,'calc wave',[450:10:650]');            % Must be a column vector
    sbjWvf = wvfComputePSF(sbjWvf);
    oi = wvf2oi(sbjWvf);
    
    % Default
%     oi = oiCreate('wvf human');
    
    %% Create the OIS 
    
    ois = oisCreate('harmonic', 'blend', stimWeights, ...
        'testParameters', hparams, 'sceneParameters', sparams,...
        'oi',oi);
    % ois.visualize('movie illuminance');
    
    %% Set the coneMosaic parameters according to the OI
    cm = coneMosaic;
    cm.integrationTime = ois.timeStep;

    % Make the cm smaller than the oi size, but never smaller than 0.2 deg
    fovDegs = max(oiGet(ois.oiFixed,'fov') - 0.2, 0.2);  % Degrees
    cm.setSizeToFOV(fovDegs);
    
    %% Calculate the absorption template for the high contrast example of the stimulus

    cm.noiseFlag = 'random';
    samplesTemp(:,:,k) = mean(squeeze(cm.compute(ois)), 3);
    trainLabels(k) = 1;
    trainContrasts(k) = scanContrast(ii);
    trainFreqs(k) = scanFreq(ii);

    % Calculate without noise
    cm.noiseFlag = 'none';
    samplesNoNoise(:,:,k) = mean(squeeze(cm.compute(ois)), 3);
    
    k = k + 1;

    %% Create a "blank" pattern without stimulus
    hparams(2).contrast  = 0.0;
    ois = oisCreate('harmonic', 'blend', stimWeights, ...
        'testParameters', hparams, 'sceneParameters', sparams,...
        'oi',oi);
    % ois.visualize('movie illuminance');

    cm.noiseFlag = 'random';
    samplesTemp(:,:,k) = mean(squeeze(cm.compute(ois)), 3);
    trainLabels(k) = 0;
    trainContrasts(k) = 0;
    trainFreqs(k) = 0;
    k = k+1;
    
end

%% Crop
trainImages = samplesTemp(1:224,1:224,:);
trainImages_NoNoise = samplesNoNoise(1:224,1:224,:);
%% Save everything

if(saveFlag)
currDate = datestr(now,'mm-dd-yy_HH_MM');    
save(sprintf('%s_%s.mat',saveName,currDate),...
    'trainImages',...
    'trainImages_NoNoise',...
    'trainLabels',...
    'trainContrasts',...
    'trainFreqs',...
    '-v7.3');
end

%% Display a random sampling of images

n = 16;
figure(1);
[dataSamp,idx] = datasample(trainImages,n,3);
labelSamp = trainLabels(idx);
contrastSamp = trainContrasts(idx);
freqSamp = trainFreqs(idx);
k = 1;
for ii = 1:n
    subplot(4,4,k);
    imagesc(dataSamp(:,:,ii));
    axis image; axis off; 
    title(sprintf('label = %i \n c = %0.4f \n | f = %0.2f',...
        labelSamp(ii),contrastSamp(ii),freqSamp(ii)))
    k = k+1;
end

