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

freqRange      = [0 30]; 
contrastRange  = [0 1]; 
nImages = 1000; % For each data point, we generate two images, one with the harmonic and one without.

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

sparams.fov = 1;

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

scanFreq = (freqRange(2)-freqRange(1)).*rand(nImages,1) + freqRange(1);
scanContrast = (contrastRange(2)-contrastRange(1)).*rand(nImages,1) + contrastRange(1);

% Aberrations
measPupilMM = 4.5; % 4.5 mm pupil size (Thibos data)
calcPupilMM = 3; % 3 mm pupil size (what we want to calculate with)
zCoeffs = VirtualEyes(nImages,measPupilMM);

%% Generate images

samplesTemp = zeros(156,156,2*nImages);
samplesNoNoise = zeros(156,156,2*nImages);
labels = zeros(2*nImages,1);
contrasts = zeros(2*nImages,1);
freqs = zeros(2*nImages,1);

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
    labels(k) = 1;
    contrasts(k) = scanContrast(ii);
    freqs(k) = scanFreq(ii);

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
    labels(k) = 0;
    contrasts(k) = 0;
    freqs(k) = 0;
    k = k+1;
    
end

%% Crop
samplesCrop = samplesTemp(1:128,1:128,:);
samplesNoNoiseCrop = samplesNoNoise(1:128,1:128,:);
%% Save everything

if(saveFlag)
currDate = datestr(now,'mm-dd-yy_HH_MM');    
save(sprintf('%s_%s.mat',saveName,currDate),...
    'samplesTemp',...
    'samplesCrop',...
    'samplesNoNoise',...
    'samplesNoNoiseCrop',...
    'labels',...
    'contrasts',...
    'freqs');
end

%% Display a random sampling of images

n = 16;
figure(1);
[dataSamp,idx] = datasample(samplesTemp,n,3);
labelSamp = labels(idx);
contrastSamp = contrasts(idx);
freqSamp = freqs(idx);
k = 1;
for ii = 1:n
    subplot(4,4,k);
    imagesc(dataSamp(:,:,ii));
    axis image; axis off; 
    title(sprintf('label = %i \n c = %0.2f | f = %0.2f',...
        labelSamp(ii),contrastSamp(ii),freqSamp(ii)))
    k = k+1;
end
