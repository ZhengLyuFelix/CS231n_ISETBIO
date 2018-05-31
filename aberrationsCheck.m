%% Compare cone output between aberrated optics and perfect optics

%% Initialize
clear; close all;
ieInit;

rng(1);

%% Set up the stimulus parameters
clear hparams
hparams(2)           = harmonicP;
hparams(2).freq      = 10;  % Set the Frequency
hparams(2).contrast  = 1;
hparams(1)           = hparams(2);
hparams(1).contrast  = 0;

sparams.fov = 1;

nTimeSteps = 20;
stimWeights = ones(1, nTimeSteps);

%% Set up cone mosaic parameters
integrationTime = 0.005;
sampleTimes = ((1:nTimeSteps) - 1) * integrationTime;   % Five ms integration time
nTrials    = 100;

%% Generate Zernike coefficients
measPupilMM = 4.5; % 4.5 mm pupil size (Thibos data)
calcPupilMM = 3; % 3 mm pupil size (what we want to calculate with)
zCoeffs = VirtualEyes(10,measPupilMM);

%% Create the oi with aberrations
    
%z = zeros(65,1);
z(1:13) = zCoeffs(1,1:13);
% z(4) = 5*z(4);

% Create the example subject
sbjWvf = wvfCreate;                                     % Initialize
sbjWvf = wvfSet(sbjWvf,'zcoeffs',z);                    % Zernike
sbjWvf = wvfSet(sbjWvf,'measured pupil',measPupilMM);   % Data
sbjWvf = wvfSet(sbjWvf,'calculated pupil',calcPupilMM); % What we calculate
sbjWvf = wvfSet(sbjWvf,'measured wavelength',550);
sbjWvf = wvfSet(sbjWvf,'calc wave',[450:10:650]');            % Must be a column vector
sbjWvf = wvfComputePSF(sbjWvf);
oiAberr = wvf2oi(sbjWvf);

%% Create default oi

oiDefault = oiCreate('wvf human');

%% Aberrated calculation

ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams,...
    'oi',oiAberr);

cm = coneMosaic;
cm.integrationTime = ois.timeStep;

% Make the cm smaller than the oi size, but never smaller than 0.2 deg
fovDegs = max(oiGet(ois.oiFixed,'fov') - 0.2, 0.2);  % Degrees
cm.setSizeToFOV(fovDegs);

cm.noiseFlag = 'none';
aberrHarmonic = mean(squeeze(cm.compute(ois)), 3);
figure;
imagesc(aberrHarmonic); colormap(gray); caxis([0 35])
title('Aberrated');
oiTmpAberr = ois.oiModulated;
oiTmpAberr = oiSet(oiTmpAberr,'name','Aberrated');
ieAddObject(oiTmpAberr);
oiWindow;

%% Default calculation

ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams,...
    'oi',oiDefault);

cm = coneMosaic;
cm.integrationTime = ois.timeStep;

% Make the cm smaller than the oi size, but never smaller than 0.2 deg
fovDegs = max(oiGet(ois.oiFixed,'fov') - 0.2, 0.2);  % Degrees
cm.setSizeToFOV(fovDegs);

cm.noiseFlag = 'none';
defaultHarmonic = mean(squeeze(cm.compute(ois)), 3);
figure;
imagesc(defaultHarmonic); colormap(gray); caxis([0 35]);

title('Default');
oiTmpDefault = ois.oiModulated;
oiTmpDefault = oiSet(oiTmpDefault,'name','Default');
ieAddObject(oiTmpDefault);
oiWindow;


%% Show difference
figure;
subplot(1,3,1);
imagesc(defaultHarmonic); colormap(gray); caxis([0 35]); colorbar;
axis image; axis off;
title('default');
subplot(1,3,2);
imagesc(aberrHarmonic); colormap(gray); caxis([0 35]); colorbar;
axis image;  axis off;
title('aberrated');
subplot(1,3,3);
imagesc(abs(aberrHarmonic - defaultHarmonic)); colorbar;
axis image;  axis off;
title('Difference');
    