% The function to calculate Poisson MLE.
% The input will be the coneMosaic absorption
% Inputs: 
%   timeInterval  - duration of each time step
%   absorption    - the coneMosaic absorption used to be classified
%   meanIsorate   - the mean isomerization rate (distribution) for each frame(time step)
%                   for the whole cones (w/ and w/o stimulus)

%% Load the absorption data (to be passed in the future)
path = [cs231nRootPath,'/data'];
cd(path)

load('sampleSet')

N = size(samplesTemp,3);

%% Initialization (temporary)

timeInterval = 0.005; %(5ms)
duration = 0.1; 
 
% what should isomerization rate be?
% finding a (constant) approximation for isomerization rate 
% with stimulus
idx = (labels == 1);
temp = double(reshape(samplesTemp(:,:,idx),[],1));
isoRate_c1 = (sum(temp)/length(temp))/timeInterval;

% without stimulus
idx = (labels == 0);
temp = double(reshape(samplesTemp(:,:,idx),[],1));
isoRate_c0 = (sum(temp)/length(temp))/timeInterval;


%% Comparing likelihoods to find if there was stimulus or not
% since the classes are balanced, class priors are equal -  we use MLE

likelihood_c0 = zeros(N,1);
likelihood_c1 = zeros(N,1);
difference = zeros(N,1);

for niter = 1:N
    % rounding to integer - it is not, but should be?
    curr_image = double(reshape(round(samplesTemp(:,:,niter)),[],1));
    
    % Finding log likelihood of cone response given no stimulus
    temp_c0 = poisspdf(curr_image,isoRate_c0*timeInterval);
    likelihood_c0(niter) = sum(log(temp_c0));
    
    % Finding log likelihood of cone response given no stimulus
    temp_c1 = poisspdf(curr_image,isoRate_c1*timeInterval);
    likelihood_c1(niter) = sum(log(temp_c1));
    
    % comparing likelihoods
    % if difference > 0 => class = stimulus
    % if difference < 0 => class = no stimulus
    difference(niter) = likelihood_c1(niter) - likelihood_c0(niter);
end

% finding prediction
predicted_labels = (difference>0);
accuracy = sum(predicted_labels == labels)/length(labels);





