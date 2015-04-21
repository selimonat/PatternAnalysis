%% will generate a 16x16 covariance matrix from simulated patterns of 8
%different faces recorded in two different phases.
w_B                 = [cos(linspace(0,2*pi,8)); ...
                       sin(linspace(0,2*pi,8))];%an example weight matrix

x = linspace(0,2*pi-2*pi/8,8);x=x-deg2rad(135);
x = make_gaussian1d(x,2,3,0);
w_T                 = [x; ...
                       0 0 0 0 0 0 0 0];%an example weight matrix

tphase                  = 2;
tcond                   = 16;
smoothness_factor       = 10;
win                     = ones(smoothness_factor,1)./smoothness_factor;
tvoxel                  = 500;
common_pattern_noiseamp = 0.1;
evoked_noiseamp         = .1;
%create a common pattern that is common to all the conditions
common_pattern          = conv(rand(tvoxel,1),win,'same');
common_pattern          = common_pattern*ones(1,tcond) +randn(tvoxel,tcond)*common_pattern_noiseamp;
% plot(common_pattern);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%different possibilities for evoked activity in the Baselin
%0. no activity just random noise added to the common activity
%1. common evoked activity for all faces (1 dimenional)
%2. a gradient of activity (1 dimensional)
%3. one specific face with the rest being similar (2 dimensional)
%4. how to generate circular covariance structure ?
%
%first create two orthogonal patterns, than fuse them according to the
%conditions above
evoked_patterns   = conv2(rand(tvoxel,2),win,'same');
%orthogonolize
evoked_patterns   = OrthogonolizeTwoVectors(evoked_patterns);
evoked_noise      = randn(tvoxel,tcond/2)*evoked_noiseamp;
evoked_activity_B = evoked_patterns*w_B + evoked_noise;
                   

%different possibilities for the Test phase as well.
%0-3. same as above
evoked_patterns   = conv2(rand(tvoxel,2),win,'same');
%orthogonolize
evoked_patterns   = OrthogonolizeTwoVectors(evoked_patterns);
evoked_noise      = randn(tvoxel,tcond/2)*evoked_noiseamp;
evoked_activity_T = evoked_patterns*w_T + evoked_noise;



%however, here we have different options to model the conditioning effect
%0. no effect of conditioning
%1. increased of variance centered on the CS+ face, more or less localized
%(two parameter: localization and strength)
%


%these are effects that are applied on the covariance matrix
%1. circular effect centered on the CS+ face, can be negative or positive
%2. localized Gaussian effect on the 4th column and row...
%3. divergent or convergent effects on the circular covariance structure
%(the number 4. above)

%%
cmat = cov([evoked_activity_B evoked_activity_T]);
imagesc(cmat);

