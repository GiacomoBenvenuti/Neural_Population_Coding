function [Phase_Locked_Amplitude, phase0_Suggested] = PhaseLokedAmp (varargin)
% Phase lock the FFT amplitude to a specific "phase0" (Flash Grating time FFT created by RunDA.m)
% Phase_Locked_Amplitude = PhaseLokedAmp --> open dialog
% Phase_Locked_Amplitude = PhaseLokedAmp (FFTmatrix, phase0)
% Phase_Locked_Amplitude = PhaseLokedAmp (FFTmatrix, []) --> use suggested
% Phase_Locked_Amplitude = PhaseLokedAmp (FFTmatrix, phase0, showImage(0/1))
% phase : the function use the phase of the stronger pixel amp
% FFTmatrix =  var 'DataTrial' from FFT file created by RunDA.m (Complex number P=Amp-iPhase)
% phase0 = phase to lock to. Should correspond to the peak of the evoked response
% Phase_Locked_Amplitude = same than 'DataTrial' but phase loked for phase0
%
% Calculate phase knwing the response peak frame
% (PeakFrame-StartFrameFFT)*FrameDuration/CycleDuration*2*pi *-1;
% - FrameDuration in sec!! but CycleDuration in ms
% 
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% EXAMPLE
%      StimulusFlshingRate = 5;
%     % Find monkey name
%     id1 = strfind(t, '/M') + 2;
%     id = str2num(t(id1:id1+1));
%     if id == 19
%         PeakFrame =  15;
%     else
%         PeakFrame =  22;
%     end
%     % Find first frame number
%     it1 = strfind(t, 'FFTS') + 4;
%     StartFrameFFT= str2num(t(it1:it1+2));
%     
%     FrameDuration =  1000/TS.Header.Imaging.FrameRate/1000 ;
%     CycleDuration =  1000/StimulusFlshingRate;
%     phase0 =  (PeakFrame-StartFrameFFT)*FrameDuration/CycleDuration*2*pi *-1;
% -------------------------------------------
% by Giacomo Benvenuti, Seidemann LAB, Autin 9/2017
% -------------------------------------------

% Get DataTrial matrix
if nargin >0
    DT = varargin{1};
else
    [Fnx  pp]=  uigetfile('','Select FFT file created by RunDA.m');
    Fn = [pp Fnx];
    load(Fn)
    DT = DataTrial;
end

% Input must be a complex number : the imaginary part is the phase
if  isreal(DT)
    error('In this matrix/File there is not phase information (Imaginary part)')
end

 % Method 1;
        tm = mean(DT,3);
        [a b ] = max(real(tm(:)));
        phase0_Suggested = imag(tm(b)) ;

% Method 2
% 10*.01/(1200/5)*2*pi *-1
% (PeakFrame-StartFrameFFT)*FrameDuration/CycleDuration*2*pi * -1 ;



% Get PHASE0
if nargin >1
    if isempty(varargin{2})
       
        phase0 =  phase0_Suggested;
    else
        phase0 = varargin{2};
    end
else
    xx = inputdlg('Phase0 (Frame response peak)', 'Phase0',1,  {num2str(phase0_Suggested)} );
    phase0 = str2num(xx{1})
end

nT= size(DT,3);
% Calculate the amplitude at phase0
for y = 1:nT % don't use index 'i'!!!
    Phase_Locked_Amplitude(:,:,y) = real(   DT(:,:,y)  *  exp( -phase0 * i )  );
end
%%
% Check before/after phase locking
if nargin>2
    showw = varargin{3} ;
else
    showw = 0;
end
if showw
figure
clear x y D
x= mean(real(Phase_Locked_Amplitude),3);
y = mean(real(DT),3);

subplot(131)
% before
imagesc(y);
c = caxis;
axis square off
sc = sum(abs(y(:)));
title(['Before: ' num2str(sc,3)])


subplot(132)
% after
imagesc(x,c);
axis square off
sc = sum(abs(x(:)));
title(['After: ' num2str(sc,3)])

subplot(133)
% differerence

D= x-y;
imagesc(D);
axis square off
sc = sum(abs(D(:)));
title(['Diff: ' num2str(sc,3)])
colormap(gray)
end
%%