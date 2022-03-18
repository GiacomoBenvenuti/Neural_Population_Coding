function B = FilterFermi2D(varargin)
% Filter A with 2D Fermi filter
% B = FilterFermi2D(A,LowCutOff,HighCutOff,SizePxl,Show,DownSampling_BinSize )
% input: 
% A is an array with at leat 2 dimentions
%  LowCutOff and HighCutOff are the low an high cut-off frequency (cyc/mm)
% SizePxl is the pixel size (mm)
% Show = 0 allows to not display the figures
% DownSampling_BinSize : Bin size for down sampling (interpolation). If
% empty NOT interpolation is carried out.
%
% Output: 
% B is an array with same size as A
%
%
% YC & GB at ES lab
% Created on Apr. 24, 2013
% Last modified on Apr. 24, 2013

% Default values
showw = 0;
BinSize=[];

in = {'A','LowCutOff','HighCutOff','SizePxl','showw','BinSize'};
for i = 1:nargin
    eval([in{i} '= varargin{i};']) ;
end
if nargin <3 
    error('not enough input arguments')
end

%% Check input/output arguements
SizeA = size(A);
Height = SizeA(1);
Width = SizeA(2);

%% Fermi filter
ParmFermiLowPass = [1,HighCutOff,0,HighCutOff*0.05];
ParmFermiHighPass = [1,LowCutOff,0,LowCutOff*0.05];

SFX = ((1:Width)-floor(Width/2)-1)/(Width-1)/SizePxl;
SFY = ((1:Height)-floor(Height/2)-1)/(Height-1)/SizePxl;
[SFXX,SFYY] = meshgrid(SFX,SFY);
SF2D = abs(SFXX+1i*SFYY);
if HighCutOff==inf
  FiltFermiLowPass = zeros(Height,Width);
else
  FiltFermiLowPass = FuncWoNFermi(ParmFermiLowPass,SF2D);
end
if LowCutOff==0
  FiltFermiHighPass = ones(Height,Width);
else
  FiltFermiHighPass = FuncWoNFermi(ParmFermiHighPass,SF2D);
end 
FiltFermi = FiltFermiHighPass-FiltFermiLowPass;

%% Show Fermi filter
if showw
thgf = figure;
ti = floor(Height/2)+1;
plot(SFX,FiltFermi(ti,:),'LineWidth',2);
axis tight;
axis([0,10,-0.1,1.1]);
line([1,1]*HighCutOff,ylim, ...
     'Color','r', ...
     'LineStyle','--', ...
     'LineWidth',2);
line([1,1]*LowCutOff,ylim, ...
     'Color','r', ...
     'LineStyle','--', ...
     'LineWidth',2);
xlabel('Spatial frequency (cycle/mm)', ...
       'FontWeight','bold', ...
       'FontSize',12);
ylabel('Amplitude', ...
       'FontWeight','bold', ...
       'FontSize',12);
title('Fermi filter (1D-slice)', ...
      'FontWeight','bold', ...
      'FontSize',12);
drawnow;pause(0.1);
end
%% Filter A
SizeA([1,2]) = 1;
B = ...
  ifft(ifft(ifftshift(ifftshift( ...
    fftshift(fftshift(fft(fft(A,[],1),[],2),1),2).* ...
    repmat(FiltFermi,SizeA),1),2),[],1),[],2);

if isreal(A)
  B = real(B);
end

%---------------
% DOWNSAMPLING
%---------------
if ~isempty(BinSize)
% calculate interp coordinates
[inH,inW,F] = size(B);
outH = round(inH / BinSize);
outW = round(inW / BinSize);
hh = linspace(0.5,inH+0.5,2*outH+1);
ww = linspace(0.5,inW+0.5,2*outW+1);
hh = hh(2:2:end);
ww = ww(2:2:end);

clear tmpx
for nn=1:F
    tmpx(:,:,nn) =  interp2(1:inW,(1:inH)',B(:,:,nn) ,ww,hh');
end
clear B
B=tmpx;
end      
        
%% Close figure
if showw
close(thgf);

end
