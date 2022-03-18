%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%   This script reproduces the majority of the
%    results reported in the paper:
% ------------------------------------
% " Scale-invariant visual capabilities explained
%  by topographic representations of luminance
%   and texture in primate V1"
% ------------------------------------
%   Authors: G. Benvenuti, Y. Chen, WS. Geisler.
%   E. Seiemann
% ------------------------------------
%   Center for Perceptual Systems, University of
%   Texas, Austin, United States
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% ------------------------------------
%  by Giacomo Benvenuti <giacomox@gmail.com>
%  8 / 1 / 2018 , Austin, Texas, USA
% modified 10/5/2018s
% ------------------------------------

% * * INSTRUCTIONS * *
%
% Start just pressing "run" on the Matlab editor, sub-functions path s
% will be added automatically
%
% To reproduce the example d'maps in Fig1C,D,E , Fig3D,F and Fig5B
% First run the section below " Analysis1: Reproduce examples d'maps Fig1C,D,E; Fig3D; Fig5B"
% to create the resutls, then goes to the sections under that called
% "Display FigureXX" to display each panel in a separate figure. You can
% also just press "Run" to display all results.
%
% Concerning the imaging data we load below ("DataTrial") :
% The image frames captured on each trial were filtered in time (ref (Chen et al., 2012);
% FFT at 5Hz for VSD and 4Hz for Calcium) resulting in a single response image for each trial (505 * 504 pixels).
% For every session, outlier  trials were removed (see ref (Chen et al.,
% 2012)) .

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% example : DataFolder = 'Desktop/Benvenuti_et_al_2018/DataFolder/'  ;
clear all
fn =  mfilename('fullpath')
[filepath,name,ext] = fileparts(fn)
DataFolder = [filepath filesep 'DataFolder' filesep]

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$s

% Just checking if you made what I asked you to do in the instructions :)
if ~exist(DataFolder)
    errordlg('Please change the "DataFolder path" line38 ','Benvenuti_et_al_2018')
end
if ~exist('header')
    errordlg('Please add "Benvenuti_et_al_2018" to Matlab path ')
end

%% Analysis1: Examples d'maps Fig1C,D,E; Fig3D; Fig5B
% Pixel-by-pixel discriminability (d’ map) across responses to two orthogonal gratings
% or two opposite phases for Row, Columnar and retinotopic signals.
disp('Calculating d'' maps. Please wait...')
%-------------------------------------
%               GENERAL PARAMETERS
% ------------------------------------
% Param Fourier spatial filtration of neural response
EnvBandPass = [0.05,0.8];  % Retinotopic responses (cycle/mm)
ColBandPass = [0.8,3];       % Columnar responses
UFBand = [0 4];                  % unfiltered response

% Downsampling images
BinSize=8;

%Common param
Ind1 =  (14:63); % (14*2:62*2);%% Crop horizontal exp1
Ind2 = (1:48); % Crop vertical exp1
SFList = {'0.5' , '2' , '8'};


for Imaging_Type = 1:4 % 1 = VSD first Order , 2 = Calcium 1st order, 3 VSD 2nd Order, 4 Calcium Retino example Fig3F
    
    % ------------------------------------
    %            LOAD AND SORT DATA
    %-------------------------------------
    % i.e. sort the imaging trials using the Meta Data of the task storerd in
    % the TS var
    clear   DataTrial TS
    
    switch Imaging_Type
        case 1 % (VSD 1st Order)
            % Exp name M15D20141223R1
            %         load([DataFolder 'ExampleVSDExp1_Imaging.mat']) % laod the neural imaging file of the example VSD experiment
             load([DataFolder 'ExampleVSDExp1_MetaData.mat']) % load corresponding metadata "TS" variable
            load([DataFolder 'ExampleVSDExp1C_Imaging.mat']) % this data file contain complex numbers! (amplitude and phase)
          
            
            % Phase Locking temporal filtration
            DataTrial= PhaseLokedAmp(DataTrial, -.0026);
            SF_num = 3;
            CondType =1;  % Order stimulus
            
        case 2 % (Calcium 1st order stimulus) Fig1E
            % Exp name M22D20171019R0
            load([DataFolder 'ExampleCalciumExp1_Imaging.mat']) % laod the neural imaging file of the example VSD experiment
            load([DataFolder 'ExampleCalciumExp1_MetaData.mat']) % load corresponding metadata "TS" variable
            SF_num = 3;
            CondType =1;
            
        case 3 % (VSD 2nd order stimuli)
            % Exp name M23D20171013R2 (LO
            load([DataFolder 'ExampleVSDExp2_Imaging.mat']) % laod the neural imaging file of the example VSD experiment
            load([DataFolder 'ExampleVSDExp2_MetaData.mat']) % load corresponding metadata "TS" variable
            SF_num = 1;
            CondType =2;
            
        case 4 % (Calcium 1st order stimulus) Fig3F
            % Exp name M23D20171013R2
            load([DataFolder 'ExampleCalciumExp2_Imaging.mat']) % laod the neural imaging file of the example VSD experiment
            load([DataFolder 'ExampleCalciumExp2_MetaData.mat']) % load corresponding metadata "TS" variable
            SF_num = 1;
            CondType =1;
    end
    
    % Get Pxl size from MetaData TS var
    PxlSize =TS.Header.Imaging.SizePxl;
    
    % Sort data
    clear BC Controlti  ControlT O P S t Ou Pu Su
    % find control condition trials
    BC = find(TS.Header.Conditions.TypeCond ==0);
    Controlti=[TS.Header.Index.iValidBLKCond{BC}];
    ControlT =DataTrial(:,:,Controlti);
    
    for Order = 1:CondType % for second order stim I need Order 1 and 2
        
        % Sort evoked trials
        O = TS.Header.Conditions.GGOrtCond;
        P= TS.Header.Conditions.GGPhsCond;
        S= TS.Header.Conditions.GGSFCond;
        
        t= find(TS.Header.Conditions.TypeCond== Order);
        Ou = unique(TS.Header.Conditions.GGOrtCond(t));
        Pu= unique(TS.Header.Conditions.GGPhsCond(t));
        Su= unique(TS.Header.Conditions.GGSFCond(t));
        
        if Pu(end) == 270;
            Pu=[90 270];
        else
            Pu=[0 180];
        end
        
        Remove_cond = setdiff(1:size(O,2) , t);
        O(Remove_cond) = nan;
        P(Remove_cond) = nan;
        S(Remove_cond) = nan;
        
        % -------------------------------------
        % SPATIAL FILTRATION AND DOWN SAMPLING
        % -------------------------------------
        clear DataTrialUF DataTrialEnv DataTrialCol
        % Filter Columanr components, remove shot noise from
        % "Unfiltered" (low pass)
        
        for  s=1:SF_num % SF
            for o=1:2 % ORT
                for p=1:2 % PHASES
                    % Sort Trials
                    clear CondNum
                    CondNum = find(O==Ou(o) & P==Pu(p) & S==Su(s))  ; % checked
                    tm=   DataTrial(:,:,TS.Header.Index.iValidBLKCond{CondNum});
                    
                    % Columnar
                    DataTrialCol{s,o,p} = ...
                        FilterFermi2D(tm,ColBandPass(1),ColBandPass(2), ...
                        TS.Header.Imaging.SizePxl,0,BinSize);
                    
                    % Retino
                    DataTrialEnv{s,o,p}= ...
                        FilterFermi2D(tm,EnvBandPass(1),EnvBandPass(2), ...
                        TS.Header.Imaging.SizePxl,0,BinSize);
                    
                    % Unfiltered (just remove high freq noise)
                    DataTrialUF{s,o,p}  = ...
                        FilterFermi2D(tm,UFBand(1),UFBand(2), ...
                        TS.Header.Imaging.SizePxl,0,BinSize);
                    
                end
            end
        end
        
        
        % Same Filtration for Blank conditions
        % Retino
        clear DataTrialEnvBlank DataTrialColBlank  DataTrialUFBlank
        DataTrialEnvBlank  = ...
            FilterFermi2D(ControlT,EnvBandPass(1),EnvBandPass(2), ...
            TS.Header.Imaging.SizePxl,0,BinSize);
        
        % Columnar
        DataTrialColBlank = ...
            FilterFermi2D(ControlT,ColBandPass(1),ColBandPass(2), ...
            TS.Header.Imaging.SizePxl,0,BinSize);
        
        % Unfiltered
        DataTrialUFBlank =  FilterFermi2D(ControlT,UFBand(1),UFBand(2), ...
            TS.Header.Imaging.SizePxl,0,BinSize);
        
        
        
        % ------------------
        % CALCUL DPRIME MAPS
        % ------------------
        r =  Imaging_Type ;
        for  s=1:SF_num % SF
            % --------------------------------------
            % Columnar responses to orthogonal orientations
            % --------------------------------------
            %     if Imaging_Type <3
            clear O1UF O2UF BL1
            % Unfilt - Combine responses to different stimulus phases
            O1UF = cat(3,DataTrialUF{s,1,1},DataTrialUF{s,1,2}) ;
            O2UF = cat(3,DataTrialUF{s,2,1},DataTrialUF{s,2,2}) ;
            BL1 =   DataTrialUFBlank;
            
            % Discriminability Dprime maps Unfitlered (Figure1C)
            DP_Ort_BL_UF{s,1,r} = CalculateDPrime(BL1,O1UF); %
            DP_Ort_BL_UF{s,2,r} = CalculateDPrime(BL1,O2UF);
            DP_O1_O2_UF{s,r} = CalculateDPrime(O1UF,O2UF);
            DP_O1_O2_UF{s,r} = DP_O1_O2_UF{s} - mean(DP_O1_O2_UF{s}(:) ); % common average subtraction
            
            % Discriminability Dprime mapsColumnar bandpass/scale (Figure1D)
            O1C=  cat(3,DataTrialCol{s,1,1},DataTrialCol{s,1,2}) ; % combine opposite phases contiotions
            O2C = cat(3,DataTrialCol{s,2,1},DataTrialCol{s,2,2});
            DP_O1_O2_Col{s,r} = CalculateDPrime(O1C,O2C);
            DP_O1_O2_Col{s,r} = DP_O1_O2_Col{s,r} - mean(DP_O1_O2_Col{s,r}(:) );
            %   end
            
            % --------------------------------------
            % Columnar responses to opposite phases
            % --------------------------------------
            
            
            % Dprime across phases of the unfiltered signal (Figure3D)
            switch  Imaging_Type
                case 3 % Second order
                    % Retinotopic signal
                    if Order == 1
                        DP_Ort_Phs_Env{s,1,Order} = CalculateDPrime(DataTrialEnv{s,1,1},DataTrialEnv{s,1,2}) ;
                        DP_Ort_Phs_Env{s,2,Order} = CalculateDPrime(DataTrialEnv{s,2,1},DataTrialEnv{s,2,2}) ;
                    else
                        DP_Ort_Phs_Env{s,1,Order} = CalculateDPrime(DataTrialEnv{s,1,2},DataTrialEnv{s,1,1}) ;
                        DP_Ort_Phs_Env{s,2,Order} = CalculateDPrime(DataTrialEnv{s,2,2},DataTrialEnv{s,2,1}) ;
                    end
                    
                    for i = 1:2
                        DP_Ort_Phs_Env{s,i,Order} = DP_Ort_Phs_Env{s,i,Order} - mean( DP_Ort_Phs_Env{s,i,Order}(:));
                    end
                    
                    
                case 1
                    DataTrialUFX = DataTrialUF; % for example lum sign Fig3A
                    DataTrialUFBlankX = DataTrialUFBlank; % for example lum sign Fig3A
                    
                    DP_Ort_Phs_UF{s,1,r} = CalculateDPrime(DataTrialUF{s,1,1},DataTrialUF{s,1,2}) ;
                    DP_Ort_Phs_UF{s,2,r} = CalculateDPrime(DataTrialUF{s,2,1},DataTrialUF{s,2,2}) ;
                    for i = 1:2
                        DP_Ort_Phs_UF{s,i,r} = DP_Ort_Phs_UF{s,i,r} - mean( DP_Ort_Phs_UF{s,i,r}(:));
                    end
                case 4
                    DP_Ort_Phs_UF{s,1,r} = CalculateDPrime(DataTrialUF{s,1,1},DataTrialUF{s,1,2}) ;
                    DP_Ort_Phs_UF{s,2,r} = CalculateDPrime(DataTrialUF{s,2,1},DataTrialUF{s,2,2}) ;
                    for i = 1:2
                        DP_Ort_Phs_UF{s,i,r} = DP_Ort_Phs_UF{s,i,r} - mean( DP_Ort_Phs_UF{s,i,r}(:));
                    end
            end
            
        end
    end
    
    
end

disp('Please wait...DONE')
%-----------------------------------
%           DISPLAY DMAPS
% -------------------------------------


%% Display Figure1C VSD unfiltered responses to orthogonal gratings stimuli
% dmap unfiltered responses to orthogonal gratings stimuli

ClipDprimeMap = [-1.2 1.2];
figure
for sf = 1:3
    clear aa qq qq1
    subplot(1,3,sf)
    qq = DP_O1_O2_UF{sf,1};
    imshow(qq(Ind2,Ind1),ClipDprimeMap)
    axis on
    GiacStyle('Color','k','Linewidth',1) % color the Dprime map frame
    title(['SF ' SFList{sf} 'cpd'])
end
header('Figure1C: d''map VSD unfiltered responses to orthogonal gratings stimuli ', 12)
h = colorbar;
h.Position = [0.9179 0.3597 0.0250 0.3120];

%% Display Figure1D VSD columnar responses to orthogonal gratings stimuli
ClipDprimeMap = [-1.2 1.2];
figure
for sf = 1:3
    clear aa qq qq1
    subplot(1,3,sf)
    qq = DP_O1_O2_Col{sf,1};
    imshow(qq(Ind2,Ind1),ClipDprimeMap)
    axis on
    GiacStyle('Color','k','Linewidth',1) % color the Dprime map frame
    title(['SF ' SFList{sf} 'cpd'])
end
header('Figure1D: d''map VSD columnar responses to orthogonal gratings stimuli ', 12)
h = colorbar  ;
h.Position = [0.9179 0.3597 0.0250 0.3120] ;

%% Display Figure1E Calcium columnar responses to orthogonal gratings stimuli
clear hh qq qq1 mm
ct=0;
for o =1:2
    ct=ct+1;
    mm(:,:,ct) =  DP_Ort_BL_UF{2,o,2}(Ind2,Ind1);
end

qq = abs(mean(mm,3));
hold on
qq1 = imgaussfilt(qq,2);
hold on
figure;      imagesc(qq1)
hold on
hh =   HandleFit(qq1,[2 2],'r');
axis square off
colormap(jet); colorbar
title('Calcium expression example exp - average detect dp maps at 2cpd')
set(gcf,'color','w')

ClipDprimeMap = [-5 5];
figure
for sf = 1:3
    clear aa qq
    subplot(1,3,sf)
    qq = DP_O1_O2_Col{sf,2};
    imshow(qq(Ind2,Ind1),ClipDprimeMap)
    axis on
    GiacStyle('Color','k','Linewidth',1) % color the Dprime map frame
    
    hold on
    hh =   HandleFit(qq1,[2 2],'r');
    
    title(['SF ' SFList{sf} 'cpd'])
end
header('Figure1E: d''map Calcium columnar responses to orthogonal gratings stimuli ', 12)
h = colorbar  ;
h.Position = [0.9179 0.3597 0.0250 0.3120] ;


%% Display Fig3A,C  Lum Signal Example vs Retinotopic Model

clear DP_1PH_ENV
for s=1:2
    for  o =1:2
        for p=1:2
              DP_1PH_ENV(:,:,s,o,p) = CalculateDPrime(DataTrialUFX{s,o,p},DataTrialUFBlankX) ;
        end
    end
end

cas =      nanmean(  nanmean(nanmean(DP_1PH_ENV(:,:,:,:,:),5) ,4),3)  ;

figure
SF = [0.5 2 8];
O=[0 90 ];
P=[90 270];
clear EX
ct =0;
for s=1
    for  o =1:2
        for p=1:2
            ct=ct+1;
            subplot(2,4,ct)
            imagesc(DP_1PH_ENV(Ind2,Ind1,s,o,p)-cas(Ind2,Ind1),[0 5])
            EX(:,:,s,o,p)  = DP_1PH_ENV(Ind2,Ind1,s,o,p)-cas(Ind2,Ind1);
            axis off square
            title(['SF ' num2str(SF(s)) ' Ort ' num2str(O(o)) ' Phs ' num2str(P(p))])
        end
    end
end

%%
%--------------------------------------
%               FIGURE 3C
%-------------------------------------

load ([DataFolder 'Retinotopy1.mat'])
% Make model low SF evoked Resp
P=[90 270 90 270]
O=[0 90 ]
clear C MO
for i = 1:2
    for y=1:2
        clear R
        R = RetinoGratingBeta(RetinotopyCartesianXValid,RetinotopyCartesianYValid,63,63,Ind2,Ind1);
        R.phase=P(i+ (2*(y-1)));
        R.theta = O(y) ;
        R.freq = 0.5 ;
        R.StimCenterXY = [-1.5 -2.2];
    
        R.makegrating
        R.CPI_Blur(1.5,0.141);
        subplot(2,4,4+i+ (2*(y-1)))
        imagesc(R.Cortical_Grating_CPI)
        axis off square
    
      
         MO(:,:,1,y,i) = R.Cortical_Grating_CPI(1:48,1:50);    
        C(i,y) = corr2(  EX(:,:,1,y,i), R.Cortical_Grating_CPI(1:48,1:50))
        title(['Corr = ' num2str(C(i,y),2) ])
    end
end
%%
header('Figure 3A,C')
colormap(gray)
ax = get(gcf,'Position');
set(gcf,'Position',[ax(1), ax(2), ax(3), ax(4)] )
%% Display Figure3D d'map unfiltered VSD responses to opposite phases stimuli 3SF
ClipDprimeRetino = [-5 5 ; -1 1; -0.5 0.5]'  ;
figure
for i = 1:3
    subplot(2,3,i)
    qq = DP_Ort_Phs_UF{i,1,1};
    imshow(   qq(Ind2,Ind1),ClipDprimeRetino(:,i))
    F3D{1,i} =  qq(Ind2,Ind1) ;
    axis off
    title(['Ort 0 ; SF ' SFList{i} ])
    
    ax =  subplot(2,3,i+3) ;
    qq = DP_Ort_Phs_UF{i,2,1};
    imshow(   qq(Ind2,Ind1),ClipDprimeRetino(:,i))
    F3D{2,i} =  qq(Ind2,Ind1) ;
    axis off
    
    % Color bar
    ax1 = ax.Position  ;
    g = colorbar( 'SouthOutside' )
    g.Position=  [ax1(1), ax1(2)-.02, ax1(3), 0.03]
    g.Ticks=ClipDprimeRetino(:,i);
    title(['Ort 90 ; SF ' SFList{i} ])
end
header('Figure3D: d''map unfiltered responses to opposite phases stimuli ', 12)
set(gcf,'color','w')


% Make model low SF evoked Resp
load ([DataFolder 'Retinotopy1.mat'])
figure
SF = [.5 2 8 ]
O=[0 90 ]
clear C
R = RetinoGratingBeta(RetinotopyCartesianXValid,RetinotopyCartesianYValid,63,63,Ind2,Ind1);
for s = 1:3
    for y=1:2
        R.phase=90;
        R.theta = O(y) ;
        R.freq = SF(s) ;
        R.StimCenterXY = [-1.5 -2.2];
        R.makegrating
        R.CPI_Blur(1.5,0.141);
        subplot(2,3,s + (3*(y-1)))
        imagesc(R.Cortical_Grating_CPI)
        axis off square
        
        C(i,y) = corr2(  F3D{y,s}, R.Cortical_Grating_CPI(1:48,1:50))
        title(['SF ' num2str(SF(s)) ' Ort ' num2str(O(y)) ' Corr = ' num2str(C(i,y),2)])
    end
end
colormap(gray)
set(gcf,'color','w')


%% Display Figure3F,G  Luminance-Retinotopic Calcium signals vs Retinotopic model
figure

subplot(221)
qq = DP_Ort_Phs_UF{1,1,4};
imshow(   qq(Ind2,Ind1),[-1 1])
axis off
title(['Ort 0 ; SF ' SFList{1} ])

ax =  subplot(223) ;
qq = DP_Ort_Phs_UF{1,2,4};
imshow(   qq(Ind2,Ind1),[-1 1])
axis off

% Color bar
ax1 = ax.Position  ;
g = colorbar( 'SouthOutside' )
g.Position=  [ax1(1), ax1(2)-.02, ax1(3), 0.03]
g.Ticks=[-1 1];
title(['Ort 90 ; SF ' SFList{i} ])


% Model
load ([DataFolder 'Retinotopy2_Calcium.mat']) % M22D20170808R0Retinotopy_Trans2.mat
CMI_sigma = 1.5; % Cortical point image blurring param

%----------------------------------------
% Find the cortical area where calcium is expressed
%----------------------------------------
clear hh qq qq1 mm
ct=0;
for o =1:2
    ct=ct+1;
    mm(:,:,ct) =  DP_Ort_BL_UF{1,o,4}(Ind2,Ind1);
end

qq = abs(mean(mm,3));
hold on
qq1 = imgaussfilt(qq,2);
hh =   HandleFit(flipud(fliplr(qq1)),[7 7],'r');
hh(:,[1 32 69]) =[];
hh(:,:) = hh(:,[1:30 67:140 31:66]);


for o =1:2
    clear R D
    R = RetinoGratingBeta(RetinotopyCartesianXValid,RetinotopyCartesianYValid,63,63,Ind2,Ind1);
    if o ==1
        R.phase=0; % 200
        R.theta = 0 ;
        R.freq =2;
    else
        R.phase=270; % 180
        R.theta = 90 ;
        R.freq =3;
    end
    
    R.StimCenterXY = [-1.5 -2.2]; % X and Y position of the center of the grating stim in this exp
    R.makegrating
    R.CPI_Blur(1.5,0.141);
  
    D =   R.Cortical_Grating_CPI;
    MM = poly2mask(hh(1,:),hh(2,:),size(D,1),size(D,2));
   D =D.*MM;
    subplot(2,2,o*2)
    imagesc(D);colormap(gray)
    set(gca,'Xtick',[],'Ytick',[]); box on
    axis square
    C = corr2(DP_Ort_Phs_UF{1,o,4}(Ind2,Ind1).*MM,D);
    title(['Corr = ' num2str(C,2)])
end

header('Figure3F Calcium Imaging', 12)

%% Display Figure 5B dmap VSD retinotopic responses to opposite phases 2nd order stimuli
OriList = {'0' '90'}
figure
for i = 1:2
    for o = 1:2
        subplot(2,2,i+(2*(o-1)))
        imshow(   DP_Ort_Phs_Env{1,o,i},[-5 5])
        set(gca,'Xtick',[],'Ytick',[]); box on
        title(['Order ' num2str(3-i) ' Ori ' OriList{o} ])
    end
end
clear  C
C(1) = corr2(DP_Ort_Phs_Env{1,1,1},DP_Ort_Phs_Env{1,1,2})
C(2) = corr2(DP_Ort_Phs_Env{1,2,1},DP_Ort_Phs_Env{1,2,2})
header(' Figure 5B dmap VSD retinotopic responses to opposite phases 2nd order stimuli', 12)


%% Analysis 2 : Overall discriminability of columar responses (Example Analysis)
A = questdlg('Are you sure to run Analysis 2? It may takes time...')
clear OD
if strcmp(A,'Yes')
    disp('Wait please....')
    
    BinSize = 1; % Imaging trials full resolution
    
    load([DataFolder 'ExampleVSDExp1_Imaging.mat']) % laod the neural imaging file of the example VSD experiment
    load([DataFolder 'ExampleVSDExp1_MetaData.mat']) % load corresponding metadata "TS" variable
    SF_num = 3;
    CondType =1;  % Order stimulus
    
    
    % Get Pxl size from MetaData TS var
    PxlSize =TS.Header.Imaging.SizePxl;
    
    % Sort data
    clear BC Controlti  ControlT O P S t Ou Pu Su
    % find control condition trials
    BC = find(TS.Header.Conditions.TypeCond ==0);
    Controlti=[TS.Header.Index.iValidBLKCond{BC}];
    ControlT =DataTrial(:,:,Controlti);
    
    % Sort evoked trials
    O = TS.Header.Conditions.GGOrtCond;
    P= TS.Header.Conditions.GGPhsCond;
    S= TS.Header.Conditions.GGSFCond;
    
    t= find(TS.Header.Conditions.TypeCond== Order);
    Ou = unique(TS.Header.Conditions.GGOrtCond(t));
    Pu= unique(TS.Header.Conditions.GGPhsCond(t));
    Su= unique(TS.Header.Conditions.GGSFCond(t));
    
    if Pu(end) == 270;
        Pu=[90 270];
    else
        Pu=[0 180];
    end
    
    Remove_cond = setdiff(1:size(O,2) , t);
    O(Remove_cond) = nan;
    P(Remove_cond) = nan;
    S(Remove_cond) = nan;
    
    % -------------------------------------
    % SPATIAL FILTRATION AND DOWN SAMPLING
    % -------------------------------------
    clear DataTrialUF DataTrialEnv DataTrialCol
    % Filter Columanr components, remove shot noise from
    % "Unfiltered" (low pass)
    
    for  s=1:SF_num % SF
        for o=1:2 % ORT
            for p=1:2 % PHASES
                % Sort Trials
                clear CondNum
                CondNum = find(O==Ou(o) & P==Pu(p) & S==Su(s))  ; % checked
                tm=   DataTrial(:,:,TS.Header.Index.iValidBLKCond{CondNum});
                
                % Columnar
                DataTrialCol{s,o,p} = ...
                    FilterFermi2D(tm,ColBandPass(1),ColBandPass(2), ...
                    TS.Header.Imaging.SizePxl,0,BinSize);
                
                % Retino
                DataTrialEnv{s,o,p}= ...
                    FilterFermi2D(tm,EnvBandPass(1),EnvBandPass(2), ...
                    TS.Header.Imaging.SizePxl,0,BinSize);
                
                % Unfiltered (just remove high freq noise)
                DataTrialUF{s,o,p}  = ...
                    FilterFermi2D(tm,UFBand(1),UFBand(2), ...
                    TS.Header.Imaging.SizePxl,0,BinSize);
                
            end
        end
        
        
        
    end
    
    
    % Same Filtration for Blank conditions
    % Retino
    clear DataTrialEnvBlank DataTrialColBlank  DataTrialUFBlank
    DataTrialEnvBlank  = ...
        FilterFermi2D(ControlT,EnvBandPass(1),EnvBandPass(2), ...
        TS.Header.Imaging.SizePxl,0,BinSize);
    
    % Columnar
    DataTrialColBlank = ...
        FilterFermi2D(ControlT,ColBandPass(1),ColBandPass(2), ...
        TS.Header.Imaging.SizePxl,0,BinSize);
    
    % Unfiltered
    DataTrialUFBlank =  FilterFermi2D(ControlT,UFBand(1),UFBand(2), ...
        TS.Header.Imaging.SizePxl,0,BinSize);
    
    %-------------------------------------
    % ORI DISCRIMINATION BASED ON COLUMNAR RESPONSES
    %-------------------------------------
    for s = 1:3 % SFs
        clear O1C O2C dp ds
        O1C= cat(3,DataTrialCol{s,1,1},DataTrialCol{s,1,2}) ;
        O2C = cat(3,DataTrialCol{s,2,1},DataTrialCol{s,2,2}) ;
        %--------------
        clear sParam
        sParam.Method =1;
        sParam.DC_Sub=1;
        sParam.Thr=1;
        %--------------
        [dp ds ] = CalculateDiscrimination (O1C,O2C,sParam) ;
        OD(1,s) = dp;
        ODS(1,s) =ds;
    end
    %-------------------------------------
    %  ORI DISCRIMINATION BASED ON RETINOTOPIC DECODING
    %-------------------------------------
    
    for s = 1:3
        clear  DecisionVar
        % Weight each trial applying Leave one out cross validatation
        for i = 1:10 %
            % Estimate the decoder
            t = setdiff(1:10,i) ; % i is the test trial number, t is the training set of trails
            clear mm ssx dp ds
            for o =1:2
                for p=1:2
                    
                    nT = size(DataTrialEnv{s,o,p},3);
                    if i<=nT
                        
                        t = setdiff(1:nT,i) ;
                        % mean
                        mm{o,p} = mean(DataTrialEnv{s,o,p}(:,:,t),3) ; % here I removed the test trial
                        mm{o,2/p} = mean(DataTrialEnv{s,o,2/p},3) ;
                        mm{2/o,p} = mean(DataTrialEnv{s,2/o,p},3) ;
                        mm{2/o,2/p} = mean(DataTrialEnv{s,2/o,2/p},3) ;
                        
                        % std
                        ssx(:,:,o,p)  =std(DataTrialEnv{s,o,p}(:,:,t),[],3) ;
                        ssx(:,:,o,2/p)  =std(DataTrialEnv{s,o,2/p},[],3) ;
                        ssx(:,:,2/o,p)  =std(DataTrialEnv{s,2/o,p},[],3) ;
                        ssx(:,:,2/o,2/p)  =std(DataTrialEnv{s,2/o,2/p},[],3) ;
                        
                        S = sqrt(mean(mean(  ssx.^2  ,3 ) , 4 )) ;
                        
                        
                        if o ==1
                            p1x = p;
                            p2x = [1 2];
                        else
                            p1x=[1 2];
                            p2x =p;
                        end
                        
                        for p1 = p1x
                            for p2= p2x
                                clear M  Weight
                                M = (mm{1,p1} - mm{1, 2/p1}) /2 -  (mm{2,p2} - mm{2, 2/p2}) /2;
                                Weight = M./S;
                                Weight = Weight - mean(Weight(:)) ; % common average subtraction
                                [h w]= size(Weight);
                                Weight(find(Weight<1 & Weight>-1)) = 0;
                                Weight = reshape(Weight,[h w]);
                                %      DecisionVar(i,p1,p2,o,p) = mean(mean(DataTrialEnv{s,o,p}(:,:,i).*Weight,1),2);
                                if o ==1
                                    DecisionVar(i,p,p2,o) = mean(mean(DataTrialEnv{s,o,p}(:,:,i).*Weight,1),2);
                                else
                                    DecisionVar(i,p1,p,o) = mean(mean(DataTrialEnv{s,o,p}(:,:,i).*Weight,1),2);
                                end
                                
                            end
                        end
                        
                    end
                end
            end
        end
        ct = 0;
        %     for p=1:2
        for p1=1:2
            for p2=1:2
                clear DP DPESSE vA vB
                vA =  DecisionVar(:,p1,p2,1);
                vB =  DecisionVar(:,p1,p2,2);
                
                DP = CalculateDPrime(vA,vB,1);
                
                % Find  Significance thashold by Resempling
                %1 Standard error by std of Bootstrap
                DPBSSE = ...
                    std(CalculateDPrime(vA(randi(10,[1000,10])), ...
                    vB(randi(10,[1000,10])),2),[],1);
                
                
                ct = ct +1;
                dp(ct) = DP;
                ds(ct) =  DPBSSE;
                
            end
        end
        %     end
        
        
        OD(2,s) = mean(dp);
        ODS(2,s) = mean(ds);
        
        
        %-------------------------------------
        %  STIMULI DETECTION BASED UNFILTERED RESPONSES
        %-------------------------------------
        BL1 =   DataTrialUFBlank;
        for s = 1:3
            clear O1UF O2UF
            O1UF = cat(3,DataTrialUF{s,1,1},DataTrialUF{s,1,2}) ;
            O2UF = cat(3,DataTrialUF{s,2,1},DataTrialUF{s,2,2}) ;
            clear dp1 dp2 ds1 ds2
            %---------------
            clear sParam
            sParam.Method =1;
            sParam.DC_Sub = 0;
            %---------------
            [dp1 ds1 ww dvA1 dvB1 ] = CalculateDiscrimination (O1UF,BL1,sParam);
            [dp2 ds2 ww dvA2 dvB2] = CalculateDiscrimination (O2UF,BL1,sParam);
            
            OD(3,s) = mean([dp1,dp2]);
            ODS(3,s) = mean([ds1,ds2]);
        end
        
        
    end
    disp('Wait please....DONE')
end
%% Display Figure 1F (Columnar decoding)
% run "Analysis 2 " section above
figure
h = bar(OD(1,:))
h.FaceColor = 'none'
h.EdgeColor = 'b'
h.LineWidth = 3
hold on
errorbar(1:3,OD(1,:),ODS(1,:),'.')
box off
set(gcf,'color','w')
set(gca,'Xticklabel',{'0.5', '2', '8' })
xlabel('Spatial Frequency (cpd)')
ylabel('d''pop')
title('Figure 1F Overall discriminability columnar responses')



%% Dispaly Figure 2 A (Detection)
% run "Analysis 2 " section above
figure
h = bar(OD(3,:))
h.FaceColor = 'none'
h.EdgeColor = 'k'
h.LineWidth = 3
hold on
errorbar(1:3,OD(3,:),ODS(3,:),'k.')
box off
set(gcf,'color','w')
set(gca,'Xticklabel',{'0.5', '2', '8' })
xlabel('Spatial Frequency (cpd)')
ylabel('d''pop')

title('Figure 2A Detectability Unfiltered signals')

%% Dispaly Figure 4 A (Retinotopic decoding)
% run "Analysis 2 " section above
figure
h = bar(OD(2,:))
h.FaceColor = 'none'
h.EdgeColor = 'c'
h.LineWidth = 3

hold on
errorbar(1:3,OD(2,:),ODS(2,:),'.')
box off
set(gcf,'color','w')
set(gca,'Xticklabel',{'0.5', '2', '8' })
xlabel('Spatial Frequency (cpd)')
ylabel('d''pop')
title('Figure 4A')
title('Figure 4A Overall discriminability retinotopic responses')

%% Fig 1G, 2B and 4B Average discriminability and detectability
clear OD ODS
 load([DataFolder 'dpop.mat']) % FlashGrat_New7
sel = find(OD(:,5,4) >2 );
col = ['kbc'];
figure

%------------------------------------------------
%   AVERAGE DISCRIMINABILITY BASED ON COLUMNAR RESP.
%------------------------------------------------
subplot(131)
t =squeeze( nanmean(OD(sel,5,[2:2:6]),1));
tse = squeeze( nanstd(OD(sel,5,[2:2:6])))  / sqrt(numel(sel)-1);

errorbar([1:3],t,tse,'.k');
hold on
h2(1) = bar( [1:3],t,.3,'Edgecolor','k','Facecolor','b','linewidth', 2);
ylim([0 12])
xlim([0.5 3.5])
box off
axis square
ylabel('d''_{pop}')
xlabel('Spatial Frequency (cpd)')
title('Retino-Luminance')


%------------------------------------------------
%   AVERAGE DISCRIMINABILITY BASED ON LUMINANCES RESP.
%------------------------------------------------
subplot(132)
t =squeeze( nanmean(OD(sel,7,[2:2:6]),1));
tse = squeeze( nanstd(OD(sel,7,[2:2:6])))  / sqrt(numel(sel)-1);

errorbar([1:3],t,tse,'.k');
hold on
h2(1) = bar( [1:3],t,.3,'Edgecolor','k','Facecolor','c','linewidth', 2);
ylim([0 10])
xlim([0.5 3.5])
box off
axis square
ylabel('d''_{pop}')
xlabel('Spatial Frequency (cpd)')

header('Average discriminability')
title('Retino-Luminance')


%------------------------------------------------
%            AVERAGE DETECTABILITY
%------------------------------------------------
subplot(133)
t =squeeze( nanmean(OD(sel,2,[2:2:6]),1));
tse = squeeze( nanstd(OD(sel,2,[2:2:6])))  / sqrt(numel(sel)-1);

errorbar([1:3],t,tse,'.k');
hold on
h2(1) = bar( [1:3],t,.3,'Edgecolor','k','Facecolor','k','linewidth', 2);
ylim([0 16])
xlim([0.5 3.5])
box off
axis square
ylabel('d''_{pop}')
xlabel('Spatial Frequency (cpd)')
title('Detectability')
header('Average discriminability')
set(gcf,'color','w')

