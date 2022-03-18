classdef RetinoGratingBeta < handle
    % Generate the retinotopic projection of a grating
    % Cortical_Grating = RetinoGrating(RetinotopyCartesianXValid, RetinotopyCartesianYValid)
    % Cortical_Grating = RetinoGrating(RetinotopyCartesianXValid, RetinotopyCartesianYValid,CortGrat_nPxl_X,CortGrat_nPxl_Y)
    % obj.DispCortGrating('CortIMG',[])
    % by Giacomo Benvenuti <giacomox@gmail.com>
    % University of Texas at Austin , USA.

    
    
    properties
        Version = '1.1.3  Updated 09/20/2018'
        Xm = []; % RetinotopyCartesianXValid
        Ym=[];  % RetinotopyCartesianYValid
        freq = 0.5; % cyc per deg
        theta = 0; % deg in visual space
        phase = 0; % deg, in visual space
        Cortical_Grating = [];
        Grating = [];
        GreenIMG = [];
        CortGrat_nPxlX = 100;
        CortGrat_nPxlY =100;
        SigmaCPI = []; % Cortical point image 2D gaussian sigmas;
        PxlSize = [];
        Cortical_Grating_CPI = [];
        cropX = [];
        cropY = [];
        StimCenterXY = [3,-3]; % in dva
        StimSizeXY = [6,6]; % in dva
        
        
    end
    
    methods
        function obj = RetinoGratingBeta(varargin)
            % Constructor:
            % RetinoGrating(RetinotopyCartesianXValid, RetinotopyCartesianYValid, ...
            %        PxlSize)
            
            if nargin >1
                obj.Xm = varargin{1};
                obj.Ym = varargin{2};
            else
                disp('Plase select a Retinotopy mat file')
                [fn pp] = uigetfile('.mat', 'Select file Retinotopy.mat');
                load([pp fn],'RetinotopyCartesianXValid','RetinotopyCartesianYValid');
                obj.Xm = RetinotopyCartesianXValid;
                obj.Ym = RetinotopyCartesianYValid;
            end
            
            if nargin>2
                obj.CortGrat_nPxlX = varargin{3};
                obj.CortGrat_nPxlY = varargin{4};
            end
            
            if nargin>4
                obj.cropX = varargin{5};
            end
            if nargin>5
                obj.cropY = varargin{6};
            end
            
            
            obj.makegrating;
            
            
        end
        
        function obj = DispCortGrating(obj,varargin)
            if nargin>0
                param = {'CortIMG','freq'};
                for i = 1: size(param,2)
                    tm = strcmp(varargin,param{i});
                    if sum(tm)>0
                        a = find(tm>0) + 1;
                        if isempty(varargin{a}) & i==1
                            [fn pp] = uigetfile('.bmp', 'Select GREEN IMAGE   ');
                            obj.GreenIMG = imread([pp fn]);
                            gm = obj.GreenIMG(:,:,[1 1 1 ]);
                            
                        else
                            
                            
                            
                            
                        end
                        
                    end
                end
            end
            
            if ~isempty(obj.GreenIMG)
                gm1 = gm-min(gm(:));
                gm2 = gm1/max(gm1(:));
                h2=  imshow(gm);
                hold on
                
            end
            
            obj.makegrating;
            
            % Giac2018
            %             cg = obj.Cortical_Grating;
            %             cg1= (cg-min(cg(:)));
            %             cg2 = cg1 / max(cg1(:));
            %             h1= imagesc(cg2)
            
            h1 =  imagesc(obj.Cortical_Grating);
            if ~isempty(obj.GreenIMG)
                set(h1, 'AlphaData',.4);
            end
            cm = jet;
            cm(1,:) = [0 0 0];
            %colormap(cm)
            colormap(parula);
            caxis([-1 1])
        end
        
        
        
        function obj = CPI_Blur(obj, SigmaCPI, PxlSize)
            % obj.CPI_Blur(obj, obj.SigmaCPI, obj.sPxlSize)
            obj.SigmaCPI = SigmaCPI;
            obj.PxlSize = PxlSize;
            g = obj.Cortical_Grating;
            Bcrop_padded = padarray(g,size(g),'both','circular'); % padding
            SiGMa = SigmaCPI /PxlSize; % scale sigmas
            A= gauss2d(ones(2*ceil(2*SiGMa)+1,2*ceil(2*SiGMa)+1), ...
                SiGMa,[ceil(2*SiGMa),ceil(2*SiGMa)]); % 2D gaussian Filter
            g2=  conv2(Bcrop_padded, A,  'same');
            [a b] = size(g2);
            obj.Cortical_Grating_CPI  = g2(a/3+1:(a/3*2),b/3+1:(b/3*2)); 
            
        end
        
        function grating = makegrating(obj)
            %----------------------------------------
            % Same algorithm to generate stimuli we presented to the monkey
            % X=0, Y=0 is the center
            % CreateStaticSinusoidalGratingCenter.m
            %----------------------------------------
            
            Xm = obj.Xm;
            Ym = obj.Ym;
            freq = obj.freq;
            phase = obj.phase;
            theta = obj.theta;
            
%             C2x = ceil(obj.StimSizeXY(1)/2);
%             C2y = ceil(obj.StimSizeXY(2)/2);
            Cx = obj.StimCenterXY(1);
            Cy = obj.StimCenterXY(2);
         
         grating = ...
             cos(2*pi*freq*sin(theta/180*pi)*(Xm-Cx)+ ...
                2*pi*freq*cos(theta/180*pi)*(Ym.*-1-Cy*-1)+ ...
                phase/180*pi) ;

                      
            grating2 = imresize(grating, 'OutputSize'  , [obj.CortGrat_nPxlX   obj.CortGrat_nPxlX]);
            
            if isempty(obj.cropX)
                obj.cropX = 1:size( grating2,1);
            end
            if isempty(obj.cropY)
                obj.cropY = 1:size( grating2,2);
            end
            
            grating2 =  grating2(obj.cropX,obj.cropY);
            
            obj.Cortical_Grating = grating2;
        end
    end % Methods
    
end

function mat = gauss2d(mat, sigma, center)
gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(R,C, sigma, center);
end

function val = gaussC(x, y, sigma, center)
xc = center(1);
yc = center(2);
exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
val       = (exp(-exponent));
end

