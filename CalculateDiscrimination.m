function [DP,DPBSSE,Weight,vA,vB] = CalculateDiscrimination(A,B,sParm)
% function [DP,DPBSSE,Weight,vA,vB] = CalculateDiscrimination(A,B,sParm)
% Calculate discrimination of A and B
% Inputs:
%   A and B must be 3D arrays with [Height,Width,nTrial]
%   Height and Width must be the same in A and B, but nTrial can be different
%   sParm.nBS define the number of bootstrap
%   sParm.Method define the model
%     1. weighted d' (Jackknife)
%     2. input sParm.Weight
%     3. input sParm.DC_Sub
%     4. input sParam.Thr : set to 0 all the Template pixels with Dprime lower than
%         this value
%
% Outputs:
%   DP is the d' between vA and vB
%   DPBSSE is the standard error of d' with bootstrape
%   vA and vB are the decision vector with Jack-knife
%
%   8/2017 Giac added DC subtraction
%
% YC & GB at ES lab
% Created on Mar. 12, 2014
% Last modified on Sept. 2017

%% Check inputs and/or outputs
[nYA,nXA,nTA] = size(A);
[nYB,nXB,nTB] = size(B);

if nYA~=nYB || nXA~=nXB || min(nTA,nTB)<=1
    error('A and B must be 3D arrays with [Height,Width,nTrial]!');
end

if ~exist('sParm','var')
    Method = 1;  % default
    nBS = 1000;  % default
else
    if ~isfield(sParm,'Method')
        Method = 1;  % default
    else
        Method = sParm.Method;
    end
    if ~isfield(sParm,'nBS')
        nBS = 1000;  % default
    else
        nBS = sParm.nBS;
    end
    
    % by GB
    if ~isfield(sParm,'DC_Sub')
        DC_Sub = 0;
    else
        DC_Sub = sParm.DC_Sub;
    end
    
    % by GB
    if isfield(sParm,'Thr')
        Thr = sParm.Thr;
    else
        Thr = 0;
    end
    
end

%% Discrimination
switch Method
    case 1  % weighted d'
        vA = zeros(1,nTA);
        for i = 1:nTA
            Weight = CalculateDPrime(A(:,:,setdiff(1:nTA,i)),B,3);

            if Thr>0
                [h w] = size(Weight);
                Weight(find(Weight<Thr & Weight>-Thr)) = 0;
                Weight = reshape(Weight,[h w]);
            end
     
            if DC_Sub==1
                Weight = Weight - nanmean(Weight(:)) ; % Giac 2017 xxxxxxxxxxxxxxxxxxxxxx
            end
          
            vA(i) = nanmean(nanmean(A(:,:,i).*Weight,1),2);
        end
        
        vB = zeros(1,nTB);
        for i = 1:nTB
            Weight = CalculateDPrime(A,B(:,:,setdiff(1:nTB,i)),3);
       
            if Thr>0
                [h w] = size(Weight);
                Weight(find(Weight<Thr & Weight>-Thr)) = 0;
                Weight = reshape(Weight,[h w]);
            end
            
          if DC_Sub==1
                Weight = Weight - nanmean(Weight(:)) ; % Giac 2017 xxxxxxxxxxxxxxxxxxxxxx
          end
            
            vB(i) = nanmean(nanmean(B(:,:,i).*Weight,1),2);
        end
        
        Weight = CalculateDPrime(A,B,3);
          if Thr>0
                [h w] = size(Weight);
                Weight(find(Weight<Thr & Weight>-Thr)) = 0;
                Weight = reshape(Weight,[h w]);
          end
       
          if DC_Sub==1
                Weight = Weight - nanmean(Weight(:)) ; % Giac 2017 xxxxxxxxxxxxxxxxxxxxxx
          end
            
    case 2  % input sParm.Weight
        if ~isfield(sParm,'Weight')
            error('Input weight for method 2!');
        end
        Weight = sParm.Weight;
     
        
        if Thr>0
            [h w] = size(Weight);
            Weight(find(Weight<Thr & Weight>-Thr)) = 0;
            Weight = reshape(Weight,[h w]);
        end
        
         if DC_Sub==1
            Weight = Weight - nanmean(Weight(:)) ; % Giac 2017 xxxxxxxxxxxxxxxxxxxx
        end
  
        vA = squeeze(nanmean(nanmean(A.*repmat(Weight,[1,1,nTA]),1),2))';
        vB = squeeze(nanmean(nanmean(B.*repmat(Weight,[1,1,nTB]),1),2))';
    otherwise
        error('Wrong method!');
end

DP = CalculateDPrime(vA,vB,2);

% Find  Significance thashold by Resempling
%1 Standard error by std of Bootstrap
DPBSSE = ...
    std(CalculateDPrime(vA(randi(nTA,[nBS,nTA])), ...
    vB(randi(nTB,[nBS,nTB])),2),[],1);

% % 2. Shuffling
% All_Trials = [vA vB];
% RandSamplingIndexes1  = randi(length(All_Trials),[length(All_Trials)/2,1])
% DPSHUFF



