function GiacStyle(varargin)
% GiacStyle('Color','r','Linewidth',2)
% by GB 2017
Param = {'Color', 'Linewidth'};
% Default values
col='k';
LW = 2;
for i = 1:size(Param,2)
    a = strcmp(varargin,Param{i});
    if sum(a)>0
        b = find(a(:)>0);
        switch i
            case 1
                col = varargin{b+1};
            case 2
                LW = varargin{b+1};
                
        end
    end
    set(gca,'Xtick',[],'Ytick',[],'Xcolor',col,'Ycolor',col,'linewidth',LW)
    
end