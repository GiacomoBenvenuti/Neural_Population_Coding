O1C= cat(3,DataTrialUF{s,1,1},DataTrialUF{s,1,2}) ;
O2C = cat(3,DataTrialUF{s,2,1},DataTrialUF{s,2,2}) ;
%----------------
clear sParam
sParam.Method =1;
sParam.DC_Sub=1;
sParam.Thr = 1;
 [dp ds Decoder vA vB ] = CalculateDiscrimination (O1C,O2C,sParam) ;
%%
mmin = min([vA vB])
vAn = vA - mmin;
vBn = vB - mmin;

mmax = max([vAn vBn])
vAn = vAn./mmax;
vBn = vBn./mmax;

vAn = vAn *2-1;
vBn = vBn *2-1;


x = -1:0.05:1
a = hist(vAn,x);
b = hist(vBn,x);

figure
bar(x,a,'r')
hold on
bar(x,b,'b')

set(gcf,'color','w')
box off

%%
%6 
for i =6
clf
tm =O1C(1:348,112:504,i) ;
tm = tm - mean(tm);
imagesc(tm)
grid
axis square
colormap(gray)
title(num2str(i))
pause(1)
end

%%
figure; imagesc(Decoder(1:348,112:504),[-5 5])

