function f=ukModifications
alpha=0.59;
%hePrepCovid19: line 47 - pr.a (default to 1)
%hePrepCovid19: line 185 - mu/occupancy threshold
%Age-matrix already 4x4 in data
%heParamsAge - Npop fed in as 4x1
%heParamsAge_oc - no mention of the over-capacity shift
%Thd, Threc - as Sri Lanka
%heSim - number of sectors manually entered
%%
%Define roadmap
%Rows are months, columns are sectors
load('xJan20toMay21.mat','xJan20toMay21')
xvalue=xJan20toMay21;
nmonths=size(xvalue,1);
%Expand from 20-sector to 63-sector configuration:
repCol=[3,1,19,1,2,1,3,5,1,4,3,1,5,4,1,1,2,2,3,1];
repColSum=[0,cumsum(repCol)];
X=zeros(nmonths,63);
for i=1:length(repCol)
    X(:,repColSum(i)+1:repColSum(i+1))=repmat(xvalue(:,i),1,repCol(i));
end
X=X'/100;
%X(X>1)=1;
tvec=[0,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31];
tvec=cumsum(tvec)+1;
%Plot X's before removing March:
plotEconX(X,alpha);
%Remove March (in orfer to fit lockdown onset):
tvec=[tvec([1,2,3:length(tvec)])];
X(:,3)=[];

[numSectors,numInt]=size(X);
data.wfhAv=[zeros(2,numSectors);repmat(ddata64.wfhAv,numInt-2,1)];
%%
%School closures:
%Summer 2021
%%
%Pingdemic - quarantine

%%
%Vaccine rollout:
%Periods: Dec/rest of 65+/over-50/19-49
%https://www.bbc.co.uk/news/health-55274833
%(10*.8+10*.65+12*.6)/32=0.6781 - 18-49, proxy for 20-49
NNage=[4064198,12192593,36577778,13005432]';
prop50=15/45;
upOver50=.9;
upUnder50=.7;%pr.upUnder50;

tpoints=[336,367,412,457,579];%1st Dec, st Jan, 15th March, 1st Apr, 1st Aug
tdiff=diff(tpoints);
tvaxed=[7e5,.9*NNage(4)-7e5,.85*prop50*NNage(3),.68*(1-prop50)*NNage(3)];
trate=tvaxed/tdiff;
dataUK.atimes=tpoints;
dataUK.arates=trate;

tpoints=[367,457,538];%1st Jan, 1st April, 1st July
trate=[3e5,2e5,1e5];%tvaxed/tdiff;
dataUK.atimes=tpoints;
dataUK.arates=trate;
dataUK.uptake=[0,0,.8,.95];
%%
%Data:
%
load('d100320to310320.mat','d100320to310320')
load('d200320to090820.mat','d200320to090820')
load('d010820to060421.mat','d010820to060421')
load('d070421to110821.mat','d070421to110821')
dataOcc=[d100320to310320(1:10),d200320to090820(1:end-9),d010820to060421,d070421to110821];
save('dataOcc100320to110821.mat','dataOcc')
%}
load('dataOcc100320to110821.mat','dataOcc')
plotData(dataOcc)
%Other:


%%
%Possibility frontiers:

end

