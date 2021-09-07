function poptim=fitEpi2(ydata,X,data)
%For fitting deltas/betamod ONLY
%%
[a,b]=size(X);
X=reshape(X',a*b,1);
xdata=336:518;%**
ystart=xdata(1)-69;
ydata=ydata(ystart:ystart+length(xdata)-1);
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

R0=2.9931;
del1=0.3836;
pfit2=[0.4092    0.4597    0.3155    0.3488    0.5203    0.7642    0.6419    0.4454    0.6498];%First one overlaps del1  -to end of 2020
pfit3=[0.6743    0.4467    0.2600    0.4910    0.3401    0.6372];%First one overlaps pfit2(end)  - to end of May 2021
tvec=[-45.3987,32,87.8544,122,153,183,214,245,275,306,336,367,398,426,457,487,518];%,548,579,610,640,671,701,731];%**

X=X(1:63*(length(tvec)-1));
numSectors=63;
numInt=length(X)/numSectors;

[pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt,R0,del1);

%%
%Deltas from May (1 to 3=Jan/Feb/Apr):
ub=2.5*ones(1,numInt-10);%-3
x0=[pfit2(end),.5*ones(1,numInt-11)];%[0.5153    0.2634    0.3994    0.4911    0.7696    0.6395    0.4485];%0.5*ones(1,numInt-3);
lb=zeros(1,numInt-10);
%%

fun=@(params,xdata)sim2fit(params,data,xdata,X,pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta,tvec,pfit2);
%{
tic
rng default;%for reproducibility
options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);%,1e7,'MaxIterations',1e5);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata,'lb',lb,'ub',ub,'options',options);
ms=MultiStart;
poptim=run(ms,problem,1);
toc
%}
%{
tic
rng default;%for reproducibility
options=optimoptions(@nlinfit,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);%,1e7,'MaxIterations',1e5);
problem=createOptimProblem('nlinfit','beta0',x0,'objective',fun,'xdata',xdata,'ydata',ydata,'options',options);
ms=MultiStart;
[poptim,R,J,CovB,~,~]=run(ms,problem,1);
toc
%}
%
tic
%[poptim,~,R,~,~,~,J]=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
%[CovB=(1/(J'*J))*((R'*R)./(numel(ydata) - numel(poptim)));
[poptim,R,J,CovB,~,~]=nlinfit(xdata,ydata,fun,x0);
conf=nlparci(poptim,R,'jacobian',J);
[Ypred,delta]=nlpredci(fun,xdata,poptim,R,'Covar',CovB);
f=poptim;
g=conf;
toc
%}
ymod=sim2fit(poptim,data,xdata,X,pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta,tvec,pfit2);
f=figure('Units','centimeters','Position',[0 0 20 20]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',15);
hold on;
bar(xdata,ydata);
plot(xdata,ymod,'linewidth',2.5,'color','red');
for i=[1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731]
    plot(i*[1,1],[0,1.25*max(ydata)],'k-','linewidth',0.01);    
end
xlim([xdata(1),xdata(end)]);
ylim([0,1.25*max(ydata)]);
axis square;
box on;
grid on;
xlabel('Time');
ylabel('Hospital Occupancy');
title('Model Fit');
%}
end


function f=sim2fit(params,data,xdata,Xfit,pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta,tvec,pin)
    Wfit=Xfit.^(1/pr.a);
    %pr.betamod(4:length(params)+3)=params;
    pr.betamod=[1,1,pin(1:end-1),params];
    [simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec,0,data);

    t=simu(:,1)';
    h=simu(:,4)';

    f=interp1(t,h,xdata); 

end