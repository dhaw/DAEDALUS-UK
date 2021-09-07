function poptim=fitEpi1(ydata,X,data)
%%
[a,b]=size(X);
X=reshape(X',a*b,1);

xdata=70:122;%245;%70:...
ydata=ydata(0+(1:length(xdata)));%0+...

ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)
%%
%Order: R0, t0, tlockdown, delta(s);
ub=[5,0,92,1];
x0=[2.8,-57,87,.45];
lb=[1,-200,61,0];
%%
fun=@(params,xdata)sim2fit(params,data,xdata,X);

tic
rng default;%for reproducibility
options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',10000000,'MaxIterations',100000);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata,'lb',lb,'ub',ub,'options',options);
ms=MultiStart;
poptim=run(ms,problem,1);
toc

% [poptim,~,R,~,~,~,J]=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
% CovB=(1/(J'*J))*((R'*R)./(numel(ydata) - numel(poptim)));
% [poptim,R,J,CovB,~,~]=nlinfit(xdata,ydata,fun,x0);
% conf=nlparci(poptim,R,'jacobian',J);
% [Ypred,delta]=nlpredci(fun,xdata,poptim,R,'Covar',CovB);
% f=poptim;
% g=conf;

%Plotting
ymod=sim2fit(poptim,data,xdata,X);
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


function f=sim2fit(params,data,xdata,Xfit)
    numSectors=63;
    numInt=length(Xfit)/numSectors;
    tvec=[params(2),32,params(3),122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731];   
    %tvec=[-41.1072,1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731];              

    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt,params(1),params(4)*ones(1,15));
    %[pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt,2.7947,[0.5241*ones(1,8),0.9009,0.4402,0.4197,0.4100,0.4084,0.3277,0.2920,0.8650,params,0.45*ones(1,6)]);
    %pr.sw=0;
    
    Wfit=Xfit.^(1/pr.a);
    
    [simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:numInt+1),0,data);

    t=simu(:,1)';
    h=simu(:,4)';

    f=interp1(t,h,xdata); 

end