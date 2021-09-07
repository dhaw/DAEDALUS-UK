function poptim=fitEpi20(ydata,X,data)
%%

%load('../Daedalus/SR35.mat','data');
[a,b]=size(X);
X=reshape(X',a*b,1);

xdata=40:245;%26:553;
%xdata=(80:244);%20th March-31st August incl.%days of year
%xdata=(26:486);%26th January-30th April incl.
%ydata=table2array(ydata);%England only, 20th March-4th February incl.
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

%%

%xdata2=26:517;
xdata2=487:553;
ydata2=ydata(xdata2-25);

del0=0.50;%*ones(1,10);

%ub=[0,      3.10,   1.50*ones(1,2)];
%x0=[-50,    2.75,   del0           ];
%lb=[-100,   2.70,   0.00*ones(1,2)];%t0, R0, deltas
% ub=[1.50*ones(1,2)];
% x0=[del0          ];
% lb=[0.00*ones(1,2)];%deltas

%Order: R0, t0, tlockdown, delta(s);
%%

fun=@(params,xdata)sim2fit(params,data,xdata2,X);

tic
rng default;%for reproducibility
options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',10000000,'MaxIterations',10000);
%options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1,'MaxIterations',1);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata2,'ydata',ydata2,'lb',lb,'ub',ub,'options',options);
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
%axis square;
box on;
grid on;
%xticks([]);    
%yticks([]);
xlabel('Time');
ylabel('Hospital Occupancy');
title('Model Fit');

end


function f=sim2fit(params,data,xdata2,Xfit)

    %addpath('../Daedalus');
    
    numInt=length(Xfit)/length(data.G);
    
    tvec=[params(1),1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731];   
    %tvec=[-41.1072,1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731];              

    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt,params(2),[params(3)*ones(1,8),params(4:end),0.45*ones(1,7)]);
    %[pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt,2.7947,[0.5241*ones(1,8),0.9009,0.4402,0.4197,0.4100,0.4084,0.3277,0.2920,0.8650,params,0.45*ones(1,6)]);
    %pr.sw=0;
    
    data.wfhAv=[zeros(3,length(data.G));repmat(data.wfhAv,numInt-2,1)];
    
    Wfit=Xfit.^(1/pr.a);
    
    [simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec,0,data);

    % if simu(1,1)>1
    %     simu=[zeros(simu(1,1)-1,1);simu(:,4)];%???
    % else
    %     simu=simu;
    % end

    t=simu(:,1)';
    h=simu(:,4)';

    f=interp1(t,h,xdata2); 

end