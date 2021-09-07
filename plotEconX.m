function f=plotEconX(X,alpha,roadmap)
%Input with/without March 2020
march=0;
%roadmap=1;

%Lasagnes
X=X.^(1/alpha);
maxx=max(max(X));
maxx=max(1,maxx);
titleString=strcat('Workforce $\alpha$=',num2str(alpha));

[numSectors,numPeriods]=size(X);

monthInt=2;%Number of months between axis labels
if march==0 && roadmap==0
    xlabels={'Jan 2020','Lockdown 1','June 2020','Aug 2020','Oct 2020','Dec 2020','Feb 2021','Apr 2020','Jun 2021'};
elseif march==0 && roadmap==1
    xlabels={'Jan 2020','Feb 2020', 'Lockdown 1','May 2020','June 2020','July 2020','Aug 2020','Sept 2020','Oct 2020','Nov 2020','Dec 2020','Jan 2021','Feb 2021','3rd Mar 2021','29th Mar 2021','12th Apr 2021','17th May 2021','18th July 2021'};
    %[1,3,6,8,10,12,14:18];
elseif march==0
    xlabels={'Jan 2020','Feb 2020', 'Lockdown 1','May 2020','June 2020','July 2020','Aug 2020','Sept 2020','Oct 2020','Nov 2020','Dec 2020','Jan 2021','Feb 2021','Mar 2021','Apr 2021','May 2021','Jun 2021','July 2021'};
else
    xlabels={'Jan 2020','March 2020','May 2020','July 2020','Sept 2020','Nov 2020','Jan 2021','March 2020','May 2021'};
end
tp=0:numPeriods;
lt=length(tp);
[tc,xc]=meshgrid(tp,(0:numSectors));

X(end+1,:)=zeros(1,numPeriods);
X(:,end+1)=zeros(numSectors+1,1);
fs=10; lw=2;

f=figure('Units','centimeters','Position',[0 0 20 12]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',fs);

yvec=cumsum([0,3,23,1,9,4,3,1,9,4]');
yvec2=yvec([1,2,3,5,6,7,9,10]);
if march==0 && roadmap==1
    xvec=tp(1:end-1);%[1,3,6,8,10,12,14:18];
else
    xvec=tp(1:monthInt:end-1);
end

ylab=num2str(yvec2+1);
hold on
h=pcolor(tc,xc,X);
for i=1:lt
    plot(tp(i)*[1,1],[-1,64],'k-','linewidth',.2)
end
for i=1:length(yvec)
    plot([tp(1),tp(end)],yvec(i)*[1,1],'k-','linewidth',.2)
end
plot([tp(1),tp(end)],[63,63],'k-','linewidth',.2)
set(h,'edgecolor','none')
set(gca,'YDir','normal')
title(titleString)
xlabel('Month')
ylabel('Sector')
set(gca,'xtick',xvec,'xticklabels',xlabels,'XTickLabelRotation',45,'ytick',yvec2,'yticklabels',ylab);%{'PRE','LD',xlab(3:end,:)} %,'fontsize',fs
axis([0,numPeriods,0,63])
caxis([0,maxx])
colorbar
grid on
box on
%imagesc(X')
end