function f=plotModelFit(data,fout,tvec,delta)%,hospData)
%Fout is output from fit - correct xdata only
%tcut=86;%data starts from tcut+1
xvec=(70:518)';%(80:80+length(fout)-1)';
data=data(1:length(xvec))*65648000/56286961;

maxy=max(max(data),max(fout+delta));
fs=10; lw=2;
cmap=lines(7);
cmap2=othercolor('BrBG4',2);
figure;
hold on

h1=bar(xvec,data,'edgecolor','k','linewidth',0.01);%,'linewidth',lw,'color',cmap(3,:));
h1.FaceColor=cmap2(1,:);%[222,184,135]/255;

%Conf intervals:
col1=cmap2(2,:);
xvec2=[xvec;flipud(xvec)];
inBetween=[fout-delta;flipud(fout+delta)];
fill(xvec2,inBetween,col1,'facealpha',.2);
plot(xvec,fout-delta,'-','linewidth',1,'color',col1);
plot(xvec,fout+delta,'-','linewidth',1,'color',col1);

maxy=max(fout+delta);

h2=plot(xvec,fout,'-','linewidth',lw,'color',cmap2(2,:));
%h3=plot(-4:length(hospData)-5,hospData,'linewidth',lw,'color',cmap(3,:));

for i=3:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxy],'--','linewidth',lw,'color',.5*[1,1,1])
end

set(gca,'fontsize',fs)
%xticks(tvec(6:end))
%xticklabels({'LD','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan'})
%xticklabels({'Aug 20','Sep 20','Oct 20','Nov 20','Dec 20','Jan 21','Feb 21','7th Mar 21'})
ax=gca;
ax.XTick=tvec;
ax.XTickLabels={'Jan 2020','Feb 2020', 'Lockdown 1','May 2020','June 2020','July 2020','Aug 2020','Sept 2020','Oct 2020','Nov 2020','Dec 2020','Jan 2021','Feb 2021','8th Mar 2021','29th Mar 2021','12th Apr 2021','17th May 2021','18th July 2021'};
ax.XTickLabelRotation=45;
xlabel('Time')
ylabel('Hospital occupancy')
axis([tvec(2),tvec(end),0,maxy])
legend([h1,h2],'Data','Fit','location','NW')
%legend([h1,h3,h2],'Data','Old data','Fit','location','NW')
grid on
%grid minor
box on
