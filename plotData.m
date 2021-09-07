function f=plotData(dataOcc)
maxMonth=20;%August2021
maxOcc=100*ceil(max(dataOcc)/100);
tvec=[0,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31];
tvec=cumsum(tvec)+1;
tvecDays=tvec(3)+9:tvec(maxMonth)+10;

xlabels={'Jan 2020','March 2020','May 2020','July 2020','Sept 2020','Nov 2020','Jan 2021','March 2020','May 2021','July 2021'};

fs=10; lw=2;
figure
cmap=lines(7);
hold on
%{
for i=1:maxMonth
    plot(tvec(i)*[1,1],[0,maxOcc],'k-','linewidth',.5);
end
%}
plot(tvec(13)*[1,1],[0,maxOcc],'k-','linewidth',2);
plot(tvecDays,dataOcc,'-','color',cmap(2,:),'linewidth',lw)
xlabel('Time')
ylabel('Hospital Occupancy (England)')
set(gca,'fontsize',fs,'xtick',tvec(1:2:end),'xticklabels',xlabels,'XTickLabelRotation',45)
axis([tvecDays(1)-9,tvecDays(end),0,maxOcc])
grid on
box on
end