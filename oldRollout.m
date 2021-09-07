%Rollout
tg=     [];%Start times of new rates
tg(1)=  ta(1);
rp=     zeros(5,4);%i=period, j=age group
rp(1,:)=[0,0,0,arate(1)];%First vax period - 65+ only

for i=1:3%Changes in rate/vax periods?
    k=5-i;%k in reverse order of age group
    j=length(tg);
        
    if (tg(j)+NNwv(k)/arate(1)) < ta(2) %no change in rate
        tg(j+1)=        tg(j)+NNwv(k)/arate(1);%Vaccinate entire age group at given rate
        rp(j+1,k-1)=    arate(1);
        
    elseif (ta(2)+(NNwv(k)-arate(1)*(ta(2)-tg(j)))/arate(2)) < ta(3) %first change in rate only
        tg(j+1)=        ta(2);
        tg(j+2)=        ta(2)+(NNwv(k)-arate(1)*(tg(j+1)-tg(j)))/arate(2);
        rp(j+1,k)=      arate(2);
        rp(j+2,k-1)=    arate(2);
    
    elseif j == i %first and second changes in rate
        tg(j+1)=        ta(2);
        tg(j+2)=        ta(3);
        tg(j+3)=        ta(3)+(NNwv(k)-arate(1)*(tg(j+1)-tg(j))-arate(2)*(tg(j+2)-tg(j+1)))/arate(3);
        rp(j+1,k)=      arate(2);   
        rp(j+2,k)=      arate(3);
        rp(j+3,k-1)=    arate(3);
        
    elseif j == i+1 %second change in rate only
        tg(j+1)=        ta(3);
        tg(j+2)=        ta(3)+(NNwv(k)-arate(2)*(tg(j+1)-tg(j)))/arate(3);
        rp(j+1,k)=      arate(3);
        rp(j+2,k-1)=    arate(3);
        
    else %no further changes in rate
        tg(j+1)=        tg(j)+NNwv(k)/arate(3);
        rp(j+1,k-1)=    arate(3);
        
    end    
end

if length(tg)<5
    tg(5:6)=Inf;
elseif length(tg)<6
    tg(6)=Inf;
end
rp(:,1)=0; 

vx.startp1= tg(1);
vx.startp2= tg(2);
vx.startp3= tg(3);
vx.startp4= tg(4);
vx.startp5= tg(5);
vx.end=     tg(6);

vx.aratep1= rp(1,:)';
vx.aratep2= rp(2,:)';
vx.aratep3= rp(3,:)';
vx.aratep4= rp(4,:)';
vx.aratep5= rp(5,:)';