%scrit file name Maxratetemp
%purpose:
%该函数用于获得各升温速率下降解速率最快时的温度
function maxtemp=Maxratetemp(ii,datain)
maxtempdata=[];
maxtempdata(:,1)=datain(:,5*ii-3);
maxtempdata(:,2)=datain(:,5*ii);
maxtempdata(find(isnan(maxtempdata(:,1))),:)=[];
maxtempdata(find(isnan(maxtempdata(:,2))),:)=[];
[~,index]=min(maxtempdata(:,2));
maxtemp=maxtempdata(index,1);