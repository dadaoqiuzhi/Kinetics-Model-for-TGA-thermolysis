%scrit file name Maxratetemp
%purpose:
%�ú������ڻ�ø����������½����������ʱ���¶�
function maxtemp=Maxratetemp(ii,datain)
maxtempdata=[];
maxtempdata(:,1)=datain(:,5*ii-3);
maxtempdata(:,2)=datain(:,5*ii);
maxtempdata(find(isnan(maxtempdata(:,1))),:)=[];
maxtempdata(find(isnan(maxtempdata(:,2))),:)=[];
[~,index]=min(maxtempdata(:,2));
maxtemp=maxtempdata(index,1);