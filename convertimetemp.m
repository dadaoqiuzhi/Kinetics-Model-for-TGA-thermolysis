%scrit file name convertimetemp
%purpose:
%�ú������ڻ�ø��¶��½���ﵽָ��ת����ʱ���õ�ʱ��
function conver_time_temp=convertimetemp(ii,jj,datain,conver,identifier)
if strcmp(identifier,'time') && mod(ii+4,5) ~= 0
    error('�����ʱ�������룬���飡����');
end
if strcmp(identifier,'temp') && mod(ii+3,5) ~= 0
    error('������¶������룬���飡����');
end

if strcmp(identifier,'time')
    converdata=[];
    converdata(:,1)=datain(:,ii-1);%ʱ��
    converdata(:,2)=datain(:,ii+2);%���ذٷֱ�
    converdata(find(isnan(converdata(:,1))),:)=[];
    converdata(find(isnan(converdata(:,2))),:)=[];
    convermax=converdata(1,2);convermin=converdata(end,2);
    deltcon=convermax-convermin;
    [raw,~]=size(converdata);converratio=conver(jj);
    for i=1:raw-1
        convera=(convermax-converdata(i,2))/deltcon;
        converb=(convermax-converdata(i+1,2))/deltcon;
        if convera < converratio && converb > converratio
            conver_time_temp=(converdata(i,1)+converdata(i+1,1))/2;
            break;
        elseif converb==converratio
            conver_time_temp=converdata(i+1,1);
            break;
        else
            continue;
        end
    end
elseif strcmp(identifier,'temp')
    converdata=[];
    converdata(:,1)=datain(:,ii);%�¶����϶�
    converdata(:,1)=converdata(:,1)+273.15;%�¶�k
    converdata(:,2)=datain(:,ii+2);%���ذٷֱ�
    converdata(find(isnan(converdata(:,1))),:)=[];
    converdata(find(isnan(converdata(:,2))),:)=[];
    convermax=converdata(1,2);convermin=converdata(end,2);
    deltcon=convermax-convermin;
    [raw,~]=size(converdata);converratio=conver(jj);
    for i=1:raw-1
        convera=(convermax-converdata(i,2))/deltcon;
        converb=(convermax-converdata(i+1,2))/deltcon;
        if convera < converratio && converb > converratio
            conver_time_temp=(converdata(i,1)+converdata(i+1,1))/2;
            break;
        elseif converb==converratio
            conver_time_temp=converdata(i+1,1);
            break;
        else
            continue;
        end
    end
else
    error('�Ƿ��ı�ʶ�����룬ֻ���ǡ�time�����ߡ�temp��������')
end

if isempty(conver_time_temp)
    fprintf('���ת����Ϊ%f\n',converb)
    error('��%d����û������Ҫ���ת����',ii+3);
end