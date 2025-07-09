%scrit file name convertimetemp
%purpose:
%该函数用于获得各温度下降解达到指定转化率时所用的时间
function conver_time_temp=convertimetemp(ii,jj,datain,conver,identifier)
if strcmp(identifier,'time') && mod(ii+4,5) ~= 0
    error('错误的时间列输入，请检查！！！');
end
if strcmp(identifier,'temp') && mod(ii+3,5) ~= 0
    error('错误的温度列输入，请检查！！！');
end

if strcmp(identifier,'time')
    converdata=[];
    converdata(:,1)=datain(:,ii-1);%时间
    converdata(:,2)=datain(:,ii+2);%热重百分比
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
    converdata(:,1)=datain(:,ii);%温度摄氏度
    converdata(:,1)=converdata(:,1)+273.15;%温度k
    converdata(:,2)=datain(:,ii+2);%热重百分比
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
    error('非法的标识符输入，只能是‘time’或者‘temp’！！！')
end

if isempty(conver_time_temp)
    fprintf('最大转化率为%f\n',converb)
    error('第%d数据没有满足要求的转化率',ii+3);
end