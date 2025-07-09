%scrit file name diffxT
%purpose:
%该函数用于获得不同升温速率聚合物降解特定转化率时不同热重曲线转化率时转化率对温度的微分,温度和时间，纵坐标升温速率，横坐标为转化率
function [dxdT,Tdata,tdata]=diffxT(datain,tempramp,conver,figans)

[~,col]=size(datain); dxdT=[];Tdata=[];
for i=1:col/5
    xTdata=[];xTdatacopy=[];
    xTdata(:,1)=datain(:,5*i-3);%温度摄氏度
    xTdata(:,1)=xTdata(:,1);%温度℃，不影响转化率对温度微分值
    xTdata(:,2)=datain(:,5*i);%质量百分比对温度微分
    xTdata(:,3)=datain(:,5*i-1);%质量百分比
    xTdata(:,4)=datain(:,5*i-4);%质量百分比
    xTdata(find(isnan(xTdata(:,1))),:)=[];
    xTdata(find(isnan(xTdata(:,2))),:)=[];
    xTdata(find(isnan(xTdata(:,3))),:)=[];
    xTdata(find(isnan(xTdata(:,4))),:)=[];
    %由质量百分比对温度微分求转化率对温度微分
    deltMper=xTdata(1,3)-xTdata(end,3);
    xTdata(:,2)=-1*xTdata(:,2)/deltMper;
    %对温度积分球转化率
    xTdatacopy(:,1)=cumtrapz(xTdata(:,1),xTdata(:,2));
    
    if strcmp(figans,'y')
        figure(i)
        [AX,H1,H2]=plotyy(xTdata(:,1),xTdata(:,2),xTdata(:,1),xTdatacopy(:,1));
        xlabel('T(k)');ylabel(AX(1),'dα/dT');ylabel(AX(2),'α');
        legend('dα/dT','α');
        set(H1,'LineStyle','-','color','k','linewidth',3);set(H2,'LineStyle','-','color','r','linewidth',3);
        titlename=strcat('升温速率',num2str(tempramp(i)),'K/min时转化率α对温度T微分图和α对温度变化图');
        title(titlename);
    end
    
    %获得每一个指定转化率时转化率对温度微分
    [row,~]=size(xTdatacopy);
    j=1;
    for k=2:row
        for ii=1:length(conver)
            if xTdatacopy(k-1,1)<conver(ii) && xTdatacopy(k,1)>=conver(ii)
                dxdT(i,j)=xTdata(k,2);%转化率
                Tdata(i,j)=xTdata(k,1)+273.15;%温度，K
                tdata(i,j)=xTdata(k,4);%时间
                j=j+1;
            end
        end
    end
    if length(dxdT(i,:))~=length(conver)
        error('\n第%d组数据有转化率对温度微分未全部找到\n',i);
    end
end



