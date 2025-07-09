%scrit file name diffxT
%purpose:
%该函数用于获得聚合物降解特定转化率时不同热重曲线转化率对温度的微分和温度
function dxdT=diffxT(datain,tempramp,conver,avenum)

[~,col]=size(datain); dxdT=[];
for i=1:col/5
    xTdata=[];xTdatacopy=[];
    xTdata(:,1)=datain(:,5*i-3);%温度摄氏度
    xTdata(:,1)=xTdata(:,1)+273.15;%温度k
    xTdata(:,2)=datain(:,5*i-1);%质量百分比
    xTdata(find(isnan(xTdata(:,1))),:)=[];
    xTdata(find(isnan(xTdata(:,2))),:)=[];
    
    control=1;%去除相邻温度一样的数据
    j=1;[row,~]=size(xTdata);
    while control
        if j>=row;
            break;
        end
        if xTdata(j+1,1)==xTdata(j,1)
            xTdata(j+1,:)=[];
        else
            j=j+1;
        end
        [row,~]=size(xTdata);
    end
    
    
    [row,~]=size(xTdata);maxconver=xTdata(1,2)-xTdata(row,2);%计算转化率
    for j=2:row
        xTdata(j,2)=(xTdata(1,2)-xTdata(j,2))/maxconver;
    end
    xTdata(1,2)=0;
    xTdatacopy=xTdata;
    
    [row,~]=size(xTdata);%计算转化率对温度微分
    for j=2:row-1
        xTdata(j,2)=(xTdata(j+1,2)-xTdata(j,2))/(xTdata(j+1,1)-xTdata(j,1));
    end
    xTdata(1,2*i)=0;
    xTdata(row,2)=xTdata(row-1,2);
    
    
    %     j=1;control=1;[row,~]=size(xTdata);
    %     while control%去除奇异点
    %         if j>row-1;
    %             break;
    %         end
    %         if abs(xTdata(j+1,2)/xTdata(j,2))>abnormal
    %             xTdata(j+1,:)=[];
    %             j=j-1;
    %             if j==0;
    %                 j=1;
    %             end
    %         else
    %             j=j+1;
    %         end
    %         [row,~]=size(xTdata);
    %     end
    
    figure(i)
    plot(xTdata(:,1),xTdata(:,2),'-k')
    if ~ischar(avenum)
        hold on;
        xTdata(:,2)=medfilt1(xTdata(:,2),avenum);%平滑
        plot(xTdata(:,1),xTdata(:,2),'*r')
    end
    xlabel('T(k)');ylabel('dα/dT');
    legendstr=strcat('β=',num2str(tempramp(i)),'K/min');
    legend(legendstr);
    titlename=strcat('不同升温速率β条件下转化率α随温度T变化图');
    title(titlename);
    
    %获得每一个指定转化率时转化率对温度微分
    [row,~]=size(xTdatacopy);[row2,~]=size(xTdata);
    for k=2:row
        for ii=1:length(conver)
            if xTdatacopy(k-1,2)<conver(ii) && xTdatacopy(k,2)>=conver(ii)
                if xTdatacopy(k,2)>=0%确保微分大于0
                    for j=1:row2
                        if xTdata(j,1)==xTdatacopy(k,1)
                            dxdT(i,ii)=xTdata(k,2);
                        end
                    end
                elseif xTdatacopy(k,2)<0
                    for jj=k+1:row
                        if xTdatacopy(jj,2)>=0
                            for j=1:row2
                                if xTdata(j,1)==xTdatacopy(jj,1)
                                    dxdT(i,ii)=xTdata(jj,2);
                                end
                            end
                        end
                    end
                else
                    error('第%d组数据，转化率为%f,未有合法的转化率对温度微分值',i,conver(ii));
                end
            end
        end
    end
end
hold off;


