%scrit file name diffxT
%purpose:
%�ú������ڻ�ò�ͬ�������ʾۺ��ｵ���ض�ת����ʱ��ͬ��������ת����ʱת���ʶ��¶ȵ�΢��,�¶Ⱥ�ʱ�䣬�������������ʣ�������Ϊת����
function [dxdT,Tdata,tdata]=diffxT(datain,tempramp,conver,figans)

[~,col]=size(datain); dxdT=[];Tdata=[];
for i=1:col/5
    xTdata=[];xTdatacopy=[];
    xTdata(:,1)=datain(:,5*i-3);%�¶����϶�
    xTdata(:,1)=xTdata(:,1);%�¶ȡ棬��Ӱ��ת���ʶ��¶�΢��ֵ
    xTdata(:,2)=datain(:,5*i);%�����ٷֱȶ��¶�΢��
    xTdata(:,3)=datain(:,5*i-1);%�����ٷֱ�
    xTdata(:,4)=datain(:,5*i-4);%�����ٷֱ�
    xTdata(find(isnan(xTdata(:,1))),:)=[];
    xTdata(find(isnan(xTdata(:,2))),:)=[];
    xTdata(find(isnan(xTdata(:,3))),:)=[];
    xTdata(find(isnan(xTdata(:,4))),:)=[];
    %�������ٷֱȶ��¶�΢����ת���ʶ��¶�΢��
    deltMper=xTdata(1,3)-xTdata(end,3);
    xTdata(:,2)=-1*xTdata(:,2)/deltMper;
    %���¶Ȼ�����ת����
    xTdatacopy(:,1)=cumtrapz(xTdata(:,1),xTdata(:,2));
    
    if strcmp(figans,'y')
        figure(i)
        [AX,H1,H2]=plotyy(xTdata(:,1),xTdata(:,2),xTdata(:,1),xTdatacopy(:,1));
        xlabel('T(k)');ylabel(AX(1),'d��/dT');ylabel(AX(2),'��');
        legend('d��/dT','��');
        set(H1,'LineStyle','-','color','k','linewidth',3);set(H2,'LineStyle','-','color','r','linewidth',3);
        titlename=strcat('��������',num2str(tempramp(i)),'K/minʱת���ʦ����¶�T΢��ͼ�ͦ����¶ȱ仯ͼ');
        title(titlename);
    end
    
    %���ÿһ��ָ��ת����ʱת���ʶ��¶�΢��
    [row,~]=size(xTdatacopy);
    j=1;
    for k=2:row
        for ii=1:length(conver)
            if xTdatacopy(k-1,1)<conver(ii) && xTdatacopy(k,1)>=conver(ii)
                dxdT(i,j)=xTdata(k,2);%ת����
                Tdata(i,j)=xTdata(k,1)+273.15;%�¶ȣ�K
                tdata(i,j)=xTdata(k,4);%ʱ��
                j=j+1;
            end
        end
    end
    if length(dxdT(i,:))~=length(conver)
        error('\n��%d��������ת���ʶ��¶�΢��δȫ���ҵ�\n',i);
    end
end



