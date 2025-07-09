%scrit file name thermo_kinetic_fit
%purpose:
%This program is used to fit models of polymer pyrolysis kinetics 
%version 1;2018.9.16
%reference literatures
%Pyrolysis kinetics and thermal decomposition behavior of polycarbonate - a TGA-FTIR study
%High-temperature pyrolysis simulation of acrylonitrile-butadiene model compound with experimental evidence
%Initiation mechanisms and kinetic analysis of the isothermal decomposition of poly(amethylstyrene): a ReaxFF molecular dynamics study
%***Isoconversional Kinetics of Thermally Stimulated Processes
%Formation of bisphenol A by thermal degradation of poly(bisphenol A carbonate)
%Thermal degradation mechanism and kinetics of polycarbonate/silica nanocomposites
%���ط��ⶨ�ۺ����Ƚ��ⷴӦ����ѧ������չ;�۱�ϩ-����̼��Ƹ��ϲ��ϵ��Ƚ��⶯��ѧ
disp('��ӭʹ�ñ�����--by ��ǿ@�Ĵ���ѧ�߷��ӿ�ѧ�빤��ѧԺ����ܽ��ڿ����飬liubinqiang@163.com')
fprintf('��������Ҫ������Ͼۺ����Ƚ��⶯��ѧģ�ͣ�����:\n**����ģ��\n11.Conventional----lnt=ln[G(��)/A]+E/RT,��<=0.05')
fprintf('\n12.����ģ��f(��) vs ������ο��ļ�����ͼƬ��ȷ��13�ֳ��÷��� @reference 4');
fprintf('\n**�ǵ���ģ��\n21.Coast-Redfern(CR)----ln[(-ln(1-��)/T^2]=ln(AR/��E)-E/RT');
fprintf('\n22.Friedman----ln(d��/dt)=ln(��d��/dT)=ln��+ln(d��/dT)=ln(Af(��))-E/RT');
fprintf('\n23.Kissinger----ln(��/Tp^2)=ln(AR/E)+ln[n(1-xp)^(n-1)]-E/RTp��ln(��/Tp^2)=-E/RTp+ln(AR/E)');
fprintf('\n24.Flynn-Wall-Ozawa (FWO)----ln(��)=ln(AE/Rg(��))-5.331-1.052E/RT��ln(��)=C-1.052*E/(RT)');
fprintf('\n25.Kissinger-Akahira-Sunose (KAS)----ln(��/T^2)=ln(AE/Rg(��))-E/RT;ln(��/T^2)=C-E/RT');
fprintf('\n26.first-order (unimolecular)ͬһ��������ͬһ��Ʒ��ͬת����ʱ----ln(d��/dT/(1-��))=lnA-E/RT');
fprintf('\n27.Freeman-Carroll�������΢�ַ���----��ln(d��/dt)/��ln(1-��)=n-E/R*��(1/T)/��ln(1-��),ͬһ����');
fprintf('\n28.Zivkovic��----ln[ln(1/(1-��))/t]=lnA-E/RT\n����ͬ����');

fprintf('##�뽫ʵ�����ݰ��¶ȸߵͻ����������ʿ�������input_data Excel�ļ��У�ÿ���������У����ΰ���ʱ�䣬�¶ȣ�����(mg)�������ٷֱȺ��������¶ȵ�һ��΢��\n')
cho=input('\n������ģ�ͱ�ţ�\n');
if exist('datain','var')
    reload=input('\ndatain�Ѵ��ڣ��Ƿ����µ����������ݣ�y/n,������:\n');
    reload=lower(reload);
    if strcmp(reload,'y');
    datain=xlsread('input_data.xlsx');
    end
else
    datain=xlsread('input_data.xlsx');
end
[~,col]=size(datain);group=col/5;
if group ~= round(group)
    error('�����������ݣ���Ҫ���ṩinput_data���ݱ�\n');
end
switch cho
    case 11
        tempramp=input('������������ݵ����¶ȣ���λ���϶�,���ÿո�����������ţ�\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);tempramp=str2double(tempramp);
        tempramp=tempramp+273.15;%�¶�Ϊk�ľ���
        if length(tempramp)~=group
            error('\n�����¶�������ʵ������������һ�£����飡����');
        end
        conver=input('��ָ��ת����,����ת�������ÿո�����������ţ�\n');
        figans=input('\n�Ƿ���ʾת���ʶ��¶�΢��ͼ��ת�������¶ȱ仯ͼ��y/n��������:\n');
        figans=lower(figans);
        disp('thermo_kinetic_fit��������11.Conversion�У���ȴ�...')
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%ת���ʾ���
        tempconver={};tempconver{1,1}='T/��';
        for i=1:length(conver)%��ʼ��ת���ʱ�ͷ
            tempconver{1,i+1}=conver(i);
        end
        [~,~,tdata]=diffxT(datain,tempramp,conver,figans);
        for i=1:group%��ʼ����һ���¶�
            tempconver{i+1,1}=tempramp(i);
            for j=1:length(conver)
                tempconver{i+1,j+1}=tdata(i,j);%�����������´ﵽ��Ӧת��������ʱ��
            end
        end
        %��Ϲ�ʽ
        xydata=[];[rawtc,coltc]=size(tempconver);%�¶ȵ�����ʱ���������
        for i=1:rawtc-1
            xydata(i,1)=1/tempconver{i+1,1};%x
            for j=1:coltc-1
                 xydata(i,j+1)=log(tempconver{i+1,j+1});%y
            end
        end
        X=[ones(rawtc-1,1) xydata(:,1)];%������Ϻ�������
        figure( 'Name', '11.Conventional' );% Plot fit with data
        legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
        xlabel( '1/T (k)' );ylabel( 'lnt ��min��' );% Label axes
        axis([min(xydata(:,1))*0.9 max(xydata(:,1))*1.2 min(xydata(:,2))*0.9 max(xydata(:,2))*1.2]);
        title('Lifetime prediction of polymer materials --11.Conversion');
        
        fittingdata={};Y=[];
        fittingdata{1,1}='Tiso/other';fittingdata{1,2}='E';fittingdata{1,3}='ln[G(��)/A]';fittingdata{1,4}='R^2';
        for i=2:coltc
            [b,bint,r,rint,stats]=regress(xydata(:,i),X);%���
            fittingdata{i,1}=conver(i-1);
            fittingdata{i,2}=b(2)*8.314;%���
            fittingdata{i,3}=b(1);%������
            fittingdata{i,4}=stats(1,1);%�ع�ϵ��
            Y(:,i-1)=b(1)+b(2)*xydata(:,1);
            plot(xydata(:,1),xydata(:,i),'ok',xydata(:,1),Y(:,i-1),'-r','linewidth',3);
            hold on;
        end
        hold off;
        fprintf('\nthermo_kinetic_fit��������11.Conversion������������û��,������ͻع�ϵ��R^2�洢��fittingdata��')
        fprintf('\nXֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��xydata��xydata��Y��,�����������ת���ʶ��¶�΢�ֺͱ仯ͼ');
        fprintf('\n11.Conversion������õ�ԭʼ���ݴ洢��xydata�����У���ϵ�Y���ݴ洢��Y�����У�ע�����𲻴�ʱͼ�ϵĵ���߻��غ�\n')
        clear bint col coltc conver E group h i j r rawtc rint stats tempramp X cho group
    
    case 12
        fprintf('\n�μ��ļ���ͼƬ������4���ݣ���ʱδд���룬�����ڴ�');
        break;
        
        
        
    case 21
        tempramp=input('��������������������ʣ���λ���϶�/min���ÿպŸ����Ҵ�����:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%��������
        if length(tempramp)~=group
            error('\n��������������ʵ������������һ�£����飡����');
        end
        conver=input('��ָ��ת����,����ת�������ÿո�����������ţ�\n');
        disp('thermo_kinetic_fit��������21.Coast-Redfern(CR)�У���ȴ�...')
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%ת���ʾ���
        tempconver={};tempconver{1,1}='��/��';
        for i=1:length(conver)%��ʼ��ת���ʱ�ͷ
            tempconver{1,i+1}=conver(i);
        end
        for i=1:group%��ʼ����һ���¶�
            tempconver{i+1,1}=tempramp(i);
            for j=1:length(conver)
                tempconver{i+1,j+1}=convertimetemp(5*i-3,j,datain,conver,'temp');%�����������´ﵽ��Ӧת����ʱ���¶�,��λk
            end
        end
        z=[];x=[];%yֵ
        for i=1:length(tempramp)
            for j=1:length(conver)
                z(j,i)=log(-log(1-conver(j))/(tempconver{i+1,j+1}*tempconver{i+1,j+1}));
                x(j,i)=1/tempconver{i+1,j+1};
            end
        end
        
        output={};subplotnum=(length(conver)+mod(length(conver),2))/2;
        fitans=input('��ѡ���������Ϸ�ʽ���,�Ƽ�2��3��1.fit����ģ����ϣ�2.ninfit��ϣ�3.lsqcurvefit���:\n');
        startpoint=input('������ģ�Ͳ����ĳ�ֵ����,�����Ƽ����3��ȡ��ʼ�㣬��ⷽ��,����Ҫ����ʷ����ֵΪ[3e+16,4e+4]��\n');
        if fitans==1
            y=[];
            for i=1:size(x,1)
                for j=1:length(tempramp)
                    y(i,j)=tempramp(j);
                end
            end
            fun='log(a/y)-e*x';
            ft=fittype(fun,'independent',{'y','x'},'dependent','z');
            opts=fitoptions('Method', 'NonlinearLeastSquares','Lower',[-Inf,-Inf],'upper',[Inf,inf]);% Set up fittype and options.
            opts.StartPoint = startpoint;%�ʵ�������ʼ�㣬�ӽ����մ�
            Z=[];output={};output{1,1}='��';output{1,2}='A';output{1,3}='E';output{1,4}='R^2';
            for i=1:length(tempramp)
                [fitresult,gof] = fit([y(:,i),x(:,i)],z(:,i),ft,opts);
                output{i+1,1}=tempramp(i);
                output{i+1,2}=fitresult.a;output{i+1,3}=fitresult.e;
                fprintf('aΪ%f,eΪ%f',fitresult.a,fitresult.e);
                output{i+1,4}=gof.rsquare;
                opts.StartPoint = cell2mat(output(i+1,2:3));
                [fitresult,gof] = fit([y(:,i),x(:,i)],z(:,i),ft,opts);%�ٴ���ϣ���С���
                Z(:,i)=log(fitresult.a./y(:,i))-fitresult.e*x(:,i);
                output{i+1,2}=fitresult.a*fitresult.e;output{i+1,3}=fitresult.e*8.314;output{i+1,4}=gof.rsquare;
                fprintf('��������Ϊ%fʱ��ת���ʵõ����Ϊ%fJ/mol��R^2Ϊ%f\n',tempramp(i),output{i+1,3},output{i+1,4})
                h=figure(1);
                set(h,'Name','21.Coast-Redfern(CR)-nolinear');
                subplot(2,subplotnum,i)
                plot( x(:,i),z(:,i),'ok',x(:,i),Z(:,i),'*-r','linewidth',3);
                legend( 'origin data','fitting data', 'Location', 'NorthEast' );
                xlabel( '1/T (k)' );ylabel( 'ln[(-ln(1-��)/T^2]' );% Label axes
                %axis([min(x(:,i))*0.9 max(x(:,i))*1.2 max(z(:,i))*1.2 min(z(:,i))*0.9]);%ע��z�����Ϊ��ֵ
                titlename=strcat('21.Coast-Redfern(CR)',' with ��=',num2str(tempramp(i)),'K/min');
                title(titlename);
            end
            fprintf('21.Coast-Redfern(CR)�������н��������E�Ͳ���A��R^2���д�����output��,xΪ�¶ȵ���1/k��zΪԭʼzֵ��ZΪ���zֵ');
            
            cho=input('\n�����������������ܺã���ֱ�Ӷ�1/T����������ϣ��Ƿ����y/n��,�������ʼ�㣬������:\n');
            cho=lower(cho);
            if strcmp(cho,'y')
                output2={};output2{1,1}='��';output2{1,2}='ln(AR/��E)';output2{1,3}='E';output2{1,4}='R^2';
                Y=[];
                for i=1:length(tempramp)
                    ft = fittype( 'poly1' );%Ҳ�ɲ���polyfit
                    [fitresult, gof]=fit(x(:,i),z(:,i),ft);%���
                    output2{i+1,1}=tempramp(i);output2{i+1,2}=fitresult.p2;output2{i+1,3}=-8.314*fitresult.p1;output2{i+1,4}=gof.rsquare;
                    fprintf('ת����Ϊ%fʱ�õ����Ϊ%fJ/mol��R^2Ϊ%f\n',tempramp(i),output2{i+1,3},output2{i+1,4})
                    Y(:,i)=fitresult.p1*x(:,i)+fitresult.p2;
                    h=figure(2);
                    set(h,'Name','21.Coast-Redfern(CR)-linear');
                    subplot(2,subplotnum,i)
                    plot( x(:,i),z(:,i),'ok',x(:,i),Y(:,i),'-r','linewidth',3);
                    legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
                    xlabel( '1/T (k)' );ylabel( 'ln[(-ln(1-��)/T^2]' );
                    titlename=strcat('21.Coast-Redfern(CR)',' with ��=',num2str(tempramp(i)),'K/min');
                    title(titlename);
                end
                fprintf('21.Coast-Redfern(CR)�������н���������ln(AR/��E)�ͻ��E��R^2���д�����output2��,xΪ�¶ȵ���1/k��zΪԭʼzֵ��YΪ���zֵ');
            end
            
            
        elseif fitans==2
            Z=[];StartPoint=startpoint;beta={};%��ʼ�����Ҫ
            beta{1,1}='��/other';beta{2,1}='A';beta{3,1}='E';
            for j=1:size(x,1)
                alpha=strcat('r',num2str(j));
                beta{3+j,1}=alpha;
            end
            for i= 1:length(tempramp)
                X=x(:,i);%yΪ�������ʵ�����k;xΪ��Ӧת�����¶ȵ�����1/k
                fz=@(beta,X)log(beta(1)/tempramp(i))-beta(2)*X(:,1);
                [b,r]=nlinfit(X,z(:,i),fz,StartPoint);
                StartPoint=b;
                [b,r]=nlinfit(X,z(:,i),fz,StartPoint);%�ٴ����
                Z(:,i)=log(b(1)/tempramp(i))-b(2)*x(:,i);%��ϼ���ֵ
                beta{1,i+1}=tempramp(i);
                beta{2,i+1}=b(1)*b(2);%�õ�A
                beta{3,i+1}=b(2)*8.314;%�õ�E
                for j=1:length(r);%��¼�в�r
                    beta{3+j,i+1}=r(j);
                end
                h=figure(1);
                set(h,'Name','21.Coast-Redfern(CR)-nolinear');
                subplot(2,subplotnum,i)
                plot( x(:,i),z(:,i),'ok',x(:,i),Z(:,i),'*-r','linewidth',3);%��Ϊ��ά�ռ䣬��Ϊ��ά
                legend( 'origin data','fitting data', 'Location', 'NorthEast' );
                xlabel( '1/T (k)' );ylabel( 'ln[(-ln(1-��)/T^2]' );% Label axes
                %axis([min(x(:,i))*0.9 max(x(:,i))*1.2 max(z(:,i))*1.2 min(z(:,i))*0.9]);%ע��z�����Ϊ��ֵ
                titlename=strcat('21.Coast-Redfern(CR)',' with ��=',num2str(tempramp(i)),'K/min');
                title(titlename);  
            end
            fprintf('21.Coast-Redfern(CR)�������н���������A�ͻ��E�������ݵ�в�r���д�����beta������\n')
       
        elseif fitans==3
            y=[];
            for i=1:size(x,1)
                for j=1:length(tempramp)
                    y(i,j)=tempramp(j);
                end
            end
            StartPoint=startpoint;beta=[];
            beta{1,1}='��/other';beta{2,1}='A';beta{3,1}='E';beta{4,1}='sum_r_res^2';
            for i= 1:length(tempramp)
                X=[y(:,i),x(:,i)];%yΪ�������ʵ�����k;xΪ��Ӧת�����¶ȵ�����1/k
                fz=@(beta,X)log(beta(1)./X(:,1))-beta(2)*X(:,2);
                [b,resnorm]=lsqcurvefit(fz,StartPoint,X,z(:,i));
                StartPoint=b;
                [b,resnorm]=lsqcurvefit(fz,StartPoint,X,z(:,i));%�ٴ����
                Z(:,i)=log(b(1)./y(:,1))-b(2)*x(:,i);%��ϼ���ֵ
                beta{1,i+1}=tempramp(i);
                beta{2,i+1}=b(1)*b(2);%�õ�A
                beta{3,i+1}=b(2)*8.314;%�õ�E
                beta{4,i+1}=resnorm(1);%�в�ƽ����
                h=figure(1);
                set(h,'Name','21.Coast-Redfern(CR)-nolinear');
                subplot(2,subplotnum,i)
                plot( x(:,i),z(:,i),'ok',x(:,i),Z(:,i),'*-r','linewidth',3);%��Ϊ��ά�ռ䣬��Ϊ��ά
                legend( 'origin data','fitting data', 'Location', 'NorthEast' );
                xlabel( '1/T (k)' );ylabel( 'ln[(-ln(1-��)/T^2]' );% Label axes
                %axis([min(x(:,i))*0.9 max(x(:,i))*1.2 max(z(:,i))*1.2 min(z(:,i))*0.9]);%ע��z�����Ϊ��ֵ
                titlename=strcat('21.Coast-Redfern(CR)',' with ��=',num2str(conver(i)),'K/min');
                title(titlename); 
            end
            fprintf('21.Coast-Redfern(CR)�������н���������A�ͻ��E���в�rƽ���Ͱ��д�����beta��\n');
            fprintf('\nXֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��x��z��Z��\n');
        end
        
        clear col conver fitresult ft fun gof h i j opts subplotnum tempramp titlename group reload cho fz fitans
         
    case 22
        disp('����regress����NonlinearLeastSquares������ln(d��/dT)=ln(Af(��))-E/RT-ln�½��ж�Ԫ�������')
        tempramp=input('��������������������ʣ���λ���϶�/min���ÿպŸ����Ҵ�����:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%��������
        if length(tempramp)~=group
            error('\n��������������ʵ������������һ�£����飡����');
        end
        conver=input('��ָ��ת����,����ת�������ÿո�����������ţ�\n');
        figans=input('\n�Ƿ���ʾת���ʶ��¶�΢��ͼ��ת�������¶ȱ仯ͼ��y/n��������:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%ת���ʾ���
        %avenum=input('������ת���ʶ��¶�΢��������ƽ����ֵ������С��5,��ƽ��������n�������ţ�\n');
        disp('thermo_kinetic_fit��������22.Friedman�У���ȴ�...')
        dxdT={};dxdT{1,1}='��/��';dxdT(1,2:length(conver)+1)=num2cell(conver);
        for i=1:length(tempramp)
            dxdT{i+1,1}=tempramp(i);
        end
        [dx_dT,Tdata,~]=diffxT(datain,tempramp,conver,figans);%�̶��������ʣ���ͬ�¶�ʱ��ת���ʶ��¶ȵ�΢��ֵ
        for i=1:length(tempramp)
            for j=1:length(conver)
                dxdT{i+1,j+1}=dx_dT(i,j);
            end
        end
        z=log(dx_dT);%z=ln(d��/dT)
        x=[];%�����������´ﵽ��Ӧת���ʵ��¶�
        for i=1:group
            for j=1:length(conver)
                x(i,j)=Tdata(i,j);%�����������´ﵽ��Ӧת����ʱ���¶�,��λk
                x(i,j)=1/x(i,j);%�¶�k����
            end
        end
        y=tempramp';
        y=log(y);
        
        fittingdata={};Y=[];Ycal=[];
        fittingdata{1,1}='other/��';fittingdata{2,1}='E';fittingdata{3,1}='ln(Af(��))';fittingdata{4,1}='R^2';
        for i=1:length(conver)
            figname=strcat('ת���ʦ�Ϊ',num2str(conver(i)));
            figure('Name',figname)
            X=[ones(length(tempramp),1) x(:,i)];%������Ϻ�������
            Y(:,i)=z(:,i)+y(:,1);%���̱�ʽΪln(d��/dT)+ln��=ln(Af(��))-E/RT
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%���
            fittingdata{1,i+1}=conver(i);
            fittingdata{2,i+1}=-b(2)*8.314;%���
            fittingdata{3,i+1}=b(1);%������
            fittingdata{4,i+1}=stats(1,1);%�ع�ϵ��
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'-r','linewidth',3);
            xlabel( '1/T (k)' );ylabel( 'ln(d��/dt)+ln��' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            title('Lifetime prediction of polymer materials --22.Friedman');
        end
        fprintf('thermo_kinetic_fit��������22.Friedman���н��������E�Ͳ���ln(Af(��))���ع�ϵ�����д�����fittingdata��\n');
        fprintf('\nXֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��x��Y��Ycal��\n');
        clear avenum control figname fitresult ft fun group i j reload cho bint r rint stats b col dx_dT 
        
    case 23
        tempramp=input('��������������������ʣ���λ���϶�/min���ÿպŸ����Ҵ�����:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%��������
        if length(tempramp)~=group
            error('\n��������������ʵ������������һ�£����飡����');
        end
        disp('thermo_kinetic_fit��������23.Kissinger�У���ȴ�...');
        Tp=[];
        for i=1:group%��󽵽������¶Ⱦ���
            Tp(i)=Maxratetemp(i,datain);
        end
        Y=[];%ln(��/Tp^2)
        X=[];%1/Tp
        for i=1:group
            Y(i)=log((tempramp(i)-273.15)/Tp(i)/Tp(i));%�µ�λΪ��/min,E������Ŷԣ������׹涨��λΪk/min����E����������������
            X(i)=1/Tp(i);
        end
        
        Y=Y';X=X';Xx=[ones(size(X)) X];y=[];
        disp('\n����regress�������\n')%������ϣ�Y=(-E/R)*x������������
        [b,bint,r,rint,stats]=regress(Y(:,1),Xx(:,:));%���
        E=-8.314*b(2);
        y(:,1)=b(2)*X(:,1)+b(1);
        plot(X(:,1),Y(:,1),'ok',X(:,1),y(:,1),'*-r','linewidth',3);
        xlabel( '1/Tp (k^-1)' );ylabel( 'ln(��/Tp^2) (��/K^2)' );% Label axes
        legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
        title('Lifetime prediction of polymer materials --23.Kissinger');
        fprintf('thermo_kinetic_fit��������23.Kissinger���н��������EΪ%.4fJ/mol���ع�ϵ��R^2Ϊ%.4f\n',E,stats(1,1));
        fprintf('\nXֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��X��Y��y��\n');
        %clear cho col group i X y Y bint r rint stats
    case 24
        tempramp=input('��������������������ʣ���λ���϶�/min���ÿպŸ����Ҵ�����:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%��������
        if length(tempramp)~=group
            error('\n��������������ʵ������������һ�£����飡����');
        end
        conver=input('��ָ��ת����,����ת�������ÿո�����������ţ�\n');
        figans=input('\n�Ƿ���ʾת���ʶ��¶�΢��ͼ��ת�������¶ȱ仯ͼ��y/n��������:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%ת���ʾ���
        disp('thermo_kinetic_fit��������24.Flynn-Wall-Ozawa (FWO)�У���ȴ�...');
        x=[];%�����������´ﵽ��Ӧת���ʵ��¶�
        [~,Tdata,~]=diffxT(datain,tempramp,conver,figans);
        for i=1:group
            for j=1:length(conver)
                x(i,j)=Tdata(i,j);%�����������´ﵽ��Ӧת����ʱ���¶�,��λk
                x(i,j)=1/x(i,j);%�¶�k����,1/T
            end
        end
        y=tempramp'-273.15;y=log(y);%ln��,��/min��E������Ŷԣ����������������ģ���Ҳ�����׹涨��λΪk/min����E����������������
        fittingdata={};Y=[];Ycal=[];X=[];
        fittingdata{1,1}='other/��';fittingdata{2,1}='E';fittingdata{3,1}='C';fittingdata{4,1}='R^2';
        for i=1:length(conver)
            figname=strcat('ת���ʦ�Ϊ',num2str(conver(i)));
            figure('Name',figname)
            X=[ones(length(tempramp),1) x(:,i)];%������Ϻ�������
            Y(:,i)=y(:,1);
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%���
            fittingdata{1,i+1}=conver(i);
            fittingdata{2,i+1}=-b(2)*8.314/1.052;%���
            fittingdata{3,i+1}=b(1);%������
            fittingdata{4,i+1}=stats(1,1);%�ع�ϵ��
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '1/T (k)' );ylabel( 'ln��(��/min)' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            title('Lifetime prediction of polymer materials --24.Flynn-Wall-Ozawa (FWO)');
        end
        fprintf('thermo_kinetic_fit��������24.Flynn-Wall-Ozawa (FWO)���н��������E�Ͳ���C���ع�ϵ�����д�����fittingdata��\n')
        fprintf('Xֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��x��Y��Ycal��\n');
        clear figname fun group i j reload cho bint r rint stats
        
    case 25
        tempramp=input('��������������������ʣ���λ���϶�/min���ÿպŸ����Ҵ�����:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%��������
        if length(tempramp)~=group
            error('\n��������������ʵ������������һ�£����飡����');
        end
        conver=input('��ָ��ת����,����ת�������ÿո�����������ţ�\n');
        figans=input('\n�Ƿ���ʾת���ʶ��¶�΢��ͼ��ת�������¶ȱ仯ͼ��y/n��������:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%ת���ʾ���
        disp('thermo_kinetic_fit��������25.Kissinger-Akahira-Sunose (KAS)�У���ȴ�...');
        x=[];%�����������´ﵽ��Ӧת���ʵ��¶�
        [~,Tdata,~]=diffxT(datain,tempramp,conver,figans);
        for i=1:group
            for j=1:length(conver)
                x(i,j)=Tdata(i,j);%�����������´ﵽ��Ӧת����ʱ���¶�,��λk
                x(i,j)=1/x(i,j);%�¶�k����,1/T
            end
        end
        y=tempramp';
        fittingdata={};Y=[];Ycal=[];X=[];
        fittingdata{1,1}='other/��';fittingdata{2,1}='E';fittingdata{3,1}='C';fittingdata{4,1}='R^2';
        for i=1:length(conver)
            figname=strcat('ת���ʦ�Ϊ',num2str(conver(i)));
            figure('Name',figname)
            X=[ones(length(tempramp),1) x(:,i)];%������Ϻ�������
            Y(:,i)=log((y(:,1)-273.15).*x(:,i).*x(:,i));%ln(��/T^2)
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%���
            fittingdata{1,i+1}=conver(i);
            fittingdata{2,i+1}=-b(2)*8.314;%���
            fittingdata{3,i+1}=b(1);%������
            fittingdata{4,i+1}=stats(1,1);%�ع�ϵ��
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '1/T (k)' );ylabel( 'ln(��/T^2) (��/K^2)' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            title('Lifetime prediction of polymer materials --25.Kissinger-Akahira-Sunose (KAS)');
        end
        fprintf('thermo_kinetic_fit��������25.Kissinger-Akahira-Sunose (KAS)���н��������E�Ͳ���C���ع�ϵ�����д�����fittingdata��\n')
        fprintf('Xֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��x��Y��Ycal��\n');
        %clear figname fun group i j reload cho bint r rint stats
        
    case 26
        tempramp=input('��������������������ʣ���λ���϶�/min���ÿպŸ����Ҵ�����:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%��������
        if length(tempramp)~=group
            error('\n��������������ʵ������������һ�£����飡����');
        end
        conver=input('������Ϊÿ������ָ������ת����,���ÿո�����������ţ�\n');
        figans=input('\n�Ƿ���ʾת���ʶ��¶�΢��ͼ��ת�������¶ȱ仯ͼ��y/n��������:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%ת���ʾ���
        disp('thermo_kinetic_fit��������26.first-order (unimolecular)�У���ȴ�...')
        dxdT={};dxdT{1,1}='��/��';dxdT(1,2:length(tempramp)+1)=num2cell(tempramp);
        for i=1:length(conver)
            dxdT{i+1,1}=conver(i);
        end
        [dxdTdata,Tdata,tdata]=diffxT(datain,tempramp,conver,figans);%�̶��������ʣ���ͬ�¶�ʱ��ת���ʶ��¶ȵ�΢��ֵ
        dx_dT=log(dxdTdata);%z=ln(d��/dT)
        for i=1:length(conver)
            for j=1:group
                dxdT{i+1,j+1}=dx_dT(j,i)-log(1-dxdT{i+1,1});%z=ln(d��/dT(1-��))
            end
        end
        
        x=[];%�����������´ﵽ��Ӧת���ʵ��¶�
        for i=1:length(conver)
            for j=1:group
                x(i,j)=Tdata(j,i);%�����������´ﵽ��Ӧת����ʱ���¶�,��λk
                x(i,j)=1/x(i,j);%�¶�k����
            end
        end
        y=cell2mat(dxdT(2:end,2:end));
        fittingdata={};Y=[];Ycal=[];
        fittingdata{1,1}='other/��';fittingdata{2,1}='E';fittingdata{3,1}='lnA';fittingdata{4,1}='R^2';
        for i=1:group
            figname=strcat('�������ʦ�Ϊ',num2str(tempramp(i)));
            figure('Name',figname)
            X=[ones(length(conver),1) x(:,i)];%������Ϻ�������
            Y(:,i)=y(:,i);
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%���
            fittingdata{1,i+1}=tempramp(i);
            fittingdata{2,i+1}=-b(2)*8.314;%���
            fittingdata{3,i+1}=b(1);%������
            fittingdata{4,i+1}=stats(1,1);%�ع�ϵ��
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '1/T (k)' );ylabel( 'ln[d��/dt/(1-��)]' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            titlename=strcat('26.first-order (unimolecular) with ��=',num2str(tempramp(i)),'K/min');
            title(titlename);
        end
        fprintf('thermo_kinetic_fit��������26.first-order (unimolecular)���н��������E�Ͳ���lnA���ع�ϵ�����д�����fittingdata��\n')
        fprintf('Xֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��x��Y��Ycal��\n');
        clear avenum figname fitresult ft fun group i j reload cho bint r rint stats
           
    case 27
        tempramp=input('��������������������ʣ���λ���϶�/min���ÿպŸ����Ҵ�����:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%��������
        if length(tempramp)~=group
            error('\n��������������ʵ������������һ�£����飡����');
        end
        conver=input('������Ϊÿ������ָ������ת����,���ÿո�����������ţ�\n');
        figans=input('\n�Ƿ���ʾת���ʶ��¶�΢��ͼ��ת�������¶ȱ仯ͼ��y/n��������:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%ת���ʾ���
        fprintf('\nthermo_kinetic_fit��������27.Freeman-Carroll�������΢�ַ���������...');
        [dxdT,Tdata,~]=diffxT(datain,tempramp,conver,figans);
        dxdt={};dxdt{1,1}='��/��';
        Tdatacopy={};Tdatacopy{1,1}='��/��';
        for i=1:length(tempramp)%��ʼ����ͷ
            dxdt{1,i+1}=tempramp(i);
            Tdatacopy{1,i+1}=tempramp(i);
        end
        for i=1:size(dxdT,1)%��dxdT��ת���ʶ��¶�΢��ת��Ϊת���ʶ�ʱ��΢��
            for j=1:size(dxdT,2)
                dxdt{j+1,1}=conver(j);
                Tdatacopy{j+1,1}=conver(j);
                dxdt{j+1,i+1}=dxdT(i,j)*(tempramp(i)-273.15);%min^-1
                Tdatacopy{j+1,i+1}=1/Tdata(i,j);%K^-1
            end
        end
        deltdt=[];deltconv=[];delttemp=[];
        for i=1:length(tempramp)
            kk=1;
            for j=1:length(conver)
                for k=1:length(conver)
                    if k>j
                        deltdt(kk,i)=log(dxdt{k+1,i+1})-log(dxdt{j+1,i+1});%��ln(d��/dt)
                        deltconv(kk,i)=log(1-dxdt{k+1,1})-log(1-dxdt{j+1,1});%��ln(1-��)
                        delttemp(kk,i)=Tdatacopy{k+1,i+1}-Tdatacopy{j+1,i+1};%��(1/T)
                        kk=kk+1;
                    end
                end
            end
        end
        y=[];%��ln(d��/dt)/��ln(1-��)
        y=deltdt./deltconv;
        x=[];%��(1/T)/��ln(1-��)
        x=delttemp./deltconv;
        
        fittingdata={};Y=[];Ycal=[];
        fittingdata{1,1}='other/��';fittingdata{2,1}='E';fittingdata{3,1}='n';fittingdata{4,1}='R^2';
        for i=1:group
            figname=strcat('�������ʦ�Ϊ',num2str(tempramp(i)));
            figure('Name',figname)
            X=[ones(length(x),1) x(:,i)];%������Ϻ�������
            Y(:,i)=y(:,i);
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%���
            fittingdata{1,i+1}=tempramp(i);
            fittingdata{2,i+1}=-b(2)*8.314;%���
            fittingdata{3,i+1}=b(1);%������
            fittingdata{4,i+1}=stats(1,1);%�ع�ϵ��
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '��(1/T)/��ln(1-��)' );ylabel( '��ln(d��/dt)/��ln(1-��)' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            titlename=strcat('27.Freeman-Carroll�����΢�ַ��� with ��=',num2str(tempramp(i)),'K/min');
            title(titlename);
        end
        fprintf('thermo_kinetic_fit��������27.Freeman-Carroll�����΢�ַ������н��������E�ͷ�Ӧ����n���ع�ϵ�����д�����fittingdata��\n')
        fprintf('Xֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��x��Y��Ycal��\n');
        clear dxdT Tdata X figname titlename

    case 28
        tempramp=input('��������������������ʣ���λ���϶�/min���ÿպŸ����Ҵ�����:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%��������
        if length(tempramp)~=group
            error('\n��������������ʵ������������һ�£����飡����');
        end
        conver=input('������Ϊÿ������ָ������ת����,���ÿո�����������ţ�\n');
        figans=input('\n�Ƿ���ʾת���ʶ��¶�΢��ͼ��ת�������¶ȱ仯ͼ��y/n��������:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%ת���ʾ���
        fprintf('\nthermo_kinetic_fit��������28.Zivkovic��������...');
        [~,Tdata,tdata]=diffxT(datain,tempramp,conver,figans);
        y=[];%ln[ln(1/(1-��))/t]
        x=[];%1/T
        for i=1:length(conver)
            for j=1:length(tempramp)
                y(j,i)=log(log(1/(1-conver(i)))/tdata(j,i));
                x(j,i)=1/Tdata(j,i);
            end
            
        end
        
        fittingdata={};Y=[];Ycal=[];
        fittingdata{1,1}='other/��';fittingdata{2,1}='E';fittingdata{3,1}='lnA';fittingdata{4,1}='R^2';
        for i=1:length(conver)
            figname=strcat('ת���ʦ�Ϊ',num2str(conver(i)));
            figure('Name',figname)
            X=[ones(length(x),1) x(:,i)];%������Ϻ�������
            Y(:,i)=y(:,i);
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%���
            fittingdata{1,i+1}=conver(i);
            fittingdata{2,i+1}=-b(2)*8.314;%���
            fittingdata{3,i+1}=b(1);%������
            fittingdata{4,i+1}=stats(1,1);%�ع�ϵ��
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '1/T (K^-1)' );ylabel( 'ln[ln(1/(1-��))/t]' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            titlename=strcat('28.Zivkovic�� with ��=',num2str(conver(i)));
            title(titlename);
        end
        fprintf('thermo_kinetic_fit��������28.Zivkovic�����н��������E�ͳ���lnA���ع�ϵ�����д�����fittingdata��\n')
        fprintf('Xֵ��Yֵ�Լ���ϵ�ֵ�ֱ�洢��x��Y��Ycal��\n');
        clear Tdata tdata figname titlename
        
    otherwise
        disp('�Ƿ����룬����������ģ�ͱ�ţ�����');
end

