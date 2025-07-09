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
%热重法测定聚合物热降解反应动力学参数进展;聚丙烯-纳米碳酸钙复合材料的热降解动力学
disp('欢迎使用本程序--by 刘强@四川大学高分子科学与工程学院李光宪教授课题组，liubinqiang@163.com')
fprintf('本程序主要用于拟合聚合物热降解动力学模型，包括:\n**等温模型\n11.Conventional----lnt=ln[G(α)/A]+E/RT,α<=0.05')
fprintf('\n12.等温模型f(α) vs α，请参考文件夹中图片，确认13种常用方程 @reference 4');
fprintf('\n**非等温模型\n21.Coast-Redfern(CR)----ln[(-ln(1-α)/T^2]=ln(AR/βE)-E/RT');
fprintf('\n22.Friedman----ln(dα/dt)=ln(βdα/dT)=lnβ+ln(dα/dT)=ln(Af(α))-E/RT');
fprintf('\n23.Kissinger----ln(β/Tp^2)=ln(AR/E)+ln[n(1-xp)^(n-1)]-E/RTp；ln(β/Tp^2)=-E/RTp+ln(AR/E)');
fprintf('\n24.Flynn-Wall-Ozawa (FWO)----ln(β)=ln(AE/Rg(α))-5.331-1.052E/RT；ln(β)=C-1.052*E/(RT)');
fprintf('\n25.Kissinger-Akahira-Sunose (KAS)----ln(β/T^2)=ln(AE/Rg(α))-E/RT;ln(β/T^2)=C-E/RT');
fprintf('\n26.first-order (unimolecular)同一升温速率同一样品不同转化率时----ln(dα/dT/(1-α))=lnA-E/RT');
fprintf('\n27.Freeman-Carroll法（差减微分法）----△ln(dα/dt)/△ln(1-α)=n-E/R*△(1/T)/△ln(1-α),同一曲线');
fprintf('\n28.Zivkovic法----ln[ln(1/(1-α))/t]=lnA-E/RT\n，不同曲线');

fprintf('##请将实验数据按温度高低或者升温速率快慢输入input_data Excel文件中，每个数据五列，依次包含时间，温度，重量(mg)，重量百分比和重量对温度的一阶微分\n')
cho=input('\n请输入模型编号：\n');
if exist('datain','var')
    reload=input('\ndatain已存在，是否重新导入输入数据？y/n,带引号:\n');
    reload=lower(reload);
    if strcmp(reload,'y');
    datain=xlsread('input_data.xlsx');
    end
else
    datain=xlsread('input_data.xlsx');
end
[~,col]=size(datain);group=col/5;
if group ~= round(group)
    error('请检查输入数据，按要求提供input_data数据表\n');
end
switch cho
    case 11
        tempramp=input('请输入各组数据等温温度，单位摄氏度,请用空格隔开，带引号：\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);tempramp=str2double(tempramp);
        tempramp=tempramp+273.15;%温度为k的矩阵
        if length(tempramp)~=group
            error('\n等温温度组数与实验数据组数不一致，请检查！！！');
        end
        conver=input('请指定转化率,多组转化率请用空格隔开，带引号：\n');
        figans=input('\n是否显示转化率对温度微分图和转化率随温度变化图？y/n，带引号:\n');
        figans=lower(figans);
        disp('thermo_kinetic_fit程序运行11.Conversion中，请等待...')
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%转化率矩阵
        tempconver={};tempconver{1,1}='T/α';
        for i=1:length(conver)%初始化转化率表头
            tempconver{1,i+1}=conver(i);
        end
        [~,~,tdata]=diffxT(datain,tempramp,conver,figans);
        for i=1:group%初始化第一列温度
            tempconver{i+1,1}=tempramp(i);
            for j=1:length(conver)
                tempconver{i+1,j+1}=tdata(i,j);%各升温速率下达到相应转化率所用时间
            end
        end
        %拟合公式
        xydata=[];[rawtc,coltc]=size(tempconver);%温度倒数和时间对数矩阵
        for i=1:rawtc-1
            xydata(i,1)=1/tempconver{i+1,1};%x
            for j=1:coltc-1
                 xydata(i,j+1)=log(tempconver{i+1,j+1});%y
            end
        end
        X=[ones(rawtc-1,1) xydata(:,1)];%线性拟合含常数项
        figure( 'Name', '11.Conventional' );% Plot fit with data
        legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
        xlabel( '1/T (k)' );ylabel( 'lnt （min）' );% Label axes
        axis([min(xydata(:,1))*0.9 max(xydata(:,1))*1.2 min(xydata(:,2))*0.9 max(xydata(:,2))*1.2]);
        title('Lifetime prediction of polymer materials --11.Conversion');
        
        fittingdata={};Y=[];
        fittingdata{1,1}='Tiso/other';fittingdata{1,2}='E';fittingdata{1,3}='ln[G(α)/A]';fittingdata{1,4}='R^2';
        for i=2:coltc
            [b,bint,r,rint,stats]=regress(xydata(:,i),X);%拟合
            fittingdata{i,1}=conver(i-1);
            fittingdata{i,2}=b(2)*8.314;%活化能
            fittingdata{i,3}=b(1);%常数项
            fittingdata{i,4}=stats(1,1);%回归系数
            Y(:,i-1)=b(1)+b(2)*xydata(:,1);
            plot(xydata(:,1),xydata(:,i),'ok',xydata(:,1),Y(:,i-1),'-r','linewidth',3);
            hold on;
        end
        hold off;
        fprintf('\nthermo_kinetic_fit程序运行11.Conversion结束，拟合所得活化能,常数项和回归系数R^2存储在fittingdata中')
        fprintf('\nX值，Y值以及拟合的值分别存储在xydata，xydata和Y中,不用理会多余的转化率对温度微分和变化图');
        fprintf('\n11.Conversion方法获得的原始数据存储在xydata矩阵中，拟合的Y数据存储在Y矩阵中，注意区别不大时图上的点和线会重合\n')
        clear bint col coltc conver E group h i j r rawtc rint stats tempramp X cho group
    
    case 12
        fprintf('\n参见文件夹图片和文献4内容，暂时未写代码，敬请期待');
        break;
        
        
        
    case 21
        tempramp=input('请输入各组数据温升速率，单位摄氏度/min，用空号隔开且带引号:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%温升矩阵
        if length(tempramp)~=group
            error('\n温升速率组数与实验数据组数不一致，请检查！！！');
        end
        conver=input('请指定转化率,多组转化率请用空格隔开，带引号：\n');
        disp('thermo_kinetic_fit程序运行21.Coast-Redfern(CR)中，请等待...')
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%转化率矩阵
        tempconver={};tempconver{1,1}='β/α';
        for i=1:length(conver)%初始化转化率表头
            tempconver{1,i+1}=conver(i);
        end
        for i=1:group%初始化第一列温度
            tempconver{i+1,1}=tempramp(i);
            for j=1:length(conver)
                tempconver{i+1,j+1}=convertimetemp(5*i-3,j,datain,conver,'temp');%各升温速率下达到相应转化率时的温度,单位k
            end
        end
        z=[];x=[];%y值
        for i=1:length(tempramp)
            for j=1:length(conver)
                z(j,i)=log(-log(1-conver(j))/(tempconver{i+1,j+1}*tempconver{i+1,j+1}));
                x(j,i)=1/tempconver{i+1,j+1};
            end
        end
        
        output={};subplotnum=(length(conver)+mod(length(conver),2))/2;
        fitans=input('请选择非线性拟合方式编号,推荐2和3，1.fit常规模型拟合，2.ninfit拟合，3.lsqcurvefit拟合:\n');
        startpoint=input('请输入模型参数的初值矩阵,优先推荐编号3获取初始点，或解方程,很重要！历史估计值为[3e+16,4e+4]：\n');
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
            opts.StartPoint = startpoint;%适当更改起始点，接近最终答案
            Z=[];output={};output{1,1}='β';output{1,2}='A';output{1,3}='E';output{1,4}='R^2';
            for i=1:length(tempramp)
                [fitresult,gof] = fit([y(:,i),x(:,i)],z(:,i),ft,opts);
                output{i+1,1}=tempramp(i);
                output{i+1,2}=fitresult.a;output{i+1,3}=fitresult.e;
                fprintf('a为%f,e为%f',fitresult.a,fitresult.e);
                output{i+1,4}=gof.rsquare;
                opts.StartPoint = cell2mat(output(i+1,2:3));
                [fitresult,gof] = fit([y(:,i),x(:,i)],z(:,i),ft,opts);%再次拟合，减小误差
                Z(:,i)=log(fitresult.a./y(:,i))-fitresult.e*x(:,i);
                output{i+1,2}=fitresult.a*fitresult.e;output{i+1,3}=fitresult.e*8.314;output{i+1,4}=gof.rsquare;
                fprintf('升温速率为%f时各转化率得到活化能为%fJ/mol，R^2为%f\n',tempramp(i),output{i+1,3},output{i+1,4})
                h=figure(1);
                set(h,'Name','21.Coast-Redfern(CR)-nolinear');
                subplot(2,subplotnum,i)
                plot( x(:,i),z(:,i),'ok',x(:,i),Z(:,i),'*-r','linewidth',3);
                legend( 'origin data','fitting data', 'Location', 'NorthEast' );
                xlabel( '1/T (k)' );ylabel( 'ln[(-ln(1-α)/T^2]' );% Label axes
                %axis([min(x(:,i))*0.9 max(x(:,i))*1.2 max(z(:,i))*1.2 min(z(:,i))*0.9]);%注意z坐标均为负值
                titlename=strcat('21.Coast-Redfern(CR)',' with β=',num2str(tempramp(i)),'K/min');
                title(titlename);
            end
            fprintf('21.Coast-Redfern(CR)方法运行结束，活化能E和参数A，R^2按列储存在output中,x为温度倒数1/k，z为原始z值，Z为拟合z值');
            
            cho=input('\n上述非线性拟合如果很好，可直接对1/T进行线性拟合，是否进行y/n？,亦可作初始点，带引号:\n');
            cho=lower(cho);
            if strcmp(cho,'y')
                output2={};output2{1,1}='β';output2{1,2}='ln(AR/βE)';output2{1,3}='E';output2{1,4}='R^2';
                Y=[];
                for i=1:length(tempramp)
                    ft = fittype( 'poly1' );%也可采用polyfit
                    [fitresult, gof]=fit(x(:,i),z(:,i),ft);%拟合
                    output2{i+1,1}=tempramp(i);output2{i+1,2}=fitresult.p2;output2{i+1,3}=-8.314*fitresult.p1;output2{i+1,4}=gof.rsquare;
                    fprintf('转化率为%f时得到活化能为%fJ/mol，R^2为%f\n',tempramp(i),output2{i+1,3},output2{i+1,4})
                    Y(:,i)=fitresult.p1*x(:,i)+fitresult.p2;
                    h=figure(2);
                    set(h,'Name','21.Coast-Redfern(CR)-linear');
                    subplot(2,subplotnum,i)
                    plot( x(:,i),z(:,i),'ok',x(:,i),Y(:,i),'-r','linewidth',3);
                    legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
                    xlabel( '1/T (k)' );ylabel( 'ln[(-ln(1-α)/T^2]' );
                    titlename=strcat('21.Coast-Redfern(CR)',' with β=',num2str(tempramp(i)),'K/min');
                    title(titlename);
                end
                fprintf('21.Coast-Redfern(CR)方法运行结束，参数ln(AR/βE)和活化能E，R^2按列储存在output2中,x为温度倒数1/k，z为原始z值，Y为拟合z值');
            end
            
            
        elseif fitans==2
            Z=[];StartPoint=startpoint;beta={};%起始点很重要
            beta{1,1}='β/other';beta{2,1}='A';beta{3,1}='E';
            for j=1:size(x,1)
                alpha=strcat('r',num2str(j));
                beta{3+j,1}=alpha;
            end
            for i= 1:length(tempramp)
                X=x(:,i);%y为升温速率倒数，k;x为相应转化率温度倒数，1/k
                fz=@(beta,X)log(beta(1)/tempramp(i))-beta(2)*X(:,1);
                [b,r]=nlinfit(X,z(:,i),fz,StartPoint);
                StartPoint=b;
                [b,r]=nlinfit(X,z(:,i),fz,StartPoint);%再次拟合
                Z(:,i)=log(b(1)/tempramp(i))-b(2)*x(:,i);%拟合计算值
                beta{1,i+1}=tempramp(i);
                beta{2,i+1}=b(1)*b(2);%得到A
                beta{3,i+1}=b(2)*8.314;%得到E
                for j=1:length(r);%记录残差r
                    beta{3+j,i+1}=r(j);
                end
                h=figure(1);
                set(h,'Name','21.Coast-Redfern(CR)-nolinear');
                subplot(2,subplotnum,i)
                plot( x(:,i),z(:,i),'ok',x(:,i),Z(:,i),'*-r','linewidth',3);%本为三维空间，简化为二维
                legend( 'origin data','fitting data', 'Location', 'NorthEast' );
                xlabel( '1/T (k)' );ylabel( 'ln[(-ln(1-α)/T^2]' );% Label axes
                %axis([min(x(:,i))*0.9 max(x(:,i))*1.2 max(z(:,i))*1.2 min(z(:,i))*0.9]);%注意z坐标均为负值
                titlename=strcat('21.Coast-Redfern(CR)',' with β=',num2str(tempramp(i)),'K/min');
                title(titlename);  
            end
            fprintf('21.Coast-Redfern(CR)方法运行结束，参数A和活化能E，各数据点残差r按列储存在beta各行中\n')
       
        elseif fitans==3
            y=[];
            for i=1:size(x,1)
                for j=1:length(tempramp)
                    y(i,j)=tempramp(j);
                end
            end
            StartPoint=startpoint;beta=[];
            beta{1,1}='β/other';beta{2,1}='A';beta{3,1}='E';beta{4,1}='sum_r_res^2';
            for i= 1:length(tempramp)
                X=[y(:,i),x(:,i)];%y为升温速率倒数，k;x为相应转化率温度倒数，1/k
                fz=@(beta,X)log(beta(1)./X(:,1))-beta(2)*X(:,2);
                [b,resnorm]=lsqcurvefit(fz,StartPoint,X,z(:,i));
                StartPoint=b;
                [b,resnorm]=lsqcurvefit(fz,StartPoint,X,z(:,i));%再次拟合
                Z(:,i)=log(b(1)./y(:,1))-b(2)*x(:,i);%拟合计算值
                beta{1,i+1}=tempramp(i);
                beta{2,i+1}=b(1)*b(2);%得到A
                beta{3,i+1}=b(2)*8.314;%得到E
                beta{4,i+1}=resnorm(1);%残差平方和
                h=figure(1);
                set(h,'Name','21.Coast-Redfern(CR)-nolinear');
                subplot(2,subplotnum,i)
                plot( x(:,i),z(:,i),'ok',x(:,i),Z(:,i),'*-r','linewidth',3);%本为三维空间，简化为二维
                legend( 'origin data','fitting data', 'Location', 'NorthEast' );
                xlabel( '1/T (k)' );ylabel( 'ln[(-ln(1-α)/T^2]' );% Label axes
                %axis([min(x(:,i))*0.9 max(x(:,i))*1.2 max(z(:,i))*1.2 min(z(:,i))*0.9]);%注意z坐标均为负值
                titlename=strcat('21.Coast-Redfern(CR)',' with β=',num2str(conver(i)),'K/min');
                title(titlename); 
            end
            fprintf('21.Coast-Redfern(CR)方法运行结束，参数A和活化能E，残差r平方和按列储存在beta中\n');
            fprintf('\nX值，Y值以及拟合的值分别存储在x，z和Z中\n');
        end
        
        clear col conver fitresult ft fun gof h i j opts subplotnum tempramp titlename group reload cho fz fitans
         
    case 22
        disp('采用regress或者NonlinearLeastSquares方法对ln(dα/dT)=ln(Af(α))-E/RT-lnβ进行二元线性拟合')
        tempramp=input('请输入各组数据温升速率，单位摄氏度/min，用空号隔开且带引号:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%温升矩阵
        if length(tempramp)~=group
            error('\n温升速率组数与实验数据组数不一致，请检查！！！');
        end
        conver=input('请指定转化率,多组转化率请用空格隔开，带引号：\n');
        figans=input('\n是否显示转化率对温度微分图和转化率随温度变化图？y/n，带引号:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%转化率矩阵
        %avenum=input('请输入转化率对温度微分正整数平滑阈值，建议小于5,不平滑请输入n并带引号：\n');
        disp('thermo_kinetic_fit程序运行22.Friedman中，请等待...')
        dxdT={};dxdT{1,1}='β/α';dxdT(1,2:length(conver)+1)=num2cell(conver);
        for i=1:length(tempramp)
            dxdT{i+1,1}=tempramp(i);
        end
        [dx_dT,Tdata,~]=diffxT(datain,tempramp,conver,figans);%固定升温速率，不同温度时的转化率对温度的微分值
        for i=1:length(tempramp)
            for j=1:length(conver)
                dxdT{i+1,j+1}=dx_dT(i,j);
            end
        end
        z=log(dx_dT);%z=ln(dα/dT)
        x=[];%各升温速率下达到相应转化率的温度
        for i=1:group
            for j=1:length(conver)
                x(i,j)=Tdata(i,j);%各升温速率下达到相应转化率时的温度,单位k
                x(i,j)=1/x(i,j);%温度k倒数
            end
        end
        y=tempramp';
        y=log(y);
        
        fittingdata={};Y=[];Ycal=[];
        fittingdata{1,1}='other/α';fittingdata{2,1}='E';fittingdata{3,1}='ln(Af(α))';fittingdata{4,1}='R^2';
        for i=1:length(conver)
            figname=strcat('转化率α为',num2str(conver(i)));
            figure('Name',figname)
            X=[ones(length(tempramp),1) x(:,i)];%线性拟合含常数项
            Y(:,i)=z(:,i)+y(:,1);%方程变式为ln(dα/dT)+lnβ=ln(Af(α))-E/RT
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%拟合
            fittingdata{1,i+1}=conver(i);
            fittingdata{2,i+1}=-b(2)*8.314;%活化能
            fittingdata{3,i+1}=b(1);%常数项
            fittingdata{4,i+1}=stats(1,1);%回归系数
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'-r','linewidth',3);
            xlabel( '1/T (k)' );ylabel( 'ln(dα/dt)+lnβ' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            title('Lifetime prediction of polymer materials --22.Friedman');
        end
        fprintf('thermo_kinetic_fit程序运行22.Friedman运行结束，活化能E和参数ln(Af(α))，回归系数按列储存在fittingdata中\n');
        fprintf('\nX值，Y值以及拟合的值分别存储在x，Y和Ycal中\n');
        clear avenum control figname fitresult ft fun group i j reload cho bint r rint stats b col dx_dT 
        
    case 23
        tempramp=input('请输入各组数据温升速率，单位摄氏度/min，用空号隔开且带引号:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%温升矩阵
        if length(tempramp)~=group
            error('\n温升速率组数与实验数据组数不一致，请检查！！！');
        end
        disp('thermo_kinetic_fit程序运行23.Kissinger中，请等待...');
        Tp=[];
        for i=1:group%最大降解速率温度矩阵
            Tp(i)=Maxratetemp(i,datain);
        end
        Y=[];%ln(β/Tp^2)
        X=[];%1/Tp
        for i=1:group
            Y(i)=log((tempramp(i)-273.15)/Tp(i)/Tp(i));%β单位为℃/min,E算出来才对，但文献规定单位为k/min，但E有两个数量级差异
            X(i)=1/Tp(i);
        end
        
        Y=Y';X=X';Xx=[ones(size(X)) X];y=[];
        disp('\n采用regress线性拟合\n')%线性拟合，Y=(-E/R)*x，不含常数项
        [b,bint,r,rint,stats]=regress(Y(:,1),Xx(:,:));%拟合
        E=-8.314*b(2);
        y(:,1)=b(2)*X(:,1)+b(1);
        plot(X(:,1),Y(:,1),'ok',X(:,1),y(:,1),'*-r','linewidth',3);
        xlabel( '1/Tp (k^-1)' );ylabel( 'ln(β/Tp^2) (℃/K^2)' );% Label axes
        legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
        title('Lifetime prediction of polymer materials --23.Kissinger');
        fprintf('thermo_kinetic_fit程序运行23.Kissinger运行结束，活化能E为%.4fJ/mol，回归系数R^2为%.4f\n',E,stats(1,1));
        fprintf('\nX值，Y值以及拟合的值分别存储在X，Y和y中\n');
        %clear cho col group i X y Y bint r rint stats
    case 24
        tempramp=input('请输入各组数据温升速率，单位摄氏度/min，用空号隔开且带引号:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%温升矩阵
        if length(tempramp)~=group
            error('\n温升速率组数与实验数据组数不一致，请检查！！！');
        end
        conver=input('请指定转化率,多组转化率请用空格隔开，带引号：\n');
        figans=input('\n是否显示转化率对温度微分图和转化率随温度变化图？y/n，带引号:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%转化率矩阵
        disp('thermo_kinetic_fit程序运行24.Flynn-Wall-Ozawa (FWO)中，请等待...');
        x=[];%各升温速率下达到相应转化率的温度
        [~,Tdata,~]=diffxT(datain,tempramp,conver,figans);
        for i=1:group
            for j=1:length(conver)
                x(i,j)=Tdata(i,j);%各升温速率下达到相应转化率时的温度,单位k
                x(i,j)=1/x(i,j);%温度k倒数,1/T
            end
        end
        y=tempramp'-273.15;y=log(y);%lnβ,℃/min，E算出来才对，文献是这样做作的，但也有文献规定单位为k/min，但E有两个数量级差异
        fittingdata={};Y=[];Ycal=[];X=[];
        fittingdata{1,1}='other/α';fittingdata{2,1}='E';fittingdata{3,1}='C';fittingdata{4,1}='R^2';
        for i=1:length(conver)
            figname=strcat('转化率α为',num2str(conver(i)));
            figure('Name',figname)
            X=[ones(length(tempramp),1) x(:,i)];%线性拟合含常数项
            Y(:,i)=y(:,1);
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%拟合
            fittingdata{1,i+1}=conver(i);
            fittingdata{2,i+1}=-b(2)*8.314/1.052;%活化能
            fittingdata{3,i+1}=b(1);%常数项
            fittingdata{4,i+1}=stats(1,1);%回归系数
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '1/T (k)' );ylabel( 'lnβ(℃/min)' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            title('Lifetime prediction of polymer materials --24.Flynn-Wall-Ozawa (FWO)');
        end
        fprintf('thermo_kinetic_fit程序运行24.Flynn-Wall-Ozawa (FWO)运行结束，活化能E和参数C，回归系数按列储存在fittingdata中\n')
        fprintf('X值，Y值以及拟合的值分别存储在x，Y和Ycal中\n');
        clear figname fun group i j reload cho bint r rint stats
        
    case 25
        tempramp=input('请输入各组数据温升速率，单位摄氏度/min，用空号隔开且带引号:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%温升矩阵
        if length(tempramp)~=group
            error('\n温升速率组数与实验数据组数不一致，请检查！！！');
        end
        conver=input('请指定转化率,多组转化率请用空格隔开，带引号：\n');
        figans=input('\n是否显示转化率对温度微分图和转化率随温度变化图？y/n，带引号:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%转化率矩阵
        disp('thermo_kinetic_fit程序运行25.Kissinger-Akahira-Sunose (KAS)中，请等待...');
        x=[];%各升温速率下达到相应转化率的温度
        [~,Tdata,~]=diffxT(datain,tempramp,conver,figans);
        for i=1:group
            for j=1:length(conver)
                x(i,j)=Tdata(i,j);%各升温速率下达到相应转化率时的温度,单位k
                x(i,j)=1/x(i,j);%温度k倒数,1/T
            end
        end
        y=tempramp';
        fittingdata={};Y=[];Ycal=[];X=[];
        fittingdata{1,1}='other/α';fittingdata{2,1}='E';fittingdata{3,1}='C';fittingdata{4,1}='R^2';
        for i=1:length(conver)
            figname=strcat('转化率α为',num2str(conver(i)));
            figure('Name',figname)
            X=[ones(length(tempramp),1) x(:,i)];%线性拟合含常数项
            Y(:,i)=log((y(:,1)-273.15).*x(:,i).*x(:,i));%ln(β/T^2)
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%拟合
            fittingdata{1,i+1}=conver(i);
            fittingdata{2,i+1}=-b(2)*8.314;%活化能
            fittingdata{3,i+1}=b(1);%常数项
            fittingdata{4,i+1}=stats(1,1);%回归系数
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '1/T (k)' );ylabel( 'ln(β/T^2) (℃/K^2)' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            title('Lifetime prediction of polymer materials --25.Kissinger-Akahira-Sunose (KAS)');
        end
        fprintf('thermo_kinetic_fit程序运行25.Kissinger-Akahira-Sunose (KAS)运行结束，活化能E和参数C，回归系数按列储存在fittingdata中\n')
        fprintf('X值，Y值以及拟合的值分别存储在x，Y和Ycal中\n');
        %clear figname fun group i j reload cho bint r rint stats
        
    case 26
        tempramp=input('请输入各组数据温升速率，单位摄氏度/min，用空号隔开且带引号:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%温升矩阵
        if length(tempramp)~=group
            error('\n温升速率组数与实验数据组数不一致，请检查！！！');
        end
        conver=input('请至少为每组数据指定三组转化率,请用空格隔开，带引号：\n');
        figans=input('\n是否显示转化率对温度微分图和转化率随温度变化图？y/n，带引号:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%转化率矩阵
        disp('thermo_kinetic_fit程序运行26.first-order (unimolecular)中，请等待...')
        dxdT={};dxdT{1,1}='α/β';dxdT(1,2:length(tempramp)+1)=num2cell(tempramp);
        for i=1:length(conver)
            dxdT{i+1,1}=conver(i);
        end
        [dxdTdata,Tdata,tdata]=diffxT(datain,tempramp,conver,figans);%固定升温速率，不同温度时的转化率对温度的微分值
        dx_dT=log(dxdTdata);%z=ln(dα/dT)
        for i=1:length(conver)
            for j=1:group
                dxdT{i+1,j+1}=dx_dT(j,i)-log(1-dxdT{i+1,1});%z=ln(dα/dT(1-α))
            end
        end
        
        x=[];%各升温速率下达到相应转化率的温度
        for i=1:length(conver)
            for j=1:group
                x(i,j)=Tdata(j,i);%各升温速率下达到相应转化率时的温度,单位k
                x(i,j)=1/x(i,j);%温度k倒数
            end
        end
        y=cell2mat(dxdT(2:end,2:end));
        fittingdata={};Y=[];Ycal=[];
        fittingdata{1,1}='other/β';fittingdata{2,1}='E';fittingdata{3,1}='lnA';fittingdata{4,1}='R^2';
        for i=1:group
            figname=strcat('升温速率β为',num2str(tempramp(i)));
            figure('Name',figname)
            X=[ones(length(conver),1) x(:,i)];%线性拟合含常数项
            Y(:,i)=y(:,i);
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%拟合
            fittingdata{1,i+1}=tempramp(i);
            fittingdata{2,i+1}=-b(2)*8.314;%活化能
            fittingdata{3,i+1}=b(1);%常数项
            fittingdata{4,i+1}=stats(1,1);%回归系数
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '1/T (k)' );ylabel( 'ln[dα/dt/(1-α)]' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            titlename=strcat('26.first-order (unimolecular) with β=',num2str(tempramp(i)),'K/min');
            title(titlename);
        end
        fprintf('thermo_kinetic_fit程序运行26.first-order (unimolecular)运行结束，活化能E和参数lnA，回归系数按列储存在fittingdata中\n')
        fprintf('X值，Y值以及拟合的值分别存储在x，Y和Ycal中\n');
        clear avenum figname fitresult ft fun group i j reload cho bint r rint stats
           
    case 27
        tempramp=input('请输入各组数据温升速率，单位摄氏度/min，用空号隔开且带引号:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%温升矩阵
        if length(tempramp)~=group
            error('\n温升速率组数与实验数据组数不一致，请检查！！！');
        end
        conver=input('请至少为每组数据指定三组转化率,请用空格隔开，带引号：\n');
        figans=input('\n是否显示转化率对温度微分图和转化率随温度变化图？y/n，带引号:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%转化率矩阵
        fprintf('\nthermo_kinetic_fit程序运行27.Freeman-Carroll法（差减微分法）运行中...');
        [dxdT,Tdata,~]=diffxT(datain,tempramp,conver,figans);
        dxdt={};dxdt{1,1}='α/β';
        Tdatacopy={};Tdatacopy{1,1}='α/β';
        for i=1:length(tempramp)%初始化表头
            dxdt{1,i+1}=tempramp(i);
            Tdatacopy{1,i+1}=tempramp(i);
        end
        for i=1:size(dxdT,1)%将dxdT中转化率对温度微分转化为转化率对时间微分
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
                        deltdt(kk,i)=log(dxdt{k+1,i+1})-log(dxdt{j+1,i+1});%△ln(dα/dt)
                        deltconv(kk,i)=log(1-dxdt{k+1,1})-log(1-dxdt{j+1,1});%△ln(1-α)
                        delttemp(kk,i)=Tdatacopy{k+1,i+1}-Tdatacopy{j+1,i+1};%△(1/T)
                        kk=kk+1;
                    end
                end
            end
        end
        y=[];%△ln(dα/dt)/△ln(1-α)
        y=deltdt./deltconv;
        x=[];%△(1/T)/△ln(1-α)
        x=delttemp./deltconv;
        
        fittingdata={};Y=[];Ycal=[];
        fittingdata{1,1}='other/β';fittingdata{2,1}='E';fittingdata{3,1}='n';fittingdata{4,1}='R^2';
        for i=1:group
            figname=strcat('升温速率β为',num2str(tempramp(i)));
            figure('Name',figname)
            X=[ones(length(x),1) x(:,i)];%线性拟合含常数项
            Y(:,i)=y(:,i);
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%拟合
            fittingdata{1,i+1}=tempramp(i);
            fittingdata{2,i+1}=-b(2)*8.314;%活化能
            fittingdata{3,i+1}=b(1);%常数项
            fittingdata{4,i+1}=stats(1,1);%回归系数
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '△(1/T)/△ln(1-α)' );ylabel( '△ln(dα/dt)/△ln(1-α)' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            titlename=strcat('27.Freeman-Carroll（差减微分法） with β=',num2str(tempramp(i)),'K/min');
            title(titlename);
        end
        fprintf('thermo_kinetic_fit程序运行27.Freeman-Carroll（差减微分法）运行结束，活化能E和反应级数n，回归系数按列储存在fittingdata中\n')
        fprintf('X值，Y值以及拟合的值分别存储在x，Y和Ycal中\n');
        clear dxdT Tdata X figname titlename

    case 28
        tempramp=input('请输入各组数据温升速率，单位摄氏度/min，用空号隔开且带引号:\n');
        tempramp=strtrim(tempramp);tempramp=strsplit(tempramp);
        tempramp=str2double(tempramp);tempramp=tempramp+273.15;%温升矩阵
        if length(tempramp)~=group
            error('\n温升速率组数与实验数据组数不一致，请检查！！！');
        end
        conver=input('请至少为每组数据指定三组转化率,请用空格隔开，带引号：\n');
        figans=input('\n是否显示转化率对温度微分图和转化率随温度变化图？y/n，带引号:\n');
        figans=lower(figans);
        conver=strtrim(conver);conver=strsplit(conver);
        conver=str2double(conver);%转化率矩阵
        fprintf('\nthermo_kinetic_fit程序运行28.Zivkovic法运行中...');
        [~,Tdata,tdata]=diffxT(datain,tempramp,conver,figans);
        y=[];%ln[ln(1/(1-α))/t]
        x=[];%1/T
        for i=1:length(conver)
            for j=1:length(tempramp)
                y(j,i)=log(log(1/(1-conver(i)))/tdata(j,i));
                x(j,i)=1/Tdata(j,i);
            end
            
        end
        
        fittingdata={};Y=[];Ycal=[];
        fittingdata{1,1}='other/α';fittingdata{2,1}='E';fittingdata{3,1}='lnA';fittingdata{4,1}='R^2';
        for i=1:length(conver)
            figname=strcat('转化率α为',num2str(conver(i)));
            figure('Name',figname)
            X=[ones(length(x),1) x(:,i)];%线性拟合含常数项
            Y(:,i)=y(:,i);
            [b,bint,r,rint,stats]=regress(Y(:,i),X(:,:));%拟合
            fittingdata{1,i+1}=conver(i);
            fittingdata{2,i+1}=-b(2)*8.314;%活化能
            fittingdata{3,i+1}=b(1);%常数项
            fittingdata{4,i+1}=stats(1,1);%回归系数
            Ycal(:,i)=b(1)+b(2)*x(:,i);
            plot(x(:,i),Y(:,i),'ok',x(:,i),Ycal(:,i),'*-r','linewidth',3);
            xlabel( '1/T (K^-1)' );ylabel( 'ln[ln(1/(1-α))/t]' );% Label axes
            legend( 'origin data','fitting curve', 'Location', 'NorthEast' );
            titlename=strcat('28.Zivkovic法 with α=',num2str(conver(i)));
            title(titlename);
        end
        fprintf('thermo_kinetic_fit程序运行28.Zivkovic法运行结束，活化能E和常数lnA，回归系数按列储存在fittingdata中\n')
        fprintf('X值，Y值以及拟合的值分别存储在x，Y和Ycal中\n');
        clear Tdata tdata figname titlename
        
    otherwise
        disp('非法输入，请重新输入模型编号！！！');
end

