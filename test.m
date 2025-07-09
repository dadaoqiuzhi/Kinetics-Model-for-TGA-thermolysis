[~,Tdata,tdata]=diffxT(datain,tempramp,conver);
y=[];%ln[ln(1/(1-¦Á))/t]
x=[];%1/T
for i=1:length(conver)
    for j=1:length(tempramp)
        y(j,i)=log(log(1/(1-conver(i)))/tdata(j,i));
        x(j,i)=1/Tdata(j,i);
    end
    
end