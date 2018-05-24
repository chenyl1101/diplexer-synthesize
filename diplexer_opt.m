%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalization for complex load Z
f1=2.620;%change L:2.478-2.568 H:2.620-2.718
f2=2.718;%change 
f0=sqrt(f1*f2);
BW=f2-f1;
omega=(f0/BW).*((freq./1e9)./f0-f0./(freq./1e9));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Z sample in passband
n=0;
for k=1:2401
    if omega(k)<=1 && omega(k)>=-1
       n=n+1;
       Z_sam(n)=Z(k);
    end
end
clear k;
sample=linspace(-1,1,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Z sample in passband
% plot(omega,real(Z));hold on;
% plot(omega,imag(Z));hold on;
% plot(sample,real(Z_sam));hold on;
% plot(sample,imag(Z_sam));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1-Ps1 2-Ps2 3 Ps3
%step1 Ps2 TZ=2.167/-2.1409
TZ=[-2.1409];%change
Ps2=poly(i*TZ);
%step2 fitting Ps1
Ps1_data=1.*polyval(Ps2,i*sample)./(2.*sqrt(real(Z_sam)));
% plot(omega,real(polyval(Ps2,i*omega)));hold on;
% plot(omega,imag(polyval(Ps2,i*omega)));hold on;
% plot(sample,real(Ps1_data));hold on;
% plot(sample,imag(Ps1_data));hold on;
Ps1=polyfit(i*sample,Ps1_data,4);
% plot(sample,real(polyval(Ps1,i*sample)));hold on;
% plot(sample,imag(polyval(Ps1,i*sample)));
%step3 synthesize Fs1 and Es1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pw1=(Ps1(1)/i).*poly(roots(Ps1)./i);
Pw1=real(Pw1);
Fw1=zeros(1,6);
Fw1(1)=1;
%derivation of Cw function
ref=[-1 -0.6 -0.2 0.2 0.6 1];
for num=1:100
    %solve the linear equations
    G=[ref(1)^4 ref(1)^3 ref(1)^2 ref(1) 1 polyval(Pw1,ref(1)),
       ref(2)^4 ref(2)^3 ref(2)^2 ref(2) 1 -polyval(Pw1,ref(2)),   
       ref(3)^4 ref(3)^3 ref(3)^2 ref(3) 1 polyval(Pw1,ref(3)),
       ref(4)^4 ref(4)^3 ref(4)^2 ref(4) 1 -polyval(Pw1,ref(4)),
       ref(5)^4 ref(5)^3 ref(5)^2 ref(5) 1 polyval(Pw1,ref(5)),
       ref(6)^4 ref(6)^3 ref(6)^2 ref(6) 1 -polyval(Pw1,ref(6))];
    V=-[ref(1)^5 ref(2)^5 ref(3)^5 ref(4)^5 ref(5)^5 ref(6)^5];
    X=inv(G)*V';
    for k=1:5
        Fw1(k+1)=X(k);
    end
    k=0;
   
    [dFw1,dPw1]=polyder(Fw1,Pw1);
    new=roots(dFw1);
    %choose inner passband extreme point
    for k=1:length(new)
        if real(new(k))>-1 && real(new(k))<1
           new(k)=new(k);
        else 
           new(k)=-2;
        end
    end
    k=0;
    %order those point
    new=real(new);%
    new=sort(new);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part is related to num of TZs
    for k=1:length(new)-3
        if new(k+3)==-2
           ref(k)=ref(k);
        else ref(k)=new(k+3);
        end
    end
    k=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
end
clear k;
% Cw after interation
Cw=poly2sym(Fw1)/(poly2sym(Pw1)*X(6));
% fplot(Cw,[-1,1],'r');hold on;
% line([-1,1],[0,0],'linestyle','--');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1=abs(polyval(Pw1,[-1])/polyval(Fw1,[-1]))/sqrt(10^2.2-1);
Fs1=poly(roots(Fw1)*i);
Ps1=Pw1(1)*poly(roots(Pw1)*i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rEs1=roots([0 Pw1]./E1-i.*Fw1);
% clear k;
% for k=1:length(rEs1)
%     if imag(rEs1(k))<0
%        rEs1(k)=conj(rEs1(k));
%     else
%        rEs1(k)=rEs1(k);
%     end
% end
% Ew1=poly(rEs1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs1_conj=[-conj(Fs1(1)),conj(Fs1(2)),-conj(Fs1(3)),conj(Fs1(4)),-conj(Fs1(5)),conj(Fs1(6))];
Ps1_conj=[conj(Ps1(1)),-conj(Ps1(2)),conj(Ps1(3)),-conj(Ps1(4)),conj(Ps1(5))];
sum1=[zeros(1,2) conv(Ps1,Ps1_conj)./E1^2]+conv(Fs1,Fs1_conj);
rsum1=roots(sum1);
rEs1=zeros(1,5);
j=1;
for k=1:length(rsum1)
    if real(rsum1(k))<0
       rEs1(j)=rsum1(k);
       j=j+1;
%     else
%        rEs1(j)=rEs1(j);
    end
end
clear j;clear k;
Es1=sqrt(-sum1(1))*poly(rEs1);%conservation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N-nfz is odd
% plot(freq,10.*log(abs(polyval(Fs1,i.*omega)./polyval(Es1,i.*omega))));
% hold on;
% plot(freq,10.*log(abs(polyval(Ps1,i.*omega)./polyval(Es1,i.*omega))./E1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step4
Fs2_data=(1+Z_sam).*polyval(Fs1,i*sample)-(1-conj(Z_sam)).*polyval(Es1,i*sample);
% figure(1);
% plot(sample,real(Fs2_data),'r');hold on;
% plot(sample,imag(Fs2_data),'b');hold on;
Fs2=polyfit(i*sample,Fs2_data,5);
% plot(sample,real(polyval(Fs2,i*sample)),'y');hold on;
% plot(sample,imag(polyval(Fs2,i*sample)),'g');
E2=E1;
%step 5
%power conservation for Es2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs2_conj=[-conj(Fs2(1)),conj(Fs2(2)),-conj(Fs2(3)),conj(Fs2(4)),-conj(Fs2(5)),conj(Fs2(6))];
sum2=[zeros(1,8) conv(Ps2,-Ps2)./E2^2]+conv(Fs2,Fs2_conj);
rsum2=roots(sum2);
rEs2=zeros(1,5);
j=1;
for k=1:length(rsum2)
    if real(rsum2(k))<0
       rEs2(j)=rsum2(k);
       j=j+1;
%     else
%        rEs2(j)=rEs2(j);
    end
end
clear j;clear k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%conservativeness condition for Es2
Es2_data=(1+conj(Z_sam)).*polyval(Es1,i*sample)-(1-Z_sam).*polyval(Fs1,i*sample);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%phase E opt
error=zeros(1,501);
j=1;
for theta_E=0*pi:0.004*pi:2*pi %501
    Es2=sqrt(-sum2(1))*exp(i*theta_E)*poly(rEs2);
    for k=1:n
        error(j)=error(j)+(abs(polyval(Es2,i*sample(k))-Es2_data(k)))^2;
    end
    error(j)=sqrt(error(j)/n);
    j=j+1;
end
clear k;clear j;
plot(linspace(0*pi,2*pi,501),error);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%phase E opt
theta_E=6.233;
Es2=sqrt(-sum2(1))*exp(i*theta_E)*poly(rEs2);%with phase optimized
% figure(3);
% plot(sample,real(Es2_data),'r');hold on;
% plot(sample,imag(Es2_data),'b');hold on;
% plot(sample,real(polyval(Es2,i*sample)),'y');hold on;
% plot(sample,imag(polyval(Es2,i*sample)),'g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ps2=i*Ps2;
% figure(4);
% plot(freq,10.*log(abs(polyval(Fs2,i.*omega)./polyval(Es2,i.*omega))));
% hold on;
% plot(freq,10.*log(abs(polyval(Ps2,i.*omega)./polyval(Es2,i.*omega))./E2));
%step 6 transmission line for Ps3 Fs3 Es3
%theta_P,theta_F,theta_E
theta_F=atan(imag(Fs2(1)/real(Fs2(1))));
Fs3=Fs2/(abs(Fs2(1))*exp(i*theta_F));%Fs3 is monic now
Es3=Es2/(abs(Es2(1))*exp(i*theta_E));%Es3 is monic now
theta_L=(theta_E-theta_F)/2;
theta_P=(theta_E+theta_F)/2;
Ps3=Ps2;
E3=E2*abs(Fs2(1));
F22s3=-[-conj(Fs3(1)),conj(Fs3(2)),-conj(Fs3(3)),conj(Fs3(4)),-conj(Fs3(5)),conj(Fs3(6))];
s22_2=exp(-i*2*theta_L)*polyval(F22s3,i.*omega)./polyval(Es3,i.*omega);
s21_2=exp(-i*theta_L)*polyval(Ps3,i.*omega)./polyval(Es3,i.*omega)./E3;
% write updated channel S-parameter
s11_2=exp(-i*2*theta_L)*polyval(Fs3,i.*omega)./polyval(Es3,i.*omega);
s12_2=exp(-i*theta_L)*polyval(Ps3,i.*omega)./polyval(Es3,i.*omega)./E3;
figure(5);
plot(freq,10.*log(abs(s22_2)),'r');hold on;
plot(freq,10.*log(abs(s21_2)),'b');
sparams=[freq'/1e9 real(s11_2)' imag(s11_2)' real(s12_2)' imag(s12_2)' real(s21_2)' imag(s21_2)' real(s22_2)' imag(s22_2)'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear n;