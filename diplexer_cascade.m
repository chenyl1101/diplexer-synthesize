ckt1=read(rfckt.passive,'junction0.s3p');
ckt2=read(rfckt.amplifier,'channel1.2.s2p');
% ckt3=read(rfckt.amplifier,'channel2.1.s2p');
freq=linspace(2.2e9,3e9,2401);
analyze(ckt1,freq);
analyze(ckt2,freq);
% analyze(ckt3,freq);
sparams_3p=ckt1.AnalyzedResult.S_Parameters;
sparams_2p_1=ckt2.AnalyzedResult.S_Parameters;
% sparams_2p_2=ckt3.AnalyzedResult.S_Parameters;

Kconn=1;
sparams_cascade_3p=cascadesparams(sparams_3p,sparams_2p_1,Kconn);%change

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=zeros(1,2401);
for k=1:2401
    Z(k)=1*(1+sparams_cascade_3p(2,2,k))./(1-sparams_cascade_3p(2,2,k));
end
clear k;
plot(freq,real(Z));hold on;
plot(freq,imag(Z));%plot Z in frequency band