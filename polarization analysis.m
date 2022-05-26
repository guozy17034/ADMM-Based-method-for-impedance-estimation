clc
clear
A=importdata('shice1.txt');
B=importdata('shice2.txt');
a=A.data;
b=B.data;
% % [timetag,data,fR,lR,fs]=load_tsn('1405B.TS5');
% % [timetag1,data1,fR1,lR1,fs1]=load_tsn('1405B.TS5');
% % 
% % y=data;
% % r=data1;
% % m=size(data,1);
% % 
% % snr=norm(y,2)/norm(y-r,2);
% % snr=log10(snr);
% % snr=20*snr;
% % 
% % n=size(y-r,1);
% % mse=(norm(y-r,2)^2)/sqrt(n);
% % 
% % 
% % fm=sum(y.*r);
% % fz=norm(y,2)*norm(r,2);
% % ncc=fm/fz;





% % a=fs*60*500;
% % b=a+40*60*fs;
% % mb=data(a:b,:);
% % figure(1)
% % H=subplot(2,1,1);
% % PPP=get(H,'pos');      %第NN张子图的当前位置PPP是一个1×4的行向量
% % PPP(2)=PPP(2 )-0.03;      %高向上方延展0.03
% % set(H,'pos',PPP)        %根据新的边界设置。
% % plot(mb(1800:3600,1),'k');
%  set(gca,'xtick',[])
%  set(gca,'ytick',[])
% % subplot(2,1,2);
% % plot(mb(1800:3600,2),'k');
% % noi=zeros(1801,1);
%rectangle%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % noi(600:650)=300;
% % % % noi(660:710)=-300;
% % % % noi(720:770)=300;
% % % % noi(780:830)=-300;
% % % % noi(840:890)=300;
% % % % noi(900:950)=-300;
% % % % noi(960:1010)=300;
% % % % noi(1020:1070)=-300;
% % % % noi(1080:1130)=300;
% % % % noi(1140:1190)=-300;
%triangle

% % my=data(996165-1800:996165,:);
% % noi(600:610)=300;
% % noi(610:650)=300*rand(1,41);
% % t=1:10:510;
% % x=300*sawtooth(t);
% % noi(660:710)=x;
% %  tt=0:0.01:4.80;
% % y=chirp(tt,0,0.1,15);
% % noi(720:1200)=300*y;


% % % noi(780:830)=-300;
% % % noi(840:890)=300;
% % % noi(900:950)=-300;
% % % noi(960:1010)=300;
% % % noi(1020:1070)=-300;
% % % noi(1080:1130)=300;
% % % noi(1140:1190)=-300;
% % % % % % % % % % % % % % % % % % % % % %% jihuaxingzhi%%%%%%%%%%%%%%%%%%%
% % uni=mb(1800:3600,:);
% % uniq(:,1)=mb(1800:3600,1)+noi;
uni=b;
uniq(:,1)=b(:,1);
fs=15;
t=0:1/fs:120;
wavename='cmor3-3';
totalscal=1024;%视窗长度
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs=cwt(uni(:,1),scals,wavename); % 求连续小波系数
cco=abs(coefs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
imagesc(t,f,abs(coefs));
set(gca,'YDir','normal')
colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('小波时频图');
for i=1:size(uniq,1)
    fu(i)=uniq(i,1)+1i*uni(i,2);
end
coefu=cwt(fu,scals,wavename); % 求连续小波系数
 
fu1=1/2*(hilbert(real(fu))+1i*hilbert(imag(fu)));
fu2=fu-fu1;
coefu1=cwt(fu1,scals,wavename); % 求连续小波系数
coefu2=cwt(fu2,scals,wavename); % 求连续小波系数

fua=abs(coefu1)+abs(coefu2);
fub=abs(abs(coefu1)-abs(coefu2));
epsilon=fub./fua;
ep=abs(epsilon);
sita=1/2*atan(abs(coefu1.*coefu2));
faia=coefu1+conj(coefu2);
faib=coefu1-conj(coefu2);
fai=atan(faia)-atan(faib)+pi/2;
faii=abs(fai);
for i=1:size(coefs,1)/8
    choucf(i,:)=coefs(8*i,:);
end
for j=1:181
    chf(:,j)=choucf(:,1+10*(j-1));
end


for i=1:size(coefs,1)/8
    chouep(i,:)=ep(8*i,:);
end
for j=1:181
    che(:,j)=chouep(:,1+10*(j-1));
end



for i=1:size(coefs,1)/8
    choust(i,:)=sita(8*i,:);
end
for j=1:181
    cht(:,j)=choust(:,1+10*(j-1));
end



for i=1:size(coefs,1)/8
    choufa(i,:)=faii(8*i,:);
end
for j=1:181
    chi(:,j)=choufa(:,1+10*(j-1));
end


for i=1:128
    for j=1:181
        if cht(i,j)<0.0001
            cht(i,j)=(cht(i-1,j)+cht(i+1,j))/2;
        end
    end
end

for i=1:128
    for j=1:181
        if cht(i,j)<0.001
            cht(i,j)=cht(i-1,j)*rand;
        end
    end
end
for i=1:100
    for j=1:181
        if cht(i,j)<0.1
            cht(i,j)=cht(i-1,j)*rand;
        end
    end
end
cht(120,:)=(cht(121,:)+cht(119,:))/2;
chtt=cht;
% % for i=1:128
% %     for j=1:181
% %         if chtt(i,j)>0.777546
% %             chtt(i,j)=(chtt(i,j)-0.777546)*100;
% %         end
% %     end
% % end
% % for i=1:107
% %     for j=1:181
% %         if chtt(i,j)<0.2
% %             chtt(i,j)=(chtt(i-1,j)+chtt(i+1,j))/2;
% %         end
% %     end
% % end
% %    for j=1:181
% %         if chtt(108,j)>0.4
% %             chtt(108,j)=chtt(108,j)-0.2*rand;
% %         end
% %    end
% %     for j=1:181
% %         if chtt(108,j)>0.4
% %             chtt(108,j)=chtt(108,j)-0.2*rand;
% %         end
% %     end