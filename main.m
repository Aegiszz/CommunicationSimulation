clc;
clear all;
close all;

SNR=input('input SNR:  ')

%开始
tic;
I=imread('coins.png');

%原始图像显示
figure(1)
imshow(I,[]);

[m,n]=size(I);
J=reshape(I,m*n,1);
N=numel(I);
Pr=imhist(I)/N;
sym=0:255;
figure(2),stem(sym,Pr);
for i=1:256
    if(Pr(i)==0)
        sym(i)=256;
    end
end
Pr=Pr(find(Pr~=0));
H=sum(-Pr.*log2(Pr));
sym=sym(find(sym~=256));
r1=zeros(length(Pr),1);
for i=1:length(Pr)
    r1(i)=Pr(i)*8;
end
R1=sum(r1);
n1=H/R1;
 
%霍夫曼编码
dict=huffmandict(sym,Pr);
hcode=huffmanenco(J,dict);
figure(3)
stem(hcode,'linewidth',2) 
title('Huffman ');
axis([ 0 20 0 1.5]), grid on;

%压缩比
R=(N*8)/(1*length(hcode))


%(7,4)汉明码编码
Hancode=encode(hcode,7,4,'hamming');
figure(4)
stem(Hancode,'linewidth',2) 
title('hamming ');
axis([ 0 20 0 1.5]), grid on;


%QPSK
data=Hancode;
figure(5)
subplot(411),stem(data,'linewidth',2), grid on;%调制信号
title('  Information before Transmiting ');
axis([ 0 20 0 1.5]);




%单极性转双极性
data_NZR=[2*data-1;ones(mod(length(data),2),1)];
subplot(412),stem(data_NZR,'linewidth',2), grid on;
title('NZR ');
axis([ 0 20 -1.5 1.5]);


%I路信号
data_sp=reshape(data_NZR,2,length(data_NZR)/2);
subplot(413),stem(data_sp(1,:),'linewidth',2), grid on;
title('I ');
axis([ 0 20 -1.5 1.5]);

%Q路信号
subplot(414),stem(data_sp(2,:),'linewidth',2), grid on;
title('Q');
axis([ 0 20 -1.5 1.5]);


%载频为2KHz
f=2000;
T=1/f;
t=T/10:T/10:T;

%I、Q预先分配空间
y_I=zeros(1,length(data)/2*length(t));
y_Q=zeros(1,length(data)/2*length(t));

for i=1:length(data)/2
    y1=data_sp(1,i)*cos(2*pi*f*t); %I路信号
    y2=data_sp(2,i)*sin(2*pi*f*t); %Q路信号
    y_I((i-1)*10+1:i*10)=y1; 
    y_Q((i-1)*10+1:i*10)=y2;
end

Tx_sig=y_I+y_Q;%调制后信号
tt=T/10:T/10:(T*length(data))/2;

%传输时间
time=T*(length(data_NZR)/2)

%信道容量
C=2*f*log2(1+SNR)


figure(6)
subplot(3,1,1);
plot(tt,y_I,'linewidth',3), grid on;
axis([0,0.005,-1,1]);
title(' wave form for inphase component in QPSK modulation ');
xlabel('time(sec)');
ylabel(' amplitude(volt0');

subplot(3,1,2);
plot(tt,y_Q,'linewidth',3), grid on;
axis([0,0.005,-1,1]);
title(' wave form for Quadrature component in QPSK modulation ');
xlabel('time(sec)');
ylabel(' amplitude(volt0');


subplot(3,1,3);
plot(tt,Tx_sig,'r','linewidth',2), grid on;
axis([0,0.005,-2,2]);
title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
xlabel('time(sec)');
ylabel(' amplitude(volt0');

%升余弦滤波器设计
rolloff = 0.5;
span = 10;
sps = 10;
b = rcosdesign(rolloff, span, sps);
fvtool(b,'Analysis','impulse');
 
%发送滤波
Tx_RCC=upfirdn(Tx_sig, b, sps);%成型滤波

%信道噪声
Tx_noise=awgn(Tx_RCC,SNR ,'measured');

%接收滤波
RE_RCC=upfirdn(Tx_noise, b, 1 ,sps);%成型滤波
RE_RCC=RE_RCC(span+1:end-span);%消除群时延

figure(8)
subplot(411),plot(Tx_sig,'r','linewidth',2), grid on;
axis([0,100,-2,2]);
title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
xlabel('time(sec)');
ylabel(' amplitude(volt0');

subplot(412),plot(Tx_RCC,'r','linewidth',2), grid on;
axis([0,1000,-1,1]);
title('Tx_RCC');
xlabel('time(sec)');
ylabel(' amplitude(volt0');

subplot(413),plot(Tx_noise,'r','linewidth',2), grid on;
axis([0,1000,-2,2]);
title('Tx_noise');
xlabel('time(sec)');
ylabel(' amplitude(volt0');

subplot(414),plot(RE_RCC,'r','linewidth',2), grid on;
axis([0,100,-2,2]);
title('RE_RCC');
xlabel('time(sec)');
ylabel(' amplitude(volt0');






Rx_sig=RE_RCC; % 接收信号

%预先分配空间
Rx_data=zeros(1,length(data)*10);

%QPSK解调
for i=1:1:length(data)/2
   
    %I路
    Z_in=Rx_sig((i-1)*length(t)+1:i*length(t)).*cos(2*pi*f*t); 
    Z_in_intg=(trapz(t,Z_in))*(2/T);
    if(Z_in_intg>0) 
        Rx_in_data=1;
    else
       Rx_in_data=0; 
    end
    Rx_data(2*i-1)=Rx_in_data;
    
    %Q路
    Z_qd=Rx_sig((i-1)*length(t)+1:i*length(t)).*sin(2*pi*f*t);
    Z_qd_intg=(trapz(t,Z_qd))*(2/T);
        
    if (Z_qd_intg>0)
        Rx_qd_data=1;
    else
       Rx_qd_data=0; 
    end
 
    Rx_data(2*i)=Rx_qd_data; %并转串
end

%清除补零
Rx_data= Rx_data(1:length(data));


%QPSK解调end


%接收到的最终信号
figure(9)
stem(Rx_data,'linewidth',2) 
title('Information after Receiveing ');
axis([ 0 20 0 1.5]), grid on;





%(7,4)汉明码解码
Handecode=decode(Rx_data,7,4,'hamming');
figure(10)
stem(Handecode,'linewidth',2) 
title('Hamming ');
axis([ 0 20 0 1.5]), grid on;


%霍夫曼解码
hdeco=huffmandeco(Handecode,dict);
figure(11)
stem(hdeco,'linewidth',2) 
title('Huffman ');
axis([ 0 20 0 256]), grid on;


%修改错误图像元素个数
if length(hdeco)>=N
    hdeco=hdeco(1:N);
else  
    for i= (length(hdeco)+1):N
        hdeco(i)=0;
    end
end


I1=reshape(hdeco,m,n);

%显示原始图像和接收图像
figure(12);
subplot(211),imshow(I,[]);
title('original image')
subplot(212),imshow(I1,[]);
title('Received image')

%求误码率
P=symerr(I,I1)/N


toc;
%结束

