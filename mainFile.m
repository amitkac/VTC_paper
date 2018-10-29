
%% the program is for the paper
% https://arxiv.org/pdf/1809.00446.pdf

clc; clear all; close all


set(0,'DefaultFigureWindowStyle','docked')
n=1e6;
lambda=2; % common rate parameter
psi=2;
N=1; % Noise variance (AWGN)

bins=100;
p_su=4;% secondary user peak power
p_pu=4; % primary user peak power
lambda_pu=2;
 
lambda_1=lambda/p_su; % Scaled rate parameter for SU

% rayleigh channels
% for primary user
Pu_ch=exprnd(1/lambda_pu,[1,n]);

% for secondary users
x1=exprnd(1/lambda,[1,n]);
x2=exprnd(1/lambda,[1,n]);
x3=exprnd(1/lambda,[1,n]);

lambda_0=lambda_pu/p_pu; % Scaled rate parameter for PU


 
 %% 1 user case: p>q
 
 num=1;
 figure (11);
 subplot(1,2,1)
 Int1=N+min(psi,p_su*x1); 
ll1= histogram(Int1,bins,'normal','pdf');
% hold on;

 plot(N+ll1.BinWidth/2:ll1.BinWidth:N+(bins-1)*ll1.BinWidth+ll1.BinWidth/2,ll1.Values,':k','LineWidth',2)
 hold on;
 
x=0:0.1:10;
x=x+N; 

y_theo=NoiseIntPdf(N,num,lambda_1,psi,x);  
 plot(x,y_theo,'ko','LineWidth',2);     
hold on;

subplot(1,2,2)
 ecdf(Int1);
 hold on;
 y_cdf=NoiseIntCdf(N,num,lambda_1,psi,x);
 plot(x,y_cdf,'ko','LineWidth',2);
 hold on;

 figure(12)
 sinr1=Pu_ch.*p_pu./Int1;
 mm1=histogram(sinr1,bins,'normal','pdf');
 plot(mm1.BinWidth/2:mm1.BinWidth:(bins-1)*mm1.BinWidth+mm1.BinWidth/2,mm1.Values,':k','LineWidth',2)
hold on;

z=0:0.1:10;
y_snr=SINR(N,num,lambda_1,lambda_0,psi,z);
plot(z,y_snr,'ko','LineWidth',2);
hold on;

figure(13)
dx=(max(z)-min(z))/length(z);
yc=cumtrapz(y_snr).*dx;
plot(z,yc,'--k'); 
hold on;
Out=OutageProb(N,num,lambda_1,lambda_0,psi,z);
plot(z,Out,'ok','LineWidth',1)
hold on;

SimCap=log(1+sinr1);
figure(14)
cc1=histogram(SimCap,bins,'normal','pdf');
plot(cc1.BinWidth/2:cc1.BinWidth:(bins-1)*cc1.BinWidth+cc1.BinWidth/2,cc1.Values,'--k','LineWidth',1)
hold on;
Out=InstCap(N,num,lambda_1,lambda_0,psi,z);
plot(z,Out,'ok','LineWidth',1)
hold on;

%% 2 users  p>q
num=2;

figure(33)
Int2=N+min(psi,p_su*x1+p_su*x2); 
ll2= histogram(Int2,bins,'normal','pdf');

figure (11);
subplot(1,2,1)
plot(N+ll2.BinWidth/2:ll2.BinWidth:N+(bins-1)*ll2.BinWidth+ll2.BinWidth/2,ll2.Values,'-k','LineWidth',2)
hold on;
 
x=0:0.1:10;
x=x+N; 

y_theo2=NoiseIntPdf(N,num,lambda_1,psi,x);  
 plot(x,y_theo2,'k>','LineWidth',2);     
hold on;

subplot(1,2,2)
ecdf(Int2);
hold on;
y_cdf2=NoiseIntCdf(N,num,lambda_1,psi,x);
plot(x,y_cdf2,'k>','LineWidth',2);
hold on;
 
figure(55)
sinr2=Pu_ch.*p_pu./Int2;
mm2=histogram(sinr2,bins,'normal','pdf');

figure(12)
plot(mm2.BinWidth/2:mm2.BinWidth:(bins-1)*mm2.BinWidth+mm2.BinWidth/2,mm2.Values,':k','LineWidth',2)
hold on;
alpha=2;
beta=lambda_1;
z=0:0.1:10;
y_snr2=SINR(N,num,lambda_1,lambda_0,psi,z);
plot(z,y_snr2,'>k','LineWidth',2);
hold on;

figure(13)
dx=(max(z)-min(z))/length(z);
yc=cumtrapz(y_snr2).*dx;
plot(z,yc,'--k'); 
hold on;
Out=OutageProb(N,num,lambda_1,lambda_0,psi,z);
plot(z,Out,'>k','LineWidth',1)
hold on;
 
SimCap2=log(1+sinr2);
figure(222)
cc1=histogram(SimCap2,bins,'normal','pdf');
figure(14)

plot(cc1.BinWidth/2:cc1.BinWidth:(bins-1)*cc1.BinWidth+cc1.BinWidth/2,cc1.Values,'--k','LineWidth',1)
hold on;
Out=InstCap(N,num,lambda_1,lambda_0,psi,z);
plot(z,Out,'>k','LineWidth',1)
hold on;

%% 3 users  p>q
num=3;
figure(66)
Int3= N + min(psi,p_su*x1+p_su*x2+p_su*x3);
ll= histogram(Int3,bins,'normal','pdf');

figure (11);
subplot(1,2,1)
plot(N+ll.BinWidth/2:ll.BinWidth:N+(bins-1)*ll.BinWidth+ll.BinWidth/2,ll.Values,'-k','LineWidth',2)
hold on;
x=0:0.1:10;
x=x+N;
xlabel('x');
ylabel('f_{NI} (x)');
xlim([0,2.9]);
ylim([0,0.5]);

y_theo3=NoiseIntPdf(N,num,lambda_1,psi,x);
plot(x,y_theo3,'ks','LineWidth',2);     
hold on;

subplot(1,2,2)
ecdf(Int3);
hold on;
y_cdf3=NoiseIntCdf(N,num,lambda_1,psi,x);
plot(x,y_cdf3,'ks','LineWidth',2);
hold on;
xlabel('x');
ylabel('F_{NI} (x)');
xlim([0,4]);
ylim([0,1]);

figure(77)
sinr3=Pu_ch.*p_pu./Int3;
mm3=histogram(sinr3,bins,'normal','pdf');

figure(12)
plot(mm3.BinWidth/2:mm3.BinWidth:(bins-1)*mm3.BinWidth+mm3.BinWidth/2,mm3.Values,':k','LineWidth',2)
hold on;

y_snr3=SINR(N,num,lambda_1,lambda_0,psi,z);
plot(z,y_snr3,'ks','LineWidth',2);

figure(13)
dx=(max(z)-min(z))/length(z);
yc=cumtrapz(y_snr3).*dx;
plot(z,yc,'--k'); 
hold on;
Out3=OutageProb(N,num,lambda_1,lambda_0,psi,z);
plot(z,Out3,'<k','LineWidth',1)
 

SimCap3=log(1+sinr3);

figure(333)
cc1=histogram(SimCap3,bins,'normal','pdf');

figure(14)
plot(cc1.BinWidth/2:cc1.BinWidth:(bins-1)*cc1.BinWidth+cc1.BinWidth/2,cc1.Values,'--k','LineWidth',1)
hold on;
Out=InstCap(N,num,lambda_1,lambda_0,psi,z);
plot(z,Out,'<k','LineWidth',1)
hold on;
xlabel('z');
ylabel('f_{z} (z)');
xlim([0,2.5]);
ylim([0,1.4]);


close(333)
close(33)
close(55)
close(222)
close(66)
close(77)
close(12)
close(13)