N=1:1:5000;
fs=1000;
delt=0.001;
t=delt*N;
x=10*sin(10*t);
xf=abs((2/5000)*fft(x));
pf=angle(fft(x));
f=1:1:2500;
ans=[xf;pf];
subplot(3,1,1)
plot(f,xf(1:2500))
subplot(3,1,2)
plot(f,pf(1:2500))
dlmwrite('mohanty.txt',ans);
autopowef=abs(((2/5000)*fft(x)).^2)






















































subplot(3,1,3)
x