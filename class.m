clear all
n=1:1:10000;
delt=0.001;
t=n*delt;
f=sin(2*pi*1.0.*t);
x=10.0*sin(2*pi*100*delt.*f.*t);
plot(t,x)