clear all; close all; clc

pkg load control
pkg load symbolic

%TP1 - PUCHETA

%Variables
R=4.7e3; L=10e-6; C=100e-9; Vin=-12; Vant=0; 
i=1; 

%Matrices
Mat_A=[-R/L  -1/L; 1/C  0]
Autovalores = eig(Mat_A)
Mat_B=[1/L  0]
Mat_C=[R  0]

%Tiempos
t95=log(0.95)/Autovalores(1)
t5=log(0.05)/Autovalores(2)

ti=t95/10 %su h
ts=t5*3 %su tf

t=[]; 
pasos=ts/ti  
aux=0.001/ti; %VER TODA

%Condiciones iniciales
vc(1)=0; il(1)=0; u(1)=12;
tant=1e-3;

X0=[0 0]'; x=[0 0]';
Vin=Vin*(-1);

while(i<=10000000+1)

  %Variables del sistema lineal
  il(i)=x(1); vc(i)=x(2);
  u(i)=Vin;
  t(i)=i*ti;

  %sistema lineal
  xp=Mat_A*(x-X0) + Mat_B*u(i);
  x=x+ti*xp;

  if(t(i)>tant) %hago variar de 12 a -12 cada 1ms
  Vin=Vin*(-1);
  tant=tant+1e-3;
  end 
i=i+1
end

subplot(3,1,1);hold on;grid on; %Grafica i, Vc y la entrada
plot(t,il,'r');title('IL , t');
subplot(3,1,2);hold on;grid on;
plot(t,vc,'r');title('Vc, t');
subplot(3,1,3);hold on;grid on;
plot(t,u,'r');title('U , t');