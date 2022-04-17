clear all; close all; clc

pkg load control
pkg load symbolic

%Variables
ti=7e-03; %Valor aproximado (necesito que sea entero) 
ts=1e-09; %Valor aproximado (necesito que sea entero)
muestras=ti/ts;
t=linspace(0,ti,muestras);
I=zeros(1,muestras);
Vc=zeros(1,muestras);
u=linspace(0,ti,muestras);
R=4.7e3; L=10e-6; C=100e-9; vin=12;

%Condiciones iniciales
I(1)=0; Vc(1)=0; u(1)=vin;

Mat_A=[-R/L  -1/L; 1/C  0];
Mat_B=[1/L  0];
Mat_C=[R  0];

Il(1)=0; Vcl(1)=0; 
x=[I(1) Vc(1)]' ; 
y(1)=0; 
Xop=[0 0]'; %punto de operacion
ii=0;

%Ciclo para generar los graficos
for i=1:muestras-1
  
  %Cambio de signo de u cada 1ms
  ii=ii+ts;
    if(ii>=1e-3)
      ii=0;
      vin=vin*-1;
    end
  %Variables de sistema lineal
  u(i)=vin;
  xp=Mat_A*(x-Xop)+Mat_B*u(i);
  x=x+xp*ts;
  Y=Mat_C*x;
  y(i+1)=Y(1);
  Il(i+1)=x(1);
  Vcl(i+1)=x(2);

  i %para poder ver por donde va la cuenta
end
figure(1)

%Grafica I, Vc y u
subplot(3,1,1);%hold on;
plot(t,Il,'b');title('Corriente , i_t');grid on;
subplot(3,1,2);%hold on;
plot(t,Vcl,'r');title('Tension Capacitor , Vc_t');grid on;
subplot(3,1,3);%hold on;
plot(t,u,'b');title(' Tension de Entrada, u_t');grid on;
