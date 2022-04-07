%RLC TAREA 1

clear all; close all; clc

%definicion de variables

R = 4700;
L = 10e-6;
C = 100e-9;
vin = 12;
i = 1;
t_anterior = 1e-3; %cada 1 ms, u cambia de signo (+12,-12)


%MATRICES
mat_A = [-R/L  -1/L ; 1/C  0]
mat_B = [1/L  0]
mat_C = [R  0]
t=[]

%Obtengo los valores propios de la matriz A 
%para obtener los tiempos dinámicos y de simulación

val=eig(mat_A)

tr=log(0.95)/val(1); %tiempo dinámica rápida
tr_2 = tr/3 %tomo 3 veces menos su valor

ts=log(0.05)/val(2); %tiempo simulación
ts_2=ts*3 %tomo al menos 3 veces su valor

muestreo = ts_2/tr_2 %(nose)

%CONDICIONES INICIALES
vc(1)=0;  %tensión capacitor
il(1)=0;  %corriente
u(1)=12;  %entrada

X0=[0 0]' 
x=[0 0]'

while(i<=muestreo+1)
il(i)=x(1); vc(i)=x(2);
u(i)=vin;
t(i)=i*tr_2;
xp=mat_A*(X0-x)+ mat_B*u(i);
x=x+ts_2*xp;
if(t(i)>t_anterior)
vin=vin*(-1);
t_anterior=t_anterior+1e-3;
end
i=i+1;
end

subplot(3,1,1)
hold on;grid on;
plot(t,il,'r')
title('Corriente,t');
%subplot(3,1,2);hold on;grid on;
%plot(t,vc,'g');title('Vcap, t');
%subplot(3,1,3);hold on;grid on;
%plot(t,u,'b');title('Vin , t');






