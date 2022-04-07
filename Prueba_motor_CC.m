pkg load control
close all; clear all;clc

%MOTOR DE CC (ejercicio pag 8)

s=tf('s');

%Valores
Laa = 366e-6;
J = 5e-9;
Ra = 55.6;
B = 0;
Ki = 6.49e-3;
Km = 6.53e-3;

%Funciones de transferencia
%FT1 = tf ([Ki], [Laa*J, (Ra*J+Laa*B), (Ra*B+Ki*Km)])
%FT = (Ki)/((s^2*Laa*J)+s*(Ra*J+Laa*B)+(Ra*B+Ki*Km));

%FT2 = (Ki/(Laa*J))/((s^2*Laa*J)/(Laa*J)+s*(Ra*J+Laa*B)/(Laa*J)+(Ra*B+Ki*Km)/(Laa*J))

%raices
%[~, d] = tfdata(FT2, 'v');
%roots(d)
%p1 = -1.5176e+05;
%p2 = -1.5260e+02;

%tiempo dinamica rapida
%tr = (log (0.95))/(p1)
%de tr tomamos el tercio de su valor
%tr2 = tr/3

%de tL tomamos el triple de su valor
%tL = (log (0.05))/(p2)
%tL2 = tL*3

%COMIENZA EL EJERCICIO
%Punto de operacion
X=-[0,0]; ii=0; t_etapa=1e-7; tF=0.001;

%CONSTANTES DEL PID
%Ejemplo 1:
Kp = 0.500;
Ki = 0.001;
Kd = 0.0001;
color_='r';

Ts = t_etapa;

A=(2*Kp*Ts+Ki*Ts^2+2*Kd)/(2*Ts)
B=(-2*Kp*Ts+Ki*Ts^2-4*Kd)/(2*Ts)
C=Kd/Ts

e=zeros(tF/t_etapa, 1); u=0; %no entiendo esta linea

function [X]=modmotor(t_etapa, xant, accion)
Va=accion;
h=1e-7;
omega= xant(1);
wp= xant(2);
for ii=1:t_etapa/h
  wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
  wp=wp+h*wpp;
  omega = omega + h*wp;
end
X=[omega,wp];

for t=0: t_etapa: tF
  ii=ii+1; k=ii+2;
  X=modmotor(t_etapa, X, u); %esta debe ser para modelar
  e(k)=wref-X(1); %ERROR
  u=u+A*e(k)+B*e(k-1)+C*e(k-2); %PID
  x1(ii)=X(1); %Omega ??
  x2(ii)=X(2); %wp ???
  acc(ii)=u; %  ???
end
t=0:e_etapa:tF;
subplot(2,1,1) ;hold on;
plot(t, x1, color_);title('salida y, \omega_t') %para mostrar la salida en radianes
subplot(2,1,2); hold on;
plot(t, acc, color_); title('Entrada u_t, v_a');
xlabel('Tiempo [seg.]');
