%CASO DEL PÉNDULO (CLASE 2)
clear all; close all; clc

pkg load symbolic
pkg load control

syms fi fi_p fi_pp p p_p p_pp M m u long Fricc g;

%disp('Para el equilibrio inestable')
%ang_inic = 0;

%Ecuaciones
%p_pp = (1/(M+m))*(u-m*long*fi_pp+m*long*fi_p^2*fi-Fricc*p_p) %Angulos pequeños
%fi_pp=(1/long)*(g*fi-p_pp) %pequeños angulos

%fi_pp=solve(fi_pp==(1/long)*(g*fi-p_pp))
%disp('fi_pp=');%pretty(simplify(fi_pp));
%p_pp=subs(p_pp, 'fi_pp', fi_pp); %la función subs reemplaza simbolos
                                %por su valor 

%Matrices
%MATRIZ A
%Mat_A=[[0 1 0 0]
%[subs(subs(subs(subs(diff(p_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0), subs(subs(subs(subs(diff(p_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0), subs(subs(subs(subs(diff(p_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0), subs(subs(subs(subs(diff(p_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]
%[0 0 0 1]
%[subs(subs(subs(subs(diff(fi_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0),subs(subs(subs(subs(diff(fi_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0),subs(subs(subs(subs(diff(fi_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0),subs(subs(subs(subs(diff(fi_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]]

%MATRIZ B
%Mat_B=[0;subs(subs(subs(subs(diff(p_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0);0;subs(subs(subs(subs(diff(fi_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0)]
%pretty(simplify(Mat_A)) %Para mi no sirven estas lineas
%pretty(simplify(Mat_B))  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Para el equilibrio estable')
ang_inic = pi;

%Ecuaciones
p_pp = (1/(M+m))*(u+m*long*fi_pp-m*long*fi_p^2*fi-Fricc*p_p)
fi_pp=(1/long)*(-g*(fi)+p_pp) %Pequeños angulos para fi~pi sin(fi)~-fi, cos(fi)=-1
%fi_pp=solve(fi_pp==(1/long)*(-g*fi+p_pp),fi_pp);
%disp('fi_pp='); pretty(simplify(fi_pp));
%p_pp=subs(p_pp,'fi_pp',fi_pp);

%Matrices
Mat_A1=[[0 1 0 0]
[subs(subs(subs(subs(diff(p_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0), subs(subs(subs(subs(diff(p_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0), subs(subs(subs(subs(diff(p_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0),subs(subs(subs(subs(diff(p_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]
[0 0 0 1]
[subs(subs(subs(subs(diff(fi_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0),subs(subs(subs(subs(diff(fi_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0),subs(subs(subs(subs(diff(fi_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0),subs(subs(subs(subs(diff(fi_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]]

Mat_B1=[0;subs(subs(subs(subs(diff(p_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0);0;subs(subs(subs(subs(diff(fi_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0)]

%pretty(simplify(Mat_A))
%pretty(simplify(Mat_B))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%VERIFICACIÓN
%Para saber si los resultados obtenidos son válidos, se tiene que
%realizar una simulación. Para ello, simulamos el comportamiento
%del péndulo invertido en su equilibrio estable.

%(lo del equilibrio inestable lo comenté porque no se usaba)

%Valores para las variables
m=.1;
Fricc=0.1; 
long=0.6;
g=9.8;
M=.5;
h=0.0001;
tiempo=(10/h);
p_pp=0;
tita_pp=0; 

t=0:h:tiempo*h;
omega(1):0; %omega(1)0
p_p=0:h:tiempo*h; 
u=linspace(0,0,tiempo+1);

%condiciones iniciales
alfa(1)=pi-0.8; 
color='b';
p(1)=0; 
p_p(1)=0; 
u(1)=0; 
p(1)=0; 
i=1;

Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0]
Mat_B=[0; 1/M; 0; 1/(long*M)]
X0=[0 0 pi 0]';x=[0 0 alfa(1) 0]';

%%%%%% CICLO WHILE %%%%
while(i<(tiempo+1))
%Variables del sistema no lineal
estado=[p(i); p_p(i); alfa(i); omega(i)];
u(i)=0;
%Sistema no lineal
p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-
Fricc*p_p(i));
tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
p_p(i+1)=p_p(i)+h*p_pp;
p(i+1)=p(i)+h*p_p(i);
omega(i+1)=omega(i)+h*tita_pp;
alfa(i+1)=alfa(i)+h*omega(i);
%Variables del sistema lineal
pl(i)=x(1); p_pl(i)=x(2);alfal(i)=x(3);omegal(i)=x(4);
%Sistema lineal
xp=Mat_A*(x-X0)+Mat_B*u(i);
x=x+h*xp;
i=i+1;
end

pl(i)=x(1); p_pl(i)=x(2);alfal(i)=x(3);omegal(i)=x(4);
figure(1);hold on;subplot(3,2,1);plot(t,omega,color);grid on; title('Velocidad Ángulo');holdon;plot(t,omegal,'k');
subplot(3,2,2);plot(t,alfa,color);hold on;plot(t,pi*ones(size(t)),'k');plot(t,alfal,'k');grid on;title('Ángulo');hold on;
subplot(3,2,3); plot(t,p,color);grid on;title('Posición carro');hold on;plot(t,pl,'k');
subplot(3,2,4);plot(t,p_p,color);grid on;title('Velocidad carro');hold on;plot(t,p_pl,'k');
subplot(3,1,3);plot(t,u,color);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;
subplot(2,2,1);plot(alfal,omegal,'k');
subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posición carro');ylabel('Velocidad carro');hold on;
subplot(2,2,2);plot(pl,p_pl,'k');







