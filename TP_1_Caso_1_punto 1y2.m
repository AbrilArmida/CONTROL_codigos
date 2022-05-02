clear all; close all; clc


% %PARTE 1 - AUTOVALORES Y TIEMPOS
% %Variables
% R = 5.6e3;
% L = 10e-6; 
% C = 100e-9; 
% vin = 12;
% 
% %Matrices
% Mat_A = [-R/L  -1/L; 1/C  0];
% Autovalores = eig(Mat_A);
% 
% %Tiempos
% t95 = log(0.95)/Autovalores(1); %dinamica rapida
% t5 = log(0.05)/Autovalores(2);  %dinamica lenta
% 
% ti = t95/0.1;                   %tiempo integracion
% ts = t5*5;                      %tiempo simulacion


%PARTE 2 - SIMULACIONES DE LA CORRIENTE Y TENSION EN EL CAPACITOR

%Tiempos y arreglos
T = 7e-3;                   %tiempo de simulacion
At = 1e-9;                  %tiempo de integracion (Delta t)
Kmax = T/At;                %pasos

t = linspace(0,T,Kmax);
I = zeros(1,Kmax);
Vc = zeros(1,Kmax);
u = linspace(0,T,Kmax);

%Variables
R = 4.7e3;                  %Caso 1
%R = 5.6e3;                 %Caso 2
L = 10e-6;
C = 100e-9;
vin = 12;

%Condiciones iniciales
I(1) = 0;                   %Corriente del circuito
Vc(1) = 0;                  %Tension en el capacitor
u(1) = vin;                 %Tension de entrada

%Matrices
A = [-R/L -1/L ; 1/C 0];
B = [1/L; 0];
E = [R 0];

%Punto de operacion
Il(1) = 0; 
Vcl(1) = 0;
x = [I(1) Vc(1)]';
y(1) = 0;
Xop = [0 0]'; 
ii = 0;

%Ciclo para graficar

for i=1:Kmax-1
    
    %Sistema real
    ii = ii+At;
    if(ii>=1e-3)            %Intercambio u de -12 a 12 cada 1ms
        ii = 0;
        vin = vin*-1;
    end
    
    u(i) = vin;
    
    %Ecuaciones del sistema real
    %Ip = (1/L)*u(i)-(1/L)*Vc(i)-(R/L)*I(i);
    %Vcp = 1/C * I(i);
    
    %I(i+1) = I(i)+Ip*At;
    %Vc(i+1) = Vc(i)+Vcp*At;
    
    %Variables del sistema lineal
    xp = A*(x-Xop)+B*u(i);
    x = x+xp*At;
    Y = E*x;
    y(i+1) = Y(1);
    Il(i+1) = x(1);
    Vcl(i+1) = x(2);
    i
end

%Graficas
figure(1)
subplot(3,1,1);%hold on;
plot(t,Il,'b');title('Corriente , i_t');grid on;
subplot(3,1,2);%hold on;
plot(t,Vcl,'r');title('Tension Capacitor , Vc_t');grid on;
subplot(3,1,3);%hold on;
plot(t,u,'b');title(' Tension de Entrada, u_t');grid on;
