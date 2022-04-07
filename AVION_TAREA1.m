%AVION TAREA

close all; clear all: clc

pkg load control
pkg load symbolic

%1) Obtener el sistema lineal en variables de estado
%% para el equilibrio x=[0 0 0 0]

syms alpha alpha_p fi fi_p fi_pp omega a b c h h_p u

%VALORES
omega = 2;
a=0.05;
b=5;
c=100;

%ECUACIONES

%alpha_p=a*(fi-alpha);
%fi_pp=-omega^2*(fi-alpha-b*u);
%h_p=c*alpha;

%MATRICES
mat_A=[(-alpha) alpha 0 0; 0 0 1 0; (omega^2) (-omega^2) 0 0; 100 0 0 0]
%mat_B=[0 0 (omega^2*b) 0 0]