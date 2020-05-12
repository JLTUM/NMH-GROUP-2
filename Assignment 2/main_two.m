close all
global sp;
sp = 1;

U0 = 10;
Gamma = 1;
points = 5;
nn = points/2;
scheme = "Upwind"; % "Central" or "Upwind or "Both"
A_D_eq(U0,Gamma,points,nn,scheme)

U0 = -10;
Gamma = 1;
points = 5;
nn = points/2;
scheme = "Upwind"; % "Central" or "Upwind or "Both"
A_D_eq(U0,Gamma,points,nn,scheme)

U0 = 10;
Gamma = 1;
points = 51;
nn = points/2;
scheme = "Upwind"; % "Central" or "Upwind or "Both"
A_D_eq(U0,Gamma,points,nn,scheme)

U0 = -10;
Gamma = 1;
points = 51;
nn = points/2;
scheme = "Upwind"; % "Central" or "Upwind or "Both"
A_D_eq(U0,Gamma,points,nn,scheme)