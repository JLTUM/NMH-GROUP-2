figure(1)


U0 = 10;
Gamma = 1;
points = 5;
nn = points/2;
i = 1
subplot(2,2,i)
scheme = "Central"; % "Central" or "Upwind"
A_D_eq(U0,Gamma,points,scheme)

U0 = -10;
Gamma = 1;
points = 5;
nn = 6;
scheme = "Central"; % "Central" or "Upwind"
i = 2
subplot(2,2,i)
A_D_eq(U0,Gamma,points,scheme)

U0 = -10;
Gamma = 1;
points = 51;
nn = 6;
scheme = "Central"; % "Central" or "Upwind"
i = 3
subplot(2,2,i)
A_D_eq(U0,Gamma,points,scheme)

U0 = 10;
Gamma = 1;
points = 51;
nn = 6;
scheme = "Central"; % "Central" or "Upwind"
i = 4
subplot(2,2,i)
A_D_eq(U0,Gamma,points,scheme)
%----------Upwind-----%
figure(2)
U0 = 10;
Gamma = 1;
points = 5;
nn = points/2;
i = 1
subplot(2,2,i)
scheme = "Upwind"; % "Upwind" or "Upwind"
A_D_eq(U0,Gamma,points,scheme)

U0 = -10;
Gamma = 1;
points = 40;
nn = 6;
scheme = "Upwind"; % "Central" or "Upwind"
i = 2
subplot(2,2,i)
A_D_eq(U0,Gamma,points,scheme)

U0 = -10;
Gamma = 1;
points = 51;
nn = 6;
scheme = "Upwind"; % "Central" or "Upwind"
i = 3
subplot(2,2,i)
A_D_eq(U0,Gamma,points,scheme)

U0 = 10;
Gamma = 1;
points = 51;
nn = 6;
scheme = "Upwind"; % "Central" or "Upwind"
i = 4
subplot(2,2,i)
A_D_eq(U0,Gamma,points,scheme)