function [R_n] = NWV_muster(Q,b,m,I,k_st)%input here zum Testen: NWV_muster(20,10,0.5,0.4/100,30)
% Q = 20;
% b = 10;
% m = 0.5;
% I = 0.401/100;
% k_st = 30;

%Q = 1.95
%kst = 30;
%b = 0.5;
%I = 0.001

 y = 0.9;
 y_n = 1;
<<<<<<< HEAD
    while abs(y - y_n) > 5e-5 %Check if they are similar enough genauigkeit beliebig Ã¤ndern
        y = y_n;
        A = b*y+m*y^2; %Allgemeine FlÃ¤chen Formel
=======
    while abs(y - y_n) > 5e-5 %Check if they are similar enough genauigkeit beliebig ändern
        y = y_n;
        A = b*y+m*y^2; %Allgemeine Flächen Formel
>>>>>>> Nick
        U = b + 2*sqrt( (m*y)^2 + y^2 ); %Allgemeine U Formel
        v = k_st * sqrt(I) * ( A / U )^(2/3); %Manning-Strickler
        
        if b ~= 0 && m == 0
            y_n = Q/(v*b); %Rechteck
        elseif b == 0 && m ~= 0
            y_n = sqrt(Q/(v*m)); %Dreieck
        else
            y_n = Q/(v*(b+m*y_n)); %Trapez  mit A = (b+m*y)*y 
        end
    end
    Fr=v/sqrt(9.81*(m*y_n+b)*y_n/(b+2*m*y_n)); %Credit: Damaris, Felix
    H_n = v^2/(2*9.81) + y_n;
    R_n = [y_n,H_n,v,Fr];
    
    %y_n soll = 1.0674
    %H_n soll = 1.2287
    %v_n soll = 1.7788
    %Fr soll = 0.5635 
end

<<<<<<< HEAD
%Julianlenz@outlook.com fÃ¼r Fragen
=======
%Julianlenz@outlook.com für Fragen
>>>>>>> Nick
