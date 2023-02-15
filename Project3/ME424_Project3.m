%% *ME424 Design Project-3*
% *Ediz Ferit Kula, 2017405087*
% 
% *Mehmet Kutay Yavuzkol, 2017405045*
% 
% _a = 7, b = 5._

clc, clear all, close all 
tic
%% Givens:

a = 7;
b = 5;
rpm = 200 + 3*b;            % frequency of applied force (rpm)
F_min = 100 + 12*a;         % minimum force (N)
F_max = 800 + 12*a;         % maximum force (N)
delta_L = 20 + 0.5*b;       % working range (mm)
M = 0.3 + 0.05*a;           % mass of the object (kg)
density = 7700;             % density of steels (kg/m^3)
G = 79;                     % shear modulus (GPa)
E = 200;                    % elastic modulus (GPa)

class_d = [0 0 0 0 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1 1.1 1.2 1.4 1.6 1.8 2 2.2 2.5 2.8 3 3.5 4 4.5 5 5.5 6 6.5 7 8 9 10 11 12 13 14;
    0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1 1.1 1.2 1.4 1.6 1.8 2 2.2 2.5 2.8 3 3.5 4 4.5 5 5.5 6 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1 1.1 1.2 1.4 1.6 1.8 2 2.2 2.5 2.8 3 3.5 4 4.5 5 5.5 6 6.5 7 8 9 10 11 12 13 14;
    0 0 0 0 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1 1.1 1.2 1.4 1.6 1.8 2 2.2 2.5 2.8 3 3.5 4 4.5 5 5.5 6 6.5 7 8 9 10 11 12 0 0;
    0 0 0 0 0 0 0 0 0 0.8 0.9 1 1.1 1.2 1.4 1.6 1.8 2 2.2 2.5 2.8 3 3.5 4 4.5 5 5.5 6 6.5 7 8 9 10 11 0 0 0]; 
% wire diameters (mm)
% 1st row: A227
% 2nd row: A228
% 3rd row: A229
% 4th row: A232
% 5th row: A401

cost_memory = 10000;          % assigning relatively high cost for Newton Raphson method    
%% Design Parameters:
% *Aim:* Minimum cost.

% d                 % wire diameter
% C                 % spring index
% material type     
%% Constraints:

% 4 < C < 12
% SF_yielding > 1
% SF_surge > 5
% SF_fatigue > 1.2
% SF_buckling > 1.5
%% Calculating Parameters Outside the Loop:

f = rpm/60*2*pi;                    % frequency of the applied force (Hz)
k = (F_max-F_min)/delta_L*1000;     % spring rate (N/m)
delta_cl = 0.1*delta_L;             % clash allowance (mm)
L_i = F_min*(10^3)/k;               % initial length (mm)
delta_max = delta_L + L_i;          % maximum deflection (mm)
F_m = (F_max+F_min)/2;              % mean force (N)
F_a = (F_max-F_min)/2;              % alternating force (N)
F_s = F_max + delta_cl*(10^(-3))*k; % solid force (N)
%% Loop starts:

for i=1:5
    for j=1:37
        if class_d(i,j) == 0
            continue;
        end 
        for C=4:0.01:12
%% Calculating Parameters Inside the Loop:

            D = C*class_d(i,j);                                                % coil diameter (mm)
            Na = (G*10^9)*(class_d(i,j)*10^(-3))/8/C^3/k;                      % number of active coils
            Nt = Na + 2;                                                       % total number of coils                                % solid length (mm)
            L_s = (Nt+1)*class_d(i,j);                                         % Solid length of spring (mm)
            m = pi*(class_d(i,j)*10^(-3))^2/4*Na*density*pi*(D*10^(-3));       % mass of the spring           
            switch i
                case 1
                    A=1753.3;
                    b=-0.1822;
                case 2
                    A=2153.5;
                    b=-0.1625;
                case 3
                    A=1831.2;
                    b=-0.1833;
                case 4
                    A=1909.9;
                    b=-0.1453;
                case 5
                    A=2059.2;
                    b=-0.0934;
            end  
            S_u = A*(class_d(i,j)^b);                % ultimate strength (MPa)
            S_sy = 0.45*S_u;             % shear yield strength w/o presetting ferrous (MPa)
            S_y = 0.75*S_u;                          % yield strength w/o presetting (MPa)            
%% Static Yielding Failure Analysis:

            K_s = 1 + 0.5/C;
            Tau_s = K_s*8*F_s*C/pi/(class_d(i,j)^2);      % max shear stress for static failure (MPa)
            SF_yielding = S_sy/Tau_s;                     % safety factor against static yielding
            if SF_yielding < 1
                continue;
            end                        
%% Fatigue Failure Analysis:

            K_w = (4*C-1)/(4*C-4) + 0.615/C;  
            Tau_a = K_w*8*F_a*C/pi/(class_d(i,j)^2);         % alternating shear stress (MPa)
            Tau_m = K_w*8*F_m*C/pi/(class_d(i,j)^2);         % mean shear stress (MPa)
            S_s = 0.375*S_u;                                      % from Figure 12.5, 10^7 life
            S_us = 0.8*S_u;                                  % ultimate shear strength (MPa)
            S_f = 0.5*S_us*S_s/(S_us-0.5*S_s);               % fatigue strength (MPa)  
            SF_fatigue = S_f*S_us/(S_us*Tau_a+S_f*Tau_m);    % safety factor against fatigue failure
            if SF_fatigue < 1.2 
                continue;
            end            
%% Buckling of Spring:

            L_f = L_s + delta_L + delta_cl + L_i;           % free length (mm)
            c_1 = E/2/(E-G);
            c_2 = 2*(pi^2)*(E-G)/(2*G+E);
            alfa=0.707;                                     % one end free the other constrained 
            lamda_effective = alfa*L_f/D;                   % effective lamda   
            delta_critical = L_f*c_1*(1-sqrt(1-c_2/(lamda_effective^2)));         % delta critical
            SF_buckling = delta_critical/delta_max;         % safety factor against buckling
            if SF_buckling < 1.5
                continue;
            end
            
%% Spring Surge:

            f = rpm/60;                                     % operating frequency (Hz)
            w = 2*pi*f;                                     % operating angular frequency (rad/s)
            x = 0.05;                                       % Newton Raphson Method
            x_old = 400;                                    % Initial values
            while abs(x_old-x) > 10^-6 && abs(x) > 10^-4
                x_old = x;
                x = x - (x - sqrt(k*m)/M*cot(x/sqrt(k/m)))/(1 + sqrt(k*m)/M/(sin(x/sqrt(k/m))^2));
            end
            w_n=abs(x);                                     % Natural angular frequency (rad/s)
            f_n = w_n/(2*pi);                               % Natural frequency (Hz)
            SF_surge = w_n/w;                               % Safety factor against spring surge
            if SF_surge < 5 
                continue;
            end           
%% Cost Analysis:

            switch i
                case 1
                    cost = 100*((class_d(i,j)*10^-1)^2)*Nt*(D*10^-1);
                case 2
                    cost = 200*((class_d(i,j)*10^-1)^2)*Nt*(D*10^-1);
                case 3
                    cost = 130*((class_d(i,j)*10^-1)^2)*Nt*(D*10^-1);
                case 4
                    cost = 250*((class_d(i,j)*10^-1)^2)*Nt*(D*10^-1);
                case 5
                    cost = 400*((class_d(i,j)*10^-1)^2)*Nt*(D*10^-1);
            end  
            if cost < cost_memory
                d_memory = class_d(i,j);
                C_memory = C;
                D_memory = D;
                Na_memory = Na;
                Nt_memory = Nt;
                L_s_memory = L_s;
                L_f_memory = L_f;
                K_s_memory = K_s;
                K_w_memory = K_w;
                slenderness_ratio_memory = L_f/D;
                deflection_ratio_memory = delta_max/L_f;
                SF_buckling_memory = SF_buckling;
                m_memory = m;
                f_n_memory = f_n; 
                SF_surge_memory = SF_surge;
                S_ut_memory = S_u;
                S_us_memory = S_us;
                S_y_memory = S_y;
                S_ys_memory = S_sy;
                Tau_s_memory = Tau_s;
                SF_yielding_memory = SF_yielding;
                Tau_a_memory = Tau_a;
                Tau_m_memory = Tau_m;
                S_s_memory = S_s;
                SF_fatigue_memory = SF_fatigue;
                S_f_memory = S_f;
                i_memory = i;
                cost_memory = cost;
            end           
        end
    end
end
toc
%% Output Parameters:

F_min
F_max
d_memory
C_memory
D_memory
delta_cl
k
Na_memory
Nt_memory
L_s_memory 
L_f_memory 
F_s
K_s_memory 
K_w_memory 
slenderness_ratio_memory 
deflection_ratio_memory 
SF_buckling_memory 
m_memory 
f_n_memory 
SF_surge_memory 
S_ut_memory 
S_us_memory 
S_y_memory 
S_ys_memory 
Tau_s_memory 
SF_yielding_memory 
Tau_a_memory 
Tau_m_memory 
S_s_memory 
SF_fatigue_memory 
% Extra output for material selection and cost
i_memory 
cost_memory
%% Printing Results:

file = fopen('output_DesignProject3_Ediz_Kutay.txt','w');
fprintf(file, "OUTPUT:\n");
fprintf(file, "\n");
fprintf(file, "Minimum force: "+F_min+" N\n");
fprintf(file, "Maximum force: "+F_max+" N\n");
fprintf(file, "Wire diameter: "+d_memory+" mm\n");
fprintf(file, "Spring index: "+C_memory+" \n");
fprintf(file, "Coil diameter: "+D_memory+" mm\n");
fprintf(file, "Clash allowance: "+delta_cl+" mm\n");
fprintf(file, "Spring rate: "+k+" N/m\n");
fprintf(file, "Number of Active Coils: "+Na_memory+" \n");
fprintf(file, "Total number of Coils: "+Nt_memory+" \n");
fprintf(file, "Solid length: "+L_s_memory+" mm\n");
fprintf(file, "Free length: "+L_f_memory+" mm\n");
fprintf(file, "The force that closes the spring: "+F_s+" N\n");
fprintf(file, "Ks: "+K_s_memory+" \n");
fprintf(file, "Kw: "+K_w_memory+" \n");
fprintf(file, "Slenderness ratio: "+slenderness_ratio_memory+" \n");
fprintf(file, "Deflection ratio: "+deflection_ratio_memory+" \n");
fprintf(file, "Safety factor against buckling : "+SF_buckling_memory+" \n");
fprintf(file, "Mass of Active Coils: "+m_memory+" kg\n");
fprintf(file, "Natural frequency: "+f_n_memory+" 1/s\n");
fprintf(file, "Safety against spring surge: "+SF_surge_memory+" \n");
fprintf(file, "Ultimate tensile strength: "+S_ut_memory+" MPa\n");
fprintf(file, "Ultimate shear strength: "+S_ys_memory+" MPa\n");
fprintf(file, "Yield strength: "+S_y_memory+" MPa\n");
fprintf(file, "Shear yield Strength: "+S_ys_memory+" MPa\n");
fprintf(file, "Stress in solid state: "+Tau_s_memory+" MPa\n");
fprintf(file, "Safety factor against yielding: "+SF_yielding_memory+" \n");
fprintf(file, "Alternating stress: "+Tau_a_memory+" MPa\n");
fprintf(file, "Mean stress: "+Tau_m_memory+" MPa\n");
fprintf(file, "Torsional stress: "+S_s_memory+" MPa\n");
fprintf(file, "Safety factor against fatigue failure: "+SF_fatigue_memory+" \n");
fprintf(file, "Material Selection: A227");
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%