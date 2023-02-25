%% *ME424 Design Project-2*
% *Ediz Ferit Kula, 2017405087*
% 
% *Mehmet Kutay Yavuzkol, 2017405045*
% 
% _a = 7, b = 5._

clc, clear all, close all 
tic
%% Givens:

a=7;
b=5;
E=207000;       % elastic modulus (MPa)
v=0.3;          % poisson ratio
t=4+0.3*a;      % cylinder thickness (mm)
D=120-4*b;      % cylinder outer diameter (mm)
L=240+12*a;     % length of the cylinder (mm)
h=20;           % thickness of the end plates (mm)
L_bolt=L+2*h;   % length of the bolt (mm)
p0=12+0.5*b;    % maximum fluid pressure (MPa)
life = 5*10^5;  % life
pm = '';        % property class

class_d=[5, 6, 7, 8, 10, 12, 14, 16, 20, 22, 24, 27, 30, 33, 36, 0, 0, 0, 0, 0;
          2,2.5,3,3.5, 4, 5, 6, 7, 8, 10, 12, 14, 16, 0, 0, 0, 0, 0, 0, 0;
          5,6,7,8,10,12,14,16,20,22,24,0,0,0,0, 0, 0, 0, 0, 0;
          16,20,22,24,27,30,33,36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
          2,2.5,3,3.5,4,5,6,7,8,10,12,14,16, 0, 0, 0, 0, 0, 0, 0;
          5,6,7,8,10,12,14,16,20,22,24,27,30,33,36, 0, 0, 0, 0, 0;
          2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 22, 24, 27, 30, 33, 36];

class_pitch=[0.8, 1, 1, 1.25, 1.5, 1.75, 2, 2, 2.5, 2.5, 3, 3, 3.5, 3.5, 4, 0, 0, 0, 0, 0;
                0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1, 1, 1.25, 1.5, 1.75, 2, 2, 0, 0, 0, 0, 0, 0, 0;
                0.8, 1, 1, 1.25, 1.5, 1.75, 2, 2, 2.5, 2.5, 3, 0, 0, 0, 0, 0, 0,0, 0, 0;
                2,2.5,2.5,3,3,3.5,3.5,4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                0.40,0.45,0.5,0.6,0.7,0.8,1,1,1.25,1.5,1.75,2,2, 0, 0, 0, 0, 0, 0, 0;
                0.8,1,1,1.25,1.5,1.75,2,2,2.5,2.5,3,3,3.5,3.5,4, 0, 0, 0, 0, 0;
                0.40,0.45,0.5,0.6,0.7,0.8,1,1,1.25,1.5,1.75,2,2,2.5,2.5,3,3,3.5,3.5,4]; 

class_S=[225, 240, 400;
            310, 340, 420;
            380, 420, 520;
            600, 660, 830;
            650, 720, 900;
            830, 940, 1040;
            970, 1100, 1220];
%% Design Parameters:
% *Aim:* minimum cost.

%d               % major diameter of bolt (mm)
%Property_class  % property class of material
%K_i             % preload factor
%% Constraints:

cost_memory = 1000;
%SF_yield > 1.6
%SF_leak > 1.6
%SF_fatigue > 1.6
%% Loop starts:

for i=1:7
    for j=1:20
        if class_d(i,j) == 0
            continue;
        end 
        for K_i=0.400:0.001:0.900
%% Parameter Assigning:

d = class_d(i,j);            % major diameter
pitch = class_pitch(i,j);    % pitch
S_p = class_S(i,1);          % proof stress
S_y = class_S(i,2);          % yield stress
S_u = class_S(i,3);          % ultimate stress
%% Calculated Parameters:   

d_p = d-3*sqrt(3)*pitch/8;      % pitch diameter (mm)
d_r = d-17*sqrt(3)*pitch/24;    % root diameter (mm)
A_t = pi/4*((d_p+d_r)/2)^2;     % tensile stress area (mm^2)
A_c = pi/4*(D^2-(D-2*t)^2);     % clamped area (mm^2)
A_b = pi/4*d^2;                 % bolt shank area (mm^2)
A_s = pi/4*((D-2*t)^2-d^2);     % surface area (mm^2)

F_i = K_i*A_t*S_p;                                  % preload (N)
n = 1/pitch/(1+F_i/E/A_b)*(F_i*L/E/A_c +...
    F_i*L_bolt/E/A_b);                         % number of turns of the nut


%==============================================================================%
delta_i = F_i*L/E/A_c;  % displacement of the cylinder after preload
r_c = (D-2*t)/2;        % cylinder inner radius

F_e = p0*A_s;                                        % maximum external seperating force (N)
F_b = ((F_i+F_e)*(L-delta_i)/(E*A_c) - (v*p0)*L/E*(1+r_c/t) + (F_i)*(L_bolt-delta_i)/(E*A_b)) /...
    ((L_bolt-delta_i)/(E*A_b) + (L-delta_i)/(E*A_c));   % maximum tension in the bolt (N)
F_c = F_b - F_e;                                        % minimum compression in the clamped member (N)

sigma_bi = F_i/A_t;     % tensile stress in the bolt after preload (MPa) 
sigma_ci = F_i/A_c;     % compressive stress in the clamped member after preload (MPa)

sigma_bf = F_b/A_t;     % axial stress in the bolt
sigma_cf = F_c/A_c;     % axial stress in the clamped member

delta_sigma_b = sigma_bf - sigma_bi;    % the increase in the axial stress of the bolt due to pressure 
delta_sigma_c = sigma_cf - sigma_ci;    % the decrease in the axial stress of the cylinder due to pressure 

%% Static Yielding Failure Analysis: 

p0_yield_sol = (S_y*A_t*((L_bolt-delta_i)/A_b+(L-delta_i)/A_c) - F_i*(L_bolt-delta_i)/A_b - F_i*(L-delta_i)/A_c) /...
    (A_s*(L-delta_i)/A_c-v*(1+r_c/t)*L);

SF_yield = p0_yield_sol/p0;       % safety factor against yielding
if SF_yield <= 1.6
    continue;
end
%% Leakage/Seperation Failure Analysis:

p0_sep_sol = (F_i*(L-delta_i)/A_c + F_i*(L_bolt-delta_i)/A_b) /...
    ((A_s*((L_bolt-delta_i)/A_b+(L-delta_i)/A_c) - A_s*(L-delta_i)/A_c + v*L*(1+r_c/t)));

SF_leak = p0_sep_sol/p0;        % safety factor against leakage
if SF_leak <= 1.6
   continue;
end
%% Fatigue Failure Analysis of the Bolt:

if i<=3                  % fatigue stress concentration factor (rolled)
    K_f = 2.2;
else
    K_f = 3;
end

F_a = (F_b-F_i)/2;          % alternating force in the bolt
F_m = (F_b+F_i)/2;          % mean force in the bolt

sigma_bif = K_f*F_i/A_t;     % initial stress in the bolt for fatigue (MPa) 
sigma_cif = K_f*F_i/A_c;     % initial compressive stress in the clamped member for fatigue (MPa)

sigma_bff = K_f*F_b/A_t;     % tensile stress in the bolt after pressure applied (MPa)
sigma_cff = K_f*F_c/A_c;     % compressive stress in the clamped member after pressure applied (MPa)

sigma_ax = (sigma_bff-sigma_bif)/2;   % axial alternating stress
sigma_ay = K_f*(p0/2);                % radial alternating stress
sigma_ea = sqrt(sigma_ax^2 + sigma_ay^2 - sigma_ax*sigma_ay);   % alternating equivalent stress

sigma_mx = (sigma_bi+sigma_bf)/2;       % axial mean stress
sigma_my = p0/2;                     % radial mean stress
sigma_em = (sigma_mx+sigma_my)/2 + abs(sigma_mx-sigma_my)/2;    % mean equivalent stress

S_n_prime = 0.5*S_u;    % standard endurance limit
C_load = 1;             % load factor (multiaxial)
if d<8                  % gradient (or size) factor
    C_size = 1;     
else
    C_size = 1.189*d^(-0.097);
end
C_surf =1;          % surface factor
C_rel = 0.753;      % reliability factor
S_n = C_load*C_size*C_surf*C_rel*S_n_prime;   % corrected endurance limit

S_f3 = 0.9*S_u;     % fatigue strength for 1000 cycles

S_f = 10^(log(S_n)+0.1003433319*(log(S_f3)-log(S_n)));   % fatigue strength 5*10^5 cycles of fatigue life

S_a = (S_u-S_y)*S_f/(S_f-S_u);      % maximum allowable alternating stress

SF_fatigue = S_a / sigma_ea;        % safety fatigue failure

if SF_fatigue <= 1.6
    continue;
end
%% Cost Calculation:

switch i
    case 1
        cost = 1.0 + 0.20*d;
        pm = '4.6';
    case 2
        cost = 1.1 + 0.22*d;      
        pm = '4.8';
    case 3
        cost = 1.3 + 0.26*d;
        pm = '5.8';
    case 4
        cost = 3 + 0.6*d;
        pm = '8.8';
    case 5
        cost = 6 + 1.2*d;
        pm = '9.8';
    case 6
        cost = 9 + 1.8*d;
        pm = '10.9';
    case 7
        cost = 11 + 2.2*d;
        pm = '12.9';
end

if cost<cost_memory    
    d_memory = d;
    pitch_memory = pitch;    
    d_p_memory = d_p;
    A_t_memory = A_t;
    pm_memory = pm;
    n_memory = n;
    F_i_memory = F_i;
    K_i_memory = K_i;
    sigma_bi_memory = sigma_bi;
    sigma_ci_memory = sigma_ci;
    F_e_memory = F_e;
    F_b_memory = F_b;
    F_c_memory = F_c;
    delta_sigma_b_memory = delta_sigma_b;
    delta_sigma_c_memory = delta_sigma_c;
    SF_yield_memory = SF_yield;
    SF_leak_memory = SF_leak;
    F_a_memory = F_a;
    F_m_memory = F_m;
    S_n_prime_memory = S_n_prime;
    C_load_memory = C_load;
    C_size_memory = C_size;
    C_surf_memory = C_surf;
    C_rel_memory = C_rel;
    S_n_memory = S_n;
    S_f_memory = S_f;
    K_f_memory = K_f;
    sigma_a_memory = sigma_ax;
    sigma_ea_memory = sigma_ea;
    S_a_memory = S_a;    
    SF_fatigue_memory = SF_fatigue;
    cost_memory = cost;
end
%% 
% 

        end
    end
end
toc
%% 
% 
%% Output Parameters:

    d_memory
    pitch_memory   
    d_p_memory
    A_t_memory
    pm_memory
    n_memory
    F_i_memory
    K_i_memory
    sigma_bi_memory
    sigma_ci_memory
    F_e_memory
    F_b_memory
    F_c_memory
    delta_sigma_b_memory
    delta_sigma_c_memory
    SF_yield_memory
    SF_leak_memory
    F_a_memory
    F_m_memory
    S_n_prime_memory
    C_load_memory
    C_size_memory
    C_surf_memory
    C_rel_memory
    S_n_memory
    S_f_memory
    K_f_memory
    sigma_a_memory
    sigma_ea_memory
    S_a_memory
    SF_fatigue_memory
    cost_memory
%% Printing Results:

file = fopen('output_DesignProject2_Ediz_Kutay.txt','w');
fprintf(file, "OUTPUT:\n");
fprintf(file, "\n");
fprintf(file, "Major diameter of the bolt: "+d_memory+" mm\n");
fprintf(file, "Pitch: "+pitch_memory+" mm\n");
fprintf(file, "Pitch diameter: "+d_p_memory+" mm\n");
fprintf(file, "Tensile stress area: "+A_t_memory+" mm^2\n");
fprintf(file, "Property class of the material: "+pm_memory+" \n");
fprintf(file, "Number of turns of the nut for preload: "+n_memory+" \n");
fprintf(file, "\n");
fprintf(file, "Initial load: "+F_i_memory+" N\n");
fprintf(file, "Preload factor: "+K_i_memory+" \n");
fprintf(file, "Tensile stress in the bolt after preload: "+sigma_bi_memory+" MPa\n");
fprintf(file, "Compressive stress in the clamped member after preload: "+sigma_ci_memory+" MPa\n");
fprintf(file, "\n");
fprintf(file, "Maximum external seperating force: "+F_e_memory+" N\n");
fprintf(file, "Maximum tension in the bolt: "+F_b_memory +" N\n");
fprintf(file, "Minimum compression in the clamped member: "+F_c_memory+" N\n");
fprintf(file, "\n");
fprintf(file, "The increase in the axial stress of the bolt due to pressure: "+delta_sigma_b_memory+" MPa\n");
fprintf(file, "The decrease in the axial stress of the cylinder due to pressure: "+delta_sigma_c_memory+" MPa\n");
fprintf(file, "\n");
fprintf(file, "Safety factor against yielding: "+SF_yield_memory+" \n");
fprintf(file, "\n");
fprintf(file, "Safety factor against leakage: "+SF_leak_memory+" \n");
fprintf(file, "\n");
fprintf(file, "Alternating force in the bolt: "+F_a_memory+" N\n");
fprintf(file, "Mean force in the bolt: "+F_m_memory+" N\n");
fprintf(file, "Standard endurance limit: "+S_n_prime_memory+" MPa\n");
fprintf(file, "Load factor: "+C_load_memory+" \n");
fprintf(file, "Gradient (size) factor: "+C_size_memory+" \n");
fprintf(file, "Surface factor: "+C_surf_memory+" \n");
fprintf(file, "Reliability factor: "+C_rel_memory+" \n");
fprintf(file, "Corrected endurance limit: "+S_n_memory+" \n");
fprintf(file, "Fatigue strength for the given life: "+S_f_memory+" MPa\n");
fprintf(file, "Fatigue stress concentration factor: "+K_f_memory+" \n");
fprintf(file, "Alternating axial stress of the bolt: "+sigma_a_memory+" MPa\n");
fprintf(file, "Alternating equivalent stress of the bolt: "+sigma_ea_memory+" MPa\n");
fprintf(file, "Maximum allowable alternating stress: "+S_a_memory+" MPa\n");
fprintf(file, "\n");
fprintf(file, "Safety factor against fatigue failure: "+SF_fatigue_memory+" \n");
fprintf(file, "\n");
fprintf(file, "Cost: "+cost_memory+" \n");
fprintf(file, "\n");