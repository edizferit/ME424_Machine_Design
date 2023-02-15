%% *ME424 Design Project-1*
% *Ediz Ferit Kula, 2017405087*
% 
% *Mehmet Kutay Yavuzkol, 2017405045*
% 
% _a = 7, b = 5._

clear all, close all, clc;
a_num = 7; b_num=5;
tic
%% Givens:

n_sun = 3000+100*a_num;     % rotation speed of the sun (rpm)
w_sun = 2*pi*n_sun/60;      % angular speed of the sun (rad/s)
teta = pi/180*20;           % pressure angle
rev_arm = 10^8;             % life of the arm (revolutions)

m_standard = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.50, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5...
    5.0, 5.5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 32, 40, 50];
P_memory =0;
%% *Design Parameters:*
% *Aim:* maximum transmitted power.

HB = 400;                                   % material hardness (HB)
% b = 14*m;                                 % face width (mm)                  
% m                                         % module (mm)
% N_sun                                     % teeth number of the sun
% N_p1                                      % teeth number of the planet 1
% N_p2                                      % teeth number of the planet 2
%% Constraints:

d_housing_max = 520;    % maximum housing inner diameter (mm)
clearance = 5;          % clearance (mm)

n_arm_max = 180+10*b_num;      % maximum rotation speed of the arm (rpm)
n_arm_min = 90+5*b_num;        % minimum rotation speed of the arm (rpm)

v_max = 25;             % maximum pitch line velocity (m/s)

HB_max = 400;           % maximum hardness 
HB_min = 200;           % minimum hardness 

SF_bending = 1.5;       % minimum safety factor against bending fatigue failure (99% reliability)
SF_surface = 1.2;       % minimum safety factor against surface fatigue failure (99% reliability)

CR_min = 1.4;           % minimum contact ratio

%ra_max>ra              % external gears interferance check
delta_N_min = 10;       % internal gears interferance check

%b_max = 14*m;          % maximum face width
%b_min = 9*m;           % minimum face width
%% *Assumptions and Decisions:*
% Assumptions:
%% 
% * Gears have minimum 12 teeth.
% * Planet 1 is rotating around the sun clockwise when the sun rotates clockwise.
%% 
% Decisions and Inferences:
%% 
% * Maximum pitch line velocity occurs between the sun and planet 1.
% * Only the sun and planet 2 will be investigated against failure. This is 
% because the sun and planet 2 have more contacts than their mating gears planet 
% 1 and ring. (It is assummed that radius of planet 2 is smaller than planet 1, 
% and radius of the ring is at least 3 times greater than planet 2's.)
% * Face width is taken as 14m. Increasing face width decreases bending and 
% surface fatigue stresses which are beneficial. It increases mounting factor 
% less.
% * HB is taken as 400. Harder the material the more strength it has if its 
% ultimate strength is less than 1400 MPa which suits our case. It decreases surface 
% factor less. 
%% Loop Starts:

for i=3:38              % To reduce the time of the process, module  ...
                        % is started from 0.4 which doesn't affect the results.
    m=m_standard(i);
%% 
% * Endurance limit calculation against bending fatigue failure:

    S_u = 3.45*HB;      % ultimate strength
    S_ns = 0.5*S_u;     % standard endurance limit
    C_L = 1;            % bending load
    if(m<=5)            % size effect
        C_G = 1;
    else
        C_G = 0.85;
    end
    C_S = -2.276e-12*HB^4 + 1.153e-09*HB^3 + ...      % surface factor
        -3.533e-07*HB^2 + -0.0004911*HB + 0.8635;
    k_r = 0.814;             % reliability factor (99%)
    k_t = 1;                 % temperature factor
    k_ms = 1.4;              % mean stress factor 
    S_n = S_ns.*C_L.*C_G.*C_S.*k_r.*k_t.*k_ms;  % ENDURANCE LIMIT
%% 
% * Calculate maximum bending stress possible:

sigma_b = S_n./SF_bending;     % Maximum bending stress
%% 
% 

    %%=====PITCH LINE VELOCITY CHECK=====%%
    upper_N_sun = 2*v_max./m./w_sun*1000;  

    if(upper_N_sun>1000)        % To reduce the time of the process, maximum teeth number ...
        upper_N_sun = 1000;     % that a gear can have is taken as 1000.
    end
    
    for N_sun=12:upper_N_sun
        
        %%=====Geometry factor calculation of the sun=====%%%
        if (N_sun<60)
            J_sun =-4.171e-05*N_sun^2 + 0.004285*N_sun + 0.1682; %0->60
        elseif (N_sun<150)
            J_sun = -8.235e-07*N_sun^2 + 0.0003829*N_sun + 0.265; %60->150
        else
            J_sun = 0.32; %150->inf
        end

        %%=====HOUSING CHECK=====%%
        upper_N_p1 = (d_housing_max-clearance*2-2*m-m*N_sun)/(2*m); 
        
        if(upper_N_p1>1000)
            upper_N_p1 = 1000;
        end
            
        for N_p1=12:upper_N_p1
            
            %%=====ARM SPEED CHECK=====%%
            lower_N_p2 = N_p1/(w_sun*N_sun/(n_arm_min*2*pi/60)/(N_sun+N_p1)-1);    
            upper_N_p2 = N_p1/(w_sun*N_sun/(n_arm_max*2*pi/60)/(N_sun+N_p1)-1);   
            
            if(upper_N_p2>1000)
                upper_N_p2 = 1000;
            end
            if(lower_N_p2>1000)
                continue;
            end
            if(lower_N_p2<12)
                lower_N_p2 = 12;
            elseif(mod(lower_N_p2,1))
                lower_N_p2 = lower_N_p2 + (1-mod(lower_N_p2,1));
            end
            if(upper_N_p2<12)
                continue;
            end
            
            for N_p2=lower_N_p2:upper_N_p2  
                
                %%=====GEOMETRIC VARIABLES=====%%
                N_ring =  N_sun+N_p1+N_p2; % teeth number of the ring
                
                d_sun = N_sun.*m;    % diameter of the sun (mm)
                d_p1 = N_p1.*m;      % diameter of the planet 1 (mm)
                d_p2 = N_p2.*m;      % diameter of the planet 2 (mm)
                d_ring = N_ring.*m;  % diameter of the ring (mm)
                d_housing = d_sun+2*d_p1+2*m; % inner diameter of the housing
                
                b = 14*m;            % face width (mm), for maximum strength
                
                r_sun = d_sun/2;     % radius of the sun (mm)
                r_p1 = d_p1/2;       % radius of the planet 1 (mm)
                r_p2 = d_p2/2;       % radius of the planet 2 (mm)
                r_ring = d_ring/2;   % radius of the ring (mm)
                c_s1 = r_sun+r_p1;   % central distance of the sun and planet 1(mm)
                c_r2 = r_ring-r_p2;  % central distance of the ring and planet 2 (mm)
                
                ra_sun = r_sun+m;    % addendum radius of the sun (mm)
                ra_p1 = r_p1+m;      % addendum radius of the planet 1 (mm)
                ra_p2 = r_p2+m;      % addendum radius of the planet 2 (mm)
                ra_ring = r_ring-m;  % addendum radius of the ring (mm)
                
                R_sun = r_sun.*cos(teta);    % base circle radius of the sun(mm)
                R_p1 = r_p1.*cos(teta);      % base circle radius of the planet 1(mm)
                R_p2 = r_p2.*cos(teta);      % base circle radius of the planet 2(mm)
                R_ring = r_ring.*cos(teta);  % base circle radiusof the ring (mm)
                
                %%=====INTERFERANCE CHECK=====%%
                ra_max_sun = sqrt(r_sun.^2*cos(teta).^2+c_s1.^2*sin(teta).^2);     % maximum addendum radius of the sun
                ra_max_p1 = sqrt(r_p1.^2*cos(teta).^2+c_s1.^2*sin(teta).^2);       % maximum addendum radius of the planet 1             
                delta_N = N_ring-N_p2;      % delta_N must be bigger than 10 
                if(ra_sun > ra_max_sun || ra_p1 > ra_max_p1 || delta_N < delta_N_min)
                    continue;
                end
                
                %%=====CONTACT RATIO CHECK=====%% 
                CR_s1 = (sqrt(ra_sun.^2-R_sun.^2) + sqrt(ra_p1.^2-R_p1.^2) - ...
                    c_s1.*sin(teta))./(pi*m*cos(teta));   % contact ratio between the sun and planet 1
                CR_r2 = (sqrt(ra_p2.^2-R_p2.^2) - sqrt(ra_ring.^2-R_ring.^2) + ...
                    c_r2.*sin(teta))./(pi*m*cos(teta)); % contact ratio between the ring and planet 2
                if(CR_s1 < CR_min || CR_r2 < CR_min)
                    continue;
                end
                
                %%=====SPEED VARIABLES=====%%
                n_arm = n_sun.*r_sun.*r_p2./(r_sun+r_p1)./(r_p1+r_p2); % rotation speed of the arm (rpm)
                
                n_p1 = n_sun*d_sun/d_p1;     % rotational velocity of the planet 1 (rpm)
                n_p2 = n_p1;                 % rotational velocity of the planet 2 (rpm)
                n_ring = 0;                  % rotational velocity of the ring (rpm)
                
                w_p1 = 2*pi*n_p1/60;         % angular speed of the planet 1 (rad/s)
                w_p2 = 2*pi*n_p2/60;         % angular speed of the planet 2 (rad/s)
                w_ring = 2*pi*n_ring/60;     % angular speed of the ring (rad/s)
                
                v_sun = w_sun.*r_sun/1000;    % pitch line speed of the sun and planet 1 (m/s)
                v_p2 = w_p2*r_p2/1000;        % average pitch line speed between planet 2 and ring (m/s)
%% Failure Analysis
% Bending Fatigue Failure
% $\sigma_b \;=\;\frac{F_t }{\textrm{mbJ}}{\;K}_v {\;K}_o {\;K}_m$ ,    $S_n 
% \;={{\;S}_n }^{\prime } {\;C}_L {\;C}_G {\;C}_S {\;k}_r {\;k}_t {\;k}_{\textrm{ms}}$  
% ,     $\textrm{SF}=\frac{S_n }{\sigma_b }$  ,    $\sigma_b =\frac{S_n }{\textrm{SF}}$  
% ,    ${\mathit{\mathbf{F}}}_{\mathit{\mathbf{t}}} =\frac{\sigma_b \;m\;b\;J}{{\;K}_v 
% {\;K}_o {\;K}_m }$
% 
% *Procedure:* Endurance limit is calculated before, safety factor is known, 
% so maximum possible bending stress is also calculated beforehand. Calculate 
% F_t for planet 1 and sun. 
%% 
% * Values for bending stress:

                %% J_sun: Geometry factor calculation of the sun is done before for optimization.
                
                %% Geometry factor calculation of the planet 2
                if (N_p2<60)
                    J_p2 =-4.171e-05*N_p2^2 + 0.004285*N_p2 + 0.1682; %0->60
                elseif (N_p2<150)
                    J_p2 = -8.235e-07*N_p2^2 + 0.0003829*N_p2 + 0.265; %60->150
                else
                    J_p2 = 0.32; %150->inf
                end
                
                K_v_sun = (1200+200*v_sun)/1200;          % velocity factor for sun
                K_v_p2 = (1200+200*v_p2)/1200;            % velocity factor for planet 2
                K_o = 1.75;                               % overload factor: uniform source, high shock
                if(b>0 && b<=50)                          % mounting factor: not high precision
                    K_m = 1.6;
                elseif(b>50 && b<=150)
                    K_m = 1.6+(b-50)/100*0.1;
                elseif(b>150 && b<=225)
                    K_m = 1.7+(b-150)/75*0.1;
                elseif(b>225 && b<=405)
                    K_m = 1.8+(b-225)/180*0.4;
                else
                    K_m = 2.2;
                end
%% 
% * Calculate maximum possible transmitted force:

                F_sun_bending = sigma_b./(K_v_sun.*K_o.*K_m).*m.*b.*J_sun; % max transmitted force at the sun (N)
                F_p2_bending = sigma_b./(K_v_p2.*K_o.*K_m).*m.*b.*J_p2; % max transmitted force at the planet 2 (N)
                F_sun_bending_from_p2 = F_p2_bending.*r_p2./r_p1;   % max transmitted force at the sun according to planet 2 (N)
                
                if(F_sun_bending>F_sun_bending_from_p2)
                    F_t_bending = F_sun_bending_from_p2;    % maximum transmitted force for the sun acc. to bending fatigue failure
                else
                    F_t_bending = F_sun_bending;
                end
% Surface Fatigue Failure
% $\sigma_H \;=C_p \sqrt{\frac{F_t }{md_p I}K_v {\;K}_o {\;K}_m }$ ,    $S_H 
% \;=S_{\textrm{fe}} {\;C}_{\textrm{Li}} {\;C}_R$ ,   $\textrm{SF}={\left(\frac{S_H 
% }{\sigma_H }\right)}^2$ ,      $\sigma_H =\frac{S_H }{\sqrt{\textrm{SF}}}$  
% ,      ${\mathit{\mathbf{F}}}_{\mathit{\mathbf{t}}} =\frac{{\sigma_H }^2 \;m\;d_p 
% \;I}{{C_p }^2 \;K_v {\;K}_o {\;K}_m }$
% 
% *Procedure:* Calculate surface fatigue strength, safety factor is known, calculate 
% surface fatigue stress, calculate F_t. (planet 2, sun)
%% 
% * Contact number calculations for required life of arm: 
% * ${\textrm{Contact}}_{\textrm{p2}} ={\textrm{Life}}_{\textrm{arm}} *\frac{N_{\textrm{ring}} 
% }{N_{\textrm{p2}} }$ ,                                           ${\textrm{Contact}}_{\textrm{sun}} 
% ={3*\textrm{Life}}_{\textrm{arm}} *\left(\frac{\left(N_{\textrm{sun}} +N_{\textrm{p1}} 
% \right)*\left(N_{\textrm{p1}} +N_{\textrm{p2}} \right)}{N_{\textrm{sun}} *N_{\textrm{p2}} 
% }-1\right)$

                contact_p2 = N_ring/N_p2*rev_arm;                                  % contact number of planet 2
                contact_sun = 3*((N_sun+N_p1)*(N_p1+N_p2)/(N_sun*N_p2)-1)*rev_arm; % contact number of sun
%% 
% * Values for surface fatigue stress (maximum contact stress) and surface fatigue 
% strength:

                C_p = 191;                                      % elastic coefficient (MPa^1/2), steel pinion and gear
                R_i = r_ring/r_p2;
                I_internal = cos(teta)*sin(teta)*R_i/(R_i-1)/2; % geometry factor internal gear (planet 2)
                R_e = r_p1/r_sun;
                I_external = cos(teta)*sin(teta)*R_e/(R_e+1)/2; % geometry factor external gear (sun)
                S_fe = 2.76*HB - 69;                            % standard surface fatigue strength (MPa)
                
                %% Life factor calculation %% a conservative choice is made
                if(contact_sun<10^7) % life factor of the sun
                    C_Li_sun = 1;
                elseif(contact_sun<10^8)
                    C_Li_sun = 0.9;
                elseif(contact_sun<10^9)
                    C_Li_sun = 0.8;
                elseif(contact_sun<10^10)
                    C_Li_sun = 0.72;
                else
                    C_Li_sun = 0.655;
                end
                if(contact_p2<10^7) % life factor of the planet 2
                    C_Li_p2 = 1;
                elseif(contact_p2<10^8)
                    C_Li_p2 = 0.9;
                elseif(contact_p2<10^9)
                    C_Li_p2 = 0.8;
                elseif(contact_p2<10^10)
                    C_Li_p2 = 0.72;
                else
                    C_Li_p2 = 0.655;
                end
                C_R = 1;                                % reliability factor (99%)
                
                S_H_p2 = S_fe.*C_Li_p2.*C_R;          % surface fatigue strength for planet 2
                S_H_sun = S_fe.*C_Li_sun.*C_R;        % surface fatigue strength for sun
%% 
% * Calculate maximum surface fstigue stress possible:

                sigma_H_p2 = S_H_p2./sqrt(SF_surface);
                sigma_H_sun = S_H_sun./sqrt(SF_surface);
%% 
% * Calculate maximum possible transmitted force:

                F_p2_surface = (sigma_H_p2/C_p)^2*b*d_p2*I_internal/(K_v_p2*K_o*K_m);
                F_sun_surface = (sigma_H_sun/C_p)^2*b*d_sun*I_external/(K_v_sun*K_o*K_m);
                F_sun_surface_from_p2 = F_p2_surface*r_p2/r_p1;   % max transmitted force at the sun according to planet 2 (N)
                
                if(F_sun_surface>F_sun_surface_from_p2)
                    F_t_surface = F_sun_surface_from_p2;   % maximum transmitted force for the sun acc. to surface fatigue failure
                else
                    F_t_surface = F_sun_surface;
                end
% Compare the Maximum Forces Possible due to Bending and Surface Fatigue Failure:

                F_t = min(F_t_bending,F_t_surface);  % Maximum transmitted force at the sun (N)
%% Power

                P = F_t.*(3*w_sun.*r_sun)./1000;     % Power (Nm/s)
                if (P>P_memory)                      % Selecting largest value of work for maximum power
                    P_memory = P;           % power transmitted
                    HB_memory = HB;         % material (HB)
                    n_sun_memory = n_sun;   % input speed (rpm)
                    n_arm_memory = n_arm;   % output speed (rpm)
                    m_memory = m;           % module (mm)
                    d_sun_memory = d_sun;   % pitch diameter of the sun (mm)
                    d_p1_memory = d_p1;     % pitch diameter of the planet 1 (mm)
                    d_p2_memory = d_p2;     % pitch diameter of the planet 2 (mm)
                    d_ring_memory = d_ring; % pitch diameter of the ring (mm)
                    N_sun_memory = N_sun;   % number of teeth of the sun
                    N_p1_memory = N_p1;     % number of teeth of the planet 1
                    N_p2_memory = N_p2;     % number of teeth of the planet 2
                    N_ring_memory = N_ring; % number of teeth of the ring
                    d_housing_memory = d_housing; % inner diameter of the housing (mm)
                    n_p1_memory = n_p1;     % speeds of the planets (rpm)
                    v_sun_memory = v_sun;   % pitch line velocity of the sun (m/s)
                    F_t_memory = F_t; % tangential force between the sun and planet 1
                    b_memory = b;           % face width (mm)
                    CR_s1_memory = CR_s1;   % contact ratio between the sun and planet 1
                    CR_r2_memory = CR_r2;   % contact ratio between the ring and planet 2                    
                    
                    %% Necessary values to specify the most critical gear:
                    F_sun_bending_memory = F_sun_bending;
                    F_sun_surface_memory = F_sun_surface;
                    F_p2_bending_memory = F_p2_bending;
                    F_p2_surface_memory = F_p2_surface;
                    F_sun_bending_from_p2_memory = F_sun_bending_from_p2;
                    F_sun_surface_from_p2_memory = F_sun_surface_from_p2;
                    sigma_b_memory = sigma_b;
                    sigma_H_sun_memory = sigma_H_sun;
                    sigma_H_p2_memory = sigma_H_p2;                        
                    
                    % Most critical gear is planet 2
                    J_cr = J_p2;
                    K_v_cr = K_v_p2;
                    K_m_cr = K_m;
                    K_o_cr = K_o;
                    sigma_b_cr = F_p2_surface_memory*K_v_cr*K_m_cr*K_o_cr/m_memory/b_memory/J_cr;
                    S_ns_cr = S_ns;
                    C_L_cr = C_L;
                    C_G_cr = C_G; 
                    C_S_cr = C_S;
                    k_r_cr = k_r;
                    k_ms_cr = k_ms;
                    S_n_cr = S_n; 
                    SF_bending_cr = S_n_cr/sigma_b_cr;
                    C_p_cr = C_p;
                    I_cr = I_internal;
                    sigma_H_cr = C_p_cr*sqrt(F_p2_surface_memory*K_v_cr*K_o_cr*K_m_cr/(b_memory*d_p2_memory*I_cr));
                    S_fe_cr = S_fe;
                    contact_cr = contact_p2;
                    C_Li_cr = C_Li_p2;
                    C_R_cr = C_R;
                    S_H_cr = S_H_p2;
                    SF_surface_cr = (S_H_cr/sigma_H_cr)^2;
                    
                end
%% 

            end
        end
    end
end
P_final = P_memory
toc
%% Calculating Parameters of the Optimum Design

% tangential force between the ring and planet 2
F_t_p2_memory = P_final*d_p1_memory/(3*w_sun*d_sun_memory/2)/d_p2_memory*1000;

% the force applied by planet 2 on the arm
F_arm_memory = (1+d_p1_memory/d_p2_memory)*P_final/(3*w_sun*d_sun_memory/2)*1000;

% interferance check   (r_ap_max-r_ap) and (delta_N-delta_N_min)
ra_max_sun = sqrt((d_sun_memory/2)^2*cos(teta)^2+(d_sun_memory/2+d_p1_memory/2)^2*sin(teta)^2);     % maximum addendum radius of the sun
ra_sun = d_sun_memory/2+m_memory;  % addendum radius of the sun
ra_max_p1 = sqrt((d_p1_memory/2)^2*cos(teta)^2+(d_sun_memory/2+d_p1_memory/2)^2*sin(teta)^2);       % maximum addendum radius of the planet 1
ra_p1 = d_p1/2+m_memory;           % addendum radius of the planet 1
int_check_sun = ra_max_sun - ra_sun;
int_check_p1 = ra_max_p1 - ra_p1;
int_check_rp2 = N_ring_memory-N_p2_memory-delta_N_min;    
%% Printing Results:

    F_sun_bending_memory
    F_sun_surface_memory
    F_p2_bending_memory
    F_p2_surface_memory
    F_sun_bending_from_p2_memory
    F_sun_surface_from_p2_memory
    sigma_b_memory
    sigma_H_sun_memory
    sigma_H_p2_memory
    
file = fopen('output_DesignProject1_Ediz_Kutay.txt','w');
fprintf(file, "OUTPUT:\n");
fprintf(file, "\n");
fprintf(file, "Power transmitted: "+P_memory+" W\n");
fprintf(file, "Material (HB): "+HB_memory+" HB\n");
fprintf(file, "Input speed: "+n_sun_memory+" rpm\n");
fprintf(file, "Output speed: "+n_arm_memory+" rpm\n");
fprintf(file, "Module: "+m_memory+" mm\n");
fprintf(file, "\n");
fprintf(file, "Diameters of the gears: \n");
fprintf(file, "Pitch diameter of the sun: "+d_sun_memory+" mm\n");
fprintf(file, "Pitch diameter of the planet 1: "+d_p1_memory+" mm\n");
fprintf(file, "Pitch diameter of the planet 2: "+d_p2_memory+" mm\n");
fprintf(file, "Pitch diameter of the ring: "+d_ring_memory+" mm\n");
fprintf(file, "\n");
fprintf(file, "Number of teeth in the gears: \n");
fprintf(file, "Number of teeth of the sun: "+N_sun_memory +" \n");
fprintf(file, "Number of teeth of the planet 1: "+N_p1_memory+" \n");
fprintf(file, "Number of teeth of the planet 2: "+N_p2_memory+" \n");
fprintf(file, "Number of teeth of the ring: "+N_ring_memory+" \n");
fprintf(file, "\n");
fprintf(file, "Inner diameter of the housing: "+d_housing_memory+" mm\n");
fprintf(file, "Speeds of the planets: "+n_p1_memory+" rpm\n");
fprintf(file, "Pitch line velocity of the sun: "+v_sun_memory+" m/s\n");
fprintf(file, "\n");
fprintf(file, "Tangential forces in the gears: \n");
fprintf(file, "Tangential force between the sun and planet 1: "+F_t_memory+" N\n");
fprintf(file, "Tangential force between the ring and planet 2: "+F_t_p2_memory+" N\n");
fprintf(file, "The force applied by planet 2 on the arm: "+F_arm_memory+" N\n");
fprintf(file, "\n");
fprintf(file, "Face width: "+b_memory+" mm\n");
fprintf(file, "\n");
fprintf(file, "Contact ratio check:\n");
fprintf(file, "Contact ratio between the sun and planet 1: "+CR_s1_memory+" \n");
fprintf(file, "Contact ratio between the ring and planet 2: "+CR_r2_memory+" \n");
fprintf(file, "\n");
fprintf(file, "Interferance check:\n");
fprintf(file, "ra_max - ra for the sun: "+int_check_sun+" mm\n");
fprintf(file, "ra_max - ra for the planet 1: "+int_check_p1+" mm\n");
fprintf(file, "delta_N - delta_N_min for internal gears: "+int_check_rp2+" \n");
fprintf(file, "\n");
fprintf(file, "Quantities for the most critical gear: Planet 2\n");
fprintf(file, "Geometry factor (bending fatigue failure): "+J_cr+" \n");
fprintf(file, "Velocity factor: "+K_v_cr+" \n");
fprintf(file, "Mounting factor: "+K_m_cr+" \n");
fprintf(file, "Overload factor: "+K_o_cr+" \n");
fprintf(file, "Bending stress: "+sigma_b_cr+" MPa\n");
fprintf(file, "Standard endurance limit: "+S_ns_cr+" MPa\n");
fprintf(file, "Load factor: "+C_L_cr+" \n");
fprintf(file, "Gradient factor: "+C_G_cr+" \n");
fprintf(file, "Surface factor: "+C_S_cr+" \n");
fprintf(file, "Reliability factor: "+k_r_cr+" \n");
fprintf(file, "Mean stress factor: "+k_ms_cr+" \n");
fprintf(file, "Endurance limit: "+S_n_cr+" MPa\n");
fprintf(file, "Safety factor against bending fatigue failure: "+SF_bending_cr+" \n");
fprintf(file, "\n");
fprintf(file, "Elastic coefficient: "+C_p_cr+" MPa^(1/2)\n");
fprintf(file, "Geometry factor (surface fatigue failure): "+I_cr+" \n");
fprintf(file, "Surface fatigue stress: "+sigma_H_cr+" MPa\n");
fprintf(file, "Strandard surface fatigue strength: "+S_fe_cr+" MPa\n");
fprintf(file, "Number of contacts: "+contact_cr+" \n");
fprintf(file, "Life factor: "+C_Li_cr+" \n");
fprintf(file, "Reliability factor: "+C_R_cr+" \n");
fprintf(file, "Surface fatigue strength: "+S_H_cr+" \n");
fprintf(file, "Safety factor against surface fatigue failure: "+SF_surface_cr+" \n");


    P_memory                % power transmitted
    HB_memory               % material (HB)
    n_sun_memory            % input speed (rpm)
    n_arm_memory            % output speed (rpm)
    m_memory                % module (mm)
    d_sun_memory            % pitch diameter of the sun (mm)
    d_p1_memory             % pitch diameter of the planet 1 (mm)
    d_p2_memory             % pitch diameter of the planet 2 (mm)
    d_ring_memory           % pitch diameter of the ring (mm)
    N_sun_memory            % number of teeth of the sun
    N_p1_memory             % number of teeth of the planet 1
    N_p2_memory             % number of teeth of the planet 2
    N_ring_memory           % number of teeth of the ring
    d_housing_memory        % inner diameter of the housing (mm)
    n_p1_memory             % speeds of the planets (rpm)
    v_sun_memory            % pitch line velocity of the sun (m/s)
    F_t_memory              % tangential force between the sun and planet 1
    F_t_p2_memory           % tangential force between the ring and planet 2
    F_arm_memory            % the force applied by planet 2 on the arm
    b_memory                % face width (mm)
    CR_s1_memory            % contact ratio between the sun and planet 1
    CR_r2_memory            % contact ratio between the ring and planet 2
    int_check_sun           % ra_max-ra for sun
    int_check_p1            % ra_max-ra for planet 1
    int_check_rp2           % delta_N-delta_N_min for ring and planet 2
        
% most critical gear is planet 2 
        J_cr 
        K_v_cr 
        K_m_cr
        K_o_cr
        sigma_b_cr 
        S_ns_cr 
        C_L_cr 
        C_G_cr 
        C_S_cr 
        k_r_cr 
        k_ms_cr 
        S_n_cr 
        SF_bending_cr 
        C_p_cr 
        I_cr 
        sigma_H_cr 
        S_fe_cr
        contact_cr 
        C_Li_cr 
        C_R_cr 
        S_H_cr 
        SF_surface_cr