% % Kalpesh Patil - 130040019
% % filter number 13
% % BSF FIR

m 				       = 13;
q_m 	               = floor(0.1*m);
r_m 	               = m - 10*q_m;
Bl_m 	               = 4 + 0.9*q_m + 2*r_m;
Bh_m 	               = Bl_m + 10;

delta_1 	           = 0.15;
delta_2 	           = 0.15;
f_sample 	           = 100*1000;

omega_p1               = (Bl_m-2)*1000;
omega_s1 			   = (Bl_m)*1000;
omega_s2               = (Bh_m)*1000;
omega_p2               = (Bh_m+2)*1000;

f_analog_array 	       = [omega_p1 omega_s1 omega_s2 omega_p2];

%normalized specs
f_digital_array        = f_analog_array.*2*pi/f_sample;

% Kaiser window parameters
del_omega1             = f_digital_array(4)-f_digital_array(3);
del_omega2             = f_digital_array(2)-f_digital_array(1);
del_omega              = min(abs(del_omega1),abs(del_omega2));
A                      = -20*log10(delta_1);

% Order
N_1 			       =  (A-8)/(2*2.285*del_omega);
N_limit 			   =  ceil(N_1);
N                      =  N_limit + 5;

if(A<21)
   	 alpha=0;
else if(A<=50)
    	alpha  		   =  0.5842*(A-21)^0.4+0.07886*(A-21);
	 else
     	alpha          =  0.1102*(A-8.7);
	 end
end
    
omega_c1			   =  (f_digital_array(2)+f_digital_array(1))*0.5;
omega_c2			   =  (f_digital_array(4)+f_digital_array(3))*0.5;

h_ideal = [];
for k = -N:N
    if(k~=0)
        h_ideal(k+N+1) = - (sin(omega_c2*k)-sin(omega_c1*k))/(pi*k);
    else
        h_ideal(k+N+1) =   0;
    end
end
    
h_ideal(N+1)		   =  1 - ((omega_c2-omega_c1)/pi);
beta				   =  alpha/N;

% Generating Kaiser window 
h_kaiser	 		   =  kaiser(2*N+1,beta);
h_org		 		   =  h_ideal.*h_kaiser';
fvtool(h_org);