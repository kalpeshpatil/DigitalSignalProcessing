% % Kalpesh Patil - 130040019
% % filter number 13
% % BPF IIR
% % Butterworth 

%finding cutoff frequencies
m 						=   13;
q_m 					=   floor(0.1*m);
r_m 					=   m - 10*q_m;
Bl_m 					=   4 + 0.7*q_m + 2*r_m;
Bh_m 					=   Bl_m + 10;

delta_1 				= 	0.15;
delta_2 				= 	0.15;
f_sample 				=   100*1000;

omega_s1 				=   (Bl_m-2)*1000;
omega_p1 				=   (Bl_m)*1000;
omega_p2 				=   (Bh_m)*1000;
omega_s2 				=   (Bh_m+2)*1000;

f_analog_array 			= 	[omega_s1 omega_p1 omega_p2 omega_s2];

%normalized specs
f_digital_array 		= 	f_analog_array.*2*pi/f_sample;
f_eqv_analog 			= 	tan(f_digital_array/2);

%analog LPF
omega_not 				= 	sqrt(f_eqv_analog(2)*f_eqv_analog(3));
B 						= 	f_eqv_analog(3) - f_eqv_analog(2);

%frequency transformation
f_eqv_analog_LPF 		= 	(f_eqv_analog.^2 - omega_not^2)./(B*f_eqv_analog);

%designing this LPF
D1                      =	(1/(1-delta_1)^2)-1;
D2                      =	(1/delta_2^2)-1;
mod_f_eqv_analog_LPF	=	abs(f_eqv_analog_LPF);
stringent_omega_s		=	min(mod_f_eqv_analog_LPF(1),mod_f_eqv_analog_LPF(4));
temp					=	log(sqrt(D2)/sqrt(D1))/log(stringent_omega_s/f_eqv_analog_LPF(3));
N 						=   ceil(temp);

%finding poles of equivalent low pass filter
omega_p 				=	f_eqv_analog_LPF(3);
omega_c 				=	((omega_p/(D1^(1/(2*N))))+(stringent_omega_s/(D2^(1/(2*N)))))/2;
thetas  				=   (2*[1:2*N] - 1)*(pi)/(2*N);

gain_LF				    = 	omega_c^N;
poles 					= 	1i*omega_c*exp(1i*thetas);
poles_LF 				= 	poles(1:N);
zeros_LF 				= 	[];
[num_LF,den_LF] 		= 	zp2tf(zeros_LF,poles_LF,gain_LF);
tf_LF 					= 	tf(num_LF,den_LF);

%converting back to band pass filter
poles_BF 				= 	zeros(1,2*N);
poles_BF(1:N) 			= 	(B/2).*(poles_LF-sqrt(poles_LF.^2 - 4*omega_not^2/(B^2)));
poles_BF(N+1:2*N) 		=	(B/2).*(poles_LF+ sqrt(poles_LF.^2 - 4*omega_not^2/(B^2)));
zeros_BF 				= 	zeros(1,N);
gain_BF 				= 	gain_LF*B^N;
[num_BF,den_BF] 		= 	zp2tf(zeros_BF',poles_BF,gain_BF);
tf_BF 					= 	tf(num_BF, den_BF);

%converting to z domain
poles_BF_z  			= 	(1 + poles_BF)./(1 - poles_BF);
temp 					= 	prod(1-poles_BF);
gain_BF_z 				= 	gain_BF/temp;
zeros_BF_z 				= 	[ones(1,N),-ones(1,N)];
[num_BF_z, den_BF_z] 	= 	zp2tf(zeros_BF_z',poles_BF_z,gain_BF_z);
DigitalFilter 			= 	tf(num_BF_z, den_BF_z,1/f_sample);
h 						= 	fvtool(num_BF_z, den_BF_z);
