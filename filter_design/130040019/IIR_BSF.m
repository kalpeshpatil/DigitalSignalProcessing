% % Kalpesh Patil - 130040019
% % filter number 13
% % BSF IIR
% % Chebyshev 

%finding cutoff frequencies
m                       =   13;
q_m                     =   floor(0.1*m);
r_m                     =   m - 10*q_m;
Bl_m                    =   4 + 0.9*q_m + 2*r_m;
Bh_m                    =   Bl_m + 10;

delta_1                 =   0.15;
delta_2                 =   0.15;
f_sample                =   100*1000;

omega_p1                =   (Bl_m-2)*1000;
omega_s1                =   (Bl_m)*1000;
omega_s2                =   (Bh_m)*1000;
omega_p2                =   (Bh_m+2)*1000;

f_analog_array          =   [omega_p1 omega_s1 omega_s2 omega_p2];

%normalized specs
f_digital_array         =   f_analog_array.*2*pi/f_sample;
f_eqv_analog            =   tan(f_digital_array/2);

%analog LPF
omega_not               =   sqrt(f_eqv_analog(1)*f_eqv_analog(4));
B                       =   f_eqv_analog(4) - f_eqv_analog(1);

%frequency transformation
f_eqv_analog_LPF        =   (B*f_eqv_analog)./(omega_not^2 - f_eqv_analog.^2);

%designing this LPF
D1                      =   (1/(1-delta_1)^2)-1;
D2                      =   (1/delta_2^2)-1;
epsilon                 =   sqrt(D1);
mod_f_eqv_analog_LPF    =   abs(f_eqv_analog_LPF);
stringent_omega_s       =   min(mod_f_eqv_analog_LPF(2),mod_f_eqv_analog_LPF(3));
temp                    =   acosh(sqrt(D2)/sqrt(D1))/acosh(stringent_omega_s/f_eqv_analog_LPF(1));
N                       =   ceil(temp);


%finding poles of equivalent low pass filter
omega_p                 =   mod_f_eqv_analog_LPF(4);
Ak                      =   ([1:2*N]*2 + 1)*pi/(2*N);
Bk                      =   asinh(1/epsilon)/N;
temp_real               =   omega_p*sin(Ak)*sinh(Bk);
temp_imag               =   omega_p*cos(Ak)*cosh(Bk);
temp                    =   temp_real + 1i*temp_imag;
poles_LF                =   zeros(1,N);
t                       =   1;
for k = 1:2*N
    if(temp_real(k) < 0)
      poles_LF(t)       =   temp(k);
      t                 =   t+1;
    end
end

if (mod(N,2) == 0)
    d                   =   1/sqrt(1+epsilon^2);
else
    d                   =   1;
end

g                       =   ((-1)^N)*prod(poles_LF);
gain_k                  =   d*g;
zeros_LF                =   [];
[num_LF,den_LF]         =   zp2tf(zeros_LF',poles_LF,gain_k);
tf_LF                   =   tf(num_LF,den_LF);

%converting back to band stop filter
poles_BF                =   zeros(1,2*N);
poles_BF(1:N)           =   (B/2).*(1./poles_LF + sqrt(1./poles_LF.^2 - 4*omega_not^2/(B^2)));
poles_BF(N+1:2*N)       =   (B/2).*(1./poles_LF - sqrt(1./poles_LF.^2 - 4*omega_not^2/(B^2)));
zeros_BF                =   [1i*omega_not*ones(1,N)' ; -1i*omega_not*ones(1,N)'];
gain_BF                 =   gain_k/g;
[num_BF,den_BF]         =   zp2tf(zeros_BF,poles_BF,gain_BF);
tf_BF                   =   tf(num_BF, den_BF);


%converting to z domain
poles_BF_z              =   (1 + poles_BF)./(1 - poles_BF);
temp                    =   prod(1-poles_BF);
gain_BF_z               =   gain_BF*((1+omega_not^2)^N)/(prod(1-poles_BF));
zeros_BF_z              =   (1 + zeros_BF)./(1 - zeros_BF);
[num_BF_z, den_BF_z]    =   zp2tf(zeros_BF_z,poles_BF_z,gain_BF_z);
DigitalFilter           =   tf(num_BF_z, den_BF_z,1/f_sample);
h                       =   fvtool(num_BF_z, den_BF_z);
