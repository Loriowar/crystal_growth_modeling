% T_amp = [0.72*10^(-2), 0.22*10^(-2), 0.72*10^(-3), 0.52*10^(-3), 0.22*10^(-3), 0.72*10^(-4), 0.52*10^(-4),0.22*10^(-4)];
% delta_x_amp = [0.566*10^(-5), 0.194*10^(-5), 0.057*10^(-5), 0.025*10^(-5), 0.011*10^(-5), 0.003*10^(-5), 0.003*10^(-5), 0.009*10^(-5)];
% Vp_amp = [0.445, 0.156, 0.042, 0.019, 0.007, 0.004, 0.0015, 0.007];
% C_amp = [0.014*10^(-4), 0.0071*10^(-4), 0.0051*10^(-4), 0.0005*10^(-4), 0.0003*10^(-4), 0.0001*10^(-4), 0.0004*10^(-4), 0.0005*10^(-4)];
% deltaC_on_C = [3.8*10^(-2), 3.7*10^(-2)	, 2.7*10^(-2), 2.6*10^(-3), 1.5^10^(-3), 5.2*10^(-4), 2.1*10^(-3), 2.6*10^(-3)];
T_amp = [0.12*10^(-3)	0.32*10^(-3)	0.52*10^(-3)	0.72*10^(-3)	0.92*10^(-3)	0.12*10^(-2)	0.32*10^(-2)	0.52*10^(-2)	1.02*10^(-2)	1.52*10^(-2)];
delta_x_amp = [0.007*10^(-5)	0.019*10^(-5)	0.029*10^(-5)	0.043*10^(-5)	0.055*10^(-5)	0.073*10^(-5)	0.193*10^(-5)	0.31*10^(-5)	0.62*10^(-5)	0.92*10^(-5)];
Vp_amp = [0.005	0.015	0.025	0.035	0.045	0.058	0.154	0.248	0.492	0.755];
C_amp = [0.00005	0.0002	0.0004	0.0007	0.0008	0.001	0.0029	0.0048	0.0096	0.014]*10^(-4);
deltaC_on_C = [0.0004	 0.0012	0.0024	 0.0042	 0.0048	 0.006	0.0174	 0.0289	 0.0578	 0.0843];

figure

subplot(3, 1, 1)
plot(T_amp, delta_x_amp)
title('Delta x')
xlabel('T, dimensionless')
ylabel('\Delta x, dimensionless')
grid

subplot(3, 1, 2)
plot(T_amp, Vp_amp)
title('Growth rate')
xlabel('T, dimensionless')
ylabel('\nu, sm/hour')
grid

subplot(3, 1, 3)
plot(T_amp, C_amp)
title('Concentration')
xlabel('T, dimensionless')
ylabel('C, dimensionless')
grid

% subplot(4, 1, 4)
figure
plot(T_amp, deltaC_on_C)
title('\delta C divide by C')
xlabel('T, dimensionless')
ylabel('\delta C/C, dimencionless')
grid

% T_amp = -log10(T_amp);
delta_x_amp = -log10(delta_x_amp);
Vp_amp = -log10(Vp_amp);
C_amp = -log10(C_amp);
deltaC_on_C = log10(deltaC_on_C);
figure

subplot(3, 1, 1)
plot(T_amp, delta_x_amp)
title('Delta x')
xlabel('T, dimensionless')
ylabel('-log10(\Delta x)')
grid

subplot(3, 1, 2)
plot(T_amp, Vp_amp)
title('Growth rate')
xlabel('T, dimensionless')
ylabel('-log10(\nu)')
grid

subplot(3, 1, 3)
plot(T_amp, C_amp)
title('Concentration')
xlabel('T, dimensionless')
ylabel('-log10(C)')
grid

% subplot(4, 1, 4)
figure
plot(T_amp, deltaC_on_C)
title('\delta C divide by C')
xlabel('T, dimensionless')
ylabel('-log10(\delta C/C), dimensionless')
grid

figure
plot(T_amp*(1414 + 273)/20, deltaC_on_C)
title('\delta C divide by C')
xlabel('T, K')
ylabel('log(\delta C/C), dimensionless')
grid

% figure
% plot(-log10(T_amp*(1414 + 273)/20), deltaC_on_C)
% title('\delta C divide by C')
% xlabel('-log10(T)')
% ylabel('-log10(\delta C/C)')
% grid

