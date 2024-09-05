%%
% HW1 Radar Imaging
% Pouya Aminaie
clc;clear; close all

%%
% Defining the parameters:
parameter.f0 = 77*10^(9);              % Carrier Frequency [Hz]
parameter.B = 1*10^(9);                % Bandwiths of system [Hz]
parameter.fs=4*10^(9);                 % Sampling Frequency [Hz]
parameter.rho = 2;                     % Angular Resolution  [deg]
parameter.num_ant = 9;                 % Number of Antenna 
parameter.C=3*10^8;                    % Speed of light [m/s]
%%
% Initialization:
lambda=parameter.C/(parameter.f0);     % Wavelength Calculation
L=lambda/(((2*pi)/180)*parameter.rho); % Length of array Calculation


% Then, the angular FOV should be limited by appling Rayleigth criterion 
% on spacing. Therefore, Delta_X < (lambda/2) ===> maximum dx = labmba/2
dx=lambda/2;                           % space between two consequential antennas
%dx=0.5;
%%
% Antenna Positioning:
Array_x=[];                            % Horizontal Position of elements
for i=1:parameter.num_ant
Array_x=[Array_x, (-floor(parameter.num_ant/2)+(i-1))*dx];
end
Array_y=zeros(1,parameter.num_ant);    % Vertical Position of element

% Targert Positioning:
Target_x= 14;                          % Horizoontal Distance[m]
Target_y= 12;                          % Vertical Distance[m]
rp=sqrt((Target_x)^2 + (Target_y)^2);  % Distance of target from the center
%%
% Distance Calculations: 
r_n=zeros(1,parameter.num_ant);  % Distant of target from element [m]

x0=Target_x;                           % Horizoontal Distance[m]
y0=Target_y;                           % Vertical Distance [m]

for i=1:parameter.num_ant
% Try to formulate geometrical distance 
r_n(i)=sqrt((Array_x(i)-x0)^2 + (Array_y(i)-y0)^2);
end

%%
% in this part, we scan the space with the resolution of 2 deg 
% in order to find the phase of received signal at each elements.

n=[-(numel(Array_x)-1)/2:1:+(numel(Array_x)-1)/2];
angles=0:parameter.rho:180-1;
angles=deg2rad(angles);



%%
%ph_global = [];
 ph = [];
syms sin_theta
for i = 1:numel(n)
   

    
        % Populate ph inside the inner loop
        ph = [ph, ((-2*pi)/lambda)*(rp - abs(n(i))*dx*sin_theta)];
    

    % Concatenate ph to ph_global after the inner loop
    
end
ph_global = ph;
clear i
Srx=exp(i*ph_global);



%%
 % % phi calculations:
 % phi=[];
 % for i=2:numel(Array_x)
 %    phi= [phi    ph_global(i)-ph_global(i-1)];
 % end

% %%
%  %DoA Calculation
% 
%  syms S(sin_theta)
%  y= asin(delta_phi*lambda/(-2*pi*dx));
%  S(sin_theta) = y;
%  DoA=rad2deg(S(sin(atan(y0/x0))));
%  disp(DoA)

%% 
% DoA estimation
syms S(sin_theta)
delta_phi=(2*pi*parameter.f0 /parameter.C)*(r_n(7)-r_n(6));
real_DoA=rad2deg(atan(y0/x0))
estimated_DoA=90+rad2deg(asin(delta_phi*lambda/(2*pi*dx)))
error=abs(abs(real_DoA)-abs(estimated_DoA))
%%
% 2.3 section 4:

% Number of Monte Carlo simulations
%num_simulations = 100;

% Range of noise powers to test
%noise_powers = logspace(-4, 0, 10);

% Initialize arrays to store results
%std_errors = zeros(size(noise_powers))

% received signal at array:
%array_response_true = exp(-2i*pi*parameter.f0*r_n/parameter.C);

% Range of noise powers to test
%noise_powers = logspace(-4, 0, 10);  % Adjust as needed

% Initialize arrays to store results
%std_errors = zeros(size(noise_powers));

%for k = 1:length(noise_powers)
  %  noise_power = noise_powers(k);
    
%    errors = zeros(1, num_simulations);
errors=[];

num_simulations = 500;
for i=1:num_simulations
        N=9;
        % repeat 100 time with 100 nosie power
        randdd=i*((randn(N, 1) + 1j * randn(N, 1)));

 %   for i = 1:num_simulations

    %    delta_phi_noisy=abs(ang((5))+abs(ang(4)));
        
        % MUSIC algorithm for DoA estimation
       % [~, DoAs] = music(Srx, M, lambda);

        % Calculate error in DoA estimation
        estimated_DoA=90+rad2deg(asin((delta_phi+randdd*lambda/(2*pi*dx))));

        errors=[errors abs(estimated_DoA-real_DoA)];
        std_errors=std(errors)
 %   end

end

%%
% plot error for 500 simulation
close all
figure
plot(std_errors/sqrt(num_simulations))

xlabel("noise power")
ylabel("std")
title("Standard Deviation for 500 simulation");
grid on
%%
% error plot for 500 simulation and a single element
figure
plot(errors(1,:)/500)
xlabel("noise power")
ylabel("error")
title("Error for 500 simulation and a single elemente");
grid on

%%
% plot std
% std_errors=[];
% x_sdt_axis=[-4 -3 -2 -1 0 1 2 3 4];
% for j=1:9
%     std_errors=[std_errors std(errors(j,:)/100) ];
% 
% 
% end

plot(x_sdt_axis, std_errors,'o--r');
title("standard deviation for all elements");
grid on
%%
        % Calculate standard deviation of the error
    %std_errors(k) = std(errors);


%end

%%

% Plot results
figure;
loglog(noise_powers, std_errors, 'o-', 'LineWidth', 2);
xlabel('Noise Power');
ylabel('Standard Deviation of DoA Estimation Error');
title('Monte Carlo Simulation: DoA Estimation with Noise');
grid on;



%%
% 2D Poositioning
% With respect to nyquist rate, we consider sampling frequency fs=4 GHz 
% in order to avoid aliasing. 

dt_n=r_n/parameter.C;                   % Corresponding delay related to each element of array

end_time= 0.5 * 10 ^ (-6);                                         
time=1/parameter.fs:1/parameter.fs:end_time;
%%

% Time of arrival calculation at rx
time_rx=[];
for i=1:numel(dt_n)
time_rx=[time_rx;time+dt_n(i)];
end

%% 
% Stx definition
syms t  g(t)

g(t)=(sinc(parameter.B*t))*(exp(2i*parameter.f0*t));
%%
% S matrix calculations

Srx=[];
Srx_demod=[];
Stx=g(time);

for i=1:numel(r_n)
Srx=[Srx; g(time_rx(i,:))];
Srx_demod=[Srx_demod;Srx(i,:).*conj(Stx)];
end

%%
% DFT Calculation


% rho = parameter.C/(2*parameter.B); % range of resolution 
% r_ax= 5:rho:10;
% Nf=length(r_ax);
% 
% 
% g=double(Srx_demod(1,:));
% [X,f]=my_dft(g,time,Nf);



%%
% 
% pow_vect=((abs(X)).^2)/Nf;
% figure(4444)
% plot(pow_vect(1,:))
% 
% 
% %%
% 
% [s_ddf,f1]=my_dft(abs(g(1,:)),time,Nf); % study the my_dft
% 
% figure(3333)
% subplot(2,1,1)
% plot(f1,abs(s_ddf(:,1)))
% xlabel("Frequency");
% ylabel("Amplitude"); 
% title "AmplitudeSpectrum"
% grid on;
% subplot(2,1,2)
% plot(f1,(angle(s_ddf(:,1))))
% xlabel("Frequency");
% ylabel("Phase"); 
% title "PhaseSpectrum"
% grid on;


%%
% our method:
%dx = lambda / 2
%phi_time=Srx_demod(6,:)-Srx_demod(5,:);
%[~, DoAs] = music(Srx_demod(6,:), 1, parameter.C / (2 * dx));


