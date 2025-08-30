%
% This code reproduces the simulations of the IK(Izhikevich) neuron network and the derived
% mean-field model descirbed in the following paper:
%
% R. Gast,S.A. Solla,& A. Kennedy, "Neural heterogeneity controls 
% computations in spiking neural networks," 
% Proc. Natl. Acad. Sci. U.S.A. 121 (3) e2311885121, 
% https://doi.org/10.1073/pnas.2311885121 (2024).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main parameters (including Lorentzian probability distribution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Regular-spiking neuron parameter
N=100; % total number of neurons 
C=100; % capacitane (pF)

v_r=-60; % resting potential (mV) 
v_0 =-100;  % reset potential (mV) %arbitrarily determined
% v_p % cutoff value; peak membrane potential(mV) %arbitrarily determined
v=v_r; % membrance potential (mV)
v_p=-40; % peak membrane potential (mV); threshold; arbitrary determined

% Applying within-population spike threshold heterogeneity/ each neuron has different one 
% We use Lorentzian probability distribution
mean_th=-40; % average threshold of all neurons (mV)
delta_v=0.5; % half-width-at-half-maximum (mV)
v_spk=mean_th+delta_v*tan(pi*(rand(N,1)-0.5)); % heterogeneity threshold vector (mV)
% Cauchy-Lorentzian PDF: 1/pi*(delta_v/((threshold-mean_threshold)^2+delta_v^2));

tau_u=33.33; % recovery variable time constant (ms)
tau_s=6; % synaptic filtering time constant (ms)
g=1; % synaptic conductance (nS)
b=-2; % scaling factor?? (nS)
kappa=10; % coupling strength?? (pA)

% Define J
J=15; % single synapse type of strength
M_J=(rand(N,N)<=0.2)*J; 
M_J=M_J/(N*0.2); % normalize

% spike & reset
v_spike = 1000;
v_reset = -1000;



k=0.7; % leakage parameter (nS/mV)
E=0; % reversal potential??(mV)
u=0; % global recovery variable % arbitrarily determined
s=0; % global post-synaptic activation term % arbitrarily determined
r=0; % average spike rate across all neurons % arbitrarily determined.. starting from 0

%% Fast-spiking neuron parameters
C=20; % capacitane (pF)
k=1; % leakage parameter (nS/mV)
v_r=-55; % resting potential (mV) 
mean_th=-40; % average threshold of all neurons (mV)
delta_v=0.2; % half-width-at-half-maximum (mV)
v_spk=mean_th+delta_v*tan(pi*(rand(N,1)-0.5)); % heterogeneity threshold vector (mV)
g=1; % synaptic conductance (nS)
E=-65; % reversal potential??(mV)
tau_u=5; % recovery variable time constant (ms)
tau_s=8; % synaptic filtering time constant (ms)
kappa=100; % coupling strength?? (pA)
b=0.025; % scaling factor?? (nS)
J=5; % single synapse type of strength

%% Low-threshold-spiking neuron parameters
C=100; % capacitane (pF)
k=1; % leakage parameter (nS/mV)
v_r=-56; % resting potential (mV) 
mean_th=-42; % average threshold of all neurons (mV)
delta_v=1; % half-width-at-half-maximum (mV)
v_spk=mean_th+delta_v*tan(pi*(rand(N,1)-0.5)); % heterogeneity threshold vector (mV)
g=1; % synaptic conductance (nS)
E=-65; % reversal potential??(mV)
tau_u=33.33; % recovery variable time constant (ms)
tau_s=8; % synaptic filtering time constant (ms)
kappa=20; % coupling strength?? (pA)
b=8; % scaling factor?? (nS)
J=5; % single synapse type of strength


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Single Population Neuron Network Simulation (SNN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a network of coupled Izhikevich (IK) neurons 

totalt=2500; % total simulation time (ms)
dt=0.01; % step size (ms)
timev=0:dt:totalt; % time vector

% input current
I=zeros(length(timev),1);
I(timev<=750|timev>=2000)=30;
I(timev>750&timev<2000)=50;


% matrix & vectors to store all values during the simultation
Matrix_v1=zeros(N,length(timev)); 
vector_r1=zeros(length(timev),1);
vector_u1=zeros(length(timev),1);
vector_s1=zeros(N,length(timev));

% set initial value
Matrix_v1(:,1)=v; vector_r1(1)=r;vector_u1(1)=u;vector_s1(:,1)=s; 
count=0; % counter

for t=timev
    if t>=totalt; break; end
    count=count+1;

    spk=Matrix_v1(:,count)>=v_spike;
    r_t=spk/dt;

    for i=1:1:N        
        vector_s1(i,count+1)=vector_s1(i,count)+(1/tau_s*(-vector_s1(i,count)+tau_s*(r_t(i))))*dt;
        % detect spikes/ need to revise
        if spk(i)==1; Matrix_v1(i,count)=v_spike; Matrix_v1(i,count+1)=v_reset; continue; end
        % Equation 1
        Matrix_v1(i,count+1)=Matrix_v1(i,count)+(1/C*(k*(Matrix_v1(i,count)-v_r)*(Matrix_v1(i,count)-v_spk(i))-vector_u1(count)+I(count)+g*((M_J(i,:)*vector_s1(:,count)))*(E-Matrix_v1(i,count))))*dt;     
    end
    
    % Equation 2-4
    % vector_r1(count+1)=sum(Matrix_v1(:,count)>=v_spk)/(N*dt);
    vector_u1(count+1)=vector_u1(count)+(1/tau_u*(-vector_u1(count)+b*(-v_r+sum(Matrix_v1(:,count),1)/N)+tau_u*kappa*(sum(r_t,1)/N)))*dt;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean-Field Equations  (M_F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We want to introduce within-population spike threshold heterogeneity.
% We assume that the spike thresholds v_threshold are distributed according
% to a Lorentzian probability distribution.


totalt=2500; % total simulation time (ms)
dt=0.01; % step size (ms)
timev=0:dt:totalt; % time vector

% input current
I=zeros(length(timev),1);
I(timev<=750|timev>=2000)=30;
I(timev>750&timev<2000)=50;

% vectors to store all values during the simultation
vector_v2=zeros(length(timev),1); 
vector_r2=zeros(length(timev),1);
vector_u2=zeros(length(timev),1);
vector_s2=zeros(length(timev),1);

% set initial value
vector_v2(1)=v;vector_r2(1)=r;vector_u2(1)=u;vector_s2(1)=s; 

count=0; % counter
spike_t2=zeros(N,2); % recent spike occurrence time vector

for t=timev
    if t>=totalt; break; end
    count=count+1;

    % define sigma-function
    if vector_v2(count)-v_r>=0; sigma_v=1; else; sigma_v=-1; end

    % Mean-field equations    
    vector_r2(count+1)=vector_r2(count)+1/C*((delta_v*k^2*sigma_v)/(pi*C)*(vector_v2(count)-v_r)+vector_r2(count)*(k*(2*vector_v2(count)-v_r-mean_th)-g*J*vector_s2(count)))*dt;
    vector_v2(count+1)=vector_v2(count)+1/C*(k*vector_v2(count)*(vector_v2(count)-v_r-mean_th)-pi*C*vector_r2(count)*(delta_v*sigma_v+pi*C*vector_r2(count)/k)+k*v_r*mean_th-vector_u2(count)+I(count)+g*J*vector_s2(count)*(E-vector_v2(count)))*dt;
    vector_u2(count+1)=vector_u2(count)+1/tau_u*(b*(vector_v2(count)-v_r)-vector_u2(count)+tau_u*kappa*vector_r2(count))*dt;
    vector_s2(count+1)=vector_s2(count)+1/tau_s*(-vector_s2(count)+tau_s*vector_r2(count))*dt;

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% external input for interaction of multiple populations
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I used the same constant value for the external input, which is "I".
% If we try to differenciate all external input by each populations, we can
% use the below variables.

n_p;% number of populations
Matrix_I=zeros(n_p,length(timev)); % input current of multiple populations
Matrix_s=zeros(n_p,length(timev)); % global post-synaptic activation term % arbitrarily determined
vector_g=zeros(n_p,1); % synaptic conductance (nS)
vector_E=zeros(n_p,1); % reversal potential(mV)

% Matrix_I = verctor_g * Matrix_s * (vector_E - Matrix_v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cauchy & Gaussian spike threshold distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cauchy Distribution
mean_th=-40; % average threshold of all neurons (mV)
delta_v=0.5; % half-width-at-half-maximum (mV)
x=mean_th-30:0.1:mean_th+30; % range of threshold % arbitraily determined
v_spk=mean_th+delta_v*tan(pi*(rand(N,1)-0.5)); % heterogeneity threshold vector (mV)
v_spk_pdf= 1/pi*(delta_v./((x-mean_th).^2+delta_v^2));

% Gaussian Distribution
mean_th=-40; % average threshold of all neurons (mV)
sigma_th=1; % arbitrary determined
v_spk_g=normrnd(mean_th,sigma_th,N,1);
v_spk_g_pdf=normpdf(x, mean_th, sigma_th); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1); clf
subplot(3,2,1);hold on;
plot(timev,mean(Matrix_v1(1:N,:)),'k','LineWidth',1)
plot(timev,vector_v2,'LineWidth',1.5);
legend('SNN', 'M-F')
xlabel('time (ms)'); ylabel('membrane potential (mV)')
title('Mean of neurons of one population')

subplot(3,2,2);
plot(timev,Matrix_v1(1,:))
xlabel('time (ms)'); ylabel('membrane potential (mV)')
title('Single neuron')

subplot(3,2,3);hold on;
plot(timev,mean(vector_s1,1)/(tau_s)*1000) 
plot(timev,vector_r2*1000,'LineWidth',1.5)
legend('SNN', 'M-F')
xlabel('time (ms)'); ylabel('r (Hz)')
title('Average Spike Rate Across All Neurons, r')

subplot(3,2,4);hold on;
plot(timev,vector_u1)
plot(timev,vector_u2,'LineWidth',1.5)
legend('SNN', 'M-F')
xlabel('time (ms)'); ylabel('u')
title('Global Recovery Variable, u')

subplot(3,2,5);hold on;
plot(timev,mean(vector_s1,1))
plot(timev,vector_s2,'LineWidth',1.5)
legend('SNN', 'M-F')
xlabel('time (ms)'); ylabel('s')
title('Global Post-Synaptic Activation Term, s')

subplot(3,2,6);
plot(timev,I)
xlabel('time (ms)'); ylabel('I')
title('Input Current, I')

%%%%%%%%%%%%%%%%%%%%%%
%% Distribution plot
%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf; 
subplot(1,2,1); hold on;
plot(x,v_spk_pdf,'r','linewidth',2);
histogram(v_spk,1000,'Normalization','pdf','FaceColor', 'b')
xlabel('threshold'); ylabel('Density');
xlim([mean_th-30,mean_th+30]); ylim([0,0.4])
grid on; title('Cauchy Distribution');


subplot(1,2,2); hold on;
plot(x,v_spk_g_pdf,'r','linewidth',2)
histogram(v_spk_g,30,'Normalization','pdf')
xlabel('threshold'); ylabel('Density');
xlim([mean_th-30,mean_th+30]); ylim([0,0.4])
grid on; title('Gaussian Distribution');


%%



