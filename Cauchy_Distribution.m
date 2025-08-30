% Cauchy distribution & PDF
%
% We create a histogram using the Cauchy distribution and put it in the 
% figure of the PDF equation of Cauchy distribution to check whether 
% it follows the theoretical PDF shape.


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main parameter & Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_0=0; % mean of all x; arbitrary setting
gamma=1; % scale-factor; arbitrary setting

x=-50:0.1:50; % x-axis vector
u=rand(length(x),1); % sample from unifrom distribution from [0,1]
 
Cauchy_D=x_0+gamma*(tan(pi*(u-1/2))); % Cauchy distribution
Cauchy_PDF=1/pi*(gamma./((x-x_0).^2+gamma^2)); % PDF of Cauchy distribution


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf; hold on;
histogram(Cauchy_D,1000,'Normalization','pdf') % 
plot(x,Cauchy_PDF,'r','linewidth',1)
xlabel('x'); ylabel('Density')
legend('Derived Equation','Theory PDF')
xlim([-150,150])
ylim([0,0.4])
grid on;hold off;








