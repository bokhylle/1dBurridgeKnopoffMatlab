clear; close all; clc
dt = 1e-3;
outputFreq = 1000;
tau_all = logspace(-3,-0.001,100);
alpha_all = logspace(-3,1.0,100);
[tau_all,alpha_all] = meshgrid(tau_all,alpha_all);
tau_all = tau_all(:); alpha_all = alpha_all(:);
tmax = Inf;
foldername = 'Results';
mkdir(foldername);

files = dir(foldername);
existing_files = [];
for filename = {files.name}
    existing_files = [existing_files,filename];
end
%%

for Nstop = 2.^[1:12]
    for runInd = 1:length(alpha_all)
        tau = tau_all(runInd)*ones(Nstop,1);
        tau_minus = tau+2;
        alpha =  alpha_all(runInd);
        filename = ['Nstop-' num2str(Nstop) '_alpha-' num2str(alpha) '_tau-' num2str(tau(1)) '.mat'];
        [x,v,a,stuck,unstickTime,t,frontTip,mean_x]=integrateSingleSystem(alpha,tau,Nstop,dt,outputFreq,tmax,'seismicMoment',tau_minus);
        save([foldername '/' filename],'x','v','a','stuck','unstickTime','t','frontTip','mean_x');
    end
end