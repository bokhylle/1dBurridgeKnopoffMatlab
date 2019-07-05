function integrateSystem(alpha_vec,tau_vec,N,filename,Nblock,dt,regularizationTerm)

[alpha_allRuns,tau_bar_allRuns]=ndgrid(alpha_vec,tau_vec);
alpha_allRuns = alpha_allRuns(:)';
tau_bar_allRuns = tau_bar_allRuns(:)';
Nsim = length(tau_bar_allRuns);

%Block level on gpu
x = zeros(N,Nblock,'double','gpuArray');
v = zeros(N,Nblock,'double','gpuArray');
a = zeros(N,Nblock,'double','gpuArray');
t = zeros(N,Nblock,'double','gpuArray');
stuck = true(N,Nblock,'gpuArray');
unstickTime = NaN(N,Nblock,'double','gpuArray');

%Save all unstick times in a separate matrix
unstickTime_allRuns = NaN(N,Nsim);
finished = zeros(1,Nblock);

simulationInd = 1:Nblock;

alpha = repmat(gpuArray(single(alpha_allRuns(simulationInd))),N,1);
tau_bar = repmat(gpuArray(single(tau_bar_allRuns(simulationInd))),N,1);

outputFreq = round(1/dt);
backupFreq = round(100/dt);

maxInd = simulationInd(end);

%Open matfile
m = matfile(filename,'Writable',true);
m.alpha_allRuns = alpha_allRuns;
m.tau_bar_allRuns = tau_bar_allRuns;
m.alpha_vec = alpha_vec;
m.tau_vec = tau_vec;
m.unstickTime_allRuns = unstickTime_allRuns;
m.N = N;
m.dt = dt;

%%
i = 0;
while true
    i=i+1;
    a = [x(1,:)+1-tau_bar(1,:);x(1:end-1,:)] -2*x-(regularizationTerm(v,tau_bar)+alpha.*v) + [x(2:end,:);x(end,:)] + tau_bar;
    unstickTime(a>=1 & stuck & isnan(unstickTime))=t(a>=1 & stuck & isnan(unstickTime));%Save first unsticking time:
    stuck(a>=1)=false;%Unstick if large enough force
    
    v(~stuck) = v(~stuck) + a(~stuck).*dt;
    stuck(v<0)=true;%Stick if negative velocity
    v(stuck)=0;%Set velocity to zero if stuck
    stuck(:,~isnan(unstickTime(end,:)))=1; %stick if all ruptured to reduce compuational cost
    x(~stuck) = x(~stuck) + v(~stuck).*dt;
    
    if(mod(i,outputFreq)==0)
        clc
        display([num2str(i) ': ' num2str(length(simulationInd)) ' simulations running, reached index ' num2str(min([maxInd,Nsim])) ' out of ' num2str(Nsim)])
    end
    
    if (mod(i,backupFreq)==0)
        
        finished( sum(stuck)==N ) = 1;
        finished( ~isnan(unstickTime(end,:)) ) = 1;
        
        ind = find(finished==1);
        if size(ind,2)>0
            display('Saving ...');
            unstickTime_allRuns(:,simulationInd(ind)) = gather(unstickTime(:,ind)); 
            m.unstickTime_allRuns = unstickTime_allRuns;
          
            unstickTime(:,ind) = NaN;
            t(:,ind) = 0;
            x(:,ind) = 0;
            v(:,ind) = 0;
            a(:,ind) = 0;
            stuck(2:end,ind) = 1;
            simulationInd(ind) = maxInd+1:maxInd+sum((finished));
            maxInd = maxInd+sum(finished);
            finished(ind)=0;
            
            removeInd = find(simulationInd>Nsim);
            
            simulationInd(removeInd)=[];
            unstickTime(:,removeInd) = [];
            t(:,removeInd) = [];
            x(:,removeInd) = [];
            v(:,removeInd) = [];
            a(:,removeInd) = [];
            stuck(:,removeInd) = [];
            finished(removeInd)=[];
            
            alpha = repmat(gpuArray(single(alpha_allRuns(simulationInd))),N,1);
            tau_bar = repmat(gpuArray(single(tau_bar_allRuns(simulationInd))),N,1);
            
        end
        if size(simulationInd,2)==0
            break;
        end
        
    end
    
    t = t+dt;
end

end

