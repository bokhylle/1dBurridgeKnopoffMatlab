clear; close all;

dt = 1e-3;
outputFreq = 1e9;
alpha_startStop = [logspace(-.5,0,1000)];
tau_startStop = NaN(size(alpha_startStop));
errorCriterion = 1e-3;
N = 1000;

ind = 0;
for alpha = alpha_startStop
    ind = ind+1
    tau_min = 1e-3;
    tau_max = 0.2;
    
    if alpha > 1
        N = 100;
    end
    if alpha>1.1
        N = 10;
    end
   
    startStop_min=detectStartStop(alpha,tau_min,N,dt,outputFreq);
    startStop_max=false;%detectStartStop(alpha,tau_max,N,dt,outputFreq);

    tau_startStop_Local = NaN;
    
    if startStop_min %Only calculate if the minimum value of tau gives a start stop
        
        while true
            tau_midpoint = (tau_max+tau_min)/2;
            
            if startStop_min ~= startStop_max %Calculate midpoint
                startStop_midpoint = detectStartStop(alpha,tau_midpoint,N,dt,outputFreq);
            end
            
            if startStop_min ~= startStop_midpoint
                tau_max = tau_midpoint;
                startStop_max = startStop_midpoint;
            else
                tau_min = tau_midpoint;
                startStop_min = startStop_midpoint;
            end
            
            error = (tau_max-tau_min)/(tau_max+tau_min);
            if error < errorCriterion
                tau_startStop_Local = (tau_max+tau_min)/2;
                break
            end
            
        end
        
    end
    
    tau_startStop(ind) = tau_startStop_Local;
end

save('startStopCurve.mat','alpha_startStop','tau_startStop');