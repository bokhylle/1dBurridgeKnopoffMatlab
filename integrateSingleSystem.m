function [x_all, v_all, a_all, stuck_all, unstickTime, t_all,frontTip,mean_x] = integrateSingleSystem(alpha,tau_bar,N,dt,outputFreq,tmax,saveMode,tau_bar_minus)
%saveMode = 'all';
%saveMode = 'seismicMoment';
%saveMode = 'arresting';

x = zeros(N,1);
v = zeros(N,1);
a = zeros(N,1);
t = zeros(N,1);
stuck = true(N,1);
unstickTime = NaN(N,1);

Nsave = 2e5;
switch saveMode
    case {'all','arresting'}
        x_all = NaN(N,Nsave);
        v_all = NaN(N,Nsave);
        a_all = NaN(N,Nsave);
        stuck_all = NaN(N,Nsave);
        t_all = NaN(1,Nsave);
        frontTip = NaN;
        mean_x = NaN;
    case 'seismicMoment'
        x_all = NaN;
        v_all = NaN;
        a_all = NaN;
        stuck_all = NaN;
        t_all = NaN(1,Nsave);
        frontTip = NaN(1,Nsave);
        mean_x = NaN(1,Nsave);
        
    otherwise
        error('savemode not given')
end

%%
i = 0;
while t<tmax
    i=i+1;
    v_old = v;
    stuck_old = stuck;
    
    tau = tau_bar;
    tau(sign(v)==-1)=tau_bar_minus(sign(v)==-1);
    
    a = [x(1,:)+1-tau_bar(1,:);x(1:end-1,:)] -2*x-alpha.*v + [x(2:end,:);x(end,:)] + tau;
    unstickTime(a>=1 & stuck & isnan(unstickTime))=t(a>=1 & stuck & isnan(unstickTime));%Save first unsticking time:
    
    stuck(a>=1)=false;%Unstick if large enough force
    stuck(a<=-1)=false;%Unstick if large enough force
    
    switch saveMode
        case {'seismicMoment','arresting'}
            stuck(end)=true;%last block always stuck for computing seismic moment
            unstickTime(end)=NaN;
    end
    
    v(~stuck) = v(~stuck) + a(~stuck).*dt;
    stuck((sign(v)~=sign(v_old) )& stuck_old==false)=true;
    v(stuck)=0;%Set velocity to zero if stuck
    x(~stuck) = x(~stuck) + v(~stuck).*dt;
    
    finished( sum(stuck)==N ) = 1;
    finished( ~isnan(unstickTime(end,:)) ) = 1;
    
    if finished
        break;
    end
    
    if(mod(i,outputFreq)==0)
        display([num2str(t(1))])
        switch saveMode
            case {'all','arresting'}
                x_all(:,i/outputFreq) = x;
                v_all(:,i/outputFreq) = v;
                a_all(:,i/outputFreq) = a;
                stuck_all(:,i/outputFreq) = stuck;
                t_all(i/outputFreq) = t(1);
                
                switch saveMode
                    case 'arresting'
                        if(~isnan(unstickTime(end-1)))
                            if(sqrt(mean(v.^2))<1e-2*tau_bar(1)/alpha)
                                finished = 1;
                            end
                        end
                end
                
            case {'seismicMoment'}
                t_all(i/outputFreq) = t(1);
                frontTip(i/outputFreq) = find(~isnan(unstickTime)==1,1,'last');
                mean_x(i/outputFreq) = mean(x(1:frontTip(i/outputFreq)));
                if(~isnan(unstickTime(end-1)))
                    if(sqrt(mean(v.^2))<1e-2*tau_bar(1)/alpha)
                        finished = 1;
                    end
                end
            otherwise
                
        end
    end 
    t = t+dt;
end

finalInd = find(isnan(t_all)==0); finalInd = finalInd(end);

switch saveMode
    case {'all','arresting'}
        x_all = x_all(:,1:finalInd);
        v_all = v_all(:,1:finalInd);
        a_all = a_all(:,1:finalInd);
        stuck_all = stuck_all(:,1:finalInd);
        t_all = t_all(1:finalInd);
    case 'seismicMoment'
        t_all = t_all(1:finalInd);
        frontTip = frontTip(1:finalInd);
        mean_x = mean_x(1:finalInd);
end

end

