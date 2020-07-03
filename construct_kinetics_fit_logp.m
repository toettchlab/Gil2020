function fcn = construct_kinetics_fit_logp(conc, t1, x1, t2, x2, toPlot)
if nargin < 6 || isempty(toPlot)
    toPlot = 0;
end

w1 = 1; % how much to weight association phase SSE
w2 = 1; % how much to weight dissociation phase SSE

% relative time to t=0 at start of the condition
t1 = t1-t1(1);
t2 = t2-t2(1);

% initialize modeled data points to a matrix of zeros
x1m = zeros(size(x1));
x2m = zeros(size(x2));

% return a function handle to the fitting function
fcn = @fit_fcn;

    function SSE = fit_fcn(p)
        % Extract named parameters
        ka  = 10^p(1); % association rate kon
        kd  = 10^p(2); % dissociation rate koff
        k2  = 10^p(3); % decay rate k2
        
        
        % Load blocks of 4 parameters for each set of data
        c = 0;
        while 1
            try
            c = c + 1;
            aon(c)  = p(4*(c-1)+4);
            bon(c)  = p(4*(c-1)+5);
            aoff(c) = p(4*(c-1)+6);
            boff(c) = p(4*(c-1)+7);
            catch, break, end
        end
        
        x1m = zeros(size(x1));
        x2m = zeros(size(x2));
        
        for i = 1:size(x1,2)
            % Model for ASSOCIATION phase
            x1m(:,i) = (aon(i)*(1-exp(-(ka*conc(i) + kd)*t1)) + bon(i)).*exp(-k2*t1);
        end
        
        for i = 1:size(x2,2)
            % Model for DISSOCIATION phase
            x2m(:,i) = (aoff(i)*exp(-kd*t2) + boff(i)).*exp(-k2*t2);
        end
        
        % Overall sum of squared error
        SSE = w1*sum((x1m(:)-x1(:)).^2) + w2*sum((x2m(:)-x2(:)).^2);
        
        if toPlot
            figure(2)
            clf
            subplot(2,1,1)
            plot(t1, x1)
            hold on
            ax = gca;
            ax.ColorOrderIndex = 1;
            plot(t1, x1m, '--')
            title(sprintf('k_{on} = %0.3g; k_{off} = %0.3g; K_D = %0.3g', 10^p(1), 10^p(2), 10^p(2)/10^p(1)))
            subplot(2,1,2)
            plot(t2, x2)
            hold on
            ax = gca;
            ax.ColorOrderIndex = 1;
            plot(t2, x2m, '--')
            SSE = 0;
            drawnow
        end
    end
end