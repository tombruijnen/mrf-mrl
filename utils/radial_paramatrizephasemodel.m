function model_parameters = radial_paramatrizephasemodel(cph,angles)
% Fitting of k0-phase (phi) to model ang/cph
% angles = 1xN double (readout-angles)
% cph = 1xN double (central-phase)

% Model
phi=@(a,theta)(a(1)*cos(theta+a(4))+a(2)*sin(theta+a(4))+a(3)); 

% Initial guess
A0=[1,1,1,0];

% Perform initial fit
[A,~,~,~,MSE0]=nlinfit(angles,cph,phi,A0);
MSE = MSE0;
A_best = A;

% If residuel is too large it is probably wrapped
% Repeated fit with addition of pi to all phases
iter=0;
while MSE > 0.1 && iter<10
    % Increment cph
    cph=cph+pi/3;
    
    % If phase > pi put it back in bounds [-pi,pi]
    for s=1:numel(cph)
        if cph(s)>pi
            cph(s)=-pi+(cph(s)-pi);
        end
    end
    
    % Repeat the fitting 
    [A,~,~,~,MSE]=nlinfit(angles,cph,phi,A0);
    if MSE < MSE0
        MSE0 = MSE;
        A_best = A;
    end
    
    % Track iteration
    iter=iter+1;
end

if iter == 10
    A = A_best;
end

%figure(7);scatter(angles,cph);hold on;plot(sort(angles),phi(A,sort(angles)));pause();hold off

% Return required model parameters
model_parameters=A(1:2);

% END
end
