function [dict,idx] = mrf_dictionary_computation(U,T1,T2,B1,W0,SP,MS,EC,fas,ga,varargin)
% Easy function to do simulations for fingerprinting
%
% Everything should be in ms.
% Init: Tom Bruijnen @ 20170807

% Check input
if nargin > 10
    julia = 1;
    
else
    julia = 0;
end
if nargin < 10
    ga = 7;
end
if nargin < 9
    fas = [];
end
if nargin < 8
    EC = 0;
end
if nargin < 7
    MS = 50; % Maximum dephased state in phase-graphs
end
if nargin < 6
    SP = 1; % (Complex) slice profile array
end
if nargin < 5 
    W0 = 0; % Off-resonance
end
if nargin < 4
    B1 = 1; 
end
if nargin < 3
    T2 = linspace(20,500,25);
end
if nargin < 2
    T1 = linspace(200,3000,25);
end

% Extract flip angles + tr
tr = U.Sequence.Repetition_time;

% Pre-pulse check (still to do)
if strcmpi(U.Sequence.PrePulse,'none')
    INV = 1;
else
    INV = -1;
end

% Check if there is a text-file present with flip angles
fas = U.Sequence.Flip_angle;

% Check if sequence is balanced or spoiled
if strcmpi(U.Sequence.Scan_method,'Gradient_spoiled_gradient_echo') % spoiled
    if julia
        [dict,idx] = julia_dictionary_computation(U,T1,T2,B1,0,fas',SP,MS,INV);
    else
        spoiled = 1;
        fas = fas(1:size(U.RawKspaceData,2));
        [dict,idx] = extendedPhaseGraphs(T1,T2,B1,SP,fas',tr',INV,MS,spoiled,EC);
        dict(1,:)=[];
    end
else % If balanced uncorrect phase cycling
    if isfield(U,'RawKspaceData')        
        try
            fas = fas(1:size(U.RawKspaceData,5));
        catch
        end
    end
    sim.INV = INV;
    sim.TR = tr;
    sim.gamma = 4258;		% Hz/G.
    sim.Trf = 0.0005;		% 100 us "hard" RF pulse.
    sim.T1 = T1; % s
    sim.T2 = T2; % s
    sim.B1 = B1;
    sim.B0 = W0; % Hz
    sim.alpha = fas; % deg
    sim.N = numel(sim.alpha); % Number of rf pulses
    sim.Tpad = (sim.TR - sim.Trf)/2;	% Sec.
    sim.t = [sim.Tpad sim.Trf sim.Tpad];
    sim.goldenangle = ga;
    [dict,idx] = bloch_mrf(sim,EC,INV);
end

% END
end