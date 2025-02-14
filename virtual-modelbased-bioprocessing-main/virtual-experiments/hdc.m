%% Double component Monod-function
%   author        : Mirko Pasquini
%
%   The function takes as input a concentration, an activation and an
%   inhibition coefficients, and returns the value of the double component
%   Monod-type kinetics function. The case of activations and inhibitions
%   only are included by putting ki = 0 or ka = 0 respectively. The neutral
%   effect is obtained by putting ki = ka = 0.
%
%   dc = hdc(c,ka,ki)
%
%   @inputs
%       c : metabolite concentration (>= 0)
%       ka : activation coefficient (>= 0)
%       ki : inhibition coefficient (>= 0)
%
%   @outputs
%       dc : double component Monod-type function i.e.
%           dc = (c/(c+ka))*(1/(1+ki*c))

function dc = hdc(c,ka,ki)
    if c == 0
        c = 1e-5; % numerical issues with concentration = 0.
    end
    if ka == 0 && ki == 0
        dc = 1;
    else
        dc = (c./(c+ka)).*(1./(1+c*ki));    
    end
end
