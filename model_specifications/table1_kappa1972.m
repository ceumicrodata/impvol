% This file describes a specific model specification for
%                       include ref here

% There are three main types of parameters that describe a model:
% 1. Model type (I-O links, Balanced trade, Labor adjustment)
% 2. Fundamental parameters (theta, eta)
% 3. Counterfactual parameters (Trade costs, Productivity shocks)



%% Model type
% io links (bool)
spec.iol = 1;

% Balanced trade (bool)
spec.bt = 0;

% Labor adjustment cost (0: infinity, int: rho)
spec.lac = 0;



%% Fundamental parameters
% theta in {2, 4, 8}
spec.th = 4;

% eta (int?)
spec.et = 4;



%% Conterfactual parameters
% trade cost (0: actual calibrated, 1: 1972, 2: free trade)
spec.tc = 1;

% productivity shock type (0: actual calibrated,
%                          1: no sectoral shocks,
%                          2: no sectoral and residual shocks)
spec.sh = 0;
