function options=optimoptions(solver,varargin)
% Octave compatibility for ORIGINAL HANKSOME scalar fsolve calls.
% This is NOT MATLAB's Levenberg-Marquardt algorithm. The original
% economic residual equation is unchanged and the scalar price root
% must be separately checked for value and branch.
if ~strcmpi(solver,'fsolve')
    error('Only original scalar fsolve options are supported');
end
if mod(numel(varargin),2)~=0
    error('Expected option name-value pairs');
end
options=optimset();
for k=1:2:numel(varargin)
    key=varargin{k};
    val=varargin{k+1};
    if strcmpi(key,'Display')
        options=optimset(options,'Display',val);
    elseif strcmpi(key,'Algorithm')
        if ~strcmpi(val,'Levenberg-Marquardt')
            error('Unsupported source scalar solver algorithm request');
        end
        % GNU Octave fsolve uses its own trust-region implementation.
        % Do not label it MATLAB's Levenberg-Marquardt or infer identical
        % root basins without comparing the resulting x numerically.
    else
        error('Unsupported original-source fsolve option: %s',key);
    end
end
end
