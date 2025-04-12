function [I] = PN_model(params, X) % params = [I0, gam, Eg, n]
    vv = X(1, :);
    tt = X(2, :);

    k = 1.380649e-23; % J/K
    Vt = 26e-3; % V
    Is = params(1) .* tt.^params(2) .* exp(-params(3)./(params(4).*k.*tt)); % A
    I = Is .* (exp(vv./(params(4).*Vt)) - 1); % A
end

volt = [-0.5, 0, 0.5, 1, 1.5, 2]; % V
temp = [300, 310, 320, 330, 340, 350]; % K

X = [volt; temp] % V, K
X(1, :)
params = [1e-12, 1.5, 1.117, 2]; % I0, gam, Eg, n


I = PN_model(params, X)

size(X)