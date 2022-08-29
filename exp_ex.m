t = 0 : .01 : 10;
tau2 = 1;
tgo = .5;
tau1 = .1;
taum = 0.4;

Ad = - 1 / taum;
bd = 1 / taum;
cd = 1;
n = length(cd);

theta  = 15 * d2r;
delta = 5 * d2r;
vm = 400;
vc = 200;
cT = cos(theta) / vc;
cm = cos(delta) / vc;
N = 3;
k = 10;

% Ad' * H + H * Ad = -4 * eye(n)
H = -4 * eye(1) / 2 / Ad;
ln = max(H);

alef = max(cm^2 * (cd * cd') + N^2 * vm^2 * ln^2 * (bd' * bd) * tau2^2 + 4 * tau2 ...
            , tau2^2 * (1 / vm^2 * (cd * cd') + k^2 * ln^2 * (bd' * bd)));
        % 5760004
        
% tf = 1 : 0.1 : 10;
% tf = 1e5 : 1e3 : 1e6;
tf = 0 : 100 : 300000;

%%
% gamma feedback asymptotic term.
%
vx = exp(-alef ./ tf);

%%
% orig no feedback asymptotic term.
%
alpha2 = 10e3;
alpha1 = 1;
% nfa = alpha2 / tau1 .* (tau2 ./ tf).^alpha1;
nfa = (tau2 ./ tf).^alpha1;

%%
% 
% 
figure
plot(tf, nfa, 'displayname', 'no fb: (\tau_2/t_f)^{\alpha_1}')
hold on
plot(tf, vx, 'displayname', '\gamma fb: e^{-\alpha/t_f}')
xlabel('t_f')
grid
legend('fontsize', 12)





% convergance test
figure
ex1 = exp(1 / tau1) * ones(size(tf));
ex2 = exp(1 / tau2) * ones(size(tf));
ex3 = exp(1 / tgo) * ones(size(tf));
plot(tf, ex1, 'displayname', '\tau_1')
plot(tf, ex2, 'displayname', '\tau_2')
plot(tf, ex3, 'displayname', 't_{go}')
grid
legend('fontsize', 12)

figure
alpha = [0.5 1 2 5 7 10 20];
tfa = meshgrid(tf, alpha);
extf = exp(-alpha' ./ tfa);
lgd=legend(strsplit(num2str(alpha)), 'fontsize', 12);
ttl = get(lgd, 'title');
set(ttl, 'string', '\alpha')
title('e^{-\alpha/t_f}', 'fontsize', 18)
xlabel('t_f')
grid







































