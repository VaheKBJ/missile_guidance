clear classes
% close all
% li = run_ppn(10, 1e3, 0);
% li.plotsim(1:5)
% 
% 
% %%
% codegen run_ppn -args {0, 0, 0}

%%
% sim 
% 
r = [50 : 50 : 950, 1000 : 100 : 5000, 5500 : 500 : 10000, 11000 : 1000 : 20000];

%%
% no feedbck original system
%
los = zeros(size(r));
wlos = zeros(5 / PPN.Ts, size(r, 2));
md = zeros(size(r));
tf = zeros(size(r));
nrm_wlos = zeros(size(r));

tau1 = 0.5;
tau2 = 1;
    
% hwait	=	waitbar(0, 'Please wait...');
hwait = waitbar(0,'Please wait...', 'windowstyle', 'modal');
frames = java.awt.Frame.getFrames();
frames(end).setAlwaysOnTop(1);

for i = 1 : length(r)
    
    waitbar(i / length(r), hwait)
    
    [li, md(i), tf(i)] = run_ppn(15, r(i), 0);
    
    t1 = tf(i) - tau2;
    t2 = tf(i) - tau1;
    
    if t1 > 0
        dt1 = li.tspan - t1;
        dt2 = li.tspan - t2;
        
        ind1 = find(dt1 < li.Ts & dt1 >= 0);
        ind2 = find(dt2 < li.Ts & dt2 >= 0);

        wlos(1 : length(li.tspan), i) = li.Xspan(3, :);
        nrm_wlos(i) = sum(abs(wlos(ind1 : ind2, i)).^2)^(1/2);
        
        am(1 : length(li.tspan), i) = li.Xspan(6, :);
        nrm_am(i) = sum(abs(am(ind1 : ind2, i)).^2)^(1/2);
    end
end
close(hwait) 

nrm_wlos_nan = nrm_wlos;
nrm_wlos_nan(nrm_wlos_nan == 0) = nan;
nrm_am_nan = nrm_am;
nrm_am_nan(nrm_am_nan == 0) = nan;

gcolor = colormap(gray(12)); 
figure
plot(tf, nrm_wlos_nan, 'DisplayName', '||\omega||', 'color', gcolor(1, :), 'linewidth', 2);
hold
plot(tf, nrm_am_nan, 'DisplayName', '||a_m||', 'color', gcolor(6, :), 'linewidth', 2);
title('||x(\tau_1, \tau_2)||', 'fontname', 'times', 'fontsize', 14)
xlabel('t_f', 'fontname', 'times', 'fontsize', 14)
set(gca, 'fontname', 'times', 'fontsize', 14);
grid
ylim([0, 1])
legend

% % figure
% % plot(wlos)
% figure('windowstyle', 'docked')
% wlosnan=wlos;
% wlosnan(wlosnan == 0) = nan;
% plot(linspace(PPN.Ts, (5 / PPN.Ts) * PPN.Ts, 5 / PPN.Ts), wlosnan)
% grid
% hold on
% % plot(tf, abs(los), 'DisplayName', 'los')
% plot(tf, abs(wlos), 'DisplayName', 'wlos')

% 
% %%
% % feedback
% %
% losfb = zeros(size(r));
% wlosfb = zeros(size(r));
% mdfb = zeros(size(r));
% tffb = zeros(size(r));
% for i = 1 : length(r)
%     [li, mdfb(i), tffb(i)] = run_ppn_mex(10, r(i), 1);
%     losfb(i) = li.Xspan(2, end);
%     wlosfb(i) = li.Xspan(3, end);
% end
% plot(tffb, abs(mdfb), 'DisplayName', 'fb md')
% % ylim([-0.2 0.2])
% hold on
% % plot(tffb, abs(losfb), 'DisplayName', 'fb los')
% plot(tffb, abs(wlosfb), 'DisplayName', 'fb wlos')
% legend('location', 'best')
% % set(gca, 'YScale', 'log')
% xlabel('flight time')
% grid
% title('x(t_f) - linear ppn')
% 
% 
% %% 
% % theoretical
% % 
% asym = tf.^(21.6e6) .* exp(-1 ./ tf);
% figure
% plot(tf, asym)
% 
% 
% 
% 





