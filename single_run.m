clear classes
close all
range = 1000; % m
initial_err = 10; % deg
li = run_ppn(initial_err, range, 0);
li.plotsim(1:5)



