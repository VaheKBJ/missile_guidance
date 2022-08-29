function [ppn, md, tf] = run_spdynamics(deltaErrDeg, rng0)
    ppn = linear_spdynamics(deltaErrDeg, rng0);
    [ppn, md, tf] = ppn.run;
end