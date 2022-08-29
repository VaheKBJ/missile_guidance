function [ppn, md, tf] = run_ppn(deltaErrDeg, rng0, k)
    ppn = linear_1order(deltaErrDeg, rng0, k);
    [ppn, md, tf] = ppn.run;
end