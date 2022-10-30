# missile_guidance
 simple proportional navigation

this project is a simple simulation of a guidance system of a missile employing pure proportional navigation.

to create a ppn object write:

ppn = linear_ideal(1000, 10)

where 1000 is range in meters, and 10 is initial heading error (azimuth) to the target. 
to run the simulation:

[ppn, md, tf] = ppn.run

md is the last miss distance, tf is the flight time.
to plot the results run:

ppn.plotsim(1 : 5)

which draws 5 basic figures of the engagement.


for single trajectory run the file single_run.m. 
the target range can be modified by 'range'. 
the initial heading error with the target can be modified by 'initial_error'.











