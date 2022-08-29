


gcolor = colormap(gray(12)); 


figure

axes1 = axes('Parent', gcf ...
            , 'YGrid', 'on' ...
            , 'XGrid', 'on' ...
            , 'fontname', 'times' ...
            , 'fontsize', 14);
box(axes1, 'on');
hold(axes1, 'all');
p1 = plot(ppn.tspan, ppn.Xspan(3,:) / d2r ... 
        , 'Parent', axes1 ...
        , 'DisplayName', '\omega' ...
        , 'color', gcolor(1, :), 'linewidth', 2);
xlabel('Time (s)')
ylabel('\omega (deg/sec)')

axes2 = axes('Parent', gcf ...
            , 'YAxisLocation', 'right' ...
            , 'YColor', [0 0 0] ...
            , 'Color', 'none' ...
            , 'fontname', 'times' ...
            , 'fontsize', 14);
hold(axes2, 'all');
p2 = plot(ppn.tspan, ppn.Xspan(6, :) / 9.8, '--' ...
       ... , 'Parent', axes2 ...
        , 'DisplayName', 'a_m' ...
        , 'color', gcolor(6, :), 'linewidth', 2);
ylabel('a_m (g)')


legend([p1, p2], '\omega', 'a_m')












































