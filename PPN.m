classdef PPN 

properties (Constant)
    Ts          	=	0.005;                  % Sampling Period  % 0.005
    vT          	=	200;                    % Target    Velocity        %m/sec
    aT          	=	0;                      % Target    Acceleration    %m/sec^2
    vM          	=	400;                    % Missile   Velocity        %m/sec
    TheTDeg     	=	15;                     % Target Velocity Angle     %Deg
    N            	=	3;                      % PN Gain Dimensionless
    AccMaxG        	=	30;                     % g
end

properties (Access = protected) 
    bf  = 1;
    rr  = 0;  
    tk  = 0;
    sinThet 
    cT 
    cM 
    tf 
    vC 
    RmtRef0
    RmtRef
    true_lambda
    lambda_tf
    Algk  
    tauM = 0.4;
end

properties 
    Xref
    Xk
    tspan
    Xspan
    Algspan
 end

methods (Abstract) 
    obj = run(obj) 
%     dx  = dX(obj, t, x)
end
    
methods 
    function delete(obj)
        warning(['PPN class is deleted ' disp( class(obj))])
    end
    function plotsim(obj, figsVector)
            
        cmf             = figure(1000);
        gcolor          = colormap(gca, gray(12)); 
        close(cmf);
        lineType.color  = gcolor(1, :);
        lineType.width  = 2;
        lineType.type1  = '-';
        lineType.type2  = '--';
        lineType.type3  = ':';
        ltitle = '';
        
        figN = 1;  %  trajectory       
        if ismember(figN, figsVector)
            figure(figN) % traj

            p.traj = plot(obj.Xspan(1, :) .* cos(obj.Xspan(2, :)), obj.Xspan(1, :) .* sin(obj.Xspan(2, :)) ...
                        , lineType.type1, 'linewidth', lineType.width, 'color', lineType.color ...
                        , 'DisplayName', ltitle);

            title('Relative Motion Trajectory', 'fontname', 'times', 'fontsize', 12)
            xlabel('X_{Rel} (m)', 'fontname', 'times', 'fontsize', 12)
            ylabel('Y_{Rel} (m)', 'fontname', 'times', 'fontsize', 12)
            set(gca, 'fontname', 'times', 'fontsize', 12)
            grid on
            box on
            hold on

            if exist('savedir', 'var') && ~isempty(savedir)
                ltitle = strrep(strrep(strrep(ltitle, '\', ''), '.', ''), ' ', '');
                saveas(1, [savedir '\traj' ltitle '.fig'])
                saveas(1, [savedir '\traj' ltitle '.jpg'])
            end

        end

        figN = 2;  %  range vs time
        if ismember(figN, figsVector)
            figure(figN) % rho

            p.rho = plot(obj.tspan, obj.Xspan(1, :) ...
                        , lineType.type1 ,'linewidth', lineType.width, 'color', lineType.color ...
                        , 'DisplayName', ltitle);

            title('Missile-Target Range', 'fontname', 'times', 'fontsize', 12)
            xlabel('Time (sec)', 'fontname', 'times', 'fontsize', 12)
            ylabel('R_{MT} (m)', 'fontname', 'times', 'fontsize', 12)
            set(gca, 'fontname', 'times', 'fontsize', 12)
            grid on
            box on
            hold on

            if exist('savedir', 'var') && ~isempty(savedir)
                ltitle = strrep(strrep(strrep(ltitle, '\', ''), '.', ''), ' ', '');
                saveas(2, [savedir '\range' ltitle '.fig'])
                saveas(2, [savedir '\range' ltitle '.jpg'])
            end
        end

        figN = 3;  %  angles vs time
        if ismember(figN, figsVector)
            figure(figN) % angles
            p.lambda = plot(obj.tspan, obj.Xspan(2,:) / d2r ...
                        , lineType.type1 ,'linewidth', lineType.width, 'color', lineType.color ...
                        , 'DisplayName', ltitle);
            hold on 
            p.gammaT = plot(obj.tspan, obj.Xspan(4,:) / d2r ...
                        , lineType.type2 ,'linewidth', lineType.width, 'color', lineType.color ...
                        , 'DisplayName', ltitle);

            p.gammaM = plot(obj.tspan, obj.Xspan(5,:) / d2r ...
                        , lineType.type3 ,'linewidth', lineType.width, 'color', lineType.color ...
                        , 'DisplayName', ltitle);

            title('System Angles', 'fontname', 'times', 'fontsize', 12)
            xlabel('Time (sec)', 'fontname', 'times', 'fontsize', 12)
            ylabel('(deg)', 'fontname', 'times', 'fontsize', 12)
            set(gca, 'fontname', 'times', 'fontsize', 12)
            legend('\lambda', '\gamma_T', '\gamma_m')
            box on
            grid on

            if exist('savedir', 'var') && ~isempty(savedir)
                ltitle = strrep(strrep(strrep(ltitle, '\', ''), '.', ''), ' ', '');
                saveas(3, [savedir '\angles' ltitle '.fig'])
                saveas(3, [savedir '\angles' ltitle '.jpg'])
            end
        end

        figN = 4;  %  los rate vs time 
        if ismember(figN, figsVector)
            figure(figN) % omega
            p.omega = plot(obj.tspan, obj.Xspan(3,:) / d2r ...
                        , lineType.type1 ,'linewidth', lineType.width, 'color', lineType.color ...
                        , 'DisplayName', ltitle);
            title('\omega', 'fontname', 'times', 'fontsize', 12)
            xlabel('Time (sec)', 'fontname', 'times', 'fontsize', 12)
            ylabel('(deg/sec)', 'fontname', 'times', 'fontsize', 12)
            set(gca, 'fontname', 'times', 'fontsize', 12)
            grid on
            box on
            hold on

            if exist('savedir', 'var') && ~isempty(savedir)
                ltitle = strrep(strrep(strrep(ltitle, '\', ''), '.', ''), ' ', '');
                saveas(4, [savedir '\omega' ltitle '.fig'])
                saveas(4, [savedir '\omega' ltitle '.jpg'])
            end
        end

        figN = 5;  %  acceleration vs time
        if ismember(figN, figsVector)
            figure(figN) % acceleration
            p.ac = plot(obj.tspan, obj.Algspan(1, :) / 9.8 ... ...
                    , lineType.type1, 'color', lineType.color, 'linewidth', lineType.width);
            hold on
            if size(obj.Xspan, 1) == 6
                p.am = plot(obj.tspan, obj.Xspan(6, :) / 9.8, lineType.type2, 'color', lineType.color, 'linewidth', lineType.width);
                legend('a_c', 'a_m')
            end    
            title('Missile Acceleration', 'fontname', 'times', 'fontsize', 12)
            xlabel('Time (sec)', 'fontname', 'times', 'fontsize', 12)
            ylabel('Acc (g)', 'fontname', 'times', 'fontsize', 12)
            grid on
            box on
            set(gca, 'fontname', 'times', 'fontsize', 12)

            if exist('savedir', 'var') && ~isempty(savedir)
                ltitle = strrep(strrep(strrep(ltitle, '\', ''), '.', ''), ' ', '');
                saveas(5, [savedir '\acc' ltitle '.fig'])
                saveas(5, [savedir '\acc' ltitle '.jpg'])
            end
        end
    end
end
end