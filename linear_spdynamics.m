classdef linear_spdynamics < PPN
    properties 
        K   % state feedback gains
        Ad  % dynamics state matrix
        bd  % dynamics input vector 
    end
    
methods 
    
    function obj = linear_spdynamics(deltaErrDeg, rng0)
        %
        %   Input Initialization
        %
        xT0       	=	rng0; 
        yT0       	=	0;
        xM0       	=	0;                                 
        yM0       	=	0;                                    
        deltaErr  	=	deltaErrDeg * d2r;

        %       Reference System
        % Initial Values of State Vector and Algebraic Variables
        %       X   = [1-Rmt, 2-Lam, 3-Omg, 4-GamT, 5-GamM, 6-deltaR, 7-q, 8-aM]
        % Algebraic Variables
        %       Alg = [1-aMc]
        obj.RmtRef0 =	sqrt((xT0 - xM0)^2 + (yT0 - yM0)^2);   
        LamRef     	=	asin((yT0 - yM0) / obj.RmtRef0);
        OmgRef     	=	0;
        ThetRef    	=	obj.TheTDeg * d2r;
        GamTRef     =	LamRef + ThetRef; 
        DelCorrect	=	asin(obj.vT / obj.vM * sin(GamTRef - LamRef));
        DelRef   	=	DelCorrect;
        GamMRef    	=	LamRef + DelRef;
        aMRef    	=	0;
        deltaR      =   0;
        q           =   0;
        obj.Xref    =	[obj.RmtRef0; LamRef; OmgRef; GamTRef; GamMRef; deltaR; q; aMRef];

        % Reference System Constants
        obj.vC      =	-(obj.vT * cos(ThetRef) - obj.vM * cos(DelRef));
        obj.tf      =	obj.RmtRef0 / obj.vC;
        obj.sinThet =   sin(ThetRef);
        obj.cT      =   cos(ThetRef) / obj.vC;
        obj.cM      =   cos(DelRef) / obj.vC;

        %       Linear System
        % Initial Values of State Vector and Algebraic Variables
        %       X   = [1-Rmt, 2-Lam, 3-Omg, 4-GamT, 5-GamM]
        % Algebraic Variables
        %       Alg = [1-aMc]
        Rmt0      	=	0;
        Lam0      	=	0;
        Omg0      	=	(obj.vT * sin(ThetRef) - obj.vM * sin(DelCorrect + deltaErr)) / obj.RmtRef0;     %X(3)
        GamT0     	=	0;
        GamM0     	=	deltaErr;
        deltaR0     =   0;
        q0          =   0;
        aM0         =   0;
        aMc0      	=	obj.N * obj.vM * Omg0;     %Alg(1)  

        if abs(aMc0 / 9.8) >= obj.AccMaxG, aMc0 = obj.AccMaxG * 9.8 * sign(aMc0); end 

        X0        	=	[Rmt0; Lam0; Omg0; GamT0; GamM0; aM0]; 
        Alg0      	=	aMc0;

        obj.RmtRef  =	obj.RmtRef0;
        obj.Xk      =	X0;
        obj.Algk    =	Alg0;
        
    end

    function [obj, miss, t_fin] = run(obj) 
        n         	=	length(obj.Xk);
        nAlg      	=	length(obj.Algk);
        Nspan_est 	=	round(2 * obj.RmtRef / obj.vC / obj.Ts);

        obj.tspan   =	zeros(1, Nspan_est);
        Xrefspan    =	zeros(n, Nspan_est); 
        Xlinspan    =	zeros(n, Nspan_est); 
        Algspan   	=	zeros(nAlg, Nspan_est);

        OptODE = odeset('AbsTol', 1e-1, 'RelTol', 1e0);
        %
        %  Simulation Cycle
        %
        kk        	=	0;
%         hwait     	=	waitbar(0, '...');
%         
        while obj.bf        
%             waitbar((obj.RmtRef0 - obj.RmtRef) / obj.RmtRef0, hwait)

            %
            kk            	=	kk + 1;
            obj.tspan(1, kk)  	=	obj.tk;
            Xrefspan(:, kk) =	obj.Xref;
            Xlinspan(:, kk) =	obj.Xk;     % actually saves Xk in k+1
            Algspan(:, kk)	=	obj.Algk;

            %
            % Reference System Solution
            % 
            obj.RmtRef      =   obj.vC * (obj.tf - obj.tk);
            obj.Xref(1)    	=	obj.RmtRef;

            %
            % Linear System Solution
            % 
            tspanode = [obj.tk, obj.tk + obj.Ts / 2, obj.tk + obj.Ts];
%             Xspank = zeros(size(tspanode))';
            xode = obj.Xk; 
            [~, Xspank] =	ode45(@(tt, xx) linear_1order.dX(tt, xx, obj), tspanode, xode, OptODE); % 
            obj.Xk      =	Xspank(end, :)'; 


            % X   = [1-Rmt, 2-Lam, 3-Omg, 4-GamT, 5-GamM, 6-aM]
            Rmt           	=	obj.Xk(1);             
%                 Lam           	=	obj.Xk(2);
            Omg           	=	obj.Xk(3);
%                 GamT          	=	obj.Xk(4);                                    
            GamM          	=	obj.Xk(5);             
            aM          	=	obj.Xk(6);             

            % Alg = [1-aMc] 
            
            aMc           	=	obj.N * obj.vM * Omg;                           %Alg(1)    
            
            %%
            %
            %
            if abs(aMc / 9.8) >= obj.AccMaxG, aMc = obj.AccMaxG * 9.8 * sign(aMc); end     
            obj.Algk        =	aMc;

            % range rate for miss distance
            rdot            =   (abs(obj.RmtRef + Rmt) - abs(Xrefspan(1, kk) + Xlinspan(1, kk))) / obj.Ts;
            if obj.rr == 0
                if rdot < 0 
                    obj.rr = 1;     % first time of negative range rate
                end
            else
                if rdot >= 0
                    obj.bf = 0;     % non-negative range rate after period of negative range rate
                end
            end
            %
            obj.tk            	=	obj.tk + obj.Ts; 
        end  %while
        
%         close(hwait) 
        kk              =	kk - 1;
        obj.tspan       =	obj.tspan(1, 1 : kk);
        t_fin        	=	obj.tspan(1, kk);

        Xrefspan     	=	Xrefspan(:, 1 : kk); 
        Xlinspan     	=	Xlinspan(:, 1 : kk);
        obj.Xspan       	=	Xrefspan + Xlinspan;
        obj.Algspan      	=	Algspan(:, 1 : kk);
        
        % miss distance

        miss          =   obj.Xspan(1, end) * sin(obj.Xspan(2, end));%;abs()
%         fprintf('\tMis dist.=%6.5g (m)   t_fin=%6.5g (sec)\n', miss, t_fin)
%         % 
% 
%         % save struct
%         obj.Xs
%         X.tspan      	=	tspan;
%         X.Xspan      	=	Xtotal;
%         X.Algspan    	=	Algspan; 
%         X.deltaErrDeg	=	deltaErrDeg; 
%         X.rho0          =   rng0;
%         X.t_fin         =   t_fin;
%         X.miss          =   miss;
%         X.lambda_tf     =   Xtotal(2, end);
%         X.true_lambda   =   get_lambda_tf(LamRef, ThetRef, deltaErr, vM, vT, N);
        
    end
    
    function plotsim(obj, figsVector)
        plotsim@PPN(obj, figsVector);
    end
end

methods (Static)
    function dx = dX(t, x, obj)     %

    %       for Linear Relative Motion Model 
    %                 with 1st order Dynamics
    %--------------------------------------------------------------------------------------------------------------------------
    % Model State Vector 
    % 		X   = [1-Rmt, 2-Lam, 3-Omg, 4-GamT, 5-GamM, 6-deltaR, 7-q, 8-aM]
    % Algebraic Variables
    % 		Alg = [1-aMc] 
    %-------------------------------------------------------------

        dx   	=	zeros(size(x));

        %
        %  X   = [1-Rmt, 2-Lam, 3-Omg, 4-GamT, 5-GamM, 6-amM]
%             rho     =   x(1);
%             lambda  =   x(2);
        omega   =   x(3);
        gammaT  =   x(4);
        gammaM  =   x(5);
%         deltaR  =   x(6);
        q       =   x(7);
        aM      =   x(8);
        %  Alg = aMc
%         k = 1;
        aMc  	=	obj.Algk(1);    % obj.N * obj.vM * Omg


        %deriv's

        dx(1)	=	-obj.vT * obj.sinThet * gammaT + obj.vT * obj.sinThet * gammaM;
        dx(2)	=	omega;

        tgo = max(obj.tf - t, 0.05);
        dx(3)	=	2 / tgo * omega + obj.cT / tgo * obj.aT - obj.cM / tgo * aM;

        if obj.vT ~= 0
            dx(4) =   obj.aT / obj.vT;
        end
        
        dx(5) =	aM / obj.vM;
        
%        
    %
    % the short period dynamics is traditionally expressed in terms of q and
    % alpha but since the guidance produces acceleration commands it is required
    % to translate alpha to am with: 
    %   α = 1/z_α ⋅ a_m - z_δ/z_α ⋅ δ
    % to translate initial conditions or equilibrium point from [delta, q, alpha] system to [delta, q, am] system            
    %       the nonsingular transformation matrix is:
    %   T = [1 0 0; 0 1 0; zd 0 za]
    % from acceleration to delta command:
    %   δ =1/(z_δ⋅v)⋅a_m - z_α/z_δ ⋅ α
    %% 
        
%         idx = time2index(0 : dt : length(dcmd) * dt, t);
%         u = dcmd(idx) - K * x(6 : 8, 1);
%           b 3x1
%           K 1x3
%           x 3x1

        alpha = 
        dcmd = - (obj.z_a / obj.z_d) * alpha + 1 / (obj.z_d * obj.vM) * aMc;
        K = [0, 0.3 * q, 0]; 
        u = dcmd - K * x(6 : 8, 1);
        
        dx(6 : 8, 1) =  Ad * x(6 : 8, 1) + bd * u;%	- aM / obj.tauM + 1 / obj.tauM * aMc;

    end
end
    
end











