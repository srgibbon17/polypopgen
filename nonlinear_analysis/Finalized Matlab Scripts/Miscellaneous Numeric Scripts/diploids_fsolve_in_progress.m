% for diploids, a nonlinear model and analysis

% determines number of values to include in 
iterations = 100;


% single point parameter values
s_val = 3.6e-8; % constant value of selection coefficient
mu_val = 5e-8; % constant value of forward mutation rate
nu_val = 1e-9; % constant value of backward mutation rate
h_val = 1; % constant value of dominance coefficient
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate

% parameter ranges 
s_val_range = logspace(-8, -6, iterations); % range of selection coefficients
mu_val_range = logspace(-9, -7, iterations); % range of forward mutation rates
nu_val_range = logspace(-8, -5, iterations); % range of backward mutation rates
h_val_range = linspace(.6, 1, iterations); % range of dominance coefficients

syms s q G0 G1 G2 g0 g1 h mu nu

specified_parameter = s;
specified_param_range = s_val_range;
specified_param_value = s_val;

variable_list = [g0, specified_parameter];
parameter_list = [s, mu, nu, h];
parameter_values = [s_val, mu_val, nu_val, h_val];
scaling_coefficient = 1e7;

s_val_list = [1e-8, 0.000000195914757, 2e-7, 5e-7, 0.000000677084801, 1e-5];

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
assume(s>=-1 & s<=1);
assume(h>=0 & h<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);

% equations to parameterize relative fitnesses
w_bar = 1 - s*(h*2*g0*g1 + g1^2);
w0 = 1/w_bar;
w1 = (1-s*h)/w_bar;
w2 = (1-s)/w_bar;

% equations for selection
sel_g0 = w1*g1*g0 + w0*g0^2; 
sel_g1 = w2*g1^2 + w1*g1*g0;

% equations for mutation
mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0;
mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1;
%removing g1 from the equations
mut_g0 = subs(mut_g0, g1, 1-g0);

%making an annotation for parameter values
s_str = strcat('s: ', string(s_val));
h_str = strcat('h: ',string(h_val)); 
mu_str = strcat('mu: ',string(mu_val));
nu_str = strcat('nu: ',string(nu_val));
mut_ratio_str = strcat(['mut-ratio ', newline, ...
      '(mu/nu): '],string(mut_ratio_val));
parameters_str_list = [s_str, mu_str, nu_str, mut_ratio_str, h_str];

%figure

% v = VideoWriter('bifn_animation_diploids.m4v', 'MPEG-4');
% v.Quality = 95;
% v.FrameRate = 20;
% open(v)
% for k = 1:3
%     for i = 1:length(s_val_range)    
%         s_val_list = s_val_range(i);
% 
%         fsolve_phase_plane_plot(mut_g0, parameter_values, parameters_str_list, s_val_list, 1);
%         frame = getframe(gcf);
%         writeVideo(v, frame)
%     end
% end
% close(v)

% figure
% movie(M, 2, 20)

%[x0_exitflags, x1_exitflags] = bifn_surface_plot_pos_sel(mu_val_range, h_val_range, nu_val, iterations);

%accurate_soln = sum(x0_exitflags == 1);
%inaccurate_soln_1 = sum(x0_exitflags == 2);
%inaccurate_soln_2 = sum(x0_exitflags == 3);


[diff_eqn_function, bifn_point_1, bifn_point_2] = fsolve_bifn_diagram_init(mut_g0, specified_parameter, parameter_list, parameter_values, g0, scaling_coefficient); 
 
fsolve_bifn_diagram_plot(diff_eqn_function, specified_parameter, parameter_list, parameters_str_list, specified_param_range, nu_val, mu_val, bifn_point_1, bifn_point_2);



% v = VideoWriter('bifn_diagram_h_animation_diploids.m4v', 'MPEG-4');
% v.Quality = 95;
% v.FrameRate = 10;
% open(v)
% figure
% for k = 1:1
%     for i = 1:length(h_val_range)
%         parameter_values(4) = h_val_range(i);
% 
%         [diff_eqn_function, bifn_point_1, bifn_point_2] = fsolve_bifn_diagram_init(mut_g0, specified_parameter, parameter_list, parameter_values, g0, scaling_coefficient); 
% 
%         fsolve_bifn_diagram_plot(diff_eqn_function, specified_parameter, parameter_list, parameters_str_list, specified_param_range, bifn_point_1, bifn_point_2);
%         frame = getframe(gcf);
%         writeVideo(v, frame)
%     end
% end
% close(v)


function [diff_eqn_value] = diff_eqn_eval_sym(mut_exp_g0, mu, mu_value, nu, nu_value, s, s_value, h, h_value, g0, g0_sub_value)

    diff_eqn_value = subs(mut_exp_g0, mu, mu_value);
    diff_eqn_value = subs(diff_eqn_value, nu, nu_value);
    diff_eqn_value = subs(diff_eqn_value, s, s_value);
    diff_eqn_value = subs(diff_eqn_value, h, h_value);
    diff_eqn_value = subs(diff_eqn_value, g0, g0_sub_value);
end

function [soln] = fixed_points_vpasolve(mut_exp_g0, mu, mu_value, nu, nu_value, s, s_value, h, h_value, g0)

    sub_eqn = subs(mut_exp_g0, mu, mu_value);
    sub_eqn = subs(sub_eqn, nu, nu_value);
    sub_eqn = subs(sub_eqn, s, s_value);
    sub_eqn = subs(sub_eqn, h, h_value);
    
    soln = vpasolve(sub_eqn, g0);
end

function [g0_soln, s_soln] = bifn_vpasolve(delta_g0, derivative_g0, mu, mu_value, nu, nu_value, h, h_value, g0, s)

    delta_g0 = subs(delta_g0, mu, mu_value);
    delta_g0 = subs(delta_g0, nu, nu_value);
    delta_g0 = subs(delta_g0, h, h_value);

    derivative_g0 = subs(derivative_g0, mu, mu_value);
    derivative_g0 = subs(derivative_g0, nu, nu_value);
    derivative_g0 = subs(derivative_g0, h, h_value);
    
    [g0_soln, s_soln] = vpasolve([delta_g0, derivative_g0], [g0, s]);
end

function F = parameterized_bifn_functions(x, h, m, n)
    a = x(2)*(1-2*h);
    b = x(2)*(3*h-2-n-h*m+h*n);
    c = -m-n+x(2)-h*x(2)+2*n*x(2)+h*m*x(2)-h*n*x(2);
    d = n*(1-x(2));

    F(1) = 1e10*(a*x(1)^3 + b*x(1)^2 + c*x(1) + d);

    F(2) = 1e10*(3*a*x(1)^2 + 2*b*x(1) + c);

end

function [x0_exit_list, x1_exit_list] = bifn_surface_plot_pos_sel(mu_val_range, h_val_range, nu_val, iterations)

    g0_max_array = zeros(iterations, iterations);
    g0_min_array = zeros(iterations, iterations);

    s_max_array = ones(iterations, iterations);
    s_min_array = zeros(iterations, iterations);

    [mu_coord, h_coord] = meshgrid(mu_val_range, h_val_range);

    x0_guess = [.5, .000001];
    x1_guess = [.01, .0000001];

    x0_exit_list = [];
    x1_exit_list = [];

    for i = 1:length(mu_val_range)
        for j = 1:length(h_val_range)
            fun_1 = @(x)parameterized_bifn_functions(x, h_coord(i,j), mu_coord(i, j), nu_val);
            [x0_soln, x0_fval, x0_exitflag] = fsolve(fun_1, x0_guess);
            [x1_soln, x1_fval, x1_exitflag] = fsolve(fun_1, x1_guess);
            x0_exit_list(end+1) = x0_exitflag;
            x1_exit_list(end+1) = x1_exitflag;

            if x0_exitflag > 0
                if x0_soln(1) > 0 && x0_soln(1) < 1 && x0_soln(2) > 0 && x0_soln(2) < 1
                    g0_max_array(i, j) = x0_soln(1); 
                    s_max_array(i, j) = x0_soln(2);
                end
            end
            if x1_exitflag > 0
                if x1_soln(1) > 0 && x1_soln(1) < 1 && x1_soln(2) > 0 && x1_soln(2) < 1
                    g0_min_array(i, j) = x1_soln(1);
                    s_min_array(i, j) = x1_soln(2);
                end
            end
        end
    end

    figure
    surf(mu_coord./nu_val, h_coord, s_min_array./s_max_array, abs(g0_min_array-g0_max_array))

    xlabel('mu/nu ratio (forward/back)')
    ylabel('h (dominance coefficient)')
    zlabel('s1/s2 (ratio of selection coefficients)')

    title("Diploid Bifurcation Surface")
    colorbar

end

function [] = fsolve_phase_plane_plot(delta_g0, param_values, param_str_list, s_val_list, subplot_index)
    
    diff_eqn = matlabFunction(delta_g0);

    s = s_val_list(subplot_index);
    mu = param_values(2);
    nu = param_values(3);
    h = param_values(4);

    diff_eqn_2 = @(g0)diff_eqn(g0, h, mu, nu, s);

    g0_values = linspace(0, 1, 1000);
    delta_g0_values = zeros(1, length(g0_values));
    
    for i = 1:length(g0_values)
        delta_g0_values(i) = diff_eqn_2(g0_values(i));
    end

    fixed_point_1 = fzero(diff_eqn_2, 1);
    fixed_point_2 = fzero(diff_eqn_2, 0);
    fixed_point_3 = fzero(diff_eqn_2, .35);

    param_str = {'Parameters:', param_str_list(1), param_str_list(2), param_str_list(3), param_str_list(4), param_str_list(5)};

    
    if length(s_val_list) == 1
        plot(g0_values, delta_g0_values, 'Color', '#0072BD', 'LineWidth', 2)
        hold on
        plot(g0_values, zeros(1, length(g0_values)), 'LineStyle','--', 'Color', 'k', 'LineWidth', 1)
        plot(fixed_point_1, 0, 'o', 'Color', "#D95319", 'LineWidth', 2, 'MarkerSize', 8)
        plot(fixed_point_3, 0, 'o', 'Color', "#77AC30", 'LineWidth', 2, 'MarkerSize', 8)
        plot(fixed_point_2, 0, 'o', 'Color', "#D95319", 'LineWidth', 2, 'MarkerSize', 8)
        hold off
        %dim = [0.5 0.5 0.3 0.3];
        %annotation('textbox', dim, 'String', param_str,'FontSize', 12, 'FitBoxToText','on')
        %legend('Δg0', '', 'stable fixed point', 'unstable fixed point', '', 'FontSize', 12)    
        ylim([-6e-8, 15e-8])
        xlabel('g0', 'FontSize', 14)
        ylabel('Δg0', 'FontSize', 14)
        sgtitle('Diploids Phase Plane', 'FontSize', 16)
        %pause(0.1)
        

        %annotation("arrow", [.3, .1], [.75, .75], 'LineWidth', 1)
        %annotation("arrow", [.3, .5], [.75, .75], 'LineWidth', 1)

        %annotation("arrow", [.3, .1], [.3, .3], 'LineWidth', 1)
        %annotation("arrow", [.3, .5], [.3, .3], 'LineWidth', 1)

    else
        subplot(2, 3, subplot_index)

        plot(g0_values, delta_g0_values, 'Color', '#0072BD', 'LineWidth', 2)
        hold on
        plot(g0_values, zeros(1, length(g0_values)), 'LineStyle','--', 'Color', 'k', 'LineWidth', 1)
        plot(fixed_point_1, 0, 'o', 'Color', "#D95319", 'LineWidth', 2, 'MarkerSize', 8)
        plot(fixed_point_3, 0, 'o', 'Color', "#77AC30", 'LineWidth', 2, 'MarkerSize', 8)
        plot(fixed_point_2, 0, 'o', 'Color', "#D95319", 'LineWidth', 2, 'MarkerSize', 8)
        dim = [0.5 0.5 0.3 0.3];
        annotation('textbox', dim, 'String', param_str,'FontSize', 12, 'FitBoxToText','on')
        legend('Δg0', '', 'stable fixed point', 'unstable fixed point', '', 'FontSize', 12)    

        xlabel('g0', 'FontSize', 14)
        ylabel('Δg0', 'FontSize', 14)
        sgtitle('Diploids Phase Plane', 'FontSize', 16)

        annotation("arrow", [.3, .1], [.75, .75], 'LineWidth', 1)
        annotation("arrow", [.3, .5], [.75, .75], 'LineWidth', 1)

        annotation("arrow", [.3, .1], [.3, .3], 'LineWidth', 1)
        annotation("arrow", [.3, .5], [.3, .3], 'LineWidth', 1)
    
    end


end


function [diff_eqn_function, varargout] = fsolve_bifn_diagram_init(diff_eqn, specified_parameter, param_list, param_values, g0, scaling_coeff) 

    %%% 
    % initializes the necessary inputs for the bifurcation diagram by
    % solving for the bifurcation points and substituting the 

    %Inputs:
        %diff_eqn
            %symbolic expression for the difference equation
        %specified_parameter
            %the specified parameter which will be plotted against in the
            %bifurcation diagram
        %param_list
            %a preset list of all parameters in their symbolic form
            %this allows subs to replace them with their respective values
        %param_values
            %a list of the parameter values on this iteration which will be
            %subbed into the difference equation and its first derivative 
        %g0
            %the symbolic variable g0 which allows subs to operate on the
            %difference equation and its derivative
        %scaling_coeff
            %a numeric value on the approximate order of s, mu, and nu
            %which acts to counteract the small scaling effects of these
            %parameters and enables fsolve to converge properly

    %Returns: 
        %diff_eqn_function
            %a Matlab function handle for the difference equation with the
            %non-specified parameter values substituted in and the
            %specified parameter as a parameter of the function
        %varargout
            %the set of all (non-trivial) bifurcation points for the given
            %parameter input values

    %%%

    syms x1 x2 y

    %simplifying the difference equation
    diff_eqn = simplify(expand(diff_eqn));

    %calculating symbolic derivative for linear stability
    %analysis/bifurcation point
    derivative_sym = diff(diff_eqn, g0);

    %substituting in new, non-symbolic variable into the
    %difference equation and its derivative
    diff_eqn_sub = subs(diff_eqn, g0, x1);
    diff_eqn_sub = subs(diff_eqn_sub, specified_parameter, y);
    diff_eqn_sub = subs(diff_eqn_sub, param_list(1), param_values(1));
    diff_eqn_sub = subs(diff_eqn_sub, param_list(2), param_values(2));
    diff_eqn_sub = subs(diff_eqn_sub, param_list(3), param_values(3));
    diff_eqn_sub = subs(diff_eqn_sub, param_list(4), param_values(4));

    diff_eqn_sub_2 = subs(diff_eqn, g0, x1);
    diff_eqn_sub_2 = subs(diff_eqn_sub_2, specified_parameter, x2);
    diff_eqn_sub_2 = subs(diff_eqn_sub_2, param_list(1), param_values(1));
    diff_eqn_sub_2 = subs(diff_eqn_sub_2, param_list(2), param_values(2));
    diff_eqn_sub_2 = subs(diff_eqn_sub_2, param_list(3), param_values(3));
    diff_eqn_sub_2 = subs(diff_eqn_sub_2, param_list(4), param_values(4));

    derivative_sub = subs(derivative_sym, g0, x1);
    derivative_sub = subs(derivative_sub, specified_parameter, x2);
    derivative_sub = subs(derivative_sub, param_list(1), param_values(1));
    derivative_sub = subs(derivative_sub, param_list(2), param_values(2));
    derivative_sub = subs(derivative_sub, param_list(3), param_values(3));
    derivative_sub = subs(derivative_sub, param_list(4), param_values(4));

    %scaling difference equation and derivative to "fix" numeric problems
    F = [scaling_coeff*diff_eqn_sub_2, scaling_coeff*derivative_sub];

    %creating a matlab function G which evaluates/solves for the symbolic  
    %equations in F
    G = matlabFunction(F);

    %modifies the matlab function G to be modG which has vectorized inputs
    %x(1) and x(2) where x(1) is g0 and x(2) is the specified parameter
    modG = @(x)G(x(1), x(2));

    diff_eqn_function = matlabFunction(scaling_coeff*diff_eqn_sub);

    if specified_parameter == param_list(1) %s

        %creates two initial conditions/guesses to pass to fsolve
        x0_guess = [.5, .000001];
        x1_guess = [.01, .0000001];
        
        %solves for the bifurcation points
        [x0_soln, x0_fval, x0_exitflag] = fsolve(modG, x0_guess)
        [x1_soln, x1_fval, x1_exitflag] = fsolve(modG, x1_guess)
        
        varargout = {[], []};
        if  x0_exitflag > 0
            varargout{1} = x0_soln;
        end
        if x1_exitflag > 0
            varargout{2} = x1_soln;
        end
        

    elseif specified_parameter == param_list(2) %mu

        %creates two initial conditions/guesses to pass to fsolve
        x0_guess = [.5, .00000001];
        x1_guess = [.01, .00001];
        
        %solves for the bifurcation points
        [x0_soln, x0_fval, x0_exitflag] = fsolve(modG, x0_guess);
        [x1_soln, x1_fval, x1_exitflag] = fsolve(modG, x1_guess);
        
        varargout = {[], []};
        if  x0_exitflag > 0
            varargout{1} = x0_soln;
        end
        if x1_exitflag > 0
            varargout{2} = x1_soln;
        end

    elseif specified_parameter == param_list(3) %nu

        %creates one initial condition/guess to pass to fsolve
        x0_guess = [.01, .0000001];
        
        %solves for the bifurcation point
        [x0_soln, x0_fval, x0_exitflag] = fsolve(modG, x0_guess);
        
        varargout = {[], []};
        if  x0_exitflag > 0
            varargout{1} = x0_soln;
        end

    elseif specified_parameter == param_list(4) %h
        %creates one initial condition/guess to pass to fsolve
        x0_guess = [.01, .6];
        
        %solves for the bifurcation point
        [x0_soln, x0_fval, x0_exitflag] = fsolve(modG, x0_guess);
        
        varargout = {[], []};
        if  x0_exitflag > 0
            varargout{1} = x0_soln;
        end
    else 
        disp('Error. Parameter not specified or specified improperly.')
    end

end

function [] = fsolve_bifn_diagram_plot(diff_eqn_function_param, specified_parameter, param_list, param_str_list, specified_param_range, nu_val, mu_val, varargin)

        stable_1_g0 = [];
        stable_2_g0 = [];
        unstable_g0 = [];

        stable_1_param = [];
        stable_2_param = [];
        unstable_param = [];

    
        if specified_parameter == param_list(1) %s
       
            bifn_1 = varargin{1};
            bifn_2 = varargin{2};

            bifn_g0_range = abs(bifn_1(1) - bifn_2(1));

            max_g0_bifn_value = max(bifn_1(1), bifn_2(1));
            min_g0_bifn_value = min(bifn_1(1), bifn_2(1));

            stable_1_guess = 1;

            stable_2_guess = 0;

            bifn_param_max = max(bifn_1(2), bifn_2(2));
            bifn_param_min = min(bifn_1(2), bifn_2(2));

            slope = (bifn_1(1) - bifn_2(1))/(bifn_1(2) - bifn_2(2));
            intercept = bifn_1(1) - slope*bifn_1(2);

            for i = 1:length(specified_param_range)
        
                y = specified_param_range(i);
        
                unstable_guess_1 = slope*y+intercept;

                diff_eqn_function = @(x1)diff_eqn_function_param(x1, y);

                if y < bifn_param_min
                    stable_2_g0_value = fzero(diff_eqn_function, stable_2_guess);
                    if stable_2_g0_value < min_g0_bifn_value
                        stable_2_g0(end+1) = stable_2_g0_value;
                        stable_2_param(end+1) = y;
                    end

                elseif y > bifn_param_min && y < bifn_param_max
                
                    stable_2_g0_value = fzero(diff_eqn_function, stable_2_guess);
                    if stable_2_g0_value < min_g0_bifn_value
                        stable_2_g0(end+1) = stable_2_g0_value;
                        stable_2_param(end+1) = y;
                    end

                    stable_1_g0_value = fzero(diff_eqn_function, stable_1_guess);
                    if stable_1_g0_value > max_g0_bifn_value
                        stable_1_g0(end+1) = stable_1_g0_value;
                        stable_1_param(end+1) = y;
                    end

                    unstable_g0_value = fzero(diff_eqn_function, unstable_guess_1);
                    if unstable_g0_value > min_g0_bifn_value && unstable_g0_value < max_g0_bifn_value
                        unstable_g0(end+1) = unstable_g0_value;
                    else
                        unstable_g0_value = fzero(diff_eqn_function, .5*(max_g0_bifn_value + min_g0_bifn_value));
                        if unstable_g0_value > min_g0_bifn_value && unstable_g0_value < max_g0_bifn_value
                            unstable_g0(end+1) = unstable_g0_value;
                        else
                            unstable_g0_value = fzero(diff_eqn_function, min_g0_bifn_value+.001);
                            if unstable_g0_value > min_g0_bifn_value && unstable_g0_value < max_g0_bifn_value
                                unstable_g0(end+1) = unstable_g0_value;
                            else 
                                unstable_g0(end+1) = fzero(diff_eqn_function, max_g0_bifn_value-.001);
                            end
                        end
                    end
                    unstable_param(end+1) = y;

                elseif y > bifn_param_max
                    stable_1_g0_value = fzero(diff_eqn_function, stable_1_guess);
                    if stable_1_g0_value > max_g0_bifn_value
                        stable_1_g0(end+1) = stable_1_g0_value;
                        stable_1_param(end+1) = y;
                    end

                end
            
            end

            param_str = {'Parameters:', param_str_list(2), param_str_list(3), param_str_list(4), param_str_list(5)};

            

            figure

            subplot(1, 2, 1)
            plot(stable_1_param, stable_1_g0, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable');
            hold on
            plot(stable_2_param, stable_2_g0, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable');
            plot(unstable_param, unstable_g0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Unstable');
            plot(bifn_1(2), bifn_1(1), '.', 'Color', '#0072BD', 'MarkerSize', 20, 'DisplayName', 'Bifurcation Point')
            plot(bifn_2(2), bifn_2(1), '.', 'Color', '#0072BD', 'MarkerSize', 20, 'DisplayName', 'Bifurcation Point')
            yline(nu_val/(nu_val+mu_val), 'Color', "#D95319", 'LineWidth', 1.5, 'LineStyle','--', 'DisplayName', 'Mutation Equilibrium')
            hold off
            xlabel('s', 'FontSize', 14)
            ylabel('g0/p', 'FontSize', 14)
            sgtitle('Diploid Bifurcation Diagram', 'FontSize', 16)  
            ylim([0, 1])
            xlim([1e-7, 1e-6])
            xscale log
            title('g0 vs s Diagram', 'FontSize', 12)
            dim = [0.2 0.2 0.3 0.3];
            annotation('textbox', dim, 'String', param_str,'FontSize', 14, 'FitBoxToText','on')
            legend('FontSize', 12)

            subplot(1, 2, 2)

            stable_1_g1 = ones(1, length(stable_1_g0)) - stable_1_g0;
            stable_2_g1 = ones(1, length(stable_2_g0)) - stable_2_g0;
            unstable_g1 = ones(1, length(unstable_g0)) - unstable_g0;

            plot(stable_1_param, stable_1_g1, 'Color', 'k', 'LineWidth', 1.5);
            hold on
            plot(stable_2_param, stable_2_g1, 'Color', 'k', 'LineWidth', 1.5);
            plot(unstable_param, unstable_g1, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
            plot(bifn_1(2), 1-bifn_1(1), '.', 'Color', '#0072BD', 'MarkerSize', 20)
            plot(bifn_2(2), 1-bifn_2(1), '.', 'Color', '#0072BD', 'MarkerSize', 20)
            yline(mu_val/(nu_val+mu_val), 'Color', "#D95319", 'LineWidth', 1.5, 'LineStyle','--', 'DisplayName', 'Mutation Equilibrium')
            hold off
            xlabel('s', 'FontSize', 14)
            ylabel('g1/q', 'FontSize', 14)
            ylim([0, 1])
            xlim([1e-7, 1e-6])
            xscale log
            title('g1 vs s Diagram', 'FontSize', 14)
            


        elseif specified_parameter == param_list(2) %mu

            bifn_1 = varargin{1};
            bifn_2 = varargin{2};

            bifn_g0_range = abs(bifn_1(1) - bifn_2(1));

            max_g0_bifn_value = max(bifn_1(1), bifn_2(1));

            unstable_guess = max_g0_bifn_value - .5*bifn_g0_range;
        
            stable_1_guess = 0;

            stable_2_guess = 1;

            percent_up = 1.04;
            percent_down = .99;

            bifn_param_max = max(bifn_1(2), bifn_2(2));
            bifn_param_min = min(bifn_1(2), bifn_2(2));

            for i = 1:length(specified_param_range)
        
                y = specified_param_range(i);
        
                diff_eqn_function = @(x1)diff_eqn_function_param(x1, y);
            
                if y < bifn_param_min*percent_down
                    
                    stable_2_g0(end+1) = fzero(diff_eqn_function, stable_2_guess);
                    stable_2_param(end+1) = y;

                elseif y > bifn_param_min*percent_up && y < bifn_param_max*percent_down
                
                    stable_2_g0(end+1) = fzero(diff_eqn_function, stable_2_guess);
                    stable_2_param(end+1) = y;

                    stable_1_g0(end+1) = fzero(diff_eqn_function, stable_1_guess);
                    stable_1_param(end+1) = y;

                    unstable_g0(end+1) = fzero(diff_eqn_function, unstable_guess);
                    unstable_param(end+1) = y;

                elseif y > bifn_param_max*percent_up

                    stable_1_g0(end+1) = fzero(diff_eqn_function, stable_2_guess);
                    stable_1_param(end+1) = y;

                end
            end

            param_str = {'Parameters:', param_str_list(1), param_str_list(3), param_str_list(5)};

            figure

            subplot(1, 2, 1)
            plot(stable_1_param, stable_1_g0, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable');
            hold on
            plot(stable_2_param, stable_2_g0, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', '');
            plot(unstable_param, unstable_g0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Unstable');
            plot(bifn_1(2), bifn_1(1), '.', 'Color', '#0072BD', 'MarkerSize', 20, 'DisplayName', 'Bifurcation Point')
            plot(bifn_2(2), bifn_2(1), '.', 'Color', '#0072BD', 'MarkerSize', 20, 'DisplayName', '')

            xlabel('mu', 'FontSize', 14)
            ylabel('g0/p', 'FontSize', 14)
            sgtitle('Diploid Bifurcation Diagram', 'FontSize', 16)  
            title('g0 vs mu Diagram', 'FontSize', 14)
            xscale log

            dim = [0.5 0.5 0.3 0.3];
            annotation('textbox', dim, 'String', param_str,'FontSize', 8, 'FitBoxToText','on')

            subplot(1, 2, 2)

            stable_1_g1 = ones(1, length(stable_1_g0)) - stable_1_g0;
            stable_2_g1 = ones(1, length(stable_2_g0)) - stable_2_g0;
            unstable_g1 = ones(1, length(unstable_g0)) - unstable_g0;

            plot(stable_1_param, stable_1_g1, 'Color', 'k', 'LineWidth', 1.5);
            hold on
            plot(stable_2_param, stable_2_g1, 'Color', 'k', 'LineWidth', 1.5);
            plot(unstable_param, unstable_g1, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
            plot(bifn_1(2), 1-bifn_1(1), '.', 'Color', '#0072BD', 'MarkerSize', 20)
            plot(bifn_2(2), 1-bifn_2(1), '.', 'Color', '#0072BD', 'MarkerSize', 20)

            xlabel('mu', 'FontSize', 14)
            ylabel('g1/q', 'FontSize', 14)
            title('g1 vs mu Diagram', 'FontSize', 14)
            xscale log
            legend('Stable', '', 'Unstable', 'Bifurcation Point', '')

        elseif specified_parameter == param_list(3) %nu
            
            bifn_1 = varargin{1};

            unstable_guess = bifn_1(1) + .1;
            
            stable_1_guess = 1;

            stable_2_guess = 0;

            percent_up = 1.01;
            percent_down = .97;

            for i = 1:length(specified_param_range)
        
                y = specified_param_range(i);
        
                diff_eqn_function = @(x1)diff_eqn_function_param(x1, y);
                
                if y < bifn_1(2)*percent_down
                    stable_2_g0(end+1) = fzero(diff_eqn_function, stable_2_guess);
                    stable_2_param(end+1) = y;

                    stable_1_g0(end+1) = fzero(diff_eqn_function, stable_1_guess);
                    stable_1_param(end+1) = y;

                    unstable_g0(end+1) = fzero(diff_eqn_function, unstable_guess);
                    unstable_param(end+1) = y;
                elseif y > bifn_1(2)*percent_up
                    stable_1_g0(end+1) = fzero(diff_eqn_function, stable_1_guess);
                    stable_1_param(end+1) = y;
                end
            end

            figure

            param_str = {'Parameters:', param_str_list(1), param_str_list(2), param_str_list(5)};

            subplot(1, 2, 1)
            plot(stable_1_param, stable_1_g0, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable');
            hold on
            plot(stable_2_param, stable_2_g0, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable');
            plot(unstable_param, unstable_g0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Unstable');
            plot(bifn_1(2), bifn_1(1), '.', 'Color', '#0072BD', 'MarkerSize', 20, 'DisplayName', 'Bifurcation Point 1')

            xlabel('nu', 'FontSize', 14)
            ylabel('g0/p', 'FontSize', 14)
            sgtitle('Diploid Bifurcation Diagram', 'FontSize', 16)  
            title('g0 vs nu Diagram', 'FontSize', 14)
            xscale log

            dim = [0.5 0.5 0.3 0.3];
            annotation('textbox', dim, 'String', param_str,'FontSize', 8, 'FitBoxToText','on')

            subplot(1, 2, 2)

            stable_1_g1 = ones(1, length(stable_1_g0)) - stable_1_g0;
            stable_2_g1 = ones(1, length(stable_2_g0)) - stable_2_g0;
            unstable_g1 = ones(1, length(unstable_g0)) - unstable_g0;

            plot(stable_1_param, stable_1_g1, 'Color', 'k', 'LineWidth', 1.5);
            hold on
            plot(stable_2_param, stable_2_g1, 'Color', 'k', 'LineWidth', 1.5);
            plot(unstable_param, unstable_g1, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
            plot(bifn_1(2), 1-bifn_1(1), '.', 'Color', '#0072BD', 'MarkerSize', 20)

            xlabel('nu', 'FontSize', 14)
            ylabel('g1/q', 'FontSize', 14)
            title('g1 vs nu Diagram', 'FontSize', 14)
            xscale log
            legend('Stable', '', 'Unstable', 'Bifurcation Point')
            

        elseif specified_parameter == param_list(4) %h
        
            bifn_1 = varargin{1};

            unstable_guess = bifn_1(1) + .1;

            stable_1_guess = 0;

            stable_2_guess = 1;

            percent_up = 1.01;
            percent_down = .97;

            for i = 1:length(specified_param_range)
        
                y = specified_param_range(i);
        
                diff_eqn_function = @(x1)diff_eqn_function_param(x1, y);

                if y < bifn_1(2)*percent_down
                    stable_2_g0(end+1) = fzero(diff_eqn_function, stable_2_guess);
                    stable_2_param(end+1) = y;
                elseif y > bifn_1(2)*percent_up
                    stable_2_g0(end+1) = fzero(diff_eqn_function, stable_2_guess);
                    stable_2_param(end+1) = y;

                    stable_1_g0(end+1) = fzero(diff_eqn_function, stable_1_guess);
                    stable_1_param(end+1) = y;

                    unstable_g0(end+1) = fzero(diff_eqn_function, unstable_guess);
                    unstable_param(end+1) = y;
                end
                
            end

            param_str = {'Parameters:', param_str_list(1), param_str_list(2), param_str_list(3), param_str_list(4)};

            figure

            subplot(1, 2, 1)
            plot(stable_1_param, stable_1_g0, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable1');
            hold on
            plot(stable_2_param, stable_2_g0, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable2');
            plot(unstable_param, unstable_g0, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Unstable');
            plot(bifn_1(2), bifn_1(1), '.', 'Color', '#0072BD', 'MarkerSize', 20, 'DisplayName', 'Bifurcation Point 1')

            xlabel('h', 'FontSize', 14)
            ylabel('g0', 'FontSize', 14)
            sgtitle('Diploid Bifurcation Diagram', 'FontSize', 16)  
            title('g0 vs h Diagram', 'FontSize', 14)

            dim = [0.5 0.5 0.3 0.3];
            annotation('textbox', dim, 'String', param_str,'FontSize', 8, 'FitBoxToText','on')

            subplot(1, 2, 2)

            stable_1_g1 = ones(1, length(stable_1_g0)) - stable_1_g0;
            stable_2_g1 = ones(1, length(stable_2_g0)) - stable_2_g0;
            unstable_g1 = ones(1, length(unstable_g0)) - unstable_g0;

            plot(stable_1_param, stable_1_g1, 'Color', 'k', 'LineWidth', 1.5);
            hold on
            plot(stable_2_param, stable_2_g1, 'Color', 'k', 'LineWidth', 1.5);
            plot(unstable_param, unstable_g1, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
            plot(bifn_1(2), 1-bifn_1(1), '.', 'Color', '#0072BD', 'MarkerSize', 20)

            xlabel('h', 'FontSize', 14)
            ylabel('g1', 'FontSize', 14)
            title('g1 vs h Diagram', 'FontSize', 14)
            legend('Stable', '', 'Unstable', 'Bifurcation Point')
            

        else 
            disp('Error. Parameter not specified or specified improperly.')
        end
    
end
