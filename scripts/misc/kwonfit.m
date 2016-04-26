function [t,fit_N_v,fit_gamma] = kwonfit(t,N_v)
    function N_v_gamma = fitfunc(gamma,t)
        function dNdt = ode(~,y)
            dNdt = -gamma(1)*y -gamma(2)*(y^2);
        end
        [~,N_v_gamma] = ode45(@ode,t,N_v(1));
    end

    plot(t,N_v,'b')
    guess = [0,0];
    [fit_gamma, ~, ~] = nlinfit(t,N_v,@fitfunc,guess);
    fit_N_v = fitfunc(fit_gamma,t);
    hold all
    plot(t,fit_N_v,'g--');
end
