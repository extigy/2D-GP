function [t,fit_N_d,fit_N_a] = kwonforcefit(t,N_v,N_d,N_a)
    N_i = N_v(1);
    dt = t(2)-t(1);
    t0 = t(1);
    function N_d_gamma = fitfunc_d(gamma_d,t)
        function dN_ddt = ode(tt,~)
            dstep = min(floor((tt-t0+dt)/dt),size(N_d,1));
            dN_ddt = gamma_d.*(N_v(dstep));
        end
        [~,N_d_gamma] = ode45(@ode,t,0);
    end
    function N_a_gamma = fitfunc_a(gamma_a,t)
        function dN_adt = ode(tt,~)
            dstep = min(floor((tt-t0+dt)/dt),size(N_d,1));
            dN_adt = gamma_a.*(N_v(dstep))^2;
        end
        [~,N_a_gamma] = ode45(@ode,t,0);
    end
    %plot(t,N_v,'b')
    hold on
    [fit_gamma_d, ~, ~] = nlinfit(t,N_d,@fitfunc_d,0)
    [fit_gamma_a, ~, ~] = nlinfit(t,N_a,@fitfunc_a,0)
    fit_N_d = fitfunc_d(fit_gamma_d,t);
    fit_N_a = fitfunc_a(fit_gamma_a,t);
    plot(t,fit_N_d,'m--');
    plot(t,fit_N_a,'r--');
end
% function [t,fit_N_v,fit_gamma] = kwonforcefit(t,N_v)
%     function N_v_gamma = fitfunc(gamma,t)
%         function dNdt = ode(~,y)
%             dNdt = -gamma(1)*y -gamma(2)*(y^2);
%         end
%         [~,N_v_gamma] = ode45(@ode,t,N_v(1));
%     end
%     fit_gamma = [0.123,0.0053];
%     fit_N_v = fitfunc(fit_gamma,t);
%     hold all
%     plot(t,fit_N_v,'m--');
% end
