function f_curr = INS_monitor(prob, prob_tight)
    q = norm(prob_tight.res);
    if q > 50
        f_curr = prob_tight.f_vec;
    else
        f_curr = zeros(length(prob.rho), 1);
    end
end