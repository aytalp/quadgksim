function [ress,err,iter]  =  quadgk_sim(fun,path_locs,tol,max_it,size_integrand)

% Adaptive Gauss-Kronrod quandrature for simultaneous integration of similar integrands

% Upgraded by Aytac Alparslan using quadgk routine of Matlab. 2021
% Based on "quadva" by Lawrence F. Shampine.
% Ref: L.F. Shampine, "Vectorized Adaptive Quadrature in Matlab",
% Journal of Computational and Applied Mathematics 211, 2008, pp.131-140.
% inputs: 
% fun: function handle including the similar integrands 
% path_locs: integration path
% tol: relative tolerance for the integrations
% max_it: maximum iteration for the routine
% size_integrand: number of similar integrals
% outputs:
% ress: the integration results
% err: errors for each of the integrations
% iter: number of iterations during the adaptive integration

pnodes = [ ...
    0.2077849550078985; 0.4058451513773972; 0.5860872354676911; ...
    0.7415311855993944; 0.8648644233597691; 0.9491079123427585; ...
    0.9914553711208126];
pwt = [ ...
    0.2044329400752989, 0.1903505780647854, 0.1690047266392679, ...
    0.1406532597155259, 0.1047900103222502, 0.06309209262997855, ...
    0.02293532201052922];
pwt7 = [0,0.3818300505051189,0,0.2797053914892767,0,0.1294849661688697,0];
nodes = [-pnodes(end:-1:1); 0; pnodes]; % nodes for the kronrod
wei = [pwt(end:-1:1), 0.2094821410847278, pwt]; % weigths for the kronrod
error_wei = wei - [pwt7(end:-1:1), 0.4179591836734694, pwt7]; %error estimator weigths

subs = [path_locs(1:end-1);path_locs(2:end)]; % sub-intervals
err = zeros(1,size_integrand); % initilize errors vector
ress = zeros(1,size_integrand); % initilize results vector
integral_OK = zeros(1,size_integrand); % initialize the completed integration flag
iter = 0; % initialize iteration counter
while true
    iter = iter + 1;
    if iter > max_it
        warning('Maximum iterations are exceeded, results may be inaccurate...');
        break;
    end
    midpt = sum(subs) / 2;   % midpoints of the subintervals
    halfh = diff(subs) / 2;  % half the lengths of the subintervals
    if min(abs(2 * halfh)) < 100 * eps(class(halfh)) * abs(path_locs(1) - path_locs(end)) %return the current values if the two nodes are very close to each other
        break;
    end
    
    kxs = nodes * halfh + midpt; %locations of the kronrod nodes along the integration path
    
    xx = reshape(kxs,1,[]);
    fvals = fun(xx); %function values at the kronrod nodes along the integration path (can be a matrix)
    fvals = reshape(fvals,numel(wei),[],size(fvals,2)); 

    ress_subs = pagemtimes(wei,fvals).* halfh; %integration result at each subinterval
    err_subs = pagemtimes(error_wei,fvals).* halfh; %error at each subinterval

    ress_subs = reshape(ress_subs,size(ress_subs,2),size(ress_subs,3));
    err_subs = reshape(err_subs,size(err_subs,2),size(err_subs,3));
    
    ind_OK3 = ones(size(ress_subs,1),1); %initialize the flag that checks if the integration is OK at the subintervals for all of the integrations
    for ii = 1:size_integrand
        if integral_OK(ii) == 0
            ind_OK2 = (abs(err_subs(:,ii)) <= (tol * abs(sum(ress_subs(:,ii))))); % find the subintervals where the integration is ok for the integral-ii
            ind_OK3 = ind_OK2.* ind_OK3; % mark the subintervals where the integration is not OK ( = 0)
            if all(ind_OK2) 
                integral_OK(ii) = 1; % if all the subinterval integrations are OK, then do not check this integration anymore.
            end
        end
    end
    
    err = err + ind_OK3' * err_subs; % add the current error to the corresponding integration's error
    ress = ress + ind_OK3' * ress_subs; % add the integral value to the corresponding integration's value

    ind_OK = find(ind_OK3 == 1); % locate the subintervals where the integration was OK
    subs(:,ind_OK) = []; % delete the subintervals where the integration was OK and handled during the current iteration
    midpt(ind_OK) = []; % delete the midpoints where the integration was OK and handled during the current iteration
    
    if isempty(midpt)
        break; % end the integation when all the subinterval integrals are OK
    end
    subs = reshape([subs(1,:); midpt; midpt; subs(2,:)],2,[]); %organize the subintervals to be used in the next iteration.

end
end

