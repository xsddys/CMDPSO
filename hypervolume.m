function v = hypervolume(P, r, N)
    % HYPERVOUME    Hypervolume indicator as a measure of Pareto front estimate.
    %   V = HYPERVOLUME(P,R,N) returns an estimation of the hypervolume (in
    %   percentage) dominated by the approximated Pareto front set P (n by d)
    %   and bounded by the reference point R (1 by d). The estimation is done
    %   through N (default is 1000) uniformly distributed random points within
    %   the bounded hyper-cuboid.
    %
    %   V = HYPERVOLUME(P,R,C) uses the test points specified in C (N by d).
    %
    % See also: paretofront, paretoGroup

 
    if nargin < 2 || nargin > 3
        error('hypervolume:InvalidInput', ...
            'Number of input arguments must be between 2 and 3.');
    end
    if nargout > 1
        error('hypervolume:InvalidOutput', ...
            'Number of output arguments must be 0 or 1.');
    end


    P = P * diag(1 ./ r);
    [n, d] = size(P);
    if nargin < 3
        N = 1000;
    end
    if ~isscalar(N)
        C = N;
        N = size(C, 1);
    else
        C = rand(N, d);
    end

    % Calculate whether dominated by Pareto front
    fDominated = false(N, 1);
    lB = min(P);
    fcheck = all(bsxfun(@gt, C, lB), 2);

    for k = 1:n
        if any(fcheck)
            f = all(bsxfun(@gt, C(fcheck, :), P(k, :)), 2);
            fDominated(fcheck) = f;
            fcheck(fcheck) = ~f;
        end
    end

    % Calculate hypervolume
    v = sum(fDominated) / N;
    
    % Plotting part
    %figure;
    %scatter(C(:, 1), C(:, 2), 10, 'b', 'filled'); % All sampling points
    %hold on;
    %scatter(C(fDominated, 1), C(fDominated, 2), 10, 'g', 'filled'); 

    %xlabel('Normalized Relative Deviation','FontSize', 14,'FontName','Times New Roman','FontWeight','bold');
    %ylabel('Normalized Activation Ratio of Gateways','FontSize', 14,'FontName','Times New Roman','FontWeight','bold');
    %title('Hypervolume of CMDPSO', 'FontSize',20,'FontWeight','bold','FontName','Times New Roman');
    %legend('Points out of Hypervolume', 'Points in Hypervolume', 'FontSize',12,'Location', 'Best','FontName','Times New Roman');
    %box on; % Add top and right scales
    %hold off;
end
