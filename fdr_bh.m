function [h, crit_p, adj_p] = fdr_bh(pvals, q, method)
    % Benjamini & Hochberg FDR controlling procedure

    if nargin < 2 || isempty(q), q = .05; end
    if nargin < 3 || isempty(method), method = 'pdep'; end

    p = pvals(:);
    V = length(p);
    [ps_sorted, sort_ids] = sort(p);

    if strcmpi(method, 'pdep')
        c_v = 1;
    elseif strcmpi(method, 'dep')
        c_v = sum(1 ./ (1:V));
    end

    % Threshold for rejection
    thresh = (1:V)' / V / c_v * q;

    % Adjusted p-values
    wtd_p = ps_sorted * V * c_v ./ (1:V)';

    % Enforce monotonicity (step-up)
    adj_p_sorted = cummin(wtd_p(end:-1:1));
    adj_p_sorted = adj_p_sorted(end:-1:1);

    % Map back to original indices
    adj_p = zeros(V, 1);
    adj_p(sort_ids) = adj_p_sorted;
    adj_p(adj_p > 1) = 1;

    rej = ps_sorted <= thresh;
    max_id = find(rej, 1, 'last');

    if isempty(max_id)
        crit_p = 0;
        h = false(V, 1);
    else
        crit_p = ps_sorted(max_id);
        h = p <= crit_p;
    end

    h = reshape(h, size(pvals));
    adj_p = reshape(adj_p, size(pvals));
end
