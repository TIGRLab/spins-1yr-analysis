%% reads in a correlation matrix CSV and uses BCT to generate graph metrics.
addpath(genpath('/projects/jdv/data/spins/1yr/outputs/bin/includes'));
sparsity = [0.01:0.01:0.20];

% degree
% regional efficiency / global efficiency
% clustering coefficient
% small-worldness
% robustness ?
% degree distribution paramaters (scale free?)

% init output datamat
data = zeros(length(sparsity), 3);
G = dlmread(filename, ',');

for sparse = sparsity;

    % binarize at some sparsity
    g = threshold_proportional(G, sparse);
    g(g > 0) = 1;
    nnodes = length(g);
    g(1:nnodes+1:nnodes*nnodes) = 0;

    % degree distribution numbers
    deg = degrees_und(g);
    degvar = var(deg);  % metric 1
    degpl = plfit(deg); % metric 2

    % efficiency
    eff = efficiency_bin(g); % metric 3

    % mean clustering coefficient
    kcoefg = mean(clustering_coef_bu(g)); % metric 4

    % small worldness
    r = randomizer_bin_und(g, 1);
    
    pthlng = 1/eff;
    pthlnr = 1/efficiency_bin(r);
    kcoefr = mean(clustering_coef_bu(r));

    smwrld = (kcoefg / kcoefr) / (pthlnr / pthlng); % metric 5

    % robustness (targeted and random)
    rbstt = zeros(nnodes, 1);
    rbstr = zeros(nnodes, 1);

    g_t = g; % targeted
    g_r = g; % random
    
    idx_r = 1:nnodes;

    for n = 1:nnodes;
        [s1, s2] = get_components(g_t); % targeted
        rbstt(n) = max(s2);             %
        [s1, s2] = get_components(g_r); % random
        rbstr(n) = max(s2);             %

        % for tageted attack, remove the (a) node with the largest cluster coef
        k = clustering_coef_bu(g_t); % clustering coefficient of each node
        idx_t = find(k == max(k));
        disp(max(k))
        if length(idx_t) > 1;
            idx_t = idx_t(1);
        end

        g_t(idx_t, :) = 0;
        g_t(:, idx_t) = 0;

        % for random attack, remove any old node that has not been removed
        idx_r = idx_r(randperm(length(idx_r)));
        g_r(idx_r(1), :) = 0;
        g_r(:, idx_r(1)) = 0;
        idx_r = idx_r(2:end);

    end

    rbstt = trapz(rbstt); % metric 6
    rbstr = trapz(rbstr); % metric 7

% % calculate characteristic path length, global efficiency
    % dist = distance_bin(G);
    % [l, e] = charpath(dist);

    % data(subj, 1) = l;
    % data(subj, 2) = e;

    % calculate average clustering coefficient
    c = mean(clustering_coef_bu(G));
    data(subj, 1) = c;

    % calculate transivity
    t = transitivity_bu(G);
    data(subj, 2) = t;

    % assortativity
    r = assortativity_bin(G, 0);
    data(subj, 3) = r;

end

for i = 1:3;
	subplot(3,1,i);

	a_ind = [1,4,7];
	b_ind = [2,5,8];
	c_ind = [3,6,9];
	a = data(a_ind, i);
	b = data(b_ind, i);
	c = data(c_ind, i);
	plot(1:3, a, 'color', 'red', 'linewidth', 2); hold all;
	plot(1:3, b, 'color', 'black', 'linewidth', 2);
	plot(1:3, c, 'color', 'green', 'linewidth', 2); hold off;
	set(gca, 'XTick', 1:3)

	if i == 1;
        title('Average Clustering Coefficient')
	end

	if i == 2;
	    title('Transitivitiy')
	end

	if i == 3;
        title('Assortativity')
        set(gca, 'XTickLabel', {'CMH', 'MRC', 'ZHH'});
	end
end

for i = 1:9;
    subplot(3,3,i);


