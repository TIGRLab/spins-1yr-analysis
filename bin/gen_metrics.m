%% reads in a correlation matrix CSV and uses BCT to generate graph metrics.
%addpath(genpath('/projects/jdv/data/spins/1yr/outputs/bin/includes'));

% degree
% regional efficiency / global efficiency
% clustering coefficient
% small-worldness
% robustness ?
% degree distribution paramaters (scale free?)

% init output datamat

function gen_metrics(filename, outputname);

sparsity = [0.01:0.01:0.20];
data = zeros(length(sparsity), 9);
G= dlmread(filename, ',');
nnodes = length(G);
G(1:nnodes+1:nnodes*nnodes) = 0;
sparsecount = 1;

for sparse = sparsity;

    % collect mean and variance of connectivity
    connmean = mean(mean(G)); % metric 1
    connvar = mean(var(G)); % metric 2

    % binarize at some sparsity
    g = threshold_proportional(G, sparse);
    g(g > 0) = 1;

    % degree distribution numbers
    deg = degrees_und(g);
    degvar = var(deg);  % metric 3
    degpl = plfit(deg); % metric 4

    % efficiency
    eff = efficiency_bin(g); % metric 5

    % mean clustering coefficient
    kcoefg = mean(clustering_coef_bu(g)); % metric 6

    % small worldness
    r = randomizer_bin_und(g, 1);

    pthlng = 1/eff;
    pthlnr = 1/efficiency_bin(r);
    kcoefr = mean(clustering_coef_bu(r));

    smwrld = (kcoefg / kcoefr) / (pthlnr / pthlng); % metric 7

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

    rbstt = trapz(rbstt); % metric 8
    rbstr = trapz(rbstr); % metric 9

    % % calculate characteristic path length, global efficiency
    % dist = distance_bin(G);
    % [l, e] = charpath(dist);

    % data(subj, 1) = l;
    % data(subj, 2) = e;

    % calculate average clustering coefficient
    %c = mean(clustering_coef_bu(G));
    %data(subj, 1) = c;

    % calculate transivity
    %t = transitivity_bu(G);
    %data(subj, 2) = t;

    % assortativity
    %r = assortativity_bin(G, 0);
    %data(subj, 3) = r;

    % load values into output array


    data(sparsecount, 1) = connmean;  
    data(sparsecount, 2) = connvar; 
    data(sparsecount, 3) = degvar; 
    data(sparsecount, 4) = degpl; 
    data(sparsecount, 5) = eff; 
    data(sparsecount, 6) = kcoefg; 
    data(sparsecount, 7) = smwrld; 
    data(sparsecount, 8) = rbstt; 
    data(sparsecount, 9) = rbstr; 

    sparsecount = sparsecount + 1;

end

% write out file
header='conn-mean,conn-var,ddist-var,ddist-plaw,eff,ddist-mean,smworld,robtarget,robrandom\n';
%fmt = repmat('%s,', 1, length(header));
fid = fopen(outputname, 'w');
%fprintf(fid, fmt, header{:});
fprintf(fid, header);
fclose(fid);

dlmwrite(outputname, data, '-append', 'delimiter', ',');

exit
end
