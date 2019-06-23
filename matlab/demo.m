clear
input_path = '../input/datasets/PBMC_4k_10X/';

%% Import expression matrix
    expression_table = readtable(fullfile(input_path, 'expression.txt'), 'delimiter', '\t');
    expression = table2array(expression_table(:, 2:end));
    gene_names = expression_table.Var1;
    cell_tags  = expression_table.Properties.VariableNames(2:end)';

%% Read sample annotations
    sample_annotations = readtable(fullfile(input_path, 'sample_annotations.txt'), 'FileType','text', 'delimiter', '\t');
    TrueLabels = sample_annotations.Labels;
    UL = unique(TrueLabels);
    [~, true_labels] = ismember(TrueLabels, UL);
    
%% Remove barely expressed genes
    counts = sum(spones(expression), 2);
    gene_filter_mask = (counts <= 10);
    expression = expression(~gene_filter_mask, :);
    gene_names = gene_names(~gene_filter_mask);

%% Normalize

if(max(expression(:)) > 100)
    expression = median(sum(expression))*normalize(expression, 'pnorm', 1);
    expression = full(spfun(@(x) log10(1+x), expression));
end


%% parameters
rows = 1:size(expression, 1);
subsample_size = 100;
subsample_no = 10000;
thread_no = 8;

%% Find differential genes between CD8/CD4
CD4_cols = find(strcmp(sample_annotations.Full_Labels, "CD4+/CD45RA+/CD25- Naive T"));
CD8_cols = find(strcmp(sample_annotations.Full_Labels, "CD8+ Cytotoxic T"));

[gene_logPvals, gene_logPvals_lower, gene_logPvals_upper] = AssessFeatures_BetweenGroups(expression, rows , CD8_cols, CD4_cols, subsample_size, subsample_no, thread_no);

[sorted_gene_pvals, perm] = sort(gene_logPvals, 'descend');
sorted_genes = gene_names(perm);
sorted_genes(1:20)

%% Find B-cell specific genes
B_cols = find(strcmp(sample_annotations.Labels, "B"));
[gene_logPvals, gene_logPvals_lower, gene_logPvals_upper] = AssessFeatures(expression, rows , B_cols, subsample_size, subsample_no, thread_no);

[sorted_gene_pvals_B, perm] = sort(gene_logPvals, 'descend');
sorted_genes = gene_names(perm);
sorted_genes(1:20)


%% Find Activity of "B module" in different cells
B_module = find(gene_logPvals > 2);

profile = ProfileModule(expression, B_module);
plot(smooth(profile, 10, 'moving'));




