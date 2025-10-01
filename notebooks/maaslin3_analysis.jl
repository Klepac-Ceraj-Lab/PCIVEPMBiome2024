#####
# Running Linear Models (a.k.a MaAsLin3 - style stuff)
#####

khula_glm_data = deepcopy(khula_pci_mbiome_data)

# khula_glm_data.InfantVisAtt = Leap.autonormalize(khula_glm_data.InfantVisAtt)
# khula_glm_data.MaternalEntropy = Leap.autonormalize(khula_glm_data.MaternalEntropy)
khula_glm_data.InfantVisAtt = Leap.rangenormalize(khula_glm_data.InfantVisAtt)
khula_glm_data.MaternalEntropy = Leap.rangenormalize(khula_glm_data.MaternalEntropy)

this_preproc = :arcsin
# this_preproc = :log2
# this_preproc = :log10
# this_preproc = :rclr
# this_preproc = :autoscale
# this_preproc = :none

## Linear models
linear_lm_results = @chain filter_prevalence(khula_glm_data, 0.05) begin
    runlms( _ , names(select(_, Not(mdata_cols))); prevalence_threshold = 0.05, abundance_threshold = 0.01, model_kind = :linear, bug_preproc = this_preproc, bugrank = false)
end

linear_lm_pvals = map(x -> subset(DataFrame(linear_lm_results[:, x]), :pvalue => x -> .!(isnan.(x))), 1:5)
map(i -> linear_lm_pvals[i].qvalue = adjust(linear_lm_pvals[i].pvalue, BenjaminiHochberg()), 1:5)
@show sort(linear_lm_pvals[1], [:pvalue, :qvalue])
@show sort(linear_lm_pvals[2], [:pvalue, :qvalue])
@show sort(linear_lm_pvals[3], [:pvalue, :qvalue])
@show sort(linear_lm_pvals[4], [:pvalue, :qvalue])
@show sort(linear_lm_pvals[5], [:pvalue, :qvalue])

## Logistic models

logistic_lm_results = @chain filter_prevalence(khula_glm_data, 0.05) begin
    runlms( _ , names(select(_, Not(mdata_cols))); prevalence_threshold = 0.05, abundance_threshold = 0.01, model_kind = :logistic, bug_preproc = this_preproc)
end

logistic_lm_pvals = map(x -> subset(DataFrame(logistic_lm_results[:, x]), :pvalue => x -> .!(isnan.(x))), 1:5)
map(i -> logistic_lm_pvals[i].qvalue = adjust(logistic_lm_pvals[i].pvalue, BenjaminiHochberg()), 1:5)
@show sort(logistic_lm_pvals[1], [:pvalue, :qvalue])
@show sort(logistic_lm_pvals[2], [:pvalue, :qvalue])
@show sort(logistic_lm_pvals[3], [:pvalue, :qvalue])
@show sort(logistic_lm_pvals[4], [:pvalue, :qvalue])
@show sort(logistic_lm_pvals[5], [:pvalue, :qvalue])
