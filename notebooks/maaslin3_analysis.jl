#####
# Running Linear Models (a.k.a MaAsLin3 - style stuff)
#####

khula_glm_data = deepcopy(khula_pci_mbiome_data)
khula_glm_data.feeding_state = Float64.(khula_glm_data.feeding_state .== "ExclBreastFed")
khula_glm_data.delivery_mode = Float64.(khula_glm_data.delivery_mode .== "Caesarean")
khula_glm_data.InfantVisAtt = Leap.rangenormalize(khula_glm_data.InfantVisAtt)
khula_glm_data.MaternalEntropy = Leap.rangenormalize(khula_glm_data.MaternalEntropy)

this_preproc = :arcsin

## Linear models
linear_lm_results = @chain filter_prevalence(khula_glm_data, 0.01) begin
    runlms( _ , names(select(_, Not(mdata_cols))); prevalence_threshold = 0.05, abundance_threshold = 0.01, model_kind = :linear, bug_preproc = this_preproc, bugrank = false)
end

linear_lm_pvals = map(x -> subset(DataFrame(linear_lm_results[:, x]), :pvalue => x -> .!(isnan.(x))), 1:5)
map(i -> linear_lm_pvals[i].qvalue = adjust(linear_lm_pvals[i].pvalue, BenjaminiHochberg()), 1:5)

println("Taxa Model 1: bug ~ sex + age + education + entropy; Stats:")
@show sort(linear_lm_pvals[1], [:pvalue, :qvalue])
println("Taxa Model 2: bug ~ sex + age + education + visual; Stats:")
@show sort(linear_lm_pvals[2], [:pvalue, :qvalue])
println("Taxa Model 3: bug ~ sex + age + education + visual*entropy; Stats for entropy:")
@show sort(linear_lm_pvals[3], [:pvalue, :qvalue])
println("Taxa Model 3: bug ~ sex + age + education + visual*entropy; Stats for visual:")
@show sort(linear_lm_pvals[4], [:pvalue, :qvalue])
println("Taxa Model 3: bug ~ sex + age + education + visual*entropy; Stats for **INTERACTION**:")
@show sort(linear_lm_pvals[5], [:pvalue, :qvalue])