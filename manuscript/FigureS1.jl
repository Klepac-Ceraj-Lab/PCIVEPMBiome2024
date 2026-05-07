Random.seed!(0)
this_preproc = :arcsin

## 1. Loading data
postexclusion_sample_set = readlines(joinpath(Base.pwd(), "data", "kept_complete_samples.txt"))
preexclusion_sample_set = readlines(joinpath(Base.pwd(), "data", "complete_samples_plus_smitis.txt"))

mdata_cols, postexclusion_khula_pci_mbiome_data, _, _, _, _ = load_pcimbiome_data(postexclusion_sample_set)
mdata_cols, preexclusion_khula_pci_mbiome_data, _, _, _, _ = load_pcimbiome_data(preexclusion_sample_set)

## Post-exclusion analysis

postexclusion_glm_data = deepcopy(postexclusion_khula_pci_mbiome_data)
postexclusion_glm_data.feeding_state = Float64.(postexclusion_glm_data.feeding_state .== "ExclBreastFed")
postexclusion_glm_data.delivery_mode = Float64.(postexclusion_glm_data.delivery_mode .== "Caesarean")
postexclusion_glm_data.InfantVisAtt = Leap.rangenormalize(postexclusion_glm_data.InfantVisAtt)
postexclusion_glm_data.MaternalEntropy = Leap.rangenormalize(postexclusion_glm_data.MaternalEntropy)

postexclusion_lm_results = @chain filter_prevalence(postexclusion_glm_data, 0.05) begin
    runlms( _ , names(select(_, Not(mdata_cols)));
    prevalence_threshold = 0.05,
    abundance_threshold = 0.001,
    model_kind = :linear,
    bug_preproc = this_preproc,
    bugrank = false,
    outdir = joinpath("manuscript", "figures", "glm_postexclusion"))
end

postexclusion_lm_pvals = map(x -> subset(DataFrame(postexclusion_lm_results[:, x]), :pvalue => x -> .!(isnan.(x))), 1:5)
map(i -> postexclusion_lm_pvals[i].qvalue = adjust(postexclusion_lm_pvals[i].pvalue, BenjaminiHochberg()), 1:5)

postexclusion_results = sort(postexclusion_lm_pvals[2], [:pvalue, :qvalue])

## pre-exclusion analysis

preexclusion_glm_data = deepcopy(preexclusion_khula_pci_mbiome_data)
preexclusion_glm_data.feeding_state = Float64.(preexclusion_glm_data.feeding_state .== "ExclBreastFed")
preexclusion_glm_data.delivery_mode = Float64.(preexclusion_glm_data.delivery_mode .== "Caesarean")
preexclusion_glm_data.InfantVisAtt = Leap.rangenormalize(preexclusion_glm_data.InfantVisAtt)
preexclusion_glm_data.MaternalEntropy = Leap.rangenormalize(preexclusion_glm_data.MaternalEntropy)

preexclusion_lm_results = @chain filter_prevalence(preexclusion_glm_data, 0.05) begin
    runlms( _ , names(select(_, Not(mdata_cols)));
    prevalence_threshold = 0.05,
    abundance_threshold = 0.001,
    model_kind = :linear,
    bug_preproc = this_preproc,
    bugrank = false,
    outdir = joinpath("manuscript", "figures", "glm_preexclusion"))
end

preexclusion_lm_pvals = map(x -> subset(DataFrame(preexclusion_lm_results[:, x]), :pvalue => x -> .!(isnan.(x))), 1:5)
map(i -> preexclusion_lm_pvals[i].qvalue = adjust(preexclusion_lm_pvals[i].pvalue, BenjaminiHochberg()), 1:5)

preexclusion_results = sort(preexclusion_lm_pvals[2], [:pvalue, :qvalue])

## loading Cook's D for outlier
cooksD_smitis = CSV.read(
    joinpath("manuscript", "figures", "glm_preexclusion", "CooksD_Streptococcus_mitis.csv"),
    DataFrame;
    stringtype = String
)

cooksD_smitis.diff_f2 = cooksD_smitis.cooksD_f2 .- mean(cooksD_smitis.cooksD_f2)
sort!(cooksD_smitis, :cooksD_f2; rev = true)

fig = Figure(size = (1000, 350))

ax_A = Axis(
    fig[1,1],
    xlabel = "rank (decreasing)",
    ylabel = "Cook's D",
    title = "Outlier exclusion analysis via Cook's D",
)
hidedecorations!(ax_A; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )

bp_A = barplot!(
    ax_A,
    1:nrow(cooksD_smitis),
    cooksD_smitis.cooksD_f2,
    color = map( x -> ((x > 3*mean(cooksD_smitis.cooksD_f2)) ? (:orange, 1.0) : (:purple, 0.5)), cooksD_smitis.cooksD_f2)
)
hlines!(mean(cooksD_smitis.cooksD_f2); color = :black, linestyle = :dash)

ax_B = Axis(
    fig[1,2],
    xlabel = "arcsin-transformed abundance",
    ylabel = "unit-range normalized VOB",
    yticks = ([0.0, 0.5, 1.0]),
    title = "Before sample exclusion",
    aspect = 1.0
)
hidedecorations!(ax_B; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )

sc_B = scatter!(
    ax_B,
    asin.(sqrt.(preexclusion_khula_pci_mbiome_data[:, :Streptococcus_mitis])),
    Leap.rangenormalize(preexclusion_khula_pci_mbiome_data.InfantVisAtt),
    color = map( x -> (x == 0 ? :red : :blue), preexclusion_khula_pci_mbiome_data[:, :Streptococcus_mitis])
)

ax_C = Axis(
    fig[1,3],
    xlabel = "arcsin-transformed abundance",
    ylabel = "unit-range normalized VOB",
    yticks = ([0.0, 0.5, 1.0]),
    title = "After sample exclusion",
    aspect = 1.0
)
hidedecorations!(ax_C; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )

sc_C = scatter!(
    ax_C,
    asin.(sqrt.(postexclusion_khula_pci_mbiome_data[:, :Streptococcus_mitis])),
    Leap.rangenormalize(postexclusion_khula_pci_mbiome_data.InfantVisAtt),
    color = map( x -> (x == 0 ? :red : :blue), postexclusion_khula_pci_mbiome_data[:, :Streptococcus_mitis])
)

Legend(## Vertical version, to the side
    fig[2,1],
    [
        PolyElement(; color = (:orange, 1.0)),
        PolyElement(; color = (:purple, 0.5)),
        LineElement(; color = :black, linestyle = :dash)
    ],
    [
        "Sample D >> 3 × mean(D)", 
        "Sample D < 3 × mean(D)", 
        "mean(D)"
    ],
    orientation = :horizontal
)

Legend(
    fig[2, 2:3],
    [
        MarkerElement(; marker=:circle, color=:red, strokewidth=1),
        MarkerElement(; marker=:circle, color=:blue, strokewidth=1),
    ],
    [
        "Absent from sample",
        "Present in sample"
    ];
    tellheight = false,
    tellwidth = false,
    nbanks = 2
)

Label(fig[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[1, 2, TopLeft()], "b", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[1, 3, TopLeft()], "c", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())

colsize!(fig.layout, 1, Relative(0.50))
colsize!(fig.layout, 2, Relative(0.25))
colsize!(fig.layout, 3, Relative(0.25))
colgap!(fig.layout, 5)

fig

save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS1.png"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS1.eps"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS1.pdf"), fig)

## Exporting table for comparison
preexclusion_select_results = @chain preexclusion_results begin
    select( [ :feature, :coef, :pvalue] )
    rename!(:coef => :pre_coef, :pvalue => :pre_pvalue)
end

postexclusion_select_results = @chain postexclusion_results begin
    select( [ :feature, :coef, :pvalue ] )
    rename!(:coef => :post_coef, :pvalue => :post_pvalue)
end

prepost_comparison_table = @chain preexclusion_select_results begin
    innerjoin(postexclusion_select_results, on = :feature)
    select([:feature, :pre_coef, :post_coef, :pre_pvalue, :post_pvalue])
end

CSV.write(
    joinpath(Base.pwd(), "manuscript", "tables", "TableS2.csv"),
    prepost_comparison_table[1:10,:]
)