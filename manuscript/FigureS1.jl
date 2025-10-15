## Supplementary Figure S1 - PC1vPC2, PC2vPC3, PC1vPC3 for taxa and genes, colro by both MatEnt and InfVisAtt
fig = Figure(; size = (1200, 1400))

ax11 = Axis(
    fig[1,1],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS2 ("*string(round(100*uni_MDS_variances[2]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)"
)
hidedecorations!(ax11; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc11 = scatter!(
    ax11,
    uni_MDS_columns[:,1],
    uni_MDS_columns[:,2];
    markersize = 10,
    color = unimdata.MaternalEntropy,
    colormap = :plasma
)

ax12 = Axis(
    fig[1,2],
    xlabel = "MDS2 ("*string(round(100*uni_MDS_variances[2]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)"
)
hidedecorations!(ax12; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc12 = scatter!(
    ax12,
    uni_MDS_columns[:,2],
    uni_MDS_columns[:,3];
    markersize = 10,
    color = unimdata.MaternalEntropy,
    colormap = :plasma
)

ax13 = Axis(
    fig[1,3],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)"
)
hidedecorations!(ax13; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc13 = scatter!(
    ax13,
    uni_MDS_columns[:,1],
    uni_MDS_columns[:,3];
    markersize = 10,
    color = unimdata.MaternalEntropy,
    colormap = :plasma
)

ax21 = Axis(
    fig[2,1],
    xlabel = "MDS1 ("*string(round(100*spe_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS2 ("*string(round(100*spe_MDS_variances[2]; digits = 2))*"%)",
    title = "Taxonomic Profiles (species)",
)
hidedecorations!(ax21; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc21 = scatter!(
    ax21,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 10,
    color = khula_pci_mbiome_data.MaternalEntropy,
    colormap = :plasma
)

ax22 = Axis(
    fig[2,2],
    xlabel = "MDS2 ("*string(round(100*spe_MDS_variances[2]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*spe_MDS_variances[3]; digits = 2))*"%)",
    title = "Taxonomic Profiles (species)"
)
hidedecorations!(ax22; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc22 = scatter!(
    ax22,
    spe_MDS_columns[:,2],
    spe_MDS_columns[:,3];
    markersize = 10,
    color = khula_pci_mbiome_data.MaternalEntropy,
    colormap = :plasma
)

ax23 = Axis(
    fig[2,3],
    xlabel = "MDS1 ("*string(round(100*spe_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*spe_MDS_variances[3]; digits = 2))*"%)",
    title = "Taxonomic Profiles (species)"
)
hidedecorations!(ax23; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc23 = scatter!(
    ax23,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,3];
    markersize = 10,
    color = khula_pci_mbiome_data.MaternalEntropy,
    colormap = :plasma
)

ax31 = Axis(
    fig[3,1],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS2 ("*string(round(100*uni_MDS_variances[2]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)"
)
hidedecorations!(ax31; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc31 = scatter!(
    ax31,
    uni_MDS_columns[:,1],
    uni_MDS_columns[:,2];
    markersize = 10,
    color = unimdata.InfantVisAtt,
    colormap = :viridis
)

ax32 = Axis(
    fig[3,2],
    xlabel = "MDS2 ("*string(round(100*uni_MDS_variances[2]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)"
)
hidedecorations!(ax32; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc32 = scatter!(
    ax32,
    uni_MDS_columns[:,2],
    uni_MDS_columns[:,3];
    markersize = 10,
    color = unimdata.InfantVisAtt,
    colormap = :viridis
)

ax33 = Axis(
    fig[3,3],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)"
)
hidedecorations!(ax33; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc33 = scatter!(
    ax33,
    uni_MDS_columns[:,1],
    uni_MDS_columns[:,3];
    markersize = 10,
    color = unimdata.InfantVisAtt,
    colormap = :viridis
)

ax41 = Axis(
    fig[4,1],
    xlabel = "MDS1 ("*string(round(100*spe_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS2 ("*string(round(100*spe_MDS_variances[2]; digits = 2))*"%)",
    title = "Taxonomic Profiles (species)",
)
hidedecorations!(ax41; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc41 = scatter!(
    ax41,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 10,
    color = khula_pci_mbiome_data.InfantVisAtt,
    colormap = :viridis
)

ax42 = Axis(
    fig[4,2],
    xlabel = "MDS2 ("*string(round(100*spe_MDS_variances[2]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*spe_MDS_variances[3]; digits = 2))*"%)",
    title = "Taxonomic Profiles (species)"
)
hidedecorations!(ax42; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc42 = scatter!(
    ax42,
    spe_MDS_columns[:,2],
    spe_MDS_columns[:,3];
    markersize = 10,
    color = khula_pci_mbiome_data.InfantVisAtt,
    colormap = :viridis
)

ax43 = Axis(
    fig[4,3],
    xlabel = "MDS1 ("*string(round(100*spe_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*spe_MDS_variances[3]; digits = 2))*"%)",
    title = "Taxonomic Profiles (species)"
)
hidedecorations!(ax43; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)
sc43 = scatter!(
    ax43,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,3];
    markersize = 10,
    color = khula_pci_mbiome_data.InfantVisAtt,
    colormap = :viridis
)

Colorbar(fig[1:2,4], sc11, vertical = true, label = "Maternal unpredictability")
Colorbar(fig[3:4,4], sc31, vertical = true, label = "Visual Orienting Behavior")

save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS1.png"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS1.eps"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS1.pdf"), fig)