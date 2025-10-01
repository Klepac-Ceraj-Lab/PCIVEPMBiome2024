## Color by Site
fig = Figure(; size = (1200, 800))

# tempcol = sort(leftjoin(khula_pci_mbiome_data, select(godhelpme, [:sample, :l1_error]), on = :sample), :sample).l1_error
# tempcol[ismissing.(tempcol)] .= 0.0
# tempcol = collect(skipmissing(tempcol))

ax1 = Axis(
    fig[1,1],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = rich("Bifidobacterium longum", font = :italic)
)
sc1 = scatter!(
    ax1,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 10,
    color = khula_pci_mbiome_data.Bifidobacterium_longum,
    # color = tempcol#,
    colormap = :PuBu
)

ax2 = Axis(
    fig[1,2],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = rich("Bifidobacterium breve", font = :italic)
)
sc2 = scatter!(
    ax2,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 10,
    color = khula_pci_mbiome_data.Bifidobacterium_breve,
    colormap = :PuBu
)

ax3 = Axis(
    fig[1,3],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = rich("Bifidobacterium bifidum", font = :italic)
)
sc3 = scatter!(
    ax3,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 10,
    color = khula_pci_mbiome_data.Bifidobacterium_bifidum,
    colormap = :PuBu
)

Colorbar(fig[2,1], sc1, vertical = false, label = rich("Bifidobacterium longum", font = :italic))
Colorbar(fig[2,2], sc2, vertical = false, label = rich("Bifidobacterium breve", font = :italic))
Colorbar(fig[2,3], sc3, vertical = false, label = rich("Bifidobacterium bifidum", font = :italic))

ax4 = Axis(
    fig[3,1],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = rich("Ruminococcus gnavus", font = :italic)
)
sc4 = scatter!(
    ax4,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 10,
    color = khula_pci_mbiome_data.Ruminococcus_gnavus,
    colormap = :PuBu
)

ax5 = Axis(
    fig[3,2],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = rich("Erysipelatoclostridium ramosum", font = :italic)
)
sc5 = scatter!(
    ax5,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 10,
    color = khula_pci_mbiome_data.Erysipelatoclostridium_ramosum,
    colormap = :PuBu
)

ax6 = Axis(
    fig[3,3],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = rich("Escherichia coli", font = :italic)
)
sc6 = scatter!(
    ax6,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 10,
    color = khula_pci_mbiome_data.Escherichia_coli,
    colormap = :PuBu
)

Colorbar(fig[4,1], sc4, vertical = false, label = rich("Ruminococcus gnavus", font = :italic))
Colorbar(fig[4,2], sc5, vertical = false, label = rich("Erysipelatoclostridium ramosum", font = :italic))
Colorbar(fig[4,3], sc6, vertical = false, label = rich("Escherichia coli", font = :italic))

Label(fig[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[1, 2, TopLeft()], "b", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[1, 3, TopLeft()], "c", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[3, 1, TopLeft()], "d", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[3, 2, TopLeft()], "e", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[3, 3, TopLeft()], "f", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())


save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS2.png"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS2.eps"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS2.pdf"), fig)