#####
# Not included in the original submission
#####

figsx = Figure(size = (1200, 500))

axSXA = Axis(
    figsx[1,1],
    xlabel = "Maternal unpredictability",
    ylabel = "Visual Orienting Behavior (VOB)",
    title = "Bifidobacterium breve"
)

axSXB = Axis(
    figsx[1,2],
    xlabel = "Maternal unpredictability",
    ylabel = "Visual Orienting Behavior (VOB)",
    title = "Bifidobacterium longum"
)

axSXC = Axis(
    figsx[1,3],
    xlabel = "Maternal unpredictability",
    ylabel = "Visual Orienting Behavior (VOB)",
    title = "Bifidobacterium class"
)

scXA = scatter!(
    axSXA,
    khula_pci_mbiome_data.MaternalEntropy,
    khula_pci_mbiome_data.InfantVisAtt;
    markersize = 12,
    color = khula_pci_mbiome_data.Bifidobacterium_breve,
    colormap = :viridis
)

scXB = scatter!(
    axSXB,
    khula_pci_mbiome_data.MaternalEntropy,
    khula_pci_mbiome_data.InfantVisAtt;
    markersize = 12,
    color = khula_pci_mbiome_data.Bifidobacterium_longum,
    colormap = :viridis
)

scXC = scatter!(
    axSXC,
    khula_pci_mbiome_data.MaternalEntropy,
    khula_pci_mbiome_data.InfantVisAtt;
    markersize = 12,
    color = [ bifido_cdict[el] for el in assign_bifido_class(khula_pci_mbiome_data) ],
    colormap = :viridis
)

Colorbar(figsx[2,1], scXA, vertical = false)
Colorbar(figsx[2,2], scXB, vertical = false)

save(joinpath(Base.pwd(),"manuscript", "figures", "FigureSZ.png"), fig)