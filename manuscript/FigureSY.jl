# Species megaplot - same plot (InfVisAttn x MaternalEntropy) colored by each of the species

## Preliminary: toy correlation analysis

cordf = DataFrame(
        :species => names(select(khula_pci_mbiome_data, Not(mdata_cols))),
        :cor => [
        cor(
            khula_pci_mbiome_data.InfantVisAtt,
            select(khula_pci_mbiome_data, Not(mdata_cols))[:,el]
        ) for el in names(select(khula_pci_mbiome_data, Not(mdata_cols))) ],
        :nzcor => [
        corspearman(
            subset(khula_pci_mbiome_data, el => x -> x .> 0.0).InfantVisAtt,
            subset(khula_pci_mbiome_data, el => x -> x .> 0.0)[:,el]
        ) for el in names(select(khula_pci_mbiome_data, Not(mdata_cols))) ],
        :spe => [
        corspearman(
            khula_pci_mbiome_data.InfantVisAtt,
            select(khula_pci_mbiome_data, Not(mdata_cols))[:,el]
        ) for el in names(select(khula_pci_mbiome_data, Not(mdata_cols))) ],
        :nzspe => [
        corspearman(
            subset(khula_pci_mbiome_data, el => x -> x .> 0.0).InfantVisAtt,
            subset(khula_pci_mbiome_data, el => x -> x .> 0.0)[:,el]
        ) for el in names(select(khula_pci_mbiome_data, Not(mdata_cols))) ],
)

subset(cordf, :cor => x -> .!(isnan.(x)))
sort!(cordf, :cor)

## Actual figure building

spectoplot = [
    "Bifidobacterium_breve",
    "Bifidobacterium_longum",
    "Streptococcus_salivarius",
    "Bifidobacterium_bifidum",
    "Ruminococcus_gnavus",
    "Escherichia_coli",
    "Enterococcus_faecalis",
    "Enterococcus_gallinarum",
    "Enterococcus_avium",
    "Bifidobacterium_kashiwanohense",
    "Veillonella_atypica",
    "Bacteroides_vulgatus",
    "Clostridium_neonatale",
    "Klebsiella_variicola",
    "Collinsella_aerofaciens",
    "Enterococcus_faecium",
    "Lactobacillus_reuteri",
    "Klebsiella_pneumoniae",
    "Bifidobacterium_animalis",
    "Veillonella_parvula",
    "Flavonifractor_plautii",
    "Erysipelatoclostridium_ramosum",
    "Parabacteroides_distasonis",
    "Bifidobacterium_pseudocatenulatum",
    "Lactococcus_lactis",
    "Clostridium_innocuum",
    "Bacteroides_fragilis",
    "Lactobacillus_gasseri",
    "Lactobacillus_vaginalis",
    "Streptococcus_thermophilus",
    "Ruminococcus_torques",
    "Faecalibacterium_prausnitzii",
    "Streptococcus_lutetiensis",
    "Megasphaera_elsdenii",
    "Bacteroides_uniformis",
    "Prevotella_copri",
    "Bacteroides_ovatus",
    "Ruthenibacterium_lactatiformans",
    "Holdemanella_biformis",
    "Klebsiella_michiganensis",
    "Parabacteroides_merdae",
    "Bacteroides_thetaiotaomicron",
    "Veillonella_sp_CAG_933",
    "Lactobacillus_oris",
    "Bifidobacterium_adolescentis",
    "Phascolarctobacterium_succinatutens",
    "Sutterella_wadsworthensis",
    "Anaerostipes_hadrus",
    "Clostridium_perfringens",
    "Bacteroides_caccae",
    "Lactobacillus_ruminis",
    "Catenibacterium_mitsuokai",
    "Fusicatenibacter_saccharivorans",
    "Collinsella_sp_CAG_289"
] 

spectoplot = sort(spectoplot)

using StatsBase

# Compute correlations safely (handles constant vectors or missing)
function safe_pearson(x, y)
    if std(x) ≈ 0 || std(y) ≈ 0
        return NaN
    else
        return cor(x, y)
    end
end

function safe_spearman(x, y)
    if length(unique(x)) < 2 || length(unique(y)) < 2
        return NaN
    else
        return corspearman(x, y)
    end
end

fig = Figure(; size = (2600, 1800))
for irow in 1:6
    for jcol in 1:9
        ispec = (irow-1)*8 + jcol
        ispec > length(spectoplot) && break

        species = spectoplot[ispec]
        x = khula_pci_mbiome_data.MaternalEntropy
        y = khula_pci_mbiome_data.InfantVisAtt
        s = khula_pci_mbiome_data[:, species]

        # Compute correlation stats
        r_p_vis = round(safe_pearson(s, y), digits=2)
        r_p_ent = round(safe_pearson(s, x), digits=2)
        r_s_vis = round(safe_spearman(s, y), digits=2)
        r_s_ent = round(safe_spearman(s, x), digits=2)

        # Text label
        label_str = """
        rₚ(Vis) = $(r_p_vis)
        rₚ(Ent) = $(r_p_ent)
        rₛ(Vis) = $(r_s_vis)
        rₛ(Ent) = $(r_s_ent)
        """

        ax = Axis(fig[irow, jcol], 
            xlabel = "Maternal Entropy", 
            ylabel = "Infant VOB", 
            title = species
        )

        scatter!(ax, x, y, color = log2.(s .+ 0.001))

        # Place label at top right of axis
        Label(fig[irow, jcol], label_str, halign = :left, valign = :top, tellwidth=false, tellheight = false)
    end
end

save(joinpath(Base.pwd(),"manuscript", "figures", "FigureSY.png"), fig)
