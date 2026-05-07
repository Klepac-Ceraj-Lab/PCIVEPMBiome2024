## PCI Topology and Bifidobacteria
## Not included in final submission

bifido_cdict = Dict(
    "dom_breve" => (:purple4, 1.0),
    "more_breve" => (:darkorchid2, 0.8),
    "null" => (:gray, 0.6),
    "more_longum" => (:orange, 0.8),
    "dom_longum" => (:darkorange2, 1.0)
)

function assign_bifido_class(tt::T where T<: Tuple)

    if (tt[1] > 0.75)
        return "dom_breve"
    elseif (tt[2] > 0.75)
        return "dom_longum"
    elseif ( (tt[1] > 0.05) & ( ( (tt[1] + 0.01) / (tt[2] + 0.01) ) > 5 ) )
        return "more_breve"
    elseif ( (tt[2] > 0.05) & ( ( (tt[2] + 0.01) / (tt[1] + 0.01) ) > 5 ) )
        return "more_longum"
    else
        return "null"
    end

end

function assign_bifido_class(df::DataFrame)

    zipobj = zip( df[:, "Bifidobacterium_breve"], df[:, "Bifidobacterium_longum"] )

    map(assign_bifido_class, zipobj)

end

fig = Figure(; size = (800, 600))

ax = Axis(
    fig[1,1],
    xlabel = "Maternal unpredictability",
    ylabel = "Visual Orienting Behavior (VOB) x10^-5",
    yticks = ([0.0, 0.00005, 0.00010], ["0.0", "5.0", "10.0"]),
    title = "Bifidobacteria across VOB x unpredictability levels",
    aspect = 1.0
)
hidedecorations!(
    ax;
    label = false,
    ticklabels = false,
    ticks = false,
    grid = true,
    minorgrid = false,
    minorticks = false
)
scD = scatter!(
    ax,
    khula_pci_mbiome_data.MaternalEntropy,
    khula_pci_mbiome_data.InfantVisAtt;
    markersize = 12,
    color = [ bifido_cdict[el] for el in assign_bifido_class(khula_pci_mbiome_data) ],
    # color = khula_pci_mbiome_data.Bifidobacterium_breve,
    colormap = :viridis
)

lgC = Legend(
    fig[2,1],
    [
        MarkerElement(; marker=:circle, color=bifido_cdict["dom_breve"], strokewidth=1),
        MarkerElement(; marker=:circle, color=bifido_cdict["more_breve"], strokewidth=1),
        MarkerElement(; marker=:circle, color=bifido_cdict["more_longum"], strokewidth=1),
        MarkerElement(; marker=:circle, color=bifido_cdict["dom_longum"], strokewidth=1),
        MarkerElement(; marker=:circle, color=bifido_cdict["null"], strokewidth=1)
    ],
    [
        rich("Dominated (>70%) by ", rich("Bifidobacterium breve", font = :italic)),
        rich(rich("B. breve", font = :italic), " above 5% and ", rich("B. breve", font = :italic), "-to-", rich("B. longum", font = :italic), " ratio > 5.0" ),
        rich(rich("B. longum", font = :italic), " above 5% and ", rich("B. longum", font = :italic), "-to-", rich("B. breve", font = :italic), " ratio > 5.0" ),
        rich("Dominated (>70%) by ", rich("Bifidobacterium longum", font = :italic)),
        rich("Not characterized by ", rich("B. breve", font = :italic), " or ", rich("B. longum", font = :italic))
    ];
    tellheight = false,
    tellwidth = false,
    nbanks = 1
)

fig

save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS3.png"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS3.eps"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS3.pdf"), fig)
