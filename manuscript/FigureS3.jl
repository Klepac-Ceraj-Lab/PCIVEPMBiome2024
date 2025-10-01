#####
# K-means clustering stuff
#####

Random.seed!(0)

function silhouette_score(X, labels)
    # Calculate the pairwise distance matrix (Euclidean distances in this example)
    D = pairwise(Euclidean(), X)

    # Initialize silhouette scores
    scores = Float64[]
    n = size(X, 1)

    for i in 1:n
        # Get distances to all other points in the same cluster (a_i)
        same_cluster = findall(labels .== labels[i])
        a_i = mean(D[i, same_cluster[2:end]])  # Exclude self in calculation

        # Get the minimum average distance to points in other clusters (b_i)
        other_clusters = unique(labels[labels .!= labels[i]])
        b_i = minimum([mean(D[i, findall(labels .== c)]) for c in other_clusters])

        # Calculate silhouette score for point i
        s_i = (b_i - a_i) / max(a_i, b_i)
        push!(scores, s_i)
    end

    # Return the average silhouette score for all points
    return mean(scores)
end

# Function to calculate Within-Cluster Sum of Squares (WCSS) for k-means
function calculate_wcss(X, k)
    result = kmeans(X, k)
    return sum(result.costs)
end

# Generate the elbow plot for a range of k values (1 to max_k)
max_k = 10  # You can adjust this based on your data

wcss = [calculate_wcss(collect(uni_MDS_results.U[:, [1,3]]'), k) for k in 2:max_k]

silhouette_scores = []
for k in 2:max_k 
    kmeans_result = kmeans(uni_MDS_results.U[:, [1,3]]', k)  # Run k-means clustering
    push!(silhouette_scores, silhouette_score(uni_MDS_results.U[:, [1,3]]', kmeans_result.assignments))
end

# Plot WCSS vs. number of clusters (elbow plot)

fig = Figure(size=(1200, 400))

AxAelbow = Axis(
    fig[1,1],
    xlabel = "Number of Clusters",
    ylabel = "Within-cluster sum of squares (WCSS)",
    title = "Clustering decision measures",
    xticks = 1:9
)

AxAsilhouette = Axis(
    fig[1,1],
    xlabel = "Number of Clusters",
    ylabel = "Silhouette score",
    xticks = 1:9,
    yaxisposition = :right
)

hidedecorations!(AxAelbow, ticklabels = false, label = false)
hidedecorations!(AxAsilhouette, ticklabels = false, label = false)

linkxaxes!(AxAelbow, AxAsilhouette)

lnAwcss = lines!(
    AxAelbow,
    collect(1:9),
    wcss,
    color = :teal,
    linewidth = 5
)

lnAsil = lines!(
    AxAsilhouette,
    collect(1:9),
    silhouette_scores,
    color = :olivedrab,
    linewidth = 5
)

axislegend(AxAelbow, [lnAwcss, lnAsil], ["WCSS", "Silhouette"], "Measure", position = :rt,
    orientation = :vertical)

## FROM THE WCSS/SILHOUETTE...
elbow_point = 4
final_clu = kmeans(collect(uni_MDS_results.U[:, [1,3]]'), elbow_point)
##

axB = Axis(
    fig[1,2],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)",
    aspect = 1.0
)
hidedecorations!(
    axB;
    label = false,
    ticklabels = false,
    ticks = false,
    grid = true,
    minorgrid = false,
    minorticks = false
)
scB = scatter!(
    axB,
    uni_MDS_columns[:,1],
    uni_MDS_columns[:,3];
    markersize = 12,
    color = final_clu.assignments,
    colormap = :berlin
)

# Calculate means and standard deviations for each cluster

clu_means = Float64[]
clu_stdns = Float64[]

for i in 1:elbow_point
    cluster_data = uniref_color_mdata.InfantVisAtt[final_clu.assignments .== i, :]
    mean_cluster = mean(cluster_data, dims=1)
    std_dev_cluster = std(cluster_data, dims=1)
    push!(clu_means, mean_cluster[1,1])
    push!(clu_stdns, std_dev_cluster[1,1])
    println("Cluster $i Mean: $mean_cluster, Standard Deviation: $std_dev_cluster")
end

axC = Axis(
    fig[1,3],
    xlabel = "Cluster",
    ylabel = "Mean VOB",
    title = "Cluster differences"
)

hidedecorations!(
    axC;
    label = false,
    ticklabels = false,
    ticks = false,
    grid = true,
    minorgrid = false,
    minorticks = false
)

bpC = barplot!(
    axC,
    1:4,
    clu_means,
    color = 1:4,
    colormap = :berlin
)

rangebars!(
    axC,
    1:4,
    clu_means .- clu_stdns/2,
    clu_means .+ clu_stdns/2,
    color = :black,
    whiskerwidth = 10
)

## Layout and Export

Label(fig[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[1, 2, TopLeft()], "b", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(fig[1, 3, TopLeft()], "c", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())

colsize!(fig.layout, 1, Relative(0.48))
colsize!(fig.layout, 2, Relative(0.35))
colsize!(fig.layout, 3, Relative(0.17))

save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS3.png"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS3.eps"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "FigureS3.pdf"), fig)
