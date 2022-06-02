#=
Visualise _Minimal_Model:
- Julia version: 1.6.0
- Author: Shauny
- Date: 2021-11-23
=#

# a ----------------

println("Visualisation of the minimal model using Escher.jl")

using Escher, CairoMakie, ColorSchemes
using COBREXA, Tulip
using Clustering

# use COBREXA to generate a flux distribution using the associated model (in .json format this time)
model = load_model("Minimal_Model.json")

fluxes_vec =  flux_balance_analysis_vec(model, Tulip.Optimizer) # FBA

logged_fluxes = log.(abs.(fluxes_vec) .+ 1e-8)

clusters = kmeans(logged_fluxes', 3)  # Colour per cluster? Cluster by size of fluxes?

centers = Dict(j=>i for (i, j) in enumerate(sortperm(clusters.centers'[:])))

order = [centers[i] for i in assignments(clusters)]

rc = Dict(rid => ColorSchemes.RdYlBu_9[10-k] # map reaction id to color
    for (rid, k) in zip(reactions(model), order))


f = Figure(resolution = (1200, 800));
ax = Axis(f[1, 1]);

escherplot!(
    ax,
    "???";
    reaction_edge_colors = rc,
    metabolite_show_text = false,
    metabolite_text_size = 9, 
    reaction_show_text = false)
hidexdecorations!(ax)
hideydecorations!(ax)
f



#save("escherplot_Minimal_Model.png", f)  # save the figure/escherplot in a .png file
