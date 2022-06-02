#=
Visualise _escherMaps:
- Julia version: 1.6.0
- Author: Shauny
- Date: 2021-10-21
=#

# a ----------------

println("Visualisation of an Escher Map using Escher.jl")

using Escher, CairoMakie, ColorSchemes
using COBREXA, Tulip
using Clustering

println("Core metabolism of E. coli with fluxes (metabolites/nodes and reactions/edges are colored/sized)")

# use COBREXA to generate a flux distribution using the associated model (in .json format this time)
model = load_model("data/iJO1366-model.json")

#=
Bin fluxes for display purposes - assigning colors to edges needs to be done
manually. The binning uses kmeans clustering on logged fluxes due to the large
differences between fluxes.
=#

fluxes_vec =  flux_balance_analysis_vec(model, Tulip.Optimizer) # FBA

logged_fluxes = log.(abs.(fluxes_vec) .+ 1e-8)

clusters = kmeans(logged_fluxes', 9)  # Colour per cluster? Cluster by size of fluxes?

centers = Dict(j=>i for (i, j) in enumerate(sortperm(clusters.centers'[:])))

order = [centers[i] for i in assignments(clusters)]

rc = Dict(rid => ColorSchemes.RdYlBu_9[10-k] # map reaction id to color
    for (rid, k) in zip(reactions(model), order))


f = Figure(resolution = (1200, 800));
ax = Axis(f[1, 1]);

escherplot!(
    ax,
    "data/iJO1366-map.json";
    reaction_edge_colors = rc,
    metabolite_show_text = true,
    metabolite_text_size = 9, 
    reaction_show_text = false)
hidexdecorations!(ax)
hideydecorations!(ax)
f


#= 

Escherplot attributes that can be used to modify the basic metabolic map figure:

metabolite_identifier = "bigg_id"                     # identifier used to extract metabolite ids
metabolite_show_text = false                          # show the metabolite identifier
metabolite_text_size = 4                              # metabolite identifier text size
metabolite_node_sizes = Dict{String, Any}()           # metabolite id => node size
metabolite_primary_node_size = 5                      # fallback primary metabolite node size
metabolite_secondary_node_size = 3                    # fallback secondary/co-factor metabolite node size
metabolite_node_colors = Dict{String, Any}()          # metabolite id => node color
metabolite_node_color = :black                        # Fallback color used for metabolite nodes
metabolite_text_color = :black                        # Color used for all metabolite text
reaction_identifier = "bigg_id"                       # identifier used to extract reaction ids
reaction_show_text = false                            # show the reaction identifier
reaction_show_name_instead_of_id = false              # show the reaction name instead of id
reaction_text_size = 4                                # reaction identifier text size
reaction_text_color = :black                          # reaction identifier text color
reaction_edge_colors = Dict{String, Any}()            # reaction id => color mapping
reaction_edge_color = :black                          # fallback color in case reaction id not present in reaction_edge_colors
reaction_edge_widths = Dict{String,Any}()             # reaction id => edge size
reaction_edge_width = 2.0                             # fallback width in case reaction id not present in reaction_edge_widths


=#

save("escherplot_3.png", f)  # save the figure/escherplot in a .png file


# b ----------------

println("Use more colours and node/edge sizes to visualise more data on the plot")

# "using ..." see above 2.1

core_model = load_model("data/core-model.json")
sol = flux_balance_analysis_dict(core_model, Tulip.Optimizer)

# Find min and max absolute fluxes for normalization
maxflux = maximum(abs.(values(sol)))
minflux = minimum(abs.(values(sol)))

# Scale width of reaction edges to fluxes
width_interp(x) = 2 + 5 * (abs(x) - minflux) / (maxflux - minflux) # widths between 2 and 5
re = Dict(k => width_interp(v) for (k, v) in sol) # map reaction id to reaction edge width

# Scale color of reaction edges to fluxes (manually binned)
color_interp(x) = begin
    normed_x = (abs(x) - minflux) / (maxflux - minflux)
    if 0 <= normed_x < 0.01
        ColorSchemes.RdYlBu_4[4]
    elseif 0.01 <= normed_x < 0.25
        ColorSchemes.RdYlBu_4[3]
    elseif 0.25 <= normed_x < 0.5
        ColorSchemes.RdYlBu_4[2]
    else
        ColorSchemes.RdYlBu_4[1]
    end
end
rc = Dict(k => color_interp(v) for (k, v) in sol) # map reaction id to reaction edge color

# metabolite node colors
mc = Dict(
    k => ColorSchemes.:Dark2_7[v] for
    (k, v) in zip(metabolites(model), rand(1:7, n_metabolites(model)))
)

# metabolite node sizes
ms = Dict(k => v for (k, v) in zip(metabolites(model), rand(3:10, n_metabolites(model))))

# Normal Makie plotting features all work (escherplot is a full recipe)
f = Figure(resolution = (1200, 800));
ax = Axis(f[1, 1]);
escherplot!(
    ax,
    "data/core-map.json";
    reaction_edge_widths = re,
    reaction_edge_colors = rc,
    metabolite_node_colors = mc,
    metabolite_node_sizes = ms,
)
hidexdecorations!(ax)
hideydecorations!(ax)
f


# c ----------------

# If incomplete reaction data (edge colors or widths) are supplied the missing reaction edges are dotted.

rc = Dict("FBA" => ColorSchemes.RdYlBu_4[4],
    "PFK" => ColorSchemes.RdYlBu_4[3],
    "PEP" => ColorSchemes.RdYlBu_4[2],
    "PYK" => ColorSchemes.RdYlBu_4[1])

escherplot!(
  ax,
  "data/core-map.json";
    reaction_edge_colors = rc,
    metabolite_show_text = false, 
    reaction_show_text = false)
    
hidexdecorations!(ax)
hideydecorations!(ax)
f

