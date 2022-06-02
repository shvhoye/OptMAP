#=
BiGG_models_trial1:
- Julia version: 1.6.0
- Author: Shauny
- Date: 2021-10-05
=#

#= 

# Downloading models from BiGG via julia

download("http://bigg.ucsd.edu/static/models/iJO1366.xml", "iJO1366.xml")

=#

# 1.1 ----------------

println("Load BiGG model iJO1366.xml into julia via COBREXA.jl")

using COBREXA   # loads the package
using Tulip     # loads the optimization solver
using DataFrames, CSV

# Open the SBML file and load the contents
SBML_model = load_model("julia/iJO1366.xml")

# What does the model contain?
metabolites(SBML_model)
n_metabolites(SBML_model) # 1805
n_genes(SBML_model)
reactions(SBML_model)
# ...

# Run a FBA (= Flux balance analysis)
fluxes = flux_balance_analysis_dict(SBML_model, Tulip.Optimizer) # The variable fluxes will now contain 
# a dictionary of the computed optimal flux of each reaction in the Model

flux_df = DataFrame(fluxes)
# CSV.write("flux_df", flux_df)

vscodedisplay(flux_df)


# 1.2 ----------------

# COBREXA

# Model variant processing

#=
The main feature of COBREXA.jl is the ability to easily 
specify and process many analyses in parallel. To 
demonstrate, let's see how the organism would perform if
some reactions were disabled independently:
=#

# convert to a model type that is efficient to modify
m = convert(StandardModel, SBML_model)

# find the model objective value if oxygen or carbon dioxide transports are disabled
screen(m, # the base model
    variants=[ # this specifies how to generate the desired model variants
        [], # one with no modifications, i.e. the base case
        [with_changed_bound("O2t", lower=0.0, upper=0.0)], # disable oxygen
        [with_changed_bound("CO2t", lower=0.0, upper=0.0)], # disable CO2
        [with_changed_bound("O2t", lower=0.0, upper=0.0),
	        with_changed_bound("CO2t", lower=0.0, upper=0.0)], # disable both
    ],
    # this specifies what to do with the model variants (received as the argument `x`)
    analysis = x ->
        flux_balance_analysis_dict(x, Tulip.Optimizer)["BIOMASS_iJO1366_w_GAM"],
)

# You should receive a result showing that missing oxygen transport makes the biomass production much harder
# Most importantly, such analyses can be easily specified by automatically generating long lists of modifications to be applied to the model, and parallelized.



#Knocking out each reaction in the model is efficiently accomplished:



#???



















# 1.3 ----------------

# Visualisation trial using iJO1366.xml

fluxes_SBML_vec = flux_balance_analysis_vec(SBML_model, Tulip.Optimizer) 

logged_fluxes_SBML = log.(abs.(fluxes_SBML_vec) .+ 1e-8)

clusters_SBML = kmeans(logged_fluxes_SBML', 9)  # Colour per cluster? Cluster by size of fluxes?

centers_SBML = Dict(j=>i for (i, j) in enumerate(sortperm(clusters_SBML.centers'[:])))

order_SBML = [centers_SBML[i] for i in assignments(clusters_SBML)]

rc_SBML = Dict(rid => ColorSchemes.RdYlBu_9[10-k] # map reaction id to color
    for (rid, k) in zip(reactions(SBML_model), order_SBML))


f_SBML = Figure(resolution = (1200, 800));
ax_SBML = Axis(f[1, 1]);
escherplot!(
    ax_SBML,
    "iJO1366-map.xml"; # doesn't work: you need a specific .json file of the map, not just the model! look at "data" directory
    reaction_edge_colors = rc_SBML,
    metabolite_show_text = false, 
    reaction_show_text = false)
hidexdecorations!(ax_SBML)
hideydecorations!(ax_SBML)
f_SBML



# 2.1 ----------------

println("Visualisation of the model using Escher.jl (plots maps of metabolic models) + ...")

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
    metabolite_show_text = false,
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


# 2.2 ----------------

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


# 2.3 ----------------

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


