# Exploring cameo

import cameo

from cameo.visualization.plotting.with_plotly import PlotlyPlotter

from cameo.visualization import plotting

#cameo.models.index_models_bigg()

model = cameo.load_model("Combining ME with GEMs/Models/e_coli_core.xml")

fba_result = cameo.fba(model)

fba_result.data_frame

fba_result.display_on_map("iJO1366.Central metabolism")

pfba_result = cameo.pfba(model)

pfba_result.objective_value

abs(fba_result.data_frame.flux).sum()

pfba_result.data_frame

plotter = PlotlyPlotter()

result = cameo.phenotypic_phase_plane(model,
                                variables=[model.reactions.BIOMASS_Ecoli_core_w_GAM],
                                objective=model.reactions.EX_ac_e)
result.plot(plotter)



from cameo.strain_design import OptGene

optgene = OptGene(model)

esult2 = optgene.run(target=model.reactions.EX_ac_e,
                     biomass=model.reactions.BIOMASS_Ecoli_core_w_GAM,
                     substrate=model.metabolites.glc__D_e,
                     max_evaluations=5000,
                     plot=False)


esult2.data_frame


#---
esult2.plot(plotter, 0)



esult2.plot(plotter, 1)

esult2.display_on_map(0, "iJO1366.Central metabolism")

esult2.data_frame


#---

import cobra

iML1515 = cobra.io.read_sbml_model("Combining ME with GEMs/Models/e_coli_flavanone.xml")

len(iML1515.metabolites) 
len(iML1515.reactions) 
len(iML1515.genes)

# 1. FBA of the original iML1515 model

model1 = iML1515.copy()

model1.objective = "4CL"

sol1 = model1.optimize()

print(sol1.fluxes["BIOMASS_Ec_iML1515_core_75p37M"])

print(model1.summary())