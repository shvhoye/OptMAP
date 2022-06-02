''' 
Minimal_metabolism:
- Author: Shauny
- Date: 2021-11-23
'''

print("Building and simulating a model based on the paper 'Minimal Pathway for the Regeneration of Redox Cofactors'")

import cobra

from cobra import Model, Reaction, Metabolite

#------------------------------------

model = Model('Minimal_metabolism_model')

## Add the 3 main reactions to the model: FDH, NADTRHD, GDR2

# 1. FDH

reaction1 = Reaction('FDH')
reaction1.name = 'Formate dehydrogenase' 
reaction1.subsystem = 'Cytosol' 
reaction1.lower_bound = 0.
reaction1.upper_bound = 1000. 

for_c = Metabolite(
    'for_c',
    formula='CH1O2',
    name='cytosolic formate',
    compartment='c') 

nad_c = Metabolite(
    'nad_c',
    formula='C21H26N7O14P2',
    name='cytosolic nad',
    compartment='c')

co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='cytosolic co2',
    compartment='c')  

nadh_c = Metabolite(
    'nadh_c',
    formula='C21H27N7O14P2',
    name='cytosolic nadh',
    compartment='c')
    
reaction1.add_metabolites({
    for_c: -1.0,
    nad_c: -1.0,
    co2_c: 1.0,
    nadh_c: 1.0

})     

reaction1.reaction # This gives a string representation of the reaction

model.add_reactions([reaction1])


# 2. NADTRHD

reaction2 = Reaction('NADTRHD')
reaction2.name = 'NAD transhydrogenase' 
reaction2.subsystem = 'Cytosol' 
reaction2.lower_bound = -1000. 
reaction2.upper_bound = 0. 

nad_c = Metabolite(
    'nad_c',
    formula='C21H26N7O14P2',
    name='cytosolic nad',
    compartment='c')

nadh_c = Metabolite(
    'nadh_c',
    formula='C21H27N7O14P2',
    name='cytosolic nadh',
    compartment='c')

nadp_c = Metabolite(
    'nadp_c',
    formula='C21H25N7O17P3',
    name='cytosolic nadp',
    compartment='c')

nadph_c = Metabolite(
    'nadph_c',
    formula='C21H26N7O17P3',
    name='cytosolic nadph',
    compartment='c')
    
reaction2.add_metabolites({
    nad_c: -1.0,
    nadph_c: -1.0,
    nadh_c: 1.0,
    nadp_c: 1.0
})     

reaction2.reaction 

model.add_reactions([reaction2])


# 3. GDR2

reaction3 = Reaction('GDR2')
reaction3.name = 'Glutathione-disulfide reductase' 
reaction3.subsystem = 'Cytosol' 
reaction3.lower_bound = 0. 
reaction3.upper_bound = 1000. 

nadp_c = Metabolite(
    'nadp_c',
    formula='C21H25N7O17P3',
    name='cytosolic nadp',
    compartment='c')

nadph_c = Metabolite(
    'nadph_c',
    formula='C21H26N7O17P3',
    name='cytosolic nadph',
    compartment='c')

h_c = Metabolite(
    'h_c',
    formula='H',
    name='cytosolic h',
    compartment='c')

gthox_c = Metabolite(
    'gthox_c',
    formula='C20H30N6O12S2',
    name='cytosolic gthox',
    compartment='c')

gthrd_c = Metabolite(
    'gthrd_c',
    formula='C10H16N3O6S',
    name='cytosolic gthrd',
    compartment='c')

reaction3.add_metabolites({
    gthox_c: -1.0,
    nadph_c: -1.0,
    h_c: -1.0,
    gthrd_c: 2.0,
    nadp_c: 1.0
})     

reaction3.reaction 

model.add_reactions([reaction3])


## Add extracellular <-> cytosol reactions for forA and co2

# Formic acid

reaction4 = Reaction('R_forA_e_to_c')
reaction4.name = 'transfer of formic acid from the extracellular space to the cytosol' 
reaction4.subsystem = 'cm' 
reaction4.lower_bound = 0. 
reaction4.upper_bound = 1000. 

forA_e = Metabolite(
    'forA_e',
    formula='CH2O2',
    name='extracellular formic acid',
    compartment='e') 

forA_c = Metabolite(
    'forA_c',
    formula='CH2O2',
    name='cytosolic formic acid',
    compartment='c') 
    
reaction4.add_metabolites({
    forA_e: -1.0,
    forA_c: 1.0
})     

reaction4.reaction 

model.add_reactions([reaction4])

# CO2

reaction5 = Reaction('R_co2_e_to_c')
reaction5.name = 'transfer of co2 from the cytosol to the extracellular space' 
reaction5.subsystem = 'cm' 
reaction5.lower_bound = 0. 
reaction5.upper_bound = 1000. 

co2_e = Metabolite(
    'co2_e',
    formula='CO2',
    name='extracellular co2',
    compartment='e') 
    
co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='cytosolic co2',
    compartment='e')       

    
reaction5.add_metabolites({
    co2_c: -1.0,
    co2_e: 1.0
})

reaction5.reaction # This gives a string representation of the reaction

model.add_reactions([reaction5])


## Add formic acid dissociation reaction

reaction6 = Reaction('R_formic_acid_dissociation_reaction')
reaction6.name = 'formic acid dissociation reaction' 
reaction6.subsystem = 'c' 
reaction6.lower_bound = 0. 
reaction6.upper_bound = 1000. 

for_c = Metabolite(
    'for_c',
    formula='CH1O2',
    name='cytosolic formate',
    compartment='c')

h_c = Metabolite(
    'h_c',
    formula='H',
    name='cytosolic h',
    compartment='c')

forA_c = Metabolite(
    'forA_c',
    formula='CH2O2',
    name='cytosolic formic acid',
    compartment='c')
    
reaction6.add_metabolites({
    forA_c: -1.0,
    for_c: 1.0,
    h_c: 1.0
})

reaction6.reaction 

model.add_reactions([reaction6])


## Add exhange reactions so formic acid can be taken up and co2 can be released into the environment. 
# Exchange (uptake) reactions 

model.add_boundary(model.metabolites.get_by_id("co2_e"), type="exchange")

model.add_boundary(model.metabolites.get_by_id("forA_e"), type="exchange")


## Demand reactions for all the metabolites that aren't formed in other reactions

# 1.0 mM NAD+
# 0.2 mM NADP 
# 2.5 mM GSSG

medium = model.medium

#gthox
model.add_boundary(model.metabolites.get_by_id("gthox_c"), type="demand")
medium["DM_gthox_c"] = 2.5

#nad_c
model.add_boundary(model.metabolites.get_by_id("nad_c"), type="demand")
medium["DM_nad_c"] = 1.0

#nadp_c
model.add_boundary(model.metabolites.get_by_id("nadp_c"), type="demand")
medium["DM_nadp_c"] = 0.2

#gthrd_c
model.add_boundary(model.metabolites.get_by_id("gthrd_c"), type="demand")
medium["DM_gthrd_c"] = 0

"""
#forA_e
model.add_boundary(model.metabolites.get_by_id("forA_e"), type="demand")
medium["DM_forA_e"] = 1000

#h_c
model.add_boundary(model.metabolites.get_by_id("h_c"), type="demand")
medium["DM_h_c"] = 1000

#forA_c
model.add_boundary(model.metabolites.get_by_id("forA_c"), type="demand")
medium["DM_forA_c"] = 1000

#for_c
model.add_boundary(model.metabolites.get_by_id("for_c"), type="demand")
medium["DM_for_c"] = 1000

#co2_c
model.add_boundary(model.metabolites.get_by_id("co2_c"), type="demand")
medium["DM_co2_c"] = 1000

#nadh_c
model.add_boundary(model.metabolites.get_by_id("nadh_c"), type="demand")
medium["DM_nadh_c"] = 1000

#nadph_c
model.add_boundary(model.metabolites.get_by_id("nadph_c"), type="demand")
medium["DM_nadph_c"] = 1000

#co2_e
model.add_boundary(model.metabolites.get_by_id("co2_e"), type="demand")
medium["DM_co2_e"] = 1000
"""

## Define the "medium" outside of the vesicle

# Add only formic acid

# Formic acid concentration from paper =  5 mM external formate 

medium["EX_forA_e"] = 5
medium["EX_co2_e"] = 0.0

model.medium = medium

model.medium

## Reaction bounds (Default values should already be ok for FDH, NADTRHD and GDR2)

#model.reactions.FDH.bounds = (0, 1000)

#model.reactions.NADTRHD.bounds = (-1000, 0)

#model.reactions.GDR2.bounds = (0, 1000)

#model.reactions.R_co2_e_to_c.bounds =  (0, 1000)

#model.reactions.R_forA_e_to_c.bounds = (-1000, 0)

#------------------------------------

## The objective is to create nicotinamide cofactors and reduce glutathione disulfide

# NADH and NADPH from an externally provided formate 

# For the objective we just want this to be the maximization of the flux in the single reaction we added

# We demonstrate that the downstream biochemical process of 
# reduction of glutathione disulfide can be driven by the transfer of 
# reducing equivalents from formate via NAD(P)H, thereby providing a versatile 
# set of electron donors for reductive metabolism.

model.objective = 'FDH'

## Run a FBA (= Flux balance analysis)

sol = model.optimize()

print(sol.fluxes)

print(model.summary())

print(model.metabolites.nadh_c.summary())

#------------------------------------

## DYNAMIC FLUX BALANCE ANALYSIS (dFBA) with the minimal metabolic model

# The model considers only basic Michaelis-Menten limited growth on formic acid 

# 1. Set up the dynamic system

import numpy as np

from tqdm import tqdm

from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt 

#----

# gthrd and forA

# Define new functions

def add_dynamic_bounds(model, y):
    """Use external concentrations to bound the uptake flux of formic acid.""" 
    gthrd, formic_acid  = y   

    # Michaelis−Menten equation: v = VMAX*[S]/ (KM + [S])
    # kCAT = Vmax/[enzyme] -> Vmax = kCAT*[enzyme] 

    # Fdh has a relatively high affinity for formate (KM = 2.15 mM)    
    # 2.0 μM Fdh, 0.08 μM SthA and 0.05 μM GorA


    #                                                           KM (mM)     kCAT (1/s)   Vmax(s)
    # Fdh -> formate:NAD+ oxidoreductase          
    #                                            NAD+           0.11         1.08         2.16
    #                                            Formate        2.15         0.87         1.74
    # SthA -> NADPH:NAD+ oxidoreductase          
    #                                            NADH           2.63         9.7          0.78
    #                                            NADP+          0.03         19.9         1.60
    # GorA -> glutathione:NADP+ oxidoreductase         
    #                                            GSSG           0.07         733.3        36.67
    #                                            NADPH          0.02         661.8        33.09
   
    
    formic_acid_max_import = -1.74 * formic_acid / (2.15 + formic_acid)  
    
    model.reactions.EX_forA_e.lower_bound = formic_acid_max_import


def dynamic_system(t, y):
    """Calculate the time derivative of external species."""
    
    gthrd, formic_acid  = y 

    # Calculate the specific exchanges fluxes at the given external concentrations.
    with model: 
        add_dynamic_bounds(model, y)
        
        cobra.util.add_lp_feasibility(model)
        feasibility = cobra.util.fix_objective_as_constraint(model)
        
        lex_constraints = cobra.util.add_lexicographic_constraints(
            model, ['DM_gthrd_c', 'EX_forA_e'], ['max', 'max'])          
    
    # Since the calculated fluxes are specific rates, we multiply them by the 
    # gthrd concentration to get the bulk exchange rates.
    fluxes = lex_constraints.values
    fluxes *= gthrd
    
    # This implementation is **not** efficient, so I display the current # simulation time using a progress bar.
    if dynamic_system.pbar is not None:
        dynamic_system.pbar.update(1) 
        dynamic_system.pbar.set_description('t = {:.3f}'.format(t))
    
    return fluxes 

dynamic_system.pbar = None


def infeasible_event(t, y): 
    """
    Determine solution feasibility.

    Avoiding infeasible solutions is handled by solve_ivp's built-in event detection.
    This function re-solves the LP to determine whether or not the solution is feasible
    (and if not, how far it is from feasibility). When the sign of this function changes
    from -epsilon to positive, we know the solution is no longer feasible. 
    
    """

    with model:
        add_dynamic_bounds(model, y)
    
        cobra.util.add_lp_feasibility(model)
        feasibility = cobra.util.fix_objective_as_constraint(model) 
    
    return feasibility - infeasible_event.epsilon


infeasible_event.epsilon = 1E-6 
infeasible_event.direction = 1 
infeasible_event.terminal = True



# 2. Run the dynamic FBA simulation

ts = np.linspace(0, 6, 100) # Desired integration resolution and interval 

y0 = [0.001, 5] 

with tqdm() as pbar: 
    dynamic_system.pbar = pbar
    
    sol = solve_ivp(
        fun=dynamic_system,
        events=[infeasible_event],
        t_span=(ts.min(), ts.max()),
        y0=y0,
        t_eval=ts,
        rtol=1e-6,
        atol=1e-8,
        method='BDF')


print(sol)


# 3. Plot the results

ax = plt.subplot(111)
ax.plot(sol.t, sol.y.T[:, 0], color='b')
ax2 = plt.twinx(ax)
ax2.plot(sol.t, sol.y.T[:, 1], color='r')
ax.set_ylabel('GSH (mM)', color='b') # GSh ~ gthrd
ax2.set_ylabel('Formic acid (mM)', color='r')

plt.show()






#------------------------------------

# Model Validation

import tempfile
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model 
with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
    write_sbml_model(model, filename=f_sbml.name)
    report = validate_sbml_model(filename=f_sbml.name)
pprint(report)

# The model is valid with no COBRA or SBML errors or warnings

#------------------------------------

# Check the model

# The objects have been added to the model

print(f'{len(model.reactions)} reactions') 
print(f'{len(model.metabolites)} metabolites') 


# Iterate through the the objects in the model
print("Reactions")

for x in model.reactions:
    print("%s : %s" % (x.id, x.reaction))

print("Metabolites") 

for x in model.metabolites:
    print('%9s : %s' % (x.id, x.formula))


#------------------------------------

## Gapfilling

from cobra.flux_analysis import gapfill

universal = cobra.Model("universal_reactions")
   
model.objective = model.add_boundary(model.metabolites.co2_c, type='demand') 

solution = gapfill(model, universal)

for reaction in solution[0]:
        print(reaction.id)

print(solution)

#------------------------------------

## Save the model

#cobra.io.write_sbml_model(model, "Minimal_Model.xml")
#cobra.io.save_json_model(model, "Minimal_Model.json")
