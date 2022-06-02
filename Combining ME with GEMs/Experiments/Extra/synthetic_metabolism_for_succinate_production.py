# Combining evolutionary algorithms with genome scale modeling for intelligent metabolic engineering - Synthetic metabolism for succinate production
# - Author: Shauny Van Hoye
# - Date: 2022-05-09

import random
import numpy as np
import cobra
import matplotlib.pyplot as plt
import seaborn as sns

## MAP-Elites functions 

def performance(x, model, objective): 
	
    # Evaluate the performance/fitness of a solution

    indices = np.where(x == 1)[0] 

    temp_model = model.copy()

    if sum(x) > 0:

        for i in indices:
            temp_model.genes[i].knock_out()

    temp_model.objective = "BIOMASS_Ec_iML1515_core_75p37M" 
    
    sol = temp_model.optimize()

    Performance =  sol.fluxes["BIOMASS_Ec_iML1515_core_75p37M"]*sol.fluxes[objective]
            
    return Performance


def niche(x):

    niche = []

    n_mutations = int(sum(x))

    if n_mutations > 10:
        niche.append(10) # if there are more than 10 mutations, they are stored in the 10th niche
    else:
        niche.append(n_mutations)

    """

    sugar = np.load("Combining ME with GEMs/Supporting scripts and data/data/sugar.npy")

    carbohydrate = np.load("Combining ME with GEMs/Supporting scripts and data/data/carbohydrate.npy")

    nucleotide = np.load("Combining ME with GEMs/Supporting scripts and data/data/nucleotide.npy")

    lipid = np.load("Combining ME with GEMs/Supporting scripts and data/data/lipid.npy")

    AA = np.load("Combining ME with GEMs/Supporting scripts and data/data/AA.npy")

    s = 0

    c = 0

    n = 0

    l = 0

    A = 0

    indices = np.where(x == 1)[0] 

    for i in indices:

        s += sugar[i]

        c += carbohydrate[i]

        n += nucleotide[i]

        l += lipid[i]

        A += AA[i]

    niche.append(int(np.argmax([s, c, n, l, A])))

    """

    niche.append(np.random.choice([0, 1, 2, 3, 4]))

    return niche
	

def random_selection(MAP_solutions, n):
	
    x1 = np.zeros(n) 

    x2 = np.zeros(n) 
    
    while (any(x1) == 0 or any(x2) == 0):
		
        x1 = random.choice(random.choice(MAP_solutions)).copy()

        x2 = random.choice(random.choice(MAP_solutions)).copy()

    return x1, x2			


def random_variation(x1, x2, n):

    p =  random.uniform(0, 1)

    if p > 0.95:

        # Crossover of the parents x1 and x2

        crossover_point = random.randrange(start = 0, stop = n , step = 1) 

        x1[crossover_point+1:] = x2[crossover_point+1:]

    else:

        # Random mutation of X

        random_mutation_site = random.randrange(start = 0, stop = n, step = 1) 

        if (x1[random_mutation_site] == 1):

            x1[random_mutation_site] = 0

        else:

            x1[random_mutation_site] = 1
        
    return x1

def initial_mutant(n, m):

    # m: the amount of knockouts you want

    # n: the amount of genes in the model

    # Generate a vector of length n with 0's and m 1's in a random plance

    x = np.zeros(n) 

    for i in range(0, m, 1):

        random_mutation_site = random.randrange(start = 0, stop = n, step = 1) 

        while (x[random_mutation_site] == 1):
            random_mutation_site = random.randrange(start = 0, stop = n, step = 1) 
        
        x[random_mutation_site] = 1            
    
    return x 


def MAP_Elites(I, G, model, objective):

    # I: the amount of iterations that will be done before the main loop is terminated

    # G: the amount of initial random solutions that will be created

    # model: the genome-scale metabolic model

    # objective: the flux to be maximized

    #----

    n = len(model.genes)

    cols = 5

    rows = 11

    MAP_solutions = [[[0] * n for i in range(cols)] for j in range(rows)] #[[0] * n for i in range(amount_of_niches)] # = X  
    MAP_performances = np.zeros((11,5)) #[0 for i in range(amount_of_niches)] # = P

    # 1. Initialize the algorithm by randomly generating G genomes and determining the performance and features of each solution + place them in the archive

    for i in range(1, G, 1): # -> OR UNTIL ALL X Y Z NICHES ARE FILLED?

        # 1.1 Generate a (random) population of vectors of 0's and 1's of length n

        # ANOTHER ALTERNATIVE IS TO START WITH ALL THE ZERO, ONE AND TWO KNOCKOUTS AS INITIAL POPULATION ASSIGN THEM TO THE ARCHIVE AND START THE NEW GENERATIONS FROM THERE

        m = np.random.choice([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) # up to 10 mutations

        #m = np.random.choice([0, 1]) # Only wildtypes (no mutations) and organisms with 1 mutation       
        
        x = initial_mutant(n, m)

        #print("sum(x)", sum(x))

        # 1.2 Evaluate the performance/fitness p of x
		
        p = performance(x, model, objective)

        print("p", p)
        
        # 1.3 Find out which niche the solution x belongs to

        niche_x_1, niche_x_2  = niche(x)

        #print("niches", niche_x_1, niche_x_2 )
        
        # 1.4 Check if the niche of the random solution is still empty or the solution is better than the current elite in its specific niche. If either one of those is true, the solution will occupy the niche.
    
        if MAP_performances[niche_x_1, niche_x_2] < p:

			# Store the performance of x in the map of elites according to its feature descriptor/niche
			
            MAP_performances[niche_x_1, niche_x_2] = p

			# Store the solution x in the map of elites according to its feature descriptor/niche

            MAP_solutions[niche_x_1][niche_x_2] = x

    #  2. Repeat the following evolutionary loop for I iterations 
    
    for i in range(1, I, 1):

        # All subsequent solutions are generated from elites in the map

		# 2.1 Randomly select 2 elites x1 and x2 from the archive MAP_solutions
		
        x1, x2 = random_selection(MAP_solutions, n)

		# 2.2 Create x_new, a randomly modified copy of x1 and x2 (via mutation and/or crossover)
			
        x_new = random_variation(x1, x2, n)

		# 2.3 Evaluate the new solution x_new and define its niche
		
		# 2.3.1 Evaluate the performance/fitness p_new of x_new
		
        p_new = performance(x_new, model, objective)

        # 2.3.2 Find out which niche the new solution x_new belongs to
 
        niche_new_1, niche_new_2  = niche(x_new)

		# 2.4 Check if the niche of the new solution is still empty or the solution is better than the current elite in its specific niche. If either one of those is true, the solution will occupy the niche.

        if MAP_performances[niche_new_1, niche_new_2] < p_new:

			# Store the performance of x in the map of elites according to its feature descriptor/niche
			
            MAP_performances[niche_new_1, niche_new_2] = p_new

			# Store the solution x in the map of elites according to its feature descriptor/niche

            MAP_solutions[niche_new_1][niche_new_2] = x_new

    return MAP_solutions, MAP_performances # Feature and performance map (P and X)
    
#---

# FINISH THE NICHES SPECIFICALY FOR THE FLAVANONE EXAMPLE -> BIOINFORMATICS_NICHE.PY = same as for Synthetic metabolism

# er is een algemeen probleem met het universeel bigg model dat nog niet opgelost is

#---

## Flavanone production

# MAP-Elites combined with the Escherichia coli str. K-12 substr. MG1655 GEM to maximize flavanone production.

# 0. Load the E.coli GEM

iML1515 = cobra.io.read_sbml_model("Combining ME with GEMs/Models/iML1515.xml")

len(iML1515.metabolites) # 1877
len(iML1515.reactions) # 2712
len(iML1515.genes) # 1516

# 1. FBA of the original iML1515 model

model1 = iML1515.copy()

model1.objective = "BIOMASS_Ec_iML1515_core_75p37M"

sol1 = model1.optimize()

print(sol1.fluxes["BIOMASS_Ec_iML1515_core_75p37M"])

print(model1.summary())

# 2. Add the metabolites, reactions and genes for the heterologous flavanone pathway. This pathway consists 
# of the reactions 4CL, CHS and CHI. 

e_coli_flavanone = iML1515

# Add all the metabolites

cma_c = cobra.Metabolite(
    'cma_c',
    formula='C9H8O3',
    name='Coumaric Acid',
    compartment='c')

cma_e = cobra.Metabolite(
    'cma_e',
    formula='C9H8O3',
    name='Coumaric Acid (Extracellular)',
    compartment='e')

cmcoa_c = cobra.Metabolite(
    'cmcoa_c',
    formula='C30H42N7O18P3S',
    name='Coumaroyl-CoA',
    compartment='c')

chal_c = cobra.Metabolite(
    'chal_c',
    formula='C15H12O5',
    name='Naringenin Chalcone',
    compartment='c')

narg_c = cobra.Metabolite(
    'narg_c',
    formula='C15H12O5',
    name='Naringenin',
    compartment='c')

narg_e = cobra.Metabolite(
    'narg_e',
    formula='C15H12O5',
    name='Naringenin (Extracellular)',
    compartment='e')

atp_c = iML1515.metabolites.get_by_id("atp_c")
coa_c = iML1515.metabolites.get_by_id("coa_c")
amp_c = iML1515.metabolites.get_by_id("amp_c")
ppi_c = iML1515.metabolites.get_by_id("ppi_c")
malcoa_c = iML1515.metabolites.get_by_id("malcoa_c")
co2_c = iML1515.metabolites.get_by_id("co2_c")

# Add 4CL

reaction1 = cobra.Reaction('4CL')
reaction1.name = '4-coumaroyl:coenzyme A (CoA) ligase'
reaction1.subsystem = 'Flavonoid biosynthesis (heterologous)'
reaction1.lower_bound = 0.  # This is the default
reaction1.upper_bound = 1000.  # This is the default

reaction1.add_metabolites({
    atp_c: -1.0,
    cma_c: -1.0,
    coa_c: -1.0,
    amp_c: 1.0,
    ppi_c: 1.0,
    cmcoa_c: 1.0
})

reaction1.reaction

reaction1.gene_reaction_rule = '4CL'
reaction1.genes

e_coli_flavanone.add_reactions([reaction1])

# Add CMAt

reaction2 = cobra.Reaction('CMAt')
reaction2.name = 'transfer of coumaric acid from the extracellular space to the cytosol'
reaction2.subsystem = 'Flavonoid biosynthesis (heterologous)'
reaction2.lower_bound = -1000.  
reaction2.upper_bound = 1000.  # This is the default

reaction2.add_metabolites({
    cma_e: -1.0,
    cma_c: 1.0
})

reaction2.reaction

reaction2.gene_reaction_rule = 'CMAt'
reaction2.genes

e_coli_flavanone.add_reactions([reaction2])

# Add CHS

reaction3 = cobra.Reaction('CHS')
reaction3.name = 'flavonoid chalcone synthase'
reaction3.subsystem = 'Flavonoid biosynthesis (heterologous)'
reaction3.lower_bound = 0.  # This is the default
reaction3.upper_bound = 1000.  # This is the default

reaction3.add_metabolites({
    malcoa_c: -3.0,
    cmcoa_c: -1.0,
    coa_c: 4.0,
    chal_c: 1.0, 
    co2_c: 3.0

})

reaction3.reaction

reaction3.gene_reaction_rule = 'CHS'
reaction3.genes

e_coli_flavanone.add_reactions([reaction3])

# Add CHI

reaction4 = cobra.Reaction('CHI')
reaction4.name = 'flavonoid chalcone isomerase'
reaction4.subsystem = 'Flavonoid biosynthesis (heterologous)'
reaction4.lower_bound = 0.  # This is the default
reaction4.upper_bound = 1000.  # This is the default

reaction4.add_metabolites({
    chal_c: -1.0, 
    narg_c: 1.0

})

reaction4.reaction

reaction4.gene_reaction_rule = 'CHI'
reaction4.genes

e_coli_flavanone.add_reactions([reaction4])

# Add NARGt

reaction5 = cobra.Reaction('NARGt')
reaction5.name = 'transfer of naringenin from the cytosol to the extracellular space'
reaction5.subsystem = 'Flavonoid biosynthesis (heterologous)'
reaction5.lower_bound = -1000.  
reaction5.upper_bound = 1000.  # This is the default

reaction5.add_metabolites({
    narg_c: -1.0,
    narg_e: 1.0
})

reaction5.reaction

reaction5.gene_reaction_rule = 'NARGt'
reaction5.genes

e_coli_flavanone.add_reactions([reaction5])

# Add cma to the medium

e_coli_flavanone.add_boundary(e_coli_flavanone.metabolites.get_by_id("cma_e"), type="exchange")

medium = e_coli_flavanone.medium

medium["EX_cma_e"] = 1000

e_coli_flavanone.medium = medium

e_coli_flavanone.medium

#---

#e_coli_flavanone.add_boundary(e_coli_flavanone.metabolites.get_by_id("cma_c"), type="demand")
#medium["DM_cma_c"] = 1000

#e_coli_flavanone.add_boundary(e_coli_flavanone.metabolites.get_by_id("coa_c"), type="demand")
#medium["DM_coa_c"] = 1000

#e_coli_flavanone.add_boundary(e_coli_flavanone.metabolites.get_by_id("atp_c"), type="demand")
#medium["DM_atp_c"] = 1000

#---

# This result in an extended model, starting with the most recent and most 
# complete Escherichia coli str. K-12 substr. MG1655 in the BIGG database and 
# the additional metabolites, reactions and genes for flavanone production.

len(e_coli_flavanone.metabolites) # 1877 -> 1883
len(e_coli_flavanone.reactions) # 2712 -> 2717
len(e_coli_flavanone.genes) # 1516 -> 1521

# 3. FBA of the extended E. coli model with the flavanone pathway

model2 = e_coli_flavanone.copy()

model2.objective = "4CL"#"BIOMASS_Ec_iML1515_core_75p37M"

sol2 = model2.optimize()

print(sol2.fluxes["BIOMASS_Ec_iML1515_core_75p37M"])

print(model2.summary())

print(sol2.fluxes["4CL"])

print(sol2.fluxes["CHS"])

print(sol2.fluxes["CHI"])

print(sol2.fluxes["CMAt"])

print(sol2.fluxes["NARGt"])

print(sol2.fluxes["ACCOAC"]) # Acetyl-CoA carboxylase: Acetyl-coa -> Malonyl-coa

print(sol2.fluxes["EX_cma_e"])

# CURRENT PROBLEM -> I DON'T GET ANY FLUX THROUGH THE NEW REACTIONS AND I DON'T KNOW WHY

#---

import cobra.test
from cobra.flux_analysis import gapfill

solution = gapfill(model2, iML1515, demand_reactions=False)

for reaction in solution[0]:
    print(reaction.id)

#---

# 4. Run MAP-Elites 

I1 = 0#5000 # Number of generations

G1 = 100#500 # Number of initial individuals generated

MAP_Elites_model = e_coli_flavanone.copy()

objective_ME = "4CL" #"NARGt" # Which reaction/metabolite should be maximized

MAP_solutions1, MAP_performances1 = MAP_Elites(I1, G1, MAP_Elites_model, objective_ME)

# 3. Analyze the results of MAP-Elites 

# The performances/fitness

np.round(MAP_performances1, 2)

# The solutions

print(MAP_solutions1)

for j in range(0, 11):
    print(j)
    for i in range(0, 5):
        print(sum(MAP_solutions1[j][i]))

# I THINK THE PROBLEM IS IN THAT ME DOESN'T FIND ANY KNOCKOUTS THAT INCREAS THE PRODUCTION OF FLAVANONES

# Heat map visualisation of the archive

ax = sns.heatmap(MAP_performances1, linewidth=0.5)

plt.show()

# Plot the performance vs the amount of knockouts/niches

#plt.bar(n_knockouts, MAP_performances)

# ...

plt.show()

# Vendiagram of (most) common genes in wellperforming organisms

# ...

# Find the elite with the highest performance and lowest amount of mutations

# ...

print(np.argmax(MAP_performances1))

# How many reactions were knocked out in this elite?

sum(MAP_solutions1[0][0])

# Which genes were knocked out?

indices = np.where(MAP_solutions1[2][3] == 1)[0] 

n = len(model1.genes)

for i in indices:
    print(model1.genes[i])

# Figure of the metabolism with names of knocked out reactions/genes in bold ~ paper flavanoids
   
# Flux profile of solution

for i in indices:
    model1.genes[i].knock_out()
   
sol_test =  model1.optimize()

print(sol_test)

print(sol_test.fluxes["BIOMASS_Ecoli_core_w_GAM"])

print(sol_test.fluxes["EX_succ_e"]) 

print(sol_test.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol_test.fluxes["EX_succ_e"])

print(model1.summary())

# Phenotypic phase plane (flux) ~ cameo

# Product envelope ~ cameo

# Table with elites, knocked out genes and reactions, target flux, biomass flux, yield, fitness ... ~ cameo

# Above was "Predict gene knockout strategies". Additionally, "Predict expression modulation targets" is of interest as well. 

# https://cameo.bio/06-predict-gene-modulation-targets.html

# + "Predict heterologous pathways"?

# + make a visual progression bar to see how far the algorithm has come in %.
















