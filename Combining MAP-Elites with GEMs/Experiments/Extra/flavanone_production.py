# Combining evolutionary algorithms with genome scale modeling for intelligent metabolic engineering - Flavanone production
# - Author: Shauny Van Hoye
# - Date: 2022-05-07

import random
import numpy as np
import cobra
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3
from cobra.flux_analysis import production_envelope

random.seed(777)

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
        
    sugar = np.load("Combining ME with GEMs/Supporting scripts and data/data/sugar_e_coli_flavanone.npy.npy")

    carbohydrate = np.load("Combining ME with GEMs/Supporting scripts and data/data/carbohydrate_e_coli_flavanone.npy")

    nucleotide = np.load("Combining ME with GEMs/Supporting scripts and data/data/nucleotide_e_coli_flavanone.npy")

    lipid = np.load("Combining ME with GEMs/Supporting scripts and data/data/lipid_e_coli_flavanone.npy")

    AA = np.load("Combining ME with GEMs/Supporting scripts and data/data/AA_e_coli_flavanone.npy")

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


def initial_mutant_v2(n, m):

    # m: the index of the gene you want to  knock out

    # n: the amount of genes in the model

    # Generate a vector of length n with 0's and m 1's in a random plance

    x = np.zeros(n) 

    x[m] = 1            
    
    return x 


def MAP_Elites(I, model, objective):

    # I: the amount of iterations that will be done before the main loop is terminated

    # model: the genome-scale metabolic model

    # objective: the flux to be maximized

    #----

    n = len(model.genes)

    cols = 5

    rows = 11

    MAP_solutions = [[[0] * n for i in range(cols)] for j in range(rows)] #[[0] * n for i in range(amount_of_niches)] # = X  
    MAP_performances = np.zeros((11,5)) #[0 for i in range(amount_of_niches)] # = P

    best_elite_per_generation = np.zeros(I)

   # 1. Initialize the algorithm by randomly generating G genomes and determining the performance and features of each solution + place them in the archive

    for i in range(0, n, 1): 

        # 1.1 Generate an initial (random) population of vectors of 0's and 1's of length n
       
        x = initial_mutant_v2(n, i)

        # 1.2 Evaluate the performance/fitness p of x
		
        p = performance(x, model, objective)
        
        # 1.3 Find out which niche the solution x belongs to

        niche_x_1, niche_x_2  = niche(x)
        
        # 1.4 Check if the niche of the random solution is still empty or the solution is better than the current elite in its specific niche. If either one of those is true, the solution will occupy the niche.
    
        if MAP_performances[niche_x_1, niche_x_2] < p:

			# Store the performance of x in the map of elites according to its feature descriptor/niche
			
            MAP_performances[niche_x_1, niche_x_2] = p

			# Store the solution x in the map of elites according to its feature descriptor/niche

            MAP_solutions[niche_x_1][niche_x_2] = x

    #  2. Repeat the following evolutionary loop for I iterations 
    
    for i in range(0, I, 1):

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

    return MAP_solutions, MAP_performances, best_elite_per_generation # Feature and performance map (P and X)
    
## Other functions

def ven_genes(performance_map, solution_map, model):

    threeBest = np.zeros(3) 

    threeBest_map = [[0] * 2 for i in range(0,3)]

    for i in range(0, 11):
        for j in range(0, 5):
            if np.round(performance_map[i, j], 2) > threeBest[0]:
                threeBest[0] = np.round(performance_map[i, j], 2)
                threeBest_map[0] = [i, j]
            elif np.round(performance_map[i, j], 2) > threeBest[1]:
                if np.round(performance_map[i, j], 2) !=  threeBest[0]:
                    threeBest[1] = np.round(performance_map[i, j], 2)
                    threeBest_map[1] = [i, j]
            elif np.round(performance_map[i, j], 2) > threeBest[2]:
                if np.round(performance_map[i, j], 2) !=  threeBest[0]:
                    if np.round(performance_map[i, j], 2) !=  threeBest[1]:
                        threeBest[2] = np.round(performance_map[i, j], 2)
                        threeBest_map[2] = [i, j]

    m1 = []

    m2 = []

    m3 = []

    for k in range(0, 3):

        #model = core_e_coli_model.copy()

        #model.objective = "BIOMASS_Ecoli_core_w_GAM"
        
        i = threeBest_map[k][0]
        
        j = threeBest_map[k][1]

        #print(i, j)
        
        x = solution_map[i][j]
        
        indices = np.where(x == 1)[0] 

        #n = len(model.genes)

        for i in indices:
            #print(model.genes[i])
            if k == 0:
                m1.append(str(model.genes[i]))
            elif k == 1:
                m2.append(str(model.genes[i]))
            else:
                m3.append(str(model.genes[i]))

        # Flux profile of solution

        #for i in indices:
            #model.genes[i].knock_out()
        
        #sol_test =  model.optimize()

        #print(sol_test.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol_test.fluxes["SUCCt3"])

    M1 = set(m1)
    M2 = set(m2)
    M3 = set(m3)

    v = venn3([M1,M2,M3], ('Mutant 1', 'Mutant 2', 'Mutant 3'))

    v.get_label_by_id('111').set_text('\n'.join(M1&M2&M3))
    v.get_label_by_id('001').set_text('\n'.join(M3-M1-M2))
    v.get_label_by_id('010').set_text('\n'.join(M2-M1-M3))
    v.get_label_by_id('011').set_text('\n'.join(M2&M3-M1))
    v.get_label_by_id('101').set_text('\n'.join(M1&M3-M2))
    v.get_label_by_id('100').set_text('\n'.join(M1-M2-M3))
    v.get_label_by_id('110').set_text('\n'.join(M1&M2-M3))

    return v


def best_elite(performance_map):

    best_elite = 0

    for i in range(0, 11):
    
        for j in range(0, 5):
    
            if np.round(performance_map[i, j], 2) > best_elite:
    
                best_elite = np.round(performance_map[i, j], 2)
    
                best_elite_position = [i, j]

    return best_elite, best_elite_position


def mapelites_heatmap(performance_map,  values):

    x_axis_labels = ["sugar", "carbohydrate", "nucleotide", "lipid", "amino acid"]

    s = sns.heatmap(performance_map, linewidth=0.5, cmap = "viridis", cbar_kws={'label': 'Performance'}, xticklabels = x_axis_labels,  annot = values)

    s.set(xlabel="Metabolism", ylabel="Number of knockouts", title = "MAP-Elites Archive")

    return s


def gene_knockouts(model, solution):

    indices = np.where(solution == 1)[0] 

    n = len(model.genes)

    knockouts = []

    for i in indices:
        #print(model.genes[i])
        knockouts.append(str(model.genes[i]))

    return knockouts


def flux_profile(model, solution):

    indices = np.where(solution == 1)[0] 

    n = len(model.genes)

    for i in indices:
        model.genes[i].knock_out()
   
    sol =  model.optimize()

    print(model.summary())

    return sol  


def PhPP(original_model, solution, objective_ME):

    core_model = e_coli_flavanone.copy()

    prod_env = production_envelope(core_model, objective_ME)

    plt.plot(prod_env.flux_maximum,prod_env.NARGt, color = 'blue')

    plt.fill_between(prod_env.flux_maximum, prod_env.NARGt, np.min(prod_env.NARGt), facecolor='blue', alpha=0.1)

    model = original_model.copy()

    indices = np.where(solution == 1)[0] 

    n = len(model.genes)

    for i in indices:
        model.genes[i].knock_out()

    prod_env = production_envelope(model, objective_ME)

    plt.plot(prod_env.flux_maximum,prod_env.NARGt, color = 'green')

    plt.fill_between(prod_env.flux_maximum, prod_env.NARGt, np.min(prod_env.NARGt), facecolor='green', alpha=0.1)

    plt.ylabel('NARGt-CHANGE to ... [mmol gDW^-1 h^-1]')
    plt.xlabel('BIOMASS-CHANGE to biomass [h^-1]')
    plt.title('Phenotypic phase plane (flux)')

    plt.legend(["Wild type", "Mutant"])

    return plt 



def genes_best_three(pef, sol, e_coli_model):

    threeBest = np.zeros(3) 

    threeBest_map = [[0] * 2 for i in range(0,3)]

    for i in range(0, 11):
        for j in range(0, 5):
            if np.round(pef[i, j], 2) > threeBest[0]:
                threeBest[0] = np.round(pef[i, j], 2)
                threeBest_map[0] = [i, j]
            elif np.round(pef[i, j], 2) > threeBest[1]:
                if np.round(pef[i, j], 2) !=  threeBest[0]:
                    threeBest[1] = np.round(pef[i, j], 2)
                    threeBest_map[1] = [i, j]
            elif np.round(pef[i, j], 2) > threeBest[2]:
                if np.round(pef[i, j], 2) !=  threeBest[0]:
                    if np.round(pef[i, j], 2) !=  threeBest[1]:
                        threeBest[2] = np.round(pef[i, j], 2)
                        threeBest_map[2] = [i, j]

    m1 = []

    m2 = []

    m3 = []

    for k in range(0, 3):

        #model = core_e_coli_model.copy()

        #model.objective = "BIOMASS_Ecoli_core_w_GAM"

        i = threeBest_map[k][0]

        j = threeBest_map[k][1]

        #print(i, j)

        x = sol[i][j]

        indices = np.where(x == 1)[0] 

        #n = len(model.genes)

        for i in indices:
            #print(model.genes[i])
            if k == 0:
                m1.append(str(e_coli_model.genes[i]))
            elif k == 1:
                m2.append(str(e_coli_model.genes[i]))
            else:
                m3.append(str(e_coli_model.genes[i]))
                
    return m1, m2, m3

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
# of the reactions CCL, CHS and CHI. 

e_coli_flavanone = iML1515.copy()

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
pi_c = iML1515.metabolites.get_by_id("pi_c")
malcoa_c = iML1515.metabolites.get_by_id("malcoa_c")
co2_c = iML1515.metabolites.get_by_id("co2_c")

# Add CCL

reaction1 = cobra.Reaction('CCL')
reaction1.name = '4-coumaroyl:coenzyme A (CoA) ligase'
reaction1.subsystem = 'Flavonoid biosynthesis (heterologous)'
reaction1.lower_bound = 0.  # This is the default
reaction1.upper_bound = 1000.  # This is the default

reaction1.add_metabolites({
    atp_c: -1.0,
    cma_c: -1.0,
    coa_c: -1.0,
    amp_c: 1.0,
    pi_c: 1.0,
    cmcoa_c: 1.0
})

reaction1.reaction

reaction1.gene_reaction_rule = 'CCL'
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

e_coli_flavanone.add_boundary(e_coli_flavanone.metabolites.get_by_id("narg_e"), type="exchange")

e_coli_flavanone.medium = medium

e_coli_flavanone.medium

# Necessary demand reaction

e_coli_flavanone.add_boundary(e_coli_flavanone.metabolites.get_by_id("chal_c"), type="demand")
medium["DM_chal_c"] = 1000

e_coli_flavanone.add_boundary(e_coli_flavanone.metabolites.get_by_id("narg_e"), type="demand")
medium["DM_narg_e"] = 1000

#---

# This result in an extended model, starting with the most recent and most 
# complete Escherichia coli str. K-12 substr. MG1655 in the BIGG database and 
# the additional metabolites, reactions and genes for flavanone production.

len(e_coli_flavanone.metabolites) # 1877 -> 1883
len(e_coli_flavanone.reactions) # 2712 -> 2721
len(e_coli_flavanone.genes) # 1516 -> 1521

# 3. FBA of the extended E. coli model with the flavanone pathway

model2 = e_coli_flavanone.copy()

model2.reactions.PDH.knock_out()

#model2.genes.b1091.knock_out()

# model2.reactions.PPC.knock_out()

# model2.reactions.DHORD5.knock_out()

# model2.reactions.DHORD2.knock_out()

# model2.reactions.SUCDi.knock_out()

# model2.genes.b3956.knock_out()

# model2.genes.b0945.knock_out()

# model2.genes.b0721.knock_out()
# model2.genes.b0722.knock_out()
# model2.genes.b0723.knock_out()
# model2.genes.b0724.knock_out()

model2.objective = "BIOMASS_Ec_iML1515_core_75p37M"

sol2 = model2.optimize()

print(sol2.fluxes["BIOMASS_Ec_iML1515_core_75p37M"])

print(model2.summary())

print(sol2.fluxes["CCL"])

print(sol2.fluxes["CHS"])

print(sol2.fluxes["CHI"])

print(sol2.fluxes["CMAt"])

print(sol2.fluxes["NARGt"])

print(sol2.fluxes["ACCOAC"]) # Acetyl-CoA carboxylase: Acetyl-coa -> Malonyl-coa (malcoa_c)

print(sol2.fluxes["EX_narg_e"])

print(sol2.fluxes["EX_glc__D_e"])

# Save the model

#cobra.io.write_sbml_model(e_coli_flavanone, "e_coli_flavanone.xml")

#---

# 4. Run MAP-Elites 

I1 = 0#5000 # Number of generations

MAP_Elites_model = e_coli_flavanone.copy()

objective_ME = "NARGt" 

MAP_solutions1, MAP_performances1, best_elite_per_generation = MAP_Elites(I1, MAP_Elites_model, objective_ME)

#np.save("MAP_solutions_flavones_25000", MAP_solutions1)

#np.save("MAP_performances_flavones_25000", MAP_performances1)

# 3. Analyze the results of MAP-Elites 

# The performances/fitness archive

np.round(MAP_performances1, 2)

# Heat map visualisation of the archive

mapelites_heatmap(MAP_performances1,  values = False)

plt.show()

# Vendiagram of (most) common genes in wellperforming mutants

ven_genes(MAP_performances1, MAP_solutions1, model1)

plt.show()

# Plot the maximum performance of the map vs the number of generations

fig3, ax3 = plt.subplots()

ax3.plot(range(0, len(best_elite_per_generation)), best_elite_per_generation)

ax3.set_ylabel('Performance')
ax3.set_xlabel('Number of generations')
ax3.set_title('Performance vs the number of generations')

plt.show()

# Plot the performance vs the number of knockouts

fig, ax = plt.subplots()

ax.plot(range(0, 11), [max(MAP_performances1[i, :]) for i in range(0, 11)])

ax.set_ylabel('Performance')
ax.set_xlabel('Number of knockouts')
ax.set_title('Performance vs the number of knockouts')

plt.show()

# Plot the performance vs the type of metabolism

fig2, ax2 = plt.subplots()

ax2.bar(range(0, 5), [np.average(MAP_performances1[:, i]) for i in range(0, 5)])

ax2.set_ylabel('Maximum performance')
ax2.set_xlabel('Type of metabolism')
ax2.set_title('Maximum performance vs the type of metabolism')
ind = np.arange(5) 
ax2.set_xticks(ind, labels=["sugar", "carbohydrate", "nucleotide", "lipid", "amino acid"])

plt.show()

# Find the elite with the highest performance and lowest amount of mutations

best_elite_performance, best_elite_position = best_elite(MAP_performances1)

print("The fitness value of the best performing mutant is ", best_elite_performance)

# How many genes were knocked out in this elite?

best_elite = MAP_solutions1[best_elite_position[0]][best_elite_position[1]]

# Which genes were knocked out?

knockouts = gene_knockouts(model1, best_elite)

# Flux profile of best_elite

sol_best_elite = flux_profile(model1, best_elite)

print(sol_best_elite.fluxes["ACCOAC"]) 

print(sol_best_elite.fluxes["BIOMASS_Ec_iML1515_core_75p37M"])

print(sol_best_elite.fluxes["BIOMASS_Ec_iML1515_core_75p37M"]*sol_best_elite.fluxes["ACCOAC"]) # = the performance

print(sol_best_elite.fluxes) # All the fluxes

# Phenotypic phase plane (flux) = Product envelope

PhPP(model1, best_elite, objective_ME)

plt.show()


best_elite = np.zeros(len(e_coli_flavanone.genes))

PhPP(model1, best_elite, objective_ME)

plt.show()

# Table with elites, knocked out genes and reactions, target flux, biomass flux, yield, fitness ... for the 3 best elites ~ cameo

# ...

# +

# Compare ME to optgene

from cameo.strain_design import OptGene
from cameo.visualization.plotting.with_plotly import PlotlyPlotter

plotter = PlotlyPlotter()

optgene = OptGene(e_coli_flavanone.copy())

result = optgene.run(target=model1.reactions.EX_ac_e,
                     biomass=model1.reactions.BIOMASS_Ecoli_core_w_GAM,
                     substrate=model1.metabolites.glc__D_e,
                     max_evaluations=5000,
                     plot=True)

result.data_frame

result.plot(plotter, 0)


# Figure of the metabolism with names of knocked out reactions/genes in bold ~ paper flavanoids -> make afterwards in Biorender
   
# ...



