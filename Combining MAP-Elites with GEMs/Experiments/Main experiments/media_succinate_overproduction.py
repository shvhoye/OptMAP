# Combining evolutionary algorithms with genome scale modeling for intelligent metabolic engineering - Succinate overproduction
# - Author: Shauny Van Hoye
# - Date: 2022-05-30

import random
import numpy as np
import pandas as pd
import cobra
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3
from cobra.flux_analysis import production_envelope
from cobra.medium import minimal_medium

#random.seed(777)

matplotlib.rcParams.update({'font.size': 14}) # python plot increase font size

## MAP-Elites functions 

def performance(x, model, objective, m, m_original): 
	
    # Evaluate the performance/fitness of a solution

    indices = np.where(x == 1)[0] 

    temp_model1 = model.copy()

    if sum(x) > 0:

        for i in indices:
            temp_model1.genes[i].knock_out()

    Performance = np.zeros(m.shape[1])

    # Media as niches

    for i in range(0, m.shape[1]):

        temp_model = temp_model1.copy()

        temp_model.medium = m[i]

        temp_model.objective = "BIOMASS_Ecoli_core_w_GAM"
    
        sol = temp_model.optimize()

        substrate = sol.fluxes["EX_fru_e"] + sol.fluxes["EX_glc__D_e"] + sol.fluxes["EX_gln__L_e"] + sol.fluxes["EX_glu__L_e"] + sol.fluxes["EX_mal__L_e"] 

        if substrate == 0:

            Performance[i] = 0
        
        else:

            Performance[i] =  sol.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol.fluxes[objective]/-substrate

        #print(sol.fluxes["BIOMASS_Ecoli_core_w_GAM"], sol.fluxes[objective])

    # temp_model = temp_model1.copy()

    # temp_model.medium = m_original

    # sol = temp_model.optimize()

    # substrate = sol.fluxes["EX_fru_e"] + sol.fluxes["EX_glc__D_e"] + sol.fluxes["EX_gln__L_e"] + sol.fluxes["EX_glu__L_e"] + sol.fluxes["EX_mal__L_e"] 

    # Performance[m.shape[1]] = sol.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol.fluxes[objective]/-substrate

    #print(sol.fluxes["BIOMASS_Ecoli_core_w_GAM"], sol.fluxes[objective])

    #print(Performance)

    return Performance


def performance2(x, model, objective): 
	
    # Evaluate the performance/fitness of a solution

    indices = np.where(x == 1)[0] 

    temp_model = model.copy()

    if sum(x) > 0:

        for i in indices:
            temp_model.genes[i].knock_out()

    temp_model.objective = "BIOMASS_Ecoli_core_w_GAM"
    
    sol = temp_model.optimize()

    objective_flux = sol.fluxes[objective] 

    biomass = sol.fluxes["BIOMASS_Ecoli_core_w_GAM"]

    Performance =  sol.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol.fluxes[objective]/-sol.fluxes["EX_glc__D_e"]
            
    return Performance, objective_flux, biomass


def performance3(x, model, objective, m, m_original): 
	
    # Evaluate the performance/fitness of a solution

    indices = np.where(x == 1)[0] 

    temp_model1 = model.copy()

    if sum(x) > 0:

        for i in indices:
            temp_model1.genes[i].knock_out()

    Performance = np.zeros(m.shape[1])

    objective_flux = np.zeros(m.shape[1])

    biomass = np.zeros(m.shape[1])

    # Media as niches

    for i in range(0, m.shape[1]):

        temp_model = temp_model1.copy()

        temp_model.medium = m[i]

        temp_model.objective = "BIOMASS_Ecoli_core_w_GAM"
    
        sol = temp_model.optimize()

        substrate = sol.fluxes["EX_fru_e"] + sol.fluxes["EX_glc__D_e"] + sol.fluxes["EX_gln__L_e"] + sol.fluxes["EX_glu__L_e"] + sol.fluxes["EX_mal__L_e"] 

        if substrate == 0:

            Performance[i] = 0
        
        else:

            Performance[i] =  sol.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol.fluxes[objective]/-substrate

            objective_flux[i] =  sol.fluxes[objective]

            biomass[i] = sol.fluxes["BIOMASS_Ecoli_core_w_GAM"]

        #print(sol.fluxes["BIOMASS_Ecoli_core_w_GAM"], sol.fluxes[objective])

    # temp_model = temp_model1.copy()

    # temp_model.medium = m_original

    # sol = temp_model.optimize()

    # substrate = sol.fluxes["EX_fru_e"] + sol.fluxes["EX_glc__D_e"] + sol.fluxes["EX_gln__L_e"] + sol.fluxes["EX_glu__L_e"] + sol.fluxes["EX_mal__L_e"] 

    # Performance[m.shape[1]] = sol.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol.fluxes[objective]/-substrate

    #print(sol.fluxes["BIOMASS_Ecoli_core_w_GAM"], sol.fluxes[objective])

    #print(Performance)

    return Performance, objective_flux, biomass


def niche(x):

    niche = []

    n_mutations = int(sum(x))

    if n_mutations > 10:
        niche.append(10) # if there are more than 10 mutations, they are stored in the 10th niche
    else:
        niche.append(n_mutations)

    # sugar = np.load("Combining ME with GEMs/Supporting scripts and data/data/sugar.npy")

    # carbohydrate = np.load("Combining ME with GEMs/Supporting scripts and data/data/carbohydrate.npy")

    # nucleotide = np.load("Combining ME with GEMs/Supporting scripts and data/data/nucleotide.npy")

    # lipid = np.load("Combining ME with GEMs/Supporting scripts and data/data/lipid.npy")

    # AA = np.load("Combining ME with GEMs/Supporting scripts and data/data/AA.npy")

    # s = 0

    # c = 0

    # n = 0

    # l = 0

    # A = 0

    # indices = np.where(x == 1)[0] 

    # for i in indices:

    #     s += sugar[i]

    #     c += carbohydrate[i]

    #     n += nucleotide[i]

    #     l += lipid[i]

    #     A += AA[i]

    # niche.append(int(np.argmax([s, c, n, l, A])))

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

        # Random mutation of x

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


def MAP_Elites(I, model, objective, m, m_original):

    # I: the number of iterations that will be done before the main loop is terminated

    # model: the genome-scale metabolic model

    # objective: the flux to be maximized

    #----

    # 0. Initialize the archive/map

    n = len(model.genes)

    cols = m.shape[1]

    rows = 11

    MAP_solutions = [[[0] * n for i in range(cols)] for j in range(rows)] # = X  
    MAP_performances = np.zeros((11,cols)) # = P

    best_elite_per_generation = np.zeros(I)

   # 1. Initialize the algorithm by randomly generating G genomes and determining the performance and features of each solution + place them in the archive

    # 1.1 Generate an initial (random) population of vectors of 0's and 1's of length n
       
    x = np.zeros(n)

    # 1.2 Evaluate the performance/fitness p of x
    
    p = performance(x, model, objective, m, m_original)
    
    # 1.3 Find out which niche the solution x belongs to

    niche_x_1  = niche(x)[0] #niche_x_1, niche_x_2  = niche(x)

    for i in range(0, m.shape[1]):

        p_i = p[i]

        niche_x_2 = i
    
        # 1.4 Check if the niche of the random solution is still empty or the solution is better than the current elite in its specific niche. If either one of those is true, the solution will occupy the niche.
    
        if MAP_performances[niche_x_1, niche_x_2] < p_i:

            # Store the performance of x in the map of elites according to its feature descriptor/niche
            
            MAP_performances[niche_x_1, niche_x_2] = p_i

            # Store the solution x in the map of elites according to its feature descriptor/niche

            MAP_solutions[niche_x_1][niche_x_2] = x

    for i in range(0, n, 1): 

        # 1.1 Generate an initial (random) population of vectors of 0's and 1's of length n
       
        x = initial_mutant_v2(n, i)

        # 1.2 Evaluate the performance/fitness p of x
		
        p = performance(x, model, objective, m, m_original)
        
        # 1.3 Find out which niche the solution x belongs to

        niche_x_1  = niche(x)[0] #niche_x_1, niche_x_2  = niche(x)

        for i in range(0, m.shape[1]):

            p_i = p[i]

            niche_x_2 = i
        
            # 1.4 Check if the niche of the random solution is still empty or the solution is better than the current elite in its specific niche. If either one of those is true, the solution will occupy the niche.
        
            if MAP_performances[niche_x_1, niche_x_2] < p_i:

                # Store the performance of x in the map of elites according to its feature descriptor/niche
                
                MAP_performances[niche_x_1, niche_x_2] = p_i

                # Store the solution x in the map of elites according to its feature descriptor/niche

                MAP_solutions[niche_x_1][niche_x_2] = x

    #  2. Repeat the following evolutionary loop for I iterations 
    
    for j in range(0, I, 1):

        # All subsequent solutions are generated from elites in the map

		# 2.1 Randomly select 2 elites x1 and x2 from the archive MAP_solutions
		
        x1, x2 = random_selection(MAP_solutions, n)

		# 2.2 Create x_new, a randomly modified copy of x1 and x2 (via mutation and/or crossover)
			
        x_new = random_variation(x1, x2, n)

		# 2.3 Evaluate the new solution x_new and define its niche
		
		# 2.3.1 Evaluate the performance/fitness p_new of x_new
		
        p_new = performance(x_new, model, objective, m, m_original)

        # 2.3.2 Find out which niche the new solution x_new belongs to
 
        niche_x_1_new  = niche(x_new)[0] #niche_x_1, niche_x_2  = niche(x)

        for i in range(0, m.shape[1]):

            p_i_new = p_new[i]

            niche_x_2_new = i
        
            # 1.4 Check if the niche of the random solution is still empty or the solution is better than the current elite in its specific niche. If either one of those is true, the solution will occupy the niche.
        
            if MAP_performances[niche_x_1_new, niche_x_2_new] < p_i_new:

                # Store the performance of x in the map of elites according to its feature descriptor/niche
                
                MAP_performances[niche_x_1_new, niche_x_2_new] = p_i_new

                # Store the solution x in the map of elites according to its feature descriptor/niche

                MAP_solutions[niche_x_1_new][niche_x_2_new] = x_new

        best_elite_per_generation[j], _ = best_elite(MAP_performances)

    return MAP_solutions, MAP_performances, best_elite_per_generation # Feature and performance map (P and X)
    
## Other functions

def ven_genes(performance_map, solution_map, model):

    threeBest = np.zeros(3) 

    threeBest_map = [[0] * 2 for i in range(0,3)]

    for i in range(0, 11):
        for j in range(0, 6):
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

    if M1&M2&M3:
        v.get_label_by_id('111').set_text('\n'.join(M1&M2&M3))
    if M3-M1-M2:
        v.get_label_by_id('001').set_text('\n'.join(M3-M1-M2))
    if M2-M1-M3: 
        v.get_label_by_id('010').set_text('\n'.join(M2-M1-M3))
    if M2&M3-M1:
        v.get_label_by_id('011').set_text('\n'.join(M2&M3-M1))
    if M1&M3-M2:
        v.get_label_by_id('101').set_text('\n'.join(M1&M3-M2))
    if M1-M2-M3:
        v.get_label_by_id('100').set_text('\n'.join(M1-M2-M3))
    if M1&M2-M3:
        v.get_label_by_id('110').set_text('\n'.join(M1&M2-M3))

    return v



def best_elite(performance_map):

    best_elite = 0

    for i in range(0, 11):
    
        for j in range(0, 6):
    
            if np.round(performance_map[i, j], 2) > best_elite:
    
                best_elite = np.round(performance_map[i, j], 2)
    
                best_elite_position = [i, j]

    return best_elite, best_elite_position


def mapelites_heatmap(performance_map,  values):

    #x_axis_labels = ["sugar", "carbohydrate", "nucleotide", "lipid", "amino acid"]

    s = sns.heatmap(performance_map, linewidth=0.5, cmap = "viridis", cbar_kws={'label': 'Performance'},  annot = values)#, xticklabels = x_axis_labels,  annot = values)

    s.set(xlabel="Medium", ylabel="Number of knockouts", title = "MAP-Elites Archive")

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


def PhPP(original_model, solution, objective_ME, m):

    core_model = core_e_coli_model.copy()

    prod_env = production_envelope(core_model, objective_ME, carbon_sources = m)

    plt.plot(prod_env.flux_maximum,prod_env.SUCCt3, color = 'blue')

    plt.fill_between(prod_env.flux_maximum, prod_env.SUCCt3, np.min(prod_env.SUCCt3), facecolor='blue', alpha=0.1)

    model = original_model.copy()

    indices = np.where(solution == 1)[0] 

    n = len(model.genes)

    for i in indices:
        model.genes[i].knock_out()

    prod_env = production_envelope(model, objective_ME, carbon_sources = m)

    plt.plot(prod_env.flux_maximum,prod_env.SUCCt3, color = 'green')

    plt.fill_between(prod_env.flux_maximum, prod_env.SUCCt3, np.min(prod_env.SUCCt3), facecolor='green', alpha=0.1)

    plt.ylabel('SUCCt3 [mmol gDW^-1 h^-1]')
    plt.xlabel('BIOMASS_Ecoli_core_w_GAM [h^-1]')
    plt.title('Phenotypic phase plane (flux)')

    plt.legend(["Wild type", "Mutant"])

    return plt   


def genes_best_three(pef, sol, e_coli_model):

    threeBest = np.zeros(3) 

    threeBest_map = [[0] * 2 for i in range(0,3)]

    for i in range(0, 11):
        for j in range(0, 6):
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
                #m1.append(i)
            elif k == 1:
                m2.append(str(e_coli_model.genes[i]))
                #m2.append(i)
            else:
                m3.append(str(e_coli_model.genes[i]))
                #m3.append(i)
                
    return m1, m2, m3




def best_three(performance_map, solution_map, model, objective, m, m_original):

    threeBest = np.zeros(3) 

    threeBest_map = [[0] * 2 for i in range(0,3)]

    for i in range(0, 11):
        for j in range(0, 6):
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


    for k in range(0, 3, 1):

        i, j = threeBest_map[k]

        x = solution_map[i][j]

        objective

        fitness, objective_flux, biomass = performance3(x, model, objective, m, m_original)

        print(np.round(fitness, 2), np.round(objective_flux, 2), np.round(biomass, 2))

    return 0
    

#---

## Succinate overproduction

# MAP-Elites combined with the core E. coli GEM to maximize succinate production.

# 0. Load the E. coli GEM

core_e_coli_model = cobra.io.read_sbml_model("Combining ME with GEMs/Models/e_coli_core.xml")

# 1. FBA of the wildtype model

model1 = core_e_coli_model.copy()

m_original = model1.medium.copy()

m1 = minimal_medium(model1, 0.8, minimize_components=8, open_exchanges=True)


for i in range(0, m1.shape[1]):

    model1.medium = m1[i]

    # for i in dir(model1): # Gets the attributes
    #     print(i)

    model1.objective = "BIOMASS_Ecoli_core_w_GAM"

    sol1 = model1.optimize()

    print(sol1.fluxes["BIOMASS_Ecoli_core_w_GAM"])

    print(sol1.fluxes["SUCCt3"])

x1 = initial_mutant(len(model1.genes), 0)

np.round(performance(x1, model1, "SUCCt3", m1, m_original), 2)

print(model1.summary())

print(sol1.fluxes["BIOMASS_Ecoli_core_w_GAM"])

print(sol1.fluxes["EX_succ_e"])

print(sol1.fluxes["SUCCt3"])

print(sol1.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol1.fluxes["SUCCt3"])

# 2. Run MAP-Elites 

core_e_coli_model = cobra.io.read_sbml_model("Combining ME with GEMs/Models/e_coli_core.xml")

I1 = 1000 # Number of generations

MAP_Elites_model = core_e_coli_model.copy()

m = minimal_medium(MAP_Elites_model, 0.8, minimize_components=8, open_exchanges=True)

objective_ME = "SUCCt3" # SUCCt3 is the final reaction to get extracellular succinate: succ_c + h_e ??? h_c + succ_e

m_original = core_e_coli_model.medium.copy()

MAP_solutions1, MAP_performances1, best_elite_per_generation = MAP_Elites(I1, MAP_Elites_model, objective_ME, m, m_original)

#m_original = {'EX_fru_e': 1000.0}

#x = np.zeros(len(MAP_Elites_model.genes)) 

# x = initial_mutant(len(MAP_Elites_model.genes), 0)

# performance(x, MAP_Elites_model, objective_ME, m, m_original)

# np.save("Combining ME with GEMs/Results and Figures/MAP_solutions_succ20000_v2", MAP_solutions1)

# np.save("Combining ME with GEMs/Results and Figures/MAP_performances_succ20000_v2", MAP_performances1)

# np.save("Combining ME with GEMs/Results and Figures/best_elite_per_generation_succ20000_v2", best_elite_per_generation)

# succ_pef = np.load("Combining ME with GEMs/Results and Figures/Succ25000/MAP_performances25000.npy")

# succ_sol = np.load("Combining ME with GEMs/Results and Figures/Succ25000/MAP_solutions25000.npy")

# MAP_performances1 = succ_pef

# MAP_solutions1 = succ_sol

# 3. Analyze the results of MAP-Elites 

# The performances/fitness archive

np.round(MAP_performances1, 2)

# Heat map visualisation of the archive

mapelites_heatmap(MAP_performances1,  values = True)

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

ax2.bar(range(0, 7), [np.average(MAP_performances1[:, i]) for i in range(0, 7)])

ax2.set_ylabel('Maximum performance')
ax2.set_xlabel('Type of metabolism')
ax2.set_title('Maximum performance vs the type of metabolism')
ind = np.arange(7) 
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

print(sol_best_elite.fluxes["SUCCt3"]) 

print(sol_best_elite.fluxes["BIOMASS_Ecoli_core_w_GAM"])

print(sol_best_elite.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol_best_elite.fluxes["SUCCt3"]) # = the performance

print(sol_best_elite.fluxes) # All the fluxes

# Phenotypic phase plane (flux) = Product envelope

z = MAP_solutions1[1][0]

performance(z, model1, "SUCCt3", m1, m_original)

carbon = ["EX_glu__L_e"]

PhPP(model1, z, objective_ME, carbon)

plt.show()

#--

z = MAP_solutions1[9][1]

performance(z, model1, "SUCCt3", m1, m_original)

carbon = ["EX_glu__L_e", "EX_glc__D_e"]

PhPP(model1, z, objective_ME, carbon)

plt.show()

#--

z = MAP_solutions1[1][0]

performance(z, model1, "SUCCt3", m1, m_original)

carbon = ["EX_glu__L_e", "EX_glc__D_e"]

PhPP(model1, z, objective_ME, carbon)

plt.show()

# Compare ME to optgene

from cameo.strain_design import OptGene
from cameo.visualization.plotting.with_plotly import PlotlyPlotter


core_e_coli_model = cobra.io.read_sbml_model("Combining ME with GEMs/Models/e_coli_core.xml")

from cameo import phenotypic_phase_plane

plotter = PlotlyPlotter()

core_e_coli_model.reactions.EX_o2_e.lower_bound = -10

result = phenotypic_phase_plane(core_e_coli_model,
                                variables=[core_e_coli_model.reactions.BIOMASS_Ecoli_core_w_GAM],
                                objective=core_e_coli_model.reactions.EX_succ_e,
                                points=10)

result.plot(plotter, 0)




optgene = OptGene(core_e_coli_model.copy())

result = optgene.run(target=core_e_coli_model.reactions.EX_succ_e,
                     biomass=core_e_coli_model.reactions.BIOMASS_Ecoli_core_w_GAM,
                     substrate=core_e_coli_model.metabolites.glc__D_e,
                     max_evaluations=5000,
                     plot=True)

result.data_frame

df = result.data_frame

df.to_csv (r'/Users/Shaun/Documents/School/2021 - 2022/Sem 2/a. Master thesis/Programming/Master Thesis/Combining ME with GEMs/Results and Figures/OptGene_results_succ.csv', index = False, header=True)

pd.set_option("display.max_rows", None, "display.max_columns", None) #display whole df

print(df)

result.plot(plotter, 0)

core_e_coli_model.medium

# Table with elites, knocked out genes and reactions, target flux, biomass flux, yield, fitness ... for the 3 best elites ~ cameo

# print(genes_best_three(MAP_performances1, MAP_solutions1, model1))

genes_best_three(MAP_performances1, MAP_solutions1, model1)

best_three(MAP_performances1, MAP_solutions1, model1, objective_ME, m1, m_original)



data = {'Solution': ['1', '2', '3'], 'reactions': [('ATPS4r', 'PFL', 'PTAr'), ('ATPS4r', ('NADTRHD', 'THD2')), ('ATPS4r', 'TPI')], 'genes': [('b3733', ('b0902', 'b3952'), ('b2297', 'b2458')), ('b3733', 'b1603'), ('b3734', 'b3919')], 'target_flux': ['128.20', '51.22', '109.34'], 'biomass_flux': ['0.37', '0.37', '0.13'], 'fitness': ['4.80', '1.92', '1.47']}

data2 = pd.DataFrame(data)

data2.to_csv (r'/Users/Shaun/Documents/School/2021 - 2022/Sem 2/a. Master thesis/Programming/Master Thesis/Combining ME with GEMs/Results and Figures/OptMAP_results_succ.csv', index = False, header=True)

# Search for the genes in the three best performing mutants using Bioreflect?

# ...

# Figure of the metabolism with names of knocked out reactions/genes in bold ~ paper flavanoids -> make afterwards in Biorender + eschermaps (julia)
   
# ...



indices_best_elite = np.where(best_elite == 1)[0] 

for i in indices_best_elite:
    print(MAP_Elites_model.genes[i])


print(MAP_Elites_model.genes[16])


# adapt code to get more explorations?


# Media as niches

model1.medium

model1.genes[2].knock_out()


minimal_medium(model1, 0.8, minimize_components=8, open_exchanges=True)

m = minimal_medium(model1, 0.8, minimize_components=8, open_exchanges=True)

m[1]

medium =  m[1]#model1.medium

# medium["EX_fru_e"] = 
# medium["EX_glc__D_e"] = 
# medium["EX_gln__L_e"] = 
# medium["EX_glu__L_e"] = 
# medium["EX_mal__L_e"] = 
# medium["EX_nh4_e"] = 
# medium["EX_o2_e"] = 
# medium["EX_pi_e"] = 

model1.medium = medium

model1.medium

cnt = 0
for i in range(0, m.shape[1]):

    cnt+=1

    print(i)#m[i])

print(cnt)


for i in range(0, 7):

    p = performance(MAP_solutions1[9][i], model1, "SUCCt3", m1, m_original)

    print(np.round(p, 2))


x = MAP_solutions1[9][2]

indices = np.where(x == 1)[0] 

#n = len(model.genes)

for i in indices:
    print(str(model1.genes[i]))



PhPP(model1, x, objective_ME)

plt.show()


indices = np.where(x == 1)[0] 

n = len(model1.genes)

for i in indices:
    model1.genes[i].knock_out()

model1.medium = m[2]

sol_best_elite =  model1.optimize()

print(sol_best_elite.fluxes["SUCCt3"]) 

print(sol_best_elite.fluxes["BIOMASS_Ecoli_core_w_GAM"])

print(sol_best_elite.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol_best_elite.fluxes["SUCCt3"]) # = the performance

substrate = sol_best_elite.fluxes["EX_fru_e"] + sol_best_elite.fluxes["EX_glc__D_e"] + sol_best_elite.fluxes["EX_gln__L_e"] + sol_best_elite.fluxes["EX_glu__L_e"] + sol_best_elite.fluxes["EX_mal__L_e"] 

sol_best_elite.fluxes["BIOMASS_Ecoli_core_w_GAM"]*sol_best_elite.fluxes["SUCCt3"]/-substrate

performance(x, model1, "SUCCt3", m1, m_original)

for i in range(0, MAP_performances1.shape[0]):

    for j in range(0, MAP_performances1.shape[1]):

        if MAP_performances1[i,j] == 0:
            
            MAP_performances1[i,j] = 0.00000000001




x = MAP_solutions1[3][5]

fitness, objective_flux, biomass = performance3(x, model1, objective_ME, m, m_original)

print(np.round(fitness, 2), np.round(objective_flux, 2), np.round(biomass, 2))

indices = np.where(x == 1)[0] 

for i in indices:

    print(str(model1.genes[i]))
