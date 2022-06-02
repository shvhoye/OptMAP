# Matching genes that are knocked out to the type of metabolism they are part of to ascribe them to a specific niche
# - Author: Shauny Van Hoye
# - Date: 2022-05-07

import re
import cobra
import numpy as np

def metabolism(gene):

    # Ecocyc

    ecocyc = open("Combining ME with GEMs/Supporting scripts and data/data/ecocyc.gaf")

    go_codes = []

    for line in ecocyc:
            
        l = line.strip()

        if re.search(gene, l):

            s = re.search('GO:([0-9]{7})', l)

            go_codes.append(s.group(1))

            
    # GO

    key = 0

    sugar = 0

    carbohydrate = 0

    nucleotide = 0

    lipid = 0

    AA = 0

    for code in go_codes:

        GO = open("Combining ME with GEMs/Supporting scripts and data/data/go.obo")

        for line in GO:
                
            l = line.strip()

            if re.match('^id:',l):

                if re.search(code, l):
                    
                    key = 1

                else:

                    key = 0

            if key == 1:

                #print(l)

                if re.search('sugar', l):

                    sugar += 1
                    
                if re.search('carbohydrate', l):

                    carbohydrate += 1

                if re.search('nucleotide', l):
                        
                    nucleotide += 1   

                if re.search('lipid', l):
                            
                    lipid += 1 

                if re.search('amino acid', l):
                        
                    AA += 1         


    return sugar, carbohydrate, nucleotide, lipid, AA


#---

# Test 1

import time

t0 = time.time()

genes = ['b1101', 'b1603', 'b3612', 'b3956', 'b3919']

for gene in genes:

    sugar, carbohydrate, nucleotide, lipid, AA = metabolism(gene)

    print(sugar, carbohydrate, nucleotide, lipid, AA)

t1 = time.time()

total = t1-t0

#---

# The next code block goes over all the genes in the core_e_coli_model and makes use of the ecocyc and 
# go database to determine which type of metabolism each gene is correlated with. This code only has to 
# be run once, afterwards, the results are saved in separate files which are in turn used to determine the 
# niches in the Map-Elites algorithm.

# Load E.coli GEM

core_e_coli_model = cobra.io.read_sbml_model("Combining ME with GEMs/Models/e_coli_core.xml")

number_of_genes = len(core_e_coli_model.genes)

sugar = np.zeros(number_of_genes) 

carbohydrate = np.zeros(number_of_genes) 

nucleotide = np.zeros(number_of_genes) 

lipid = np.zeros(number_of_genes) 

AA = np.zeros(number_of_genes) 

for i, gene in enumerate(core_e_coli_model.genes):
    #print(type(str(gene)))
    #print(i)
    #print(AA[i])
    sugar[i], carbohydrate[i], nucleotide[i], lipid[i], AA[i] = metabolism(str(gene))

    if i > 2:
        break

np.save("sugar", sugar)

np.save("carbohydrate", carbohydrate)

np.save("nucleotide", nucleotide)

np.save("lipid", lipid)

np.save("AA", AA)

#---

s = np.load("sugar.npy")

c = np.load("carbohydrate.npy")

n = np.load("nucleotide.npy")

l = np.load("lipid.npy")

A = np.load("AA.npy")

# Calculate averages, variance ... to determine how the niches will be structured

np.average(s)

np.average(c)

np.average(n)

np.average(l)

np.average(A)

np.var(s)

np.var(c)

np.var(n)

np.var(l)

np.var(A)

#---

# The next code block goes over all the genes in the iML1515 model and makes use of the ecocyc and 
# go database to determine which type of metabolism each gene is correlated with. This code only has to 
# be run once, afterwards, the results are saved in separate files which are in turn used to determine the 
# niches in the Map-Elites algorithm.

# Load the iML1515 E. coli GEM

iML1515 = cobra.io.read_sbml_model("Combining ME with GEMs/Models/iML1515.xml")

number_of_genes = len(iML1515.genes)

sugar = np.zeros(number_of_genes) 

carbohydrate = np.zeros(number_of_genes) 

nucleotide = np.zeros(number_of_genes) 

lipid = np.zeros(number_of_genes) 

AA = np.zeros(number_of_genes) 

for i, gene in enumerate(iML1515.genes):
    #print(type(str(gene)))
    #print(i)
    #print(AA[i])
    sugar[i], carbohydrate[i], nucleotide[i], lipid[i], AA[i] = metabolism(str(gene))
    if i > 0:
        break


sugar
carbohydrate
nucleotide
lipid
AA



np.save("sugar", sugar)

np.save("carbohydrate", carbohydrate)

np.save("nucleotide", nucleotide)

np.save("lipid", lipid)

np.save("AA", AA)


np.save("Combining ME with GEMs/Supporting scripts and data/data/sugar_f", sugar)

np.save("Combining ME with GEMs/Supporting scripts and data/data/sugar_carbohydrate", carbohydrate)

np.save("Combining ME with GEMs/Supporting scripts and data/data/sugar_nucleotide", nucleotide)

np.save("Combining ME with GEMs/Supporting scripts and data/data/sugar_lipid", lipid)

np.save("Combining ME with GEMs/Supporting scripts and data/data/sugar_AA", AA)


















#---

# Get fluxes for 500-1000 random knockouts to see the distributions. All slingle knockouts + all double + randoms?





