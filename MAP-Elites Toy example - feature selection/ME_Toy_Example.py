# MAP-Elites Toy Example - feature selection:
# - Author: Shauny Van Hoye
# - Date: 2022-02-17

import random
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn import linear_model
import matplotlib
import matplotlib.pyplot as plt

random.seed(777)

matplotlib.rcParams.update({'font.size': 14})

## Preparation for Naive Bayes

gnb = GaussianNB()

## MAP-Elites functions 

def performance(x, X_train, y_train, X_test, y_test, application):
	
    # Evaluate the performance/fitness of a solution

    # Evaluation of the solutions is done by using Bayes Classifier

    # Naive Bayes

    feature_indices = np.where(x == 1)[0] 
    
    if feature_indices.size == 0:

         Performance = 0

    else:

        if application == 'regression':

            # Create linear regression object
            regr = linear_model.LinearRegression()

            # Train the model using the training sets
            regr.fit(X_train[:,feature_indices], y_train)

            # Make predictions using the testing set
            y_pred = regr.predict(X_test[:,feature_indices])

            # MSE

            Performance = 1/sum(np.square(np.subtract(y_test, y_pred)))
        
        else:

            y_pred = gnb.fit(X_train[:,feature_indices], y_train).predict(X_test[:,feature_indices])
            
            # Accuracy

            Performance = (y_test == y_pred).sum()/X_test[:,feature_indices].shape[0]

    return Performance


def niche(x):

    niche = int(sum(x))

    return niche
	

def random_selection(MAP_solutions, n_att):
	
    x1 = np.zeros(n_att) 
    
    x2 = np.zeros(n_att) 
    
    while (any(x1) == 0 or any(x2) == 0):
		
        x1 = random.choice(MAP_solutions)

        x2 = random.choice(MAP_solutions)

    return x1, x2			


def random_variation(x1, x2, n_att):

    x = np.zeros(n_att)

    # Crossover of the parents x1 and x2

    #crossover_point = random.randrange(start = 0, stop = n_att , step = 1) 

    #x[0:crossover_point] = x1[0:crossover_point]

    #x[crossover_point+1:] = x2[crossover_point+1:]

    # Random mutation of X

    random_mutation_site = random.randrange(start = 0, stop = n_att , step = 1) 

    if (x[random_mutation_site] == 1):

        x[random_mutation_site] = 0

    else:

        x[random_mutation_site] = 1
    
    return x 

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


def MAP_Elites(X_train, y_train, X_test, y_test, I, G, application):

    # X_train, y_train, X_test, y_test: the training and test data

    # I: the amount of iterations that will be done before the main loop is terminated

    # G: the amount of initial random solutions that will be created

    n_att = len(X_train[1,:])

    amount_of_niches = n_att

    MAP_solutions = [np.zeros(n_att) for i in range(amount_of_niches)] # = X 
    MAP_performances = [0 for i in range(amount_of_niches)] # = P

    # 1. Initialize the algorithm by randomly generating G genomes and determining the performance and features of each solution + place them in the archive

    for i in range(1, G, 1):

        # 1.1 Generate a random vector of 0's and 1's of length n_att

        #x = np.random.choice([0, 1], size=(n_att))

        m = random.randrange(start = 0, stop = amount_of_niches , step = 1) 
        
        x = initial_mutant(n_att, m)

        # 1.2 Evaluate the performance/fitness p of x
		
        p = performance(x, X_train, y_train, X_test, y_test, application)

        # 1.3 Find out which niche the solution x belongs to

        niche_x = niche(x)

        # 1.4 Check if the niche of the random solution is still empty or the solution is better than the current elite in its specific niche. If either one of those is true, the solution will occupy the niche.
    
        if MAP_performances[niche_x] < p:

			# Store the performance of x in the map of elites according to its feature descriptor/niche
			
            MAP_performances[niche_x] = p

			# Store the solution x in the map of elites according to its feature descriptor/niche

            MAP_solutions[niche_x] = x

    #  2. Repeat the following evolutionary loop for I iterations 

    for i in range(1, I, 1):

        # All subsequent solutions are generated from elites in the map

		# 2.1 Randomly select 2 elites x1 and x2 from the archive MAP_solutions
		
        x1, x2 = random_selection(MAP_solutions, n_att)

		# 2.2 Create x_new, a randomly modified copy of x1 and x2 (via mutation and/or crossover)
			
        x_new = random_variation(x1, x2, n_att)

		# 2.3 Evaluate the new solution x_new and define its niche
		
		# 2.3.1 Evaluate the performance/fitness p_new of x_new
		
        p_new = performance(x_new, X_train, y_train, X_test, y_test, application)

        # 2.3.2 Find out which niche the new solution x_new belongs to
 
        niche_new = niche(x_new)

		# 2.4 Check if the niche of the new solution is still empty or the solution is better than the current elite in its specific niche. If either one of those is true, the solution will occupy the niche.
        
        if MAP_performances[niche_new] < p_new:

			# Store the performance of x in the map of elites according to its feature descriptor/niche
			
            MAP_performances[niche_new] = p_new

			# Store the solution x in the map of elites according to its feature descriptor/niche

            MAP_solutions[niche_new] = x_new

    return MAP_solutions, MAP_performances # Feature and performance map (P and X)

#----------

## Preparation of data for Naive Bayes (has to be done for each dataset separately, because they are all constructed differently)

## Ionosphere (classification)

# 1. Load the Ionosphere data file (had to use absolute path, because relative path did not work in VSCODE for me)

Ionosphere_data = open("/Users/Shaun/Documents/School/2021 - 2022/Sem 2/a. Master thesis/Programming/Master Thesis/MAP-Elites Toy example/Data/ionosphere.data")

# 2. y = and X = from Ionosphere data

X_Ionosphere = np.zeros((351, 34))

y_Ionosphere = []

cnt = 0

for l in Ionosphere_data: 

    lines = l.strip()

    for i in range(0, 34, 1):

        X_Ionosphere[cnt, i] = float(lines[0:].split(',')[i])

    y_Ionosphere.append(lines[-1])

    cnt += 1 

y_final_Ionosphere = np.zeros((351))

cnt2 = 0

for i in y_Ionosphere:

    if i == 'g':
       
       y_final_Ionosphere[cnt2] = 1
       
       cnt2 += 1 

    else:
         
         y_final_Ionosphere[cnt2] = 0

         cnt2 += 1 

# 3. Split data in train and test set

X_train_Ionosphere, X_test_Ionosphere, y_train_Ionosphere, y_test_Ionosphere = train_test_split(X_Ionosphere, y_final_Ionosphere, test_size=0.5, random_state=0)

# MAP-Elites for Ionosphere feature selection

I = 5000 # Number of generations

G = 500 # Number of initial individuals generated

application_Ionosphere = 'classification'

MAP_solutions_Ionosphere, MAP_performances_Ionosphere = MAP_Elites(X_train_Ionosphere, y_train_Ionosphere, X_test_Ionosphere, y_test_Ionosphere, I, G, application_Ionosphere)

print(MAP_solutions_Ionosphere)
np.round(MAP_performances_Ionosphere, 2)

plt.bar(range(0, 34, 1), MAP_performances_Ionosphere)

plt.plot(range(0, 34, 1), MAP_performances_Ionosphere)

plt.show()

sum(MAP_solutions_Ionosphere[np.argmax(MAP_performances_Ionosphere)])

np.max(MAP_performances_Ionosphere)

for i in MAP_solutions_Ionosphere:
    print(sum(i))

# Map-Elites Algorithm for Features Selection Problem paper reported a fitness of ~ 92 (corresponds to accuracy in % I thing, because it isn't explicitly stated) with ~ 13 atributes

# The results here are similar or just slightly worse. I think that in the paper, they ran the algorithm multiple times and took averages of fitness and amount of atributes.

performances_Ionosphere = np.zeros((21))

number_of_attributes_Ionosphere = np.zeros((21))

for i in range(0,21):
   
    map_sol, map_perf = MAP_Elites(X_train_Ionosphere, y_train_Ionosphere, X_test_Ionosphere, y_test_Ionosphere, I, G, application_Ionosphere)

    id = np.argmax(map_perf)

    performances_Ionosphere[i] = map_perf[id]  

    number_of_attributes_Ionosphere[i]  = sum(map_sol[id])
    
# average_performance_Ionosphere

sum(performances_Ionosphere)/len(performances_Ionosphere)

# average_number_of_attributes_Ionosphere

sum(number_of_attributes_Ionosphere)/len(number_of_attributes_Ionosphere)

plt.bar(number_of_attributes_Ionosphere, performances_Ionosphere)

plt.show()

## Glass (classification)

# 1. Load the Glass data file (had to use absolute path, because relative path did not work in VSCODE for me)

Glass_data = open("/Users/Shaun/Documents/School/2021 - 2022/Sem 2/a. Master thesis/Programming/Master Thesis/MAP-Elites Toy example/Data/glass.data")

# 2. y = and X = from Glass data

X_Glass = np.zeros((214, 9))

y_Glass = []

cnt = 0

for l in Glass_data: 

    lines = l.strip()

    for i in range(1, 10, 1):

        X_Glass[cnt, i-1] = float(lines[1:].split(',')[i])

    y_Glass.append(int(lines[-1]))

    cnt += 1 


# 3. Split data in train and test set

X_train_Glass, X_test_Glass, y_train_Glass, y_test_Glass = train_test_split(X_Glass, y_Glass, test_size=0.5, random_state=0)

# MAP-Elites for Glass feature selection

I = 5000 # Number of generations

G = 500 # Number of initial individuals generated

application_Glass = 'classification'

MAP_solutions_Glass, MAP_performances_Glass = MAP_Elites(X_train_Glass, y_train_Glass, X_test_Glass, y_test_Glass, I, G, application_Glass)

print(MAP_solutions_Glass)
print(MAP_performances_Glass)

plt.bar(range(0, 9, 1), MAP_performances_Glass)


plt.plot(range(0, 9, 1), MAP_performances_Glass)

plt.show()

# Averages

performances_Glass = np.zeros((21))

number_of_attributes_Glass = np.zeros((21))

for i in range(0,21):
   
    map_sol_Glass, map_perf_Glass = MAP_Elites(X_train_Glass, y_train_Glass, X_test_Glass, y_test_Glass, I, G, application_Glass)

    id_Glass = np.argmax(map_perf_Glass)

    performances_Glass[i] = sum(map_sol_Glass[id_Glass])

    number_of_attributes_Glass[i]  = map_perf_Glass[id_Glass]
        
# average_performance_Glass

sum(performances_Glass)/len(performances_Glass)

# average_number_of_attributes_Glass

sum(number_of_attributes_Glass)/len(number_of_attributes_Glass)

# These results also roughly correspond to what's in the paper.
# Is it possible that there might be a difference (between these result and the paper) because of the exact implementation of Naive Bayes used? (+ no 2 fold cross validation)


## Diagnostic Breast Cancer (classification)

# 1. Load the Cancer data file (had to use absolute path, because relative path did not work in VSCODE for me)

Cancer_data = open("/Users/Shaun/Documents/School/2021 - 2022/Sem 2/a. Master thesis/Programming/Master Thesis/MAP-Elites Toy example/Data/wdbc.data")

# 2. y = and X = from Cancer data

X_Cancer = np.zeros((569, 30))

y_Cancer = []

cnt = 0

for l in Cancer_data: 

    lines = l.strip()

    for i in range(1, 31, 1):

        X_Cancer[cnt, i-1] = float(lines.split(',')[i+1])

    y_Cancer.append(lines.split(',')[1])

    cnt += 1 

cnt2 = 0

y_final_Cancer = np.zeros((569))

for i in y_Cancer:

    if i == 'B':
       
       y_final_Cancer[cnt2] = 1
       
       cnt2 += 1 

    else:
         
         y_final_Cancer[cnt2] = 0

         cnt2 += 1 

# 3. Split data in train and test set

X_train_Cancer, X_test_Cancer, y_train_Cancer, y_test_Cancer = train_test_split(X_Cancer, y_final_Cancer, test_size=0.5, random_state=0)

# MAP-Elites for Cancer feature selection

I = 5000 # Number of generations

G = 500 # Number of initial individuals generated

application_Cancer = 'classification'

MAP_solutions_Cancer, MAP_performances_Cancer = MAP_Elites(X_train_Cancer, y_train_Cancer, X_test_Cancer, y_test_Cancer, I, G, application_Cancer)

print(MAP_solutions_Cancer)
print(MAP_performances_Cancer)

plt.bar(range(0, 30, 1), MAP_performances_Cancer)

plt.plot(range(0, 30, 1), MAP_performances_Cancer)

plt.show()

# Averages

performances_Cancer= np.zeros((21))

number_of_attributes_Cancer = np.zeros((21))

for i in range(0,21):
   
    map_sol_Cancer, map_perf_Cancer =  MAP_Elites(X_train_Cancer, y_train_Cancer, X_test_Cancer, y_test_Cancer, I, G, application_Cancer)

    id_Cancer = np.argmax(map_perf_Cancer)

    performances_Cancer[i] = sum(map_sol_Cancer[id_Cancer])

    number_of_attributes_Cancer[i]  = map_perf_Cancer[id_Cancer]
        
# average_performance_Cancer

sum(performances_Cancer)/len(performances_Cancer)

# average_number_of_attributes_Cancer

sum(number_of_attributes_Cancer)/len(number_of_attributes_Cancer)


## Parkinson's disease (regression)

# 1. Load the Parkinsons data file (had to use absolute path, because relative path did not work in VSCODE for me)

Parkinsons_data = open("/Users/Shaun/Documents/School/2021 - 2022/Sem 2/a. Master thesis/Programming/Master Thesis/MAP-Elites Toy example/Data/parkinsons_updrs.data")

# 2. y = and X = from Parkinsons data

X_Parkinsons = np.zeros((5875, 21))

y_Parkinsons = []

cnt = 0

for l in Parkinsons_data: 

    if cnt > 0:

        lines = l.strip()

        cnt_extra = 0

        for i in range(0, 22, 1):

            if i != 5:

                X_Parkinsons[cnt-1, cnt_extra] = float(lines.split(',')[i])

                cnt_extra += 1
        
        y_Parkinsons.append(float(lines.split(',')[5]))

    cnt += 1 



# 3. Split data in train and test set

X_train_Parkinsons, X_test_Parkinsons, y_train_Parkinsons, y_test_Parkinsons = train_test_split(X_Parkinsons, y_Parkinsons, test_size=0.5, random_state=0)

# MAP-Elites for Parkinsons feature selection

I = 5000 # Number of generations

G = 500 # Number of initial individuals generated

application_Parkinsons = 'regression'

MAP_solutions_Parkinsons, MAP_performances_Parkinsons = MAP_Elites(X_train_Parkinsons, y_train_Parkinsons, X_test_Parkinsons, y_test_Parkinsons, I, G, application_Parkinsons)

print(MAP_solutions_Parkinsons)
print(MAP_performances_Parkinsons)

for i in MAP_performances_Parkinsons:
    print(i)

plt.bar(range(0, 21, 1), MAP_performances_Parkinsons)

park =[]

for i in MAP_performances_Parkinsons:
    if i > 0:
        park.append(1/i)


plt.plot(range(1, 21, 1), park, 'o')

plt.ylabel('MSE')
plt.xlabel('Complexity')
plt.title('Performance vs Complexity for the Parkinsons results')
plt.xticks(range(1, 21, 2))


plt.show()

1/MAP_performances_Parkinsons[-1]

# Averages

performances_Parkinsons = np.zeros((21))

number_of_attributes_Parkinsons = np.zeros((21))

for i in range(0,21):
   
    map_sol_Parkinsons, map_perf_Parkinsons = MAP_Elites(X_train_Parkinsons, y_train_Parkinsons, X_test_Parkinsons, y_test_Parkinsons, I, G, application_Parkinsons)

    id_Parkinsons = np.argmax(map_perf_Parkinsons)

    performances_Parkinsons[i] = sum(map_sol_Parkinsons[id_Parkinsons])

    number_of_attributes_Parkinsons[i]  = map_perf_Parkinsons[id_Parkinsons]
        
# average_performance_Parkinsons

sum(performances_Parkinsons)/len(performances_Parkinsons)

# average_number_of_attributes_Parkinsons

1/(sum(number_of_attributes_Parkinsons)/len(number_of_attributes_Parkinsons))
