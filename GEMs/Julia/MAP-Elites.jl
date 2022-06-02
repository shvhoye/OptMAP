#=
MAP-Elites:
- Author: Shauny Van Hoye
- Date: 2022-02-16
=#

using Random, Distributions

## Functions and variables

dimentions = 2

A = 10

function random_solution()

	d = Uniform(-2*pi, 2*pi)

	return rand(d , dimentions)

end


function random_selection(X)
	
	selection = [100.0, 100.0]

	map = copy(X)
	
	while selection == [100.0, 100.0]
		
			selection = rand(map)
	
	end

    return selection		

end


function random_variation(z; scale=1)

    x = copy(z)

    d = Normal(0.0, scale)

    variation1 = rand(d)
    
    variation2 = rand(d)

    x[1] = x[1] + variation1

    x[2] = x[2] + variation2
    
    return x 
end


function performance(x)
	
	# The rastrigin function
	    
	# Since MAP_elites maximizes a function: performance(x) = -rastrigin
	
	# Maximum = 0
	
	n = length(x)
	
	return -(A*n + sum(x.^2 - A*cos.(2*pi*x)))

end


function niche(x)

	if x[1] > 0

		i = 1
		
	else
		
		i = 2

	end
	
	if x[2] > 0

		j = 1
		
	else
		
		j = 2

	end

	
	return i, j
	
end


function MAP_Elites(random_solution, random_selection, random_variation, performance, niche, N; max_iteration = 10, a = [100.0, 100.0], b=100.0)
	
	# 0. Create an empty, N-dimensional map of elites: {solutions X and their performances P}
	
	MAP_solutions = fill(a, N, N) # = X
	MAP_performances = fill(b, N, N) # = P

	#  Repeat for I iterations 

	iterations = 0
		
	while iterations < max_iteration
		
		if iterations < 1 

			# 0. Initialize by generating G random solutions
			
			x′ = random_solution() 
		
		else
			
			# All subsequent solutions are generated from elites in the map

			# 1. Randomly select an elite x from the map X (MAP_solutions)
			
			x = random_selection(MAP_solutions) 

			# 2. Create x′, a randomly modified copy of x (via mutation and/or crossover)
			x′ = random_variation(x)
			
		end

		# 3. Score + define nich of new individual
	
			# SCORE
			
			# current_score = ...?
		
		# 3. Simulate the candidate solution x′ and record its feature descriptor b′

		# b′ = feature_descriptor(x′)
		
		# NICHE
		
		i_current, j_current = niche(x′)

		# Record the performance p′ of x′
		
		p′ = performance(x′)



		# 4. Check if the new individual is better than the current elite in its specific niche. If the appropriate cell is empty or its occupants’s performance is ≤ p′, then:

		if MAP_performances[i_current, j_current] == 100.0 || MAP_performances[i_current, j_current] < p′ 

			# store the performance of x′ in the map of elites according to its feature descriptor b′
			
			MAP_performances[i_current, j_current] = p′

			# store the solution x′ in the map of elites according to its feature descriptor b′

			MAP_solutions[i_current, j_current] = x′

		end


		# 5. Update iterations
	
			iterations += 1

	end

	
	return MAP_solutions, MAP_performances # feature-performance map (P and X)

end


function ackley(x; a=20, b=0.2, c=2π)
   
    d = length(x)
    
    return -(-a * exp(-b*sqrt(sum(x.^2)/d)) - exp(sum(cos.(c .* x))/d))

end

	
function niche_ackley(x)

	if x[1] < -0.5

		i = 1
		
	elseif -0.5 < x[1] < 0.5
		
		i = 2

	else
		
		i = 3
		
	end
	
	if x[2] < -0.5

		j = 1
		
	elseif -0.5 < x[2] < 0.5
		
		j = 2

	else
		
		j = 3
		
	end

	
	return i, j
	
end

# Tests

MAP_solutions, MAP_performances = MAP_Elites(random_solution, random_selection, random_variation, ackley, niche_ackley, 3;  max_iteration = 10000000)

print(MAP_solutions, MAP_performances)