sub_select is a tool set for finding "counterfactual" subsequences in a time-series.

Contents:
	Introduction
	Usage
	Related Work
	Benchmarks and Analysis
	Future Work
	Acknowledgements
	Disclaimers
	References
	
Introduction:
	Consider a stochastic process. If you want to discover how this process might change over time from a given state, you might simulate the process. You could then start the simulation at the given state and, through the power of Monte-Carlo, develop a probability distribution of states the process visits after the given state. But what if you don't possess a simulation for the process? How can you discover how the process might behave from the same given state - how can you discover its counterfactual behavior? 
	A stochastic system that is left to its own devices will, eventually, pass through each of its possible states multiple times. By looking at how such a system behaves after each visit to a given state provides researchers an understanding of the capacity that system has for change. Each time the system revisits the given state, it is like you are able to record its counterfactual behavior.
	This counterfactual behavior is central to causal analysis, and is the core of experimental methods [1]. However, there are few tools available for discovering this counterfactual data in time-series of stochastic processes. Further, finding the maximum number of counterfactual subsequences in a time-series is not a trivial task, as the selection of one subsequence may prevent the selection of others. This process is further complicated if counterfactuals relating to more than one given state are required.
	Here I present a heuristic approach to finding the maximum number of counterfactual subsequences in a time-series corresponding to two given states, such that no subsequence overlaps another. This approach is consistently better than random selection, often close to optimal, and much quicker than solving a linear program.
	
Usage:
	The most common use case scenario of this software follows the following protocol:
	# Import your time-series
	> data = numpy.load('data.npy')
	# Select indeces to start subsequences
	# We're interested in initial conditions 2 and 5, and want each subsequence to be 50 steps long
	> indeces = score_and_start(data, 2, 5, 50)
	# Display each subsequence corresponding to initial condition of 2
	> select_sequences(indeces[0], data, 50)
	
	Further documentation is available within python.
	
Related Work:
	Much work has been done in the past on finding subsequences that exhibit similar patterns. This includes the matrix profile of the time series [2], and similarity search [3]. However, these algorithms would purposefully omit large ranges of dynamics because they are dissimilar to previous behavior, and do not attempt to find an optimal number of subsequences - the core of the problem this software sets out to solve.
	The problem I have outlined can be expressed as a linear program [4], and can be solved as one. This approach, as opposed to our own, does provide a guarantee of optimality, however, due to the structure of the problem, it takes far too long to complete. For this reason, I make no comparison to this approach in our "Benchmarks" section, but do suggest that the margin between optimal and our solution is much smaller than the margin between our solution and random selection.
	
Benchmarks and Analysis:
	The image in "performance.pdf" presents the margin between random selection and our solution, as performed on sets of random data. For small data sets, there is not a very large difference between the performance of random subsequence selection, and our algorithm. However, the margine grows greatly with the data, but is tamped by the window size. This is to be expected, as smaller data sets have a small number of possible subsequences, and larger window sizes also leads to fewer possible subsequences.
	My personal work in causal analysis suggests that the number of counterfactual examples one can produce is often more important than their size. These conditions are where the given solution shines - providing great benefit over random selection.
	It should be noted that the given approach does not always perform better than random for some pathological cases. However, in such cases, the margin between random and this method is quite small.
	
Future Work:
	The first avenue of future work pertains to the development of a more nuanced optimization method for this problem. I have begun looking into Lagrangian Relaxation methods [5], where we take the linear program problem, and solve simpler versions of it using dual variables. Also under consideration is the submodular approach [6], as this problem may be able to be expressed directly as a selection of subsets from two sets. In either case, speed of the approach will be taken into consideration as the solution to this problem does not gain a siginificant advantage from being exactly optimal - a fast approach that is close to optimal is greatly prefferable.
	The second avenue of future work pertains to the expansion of utility that the algorithm provides. I would like to have the algorithm be able to select subsequences that do not begin with exactly the same initial condition. This must be considered carefully, however, as my personal experience suggests that the precision of the distribution of selected initial conditions is more important still than the number of subsequences, or their length. I would also like to introduce the ability to analyze higher-dimensional data, allowing for an enormous increase in the applicability of this algorithm.
	
Acknowledgements:
	I would like to thank Dr. Bert Huang for his guidance and suggestions in the topics of optimization, Dr. Benjamin Jantzen for helping to develop the problem I set out to solve, and Dr. Subhradeep Roy for the pathological data and direct application of this algorithm.
	
Disclaimers:
	This work is distributed under the MIT license, and is provided as-is, and without warranty. If you have suggestions, or comments, please contact me at c0lin@vt.edu. 

	
References:
	[1] J. Pearl, "Probabilities of Causation: Three Counterfactual Interpretations and Their Identification," Synthese, vol. 121, (1/2), pp. 93-149, 1999.
	
	[2] C. M. Yeh et al, "Time series joins, motifs, discords and shapelets: a unifying view that exploits the matrix profile," Data Mining and Knowledge Discovery, vol. 32, (1), pp. 83-123, 2018.
	
	[3] T. Fu, "A review on time series data mining," Engineering Applications of Artificial Intelligence, vol. 24, (1), pp. 164-181, 2011.
	
	[4] D. Bertsimas and J. N. Tsitsiklis, Introduction to Linear Optimization. Belmont, Mass: Athena Scientific, 1997.
	
	[5] K. C. Kiwiel et al, "Lagrangian Relaxation via Ballstep Subgradient Methods," Mathematics of Operations Research, vol. 32, (3), pp. 669-686, 2007.
	
	[6] F. Bach, "Learning with submodular functions: a convex optimization perspective," 2013.