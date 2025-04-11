# paternal (maternal) offspring chromosome is mozaic of its 2 sire's (dam's) chromosomes
# viterbi algoritm to find the mozaic (called best_path in viterby terms) 
function viterbi(OBS,  recomb_rate, error_rate)
    # Initialize the Viterbi trellis and backpointers
        T = size(OBS,1)
	N = 2  #2 states
	V = zeros(N, T)  #Viterbi Path-probability (state(i),time(j))
	B = zeros(Int, N, T) #trace of maximum path probs
 trans_p= [(1-recomb_rate) recomb_rate       # state transition from position t to t+1
           recomb_rate  (1-recomb_rate)]    # recombination prob is recomb_rate

    # set up emit_prob
    obs=deepcopy(OBS)
    obs[obs.==0].=3   #missing
	emit_p=[(1-error_rate) error_rate 1.0
	        error_rate (1-error_rate) 1.0] #emission probabilities; 1.0 in row&column 3 denotes no info (missing data)
	  
    # Initialize the first column of the Viterbi-probs  (prior_prob=0.5)
       for i=1:N
          V[i,1]=0.5*emit_p[i,obs[1]]
	end  
 
    # Perform the Viterbi algorithm
        for t=2:T
	        for i=1:N
		      (max_prob,max_state)=findmax(V[:,t-1].*trans_p[:,i])
            	      V[i, t] = max_prob*emit_p[i,obs[t]]
	              B[i, t] = max_state
	        end
    	end

    # Find the best path by backtracking
        best_path = zeros(Int, T)
	(best_prob,best_state)=findmax(V[:,T])
     best_path[T] = best_state
    for t in T-1:-1:1
        best_path[t] = B[best_path[t+1], t+1]
    end
    return best_path
    end

