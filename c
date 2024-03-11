# are there degeneracires in fits between ebv and beta? does curve fit always 
# recover correct input values?

# reduced xi squared = xi squared over dof (N - no of free params)

# plotting input ebv against resulting beta, with reduced xi^2 color coding
# plotting input beta against resulting ebv, with reduced xi^2 color coding

# how large can simulated uncertainties be before degeneracies between beta and ebv

# how does xi^2 vary if fit simulated extinguished afterglow with different
# extinction laws? should give worse results for wrong extinction law

...

# introduce x ray data: dominated by absorption from metal ions, not dust
# turnover at lower frequency end of spectrum
# photoelectric absorption 
# extend simulation to x rays, with 50 data points (.3 - 10 KeV)
# paper on absorption of x rays in ism: observed intensity as a fn of energy
# Iobs(E) = I0(E) exp[-sigma(E) NH]
# NH = column density - how much mateial along line of sight
# amount of H atoms
# use sigma function from code given then use np where to find values
