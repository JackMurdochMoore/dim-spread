# dim-spread
Code implementing dimensional spreading SIR model and comparing it with ground truth and other models

__Associated with the paper__\
"Network spreading from network dimension"\
by\
Jack Murdoch Moore, Michael Small, Gang Yan, Huijie Yang, Changgui Gu, and Haiying Wang.

__Correlation dimension__ D is defined by the power-law†\
c(s) ∝ s^(D-1),\
where D is correlation dimension, and c(s) is correlation (fraction of distinct nodes at distance s).

__To produce some figures__ from the manuscript, please run:\
run_model_comparison\
run_alt_model_comparison

__Folders:__\
_networks_: Data defining empirical networks.\
_results-time-series-new_: Folder in which results are saved.

__Functions and scripts:__\
_BA_mod_2.m_: Generate BA scale-free network.\
_count_distances.m_: Return vector of network distances and number of pairs of distinct nodes at each network distance.\
_est_corr_dim_new_1.m_: Estimate correlation dimension† and scaling interval of a networks using different methods and model c(s) ∝ s^(D-1).\
_inclusivity.m_: Generate correlated version\# of BA scale-free network.\
_load_network.m_: Load an empirical network from data in folder "networks".\
_run_alt_model_comparison.m_: Compare Monte Carlo (ground truth) with dimensional spreading, reduced effective degree, homogeneous pair approximation, and pair-based models. Also return mean Euclidean error for each model.\
_run_model_comparison.m_: Compare Monte Carlo (ground truth) with dimensional spreading, homogeneous mean field, heterogeneous mean field and PDMC models. Also return mean Euclidean error and R_0 for each model.\
_run_sir_0_mod_2.m_: Run ground truth discrete SIR model.\
_run_sir_dim_new_1.m_: Run dimensional spreading SIR model.\
_run_sir_pair_based_model.m_: Run pair-based SIR model§.\
_small_world_manhattan.m_: Generate lattice\* or small world network\*.\
_small_world_manhattan_lcc.m_: Generate lattice\* or small world network* and retain only its largest connected component.\

\* Lattices/small world networks are/are derived from regular d-dimensional toroidal lattices defined using a periodic version of the city block (or Manhattan or taxi cab) metric mentioned but not explored in "Epidemic dynamics on higher-dimensional small world networks", Applied Mathematics and Computation 421, 126911, by H. Wang, J. M. Moore, M. Small, J. Wang, H. Yang and C. Gu (2022) (associated code at https://github.com/JackMurdochMoore/small-world). 

\# Correlated versions of BA networks are inclusivity model networks introduced in "Inclusivity enhances robustness and efficiency of social networks”, Physica A 563, 125490, by J.M. Moore, M. Small, and G. Yan (2021).

† Correlation dimension is estimated following "Correlation dimension in empirical networks”, "Correlation dimension in empirical networks" by J.M. Moore, H. Wang, M. Small, G. Yan, H. Yang, and C. Gu (associated code at https://github.com/JackMurdochMoore/net-corr-dim).

§ Code for the pair-based model is adapted from pair_based_model.m by K. J. Sharkey (2010), which is available in the supplementary material for “Deterministic epidemic models on contact networks: correlations and unbiological terms”, Theoretical Population Biology, Volume 79, Issue 4, % June 2011, Pages 115-129.


