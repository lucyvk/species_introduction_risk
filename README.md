# species_introduction_risk

This repository consists of the scripts used for analysis for the paper "A predictive framework for species introductions under uncertainty and multiple management objectives"

ATN.py
- Script with functions used for running dynamic simulations of the allometric trophic network model

niche.py
- Script with functions used for generating synthetic networks with the niche model

paper_webs_0624.zip
- A zipped folder containing the structure of the 37 networks used for analysis in the paper
- The intermediate files for the structure of the network after every possible introduction is included for one example web for demonstration though otherwise not included due to size; the rest of these intermediate files can be produced using the scripts if desired
- The files for the final sets of metrics used in the visualization notebook are also included 

The following notebooks are run in order:

1_species_risk_paper_generate_webs
- Jupyter notebook used to generate a plausible set of food webs to be used for analysis

2_species_risk_paper_generate_invaders
- Jupyter notebook used to generate niche model parameters for introduced species

3_species_risk_paper_generate_invaded_structure
- Generate and save the structure of each of the networks after species are introduced

4_species_risk_paper_run_simulations
- Run the dynamics for all of the networks after introduction (and in counterfactual case without any introduction)

5_species_risk_paper_generate_metrics
- Generate metrics related to management outcomes from the final biomasses of species within the networks after the dynamics are run

6_species_risk_paper_viz
- Produce plots based on the final metrics including tradeoff plots

7_extra_si_figures
- Additional plotting code for si figures (exploring realized trophic level and in and out degree of introduced species)
