# Penguins & Oceanography

Code supporting the manuscript **Mesoscale eddies in the Southern Ocean influence foraging trip behaviour of multiple penguin species** by Wilson et al. (2026).

---

## Overview

This repository contains scripts used to relate penguin foraging trajectories to mesoscale eddies, by fitting Hidden Markov Models (HMMs) to isolate Area-Restricted-Search behaviour and then using Generalised Additive Mixed Models (GAMMs) to relate ARS incidence to eddies.

---

## Repository Structure

```text
├── functions/                                  # Functions used in other scripts
│
├── 01a_hidden_markov_models.R                  # Fit Hidden Markov Models to tracks
├── 01b_hmm_quality_control.R                   # Inspect and remove trips with poor HMM assignment
│
├── 02a_background_central_place_foraging.R     # Create background samples for central place foraging case studies
├── 02b_backgorund_free_ranging.R               # Create background samples for other case studies
│
├── 03a_environmental_extraction.R              # Extract environmental variables to ARS locations and background samples
├── 03b_thinning_tracks.R                       # Decorrelate tracks and prepare data for modelling
│
├── 04a_gamm_meta_king_macaroni.R               # Fit GAMMs for king and macaroni penguin case studies (no sea ice) using the Mesoscale Eddy Trajectory Atlas
├── 04b_gamm_meta_chinstrap_adelie_emperor.R    # Fit GAMMs for chinstrap, Adelie, and emperor penguin case studies (with sea ice) using the Mesoscale Eddy Trajectory Atlas
├── 04c_gamm_sea_ice_eddies.R                   # Fit GAMMs for chinstrap, Adelie, and emperor penguin case studies (with sea ice) using a specialised sea ice eddy detection algorithm
│
├── 05_summarise_specialisation_by_group.R      # Compare how regional, interspecific, and behavioural variation affect the proportion of case studies that are specialised
│
├── 06_*.R                                      # Investigate five comprehensive case studies in more detail
│
├── 07_*.R                                      # Match and test how eddy maturity, amplitude, and intensity vary in the five case studies
│
├── 99_create_eddy_raster.R                     # Create a daily circumpolar raster using the Mesoscale Eddy Trajectory Atlas
├── 99_create_eddy_sea_ice_raster.R             # Create a daily circumpolar raster using a specialised sea ice eddy detection algorithm
├── 99_plot_all_case_studies.R                  # Plot each case study for supplementary information
├── 99_plot_all_case_studies_cross_dateline.R   # Plot each case study for supplementary information (for case studies that cross the international dateline)
├── 99_schematic_plots.R                        # Produce plots for the schematic (Fig. 2)

```

---

## Requirements

The code is written entirely in **R**

---

## Citation

If you use this repository in your research, please cite the associated publication:

Joshua C. Wilson, Philip N. Trathan, Hugh J. Venables, David G. Ainley, Matthis Auger, Alistair M.M. Baylis, Charles-André Bost, Louise Emmerson, Kimberley T. Goetz, Mark A. Hindell, Jefferson T. Hinke, Catharine Horswill, Aymeric Houstin, Akiko Kato, Nobuo Kokubun, Gerald L. Kooyman, Małgorzata Korczak-Abshire, Sara Labrousse, Céline Le Bohec, Andrew D. Lowther, Phil O’B Lyver, Ana Laura Machado-Gaye, Azwianewi B Makhado, Silvia Olmastroni, Pierre A. Pistorius, Klemens Pütz, Norman Ratcliffe, Yan Ropert-Coudert, Peter G. Ryan, Jean-Baptiste Sallée, Mercedes Santos, Richard B. Sherley, Alvaro Soutullo, Akinori Takahashi, Barbara Wienecke, Daniel P. Zitterbart, Ryan R. Reisinger (2026). Mesoscale eddies in the Southern Ocean influence foraging trip behaviour of multiple penguin species, Progress in Oceanography, 247, 103798, https://doi.org/10.1016/j.pocean.2026.103798.


---

## Contact

**Author:** Joshua Wilson
