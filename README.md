# Simulation code for *Savanna canopy trees under fire: Long-term persistence and transient dynamics from a stage-based matrix population model*
## Patricia A Werner and Stephanie J Peacock*
## Accepted to Ecosphere March 11, 2019
*All questions regarding this code should be directed to stephanie at <stephanie.j.peacock@gmail.com>

## Contents

```1-baseline_matrices.R``` : stable stage distributions and sensitivity to the baseline matrices given in Table 3. 

```2a-scenarios-long-term-dynamics.R```: Calculation of lambda (deterministic and stochastic) for different fire scenarios given in Table 5, with stochastic simulations run for a large number of timesteps.

```2b-scenarios-transient-dynamics.R```: Calculation of lambda (deterministic and stochastic) for different fire scenarios given in Table 5 over short-term (100-year "transient") dynamics, repeating stochastic simulations for 1000 MC trials.

```3a-recruitment_and_seedling_survival.R```: sensitivity of lambda to changing recruitment and seedling survival (*results not included in Ecosphere paper*)

```3b-patchiness.R```: sensitivity to changing fire patchiness for early dry season fires (*results not included in Ecosphere paper*)

```4-Kapalga-simulations.R```: Code to run the simulations from 1982-2002 for Kapalga National Park.

```functions.R```: functions called on by other .R files

```matrices.csv```: the different transition matrices used in the model simulations (Table 3)
