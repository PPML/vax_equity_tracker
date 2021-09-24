# Tracking Racial/Ethnic Disparities in COVID-19 Vaccination

## *Projections Last Updated August 16, 2021*

### Summary: Tracking racial/ethnic disparities in COVID-19 vaccination coverage across US states is essential to advancing vaccine equity. State-level data needed to monitor levels and trends in vaccination disparities are incomplete and diverse. We applied consistent rules to reconcile data, estimate historical coverage, and project future coverage of COVID-19 vaccinations.

#### This repository contains additional details on methods, analytic code, input data, and coverage estimates in table format. If you have questions or comments, please contact us at <vax.equity.tracker@gmail.com>.

#### Coverage Estimates

You can access our coverage estimates in table format [here](https://github.com/PPML/vax_equity_tracker/tree/main/estimates).


#### Methods

We use historical state-reported vaccination data by race/ethnicity to project vaccine coverage, by state and nationally, among people ages 12 and older for four mutually exclusive racial/ethnic groups (White, Black, Hispanic, and Asian). Specifically, we used data on distribution of vaccines administered by race/ethnicity extracted from state reporting dashboards by the Kaiser Family Foundation, the percentage of people ages 12 and older who have received at least one dose from the Centers for Disease Control and Prevention, and population data from the 2019 American Community Survey to estimate the share of people ages 12 years and older receiving one or more COVID-19 vaccination doses, by state and race/ethnicity. We project coverage rates through December 31, 2021 based on the average vaccination rate between July 19 and August 16 for each racial/ethnic group, by state.

Data on vaccination coverage by race/ethnicity vary by state in terms of reporting groups and completeness. We applied the following data processing steps to produce comparable estimates. We assumed vaccinations reported as “unknown” race/ethnicity were distributed proportional to shares of vaccinations with known race/ethnicity in each state. Examining vaccinations reported as “other” race/ethnicity, we found that in most states, the shares attributed to “other” greatly exceeded population shares (implying coverage >100%). We therefore adjusted shares by assuming “other” were vaccinated proportional to eligible population, and proportionally redistributed remaining vaccinations among specified racial/ethnic groups. We adjusted shares to avoid double-counting in states that report shares by race separate from shares by ethnicity. For racial/ethnic groups not reported by specific states, we assumed these groups were vaccinated proportional to population size and scaled down shares of vaccines to reported groups accordingly. We capped coverage among any racial/ethnic group at 100% of the eligible population, and in cases where implied coverage exceeded 100%, we proportionally redistributed the excess across other groups. Finally, we enforced a monotonicity constraint on cumulative coverage, removing and interpolating noisy data points that would have implied decreasing cumulative coverage over time within state and race/ethnic group (3% of data points).

Limitations of this analysis include reliance on several assumptions to address incomplete and heterogeneous reporting of vaccination data by race/ethnicity across states. Previous reporting on racial/ethnic disparities in vaccination through the CDC and other sources has not adjusted for these data discrepancies, resulting in reported coverage levels that likely underestimate actual population coverage. Although we have adopted a standard set of definitions and rules for reconciling unknown or discrepant data elements to enable transparent and comparable estimation of coverage over time and place, results must be interpreted as approximations in the context of missing and sometimes noisy data.


#### Input data

We incorporate demographic data from the American Community Survey, overall vaccination rates from the Centers for Disease Control and Prevention, and vaccination breakdowns by race/ethnicity from the Kaiser Family Foundation. You can view the input data [here](https://github.com/PPML/vax_equity_tracker/tree/main/input). 


#### Analytic Code

You can access our analytic code [here](https://github.com/PPML/vax_equity_tracker/tree/main/code).


#### Acknowledgements

The research team are members of [Stanford University’s Department of Health Policy](https://healthpolicy.fsi.stanford.edu/). We are part of the [Prevention Policy Modeling Lab (PPML)](https://ppml.stanford.edu/) and the [Stanford-CIDE Coronavirus Simulation Modeling (SC-COSMO) Consortium](https://www.sc-cosmo.org/).

Our work is supported in part by funding from the Centers for Disease Control and Prevention and the Council of State and Territorial Epidemiologists, the National Institute on Drug Abuse, the State of California, Stanford's Clinical and Translational Science Award (CTSA) to Spectrum, and Stanford's COVID-19 Research Fund. The content is solely the responsibility of the authors and does not necessarily represent the official views of funders.

We are grateful for the tireless work of analysts at the [Kaiser Family Foundation](https://www.kff.org/other/state-indicator/covid-19-vaccinations-by-race-ethnicity), who routinely extract and compile data from state dashboards to provide breakdowns of vaccinations by race/ethnicity.
 

