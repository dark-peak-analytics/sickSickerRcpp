# sickSickerRcpp | DPA <img src="https://github.com/RobertASmith/darkpeak/blob/main/man/figures/logo_concise.PNG" align="right" width="120" />

## Background

Microsimulation allows for the incorporation of indivividual level hetrogeneous effects into models. However, they are often slow to run relative to cohort models. We have written this example to demonstrate how a model can be sped up using vectorisation and C++. These methods can be integrated with existing statistical analysis and reporting pipelines, including automated (living) HTA updates and web-based applications for models.

We uses the sick-sicker model published by [Krijkamp et al. 2018.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6349385/) as a case study model to demonstract functionality.

We vectorise the model in C++ and use the Rcpp package to run the model from R. This runs MUCH faster because:  
1. all individuals run through the model simultaneously (vectorised)  
2. C++ is faster than R for these operations  
Of these two things (1) is probably more important!


## Instructions

Please ensure that you open this project by launching the sickSickerRcpp.Rproj file.

If this is not possible then:
1. set the working directory to the project folder (sickSickerRcpp/)
e.g. `setwd("C:/Users/r_a_s/Documents/Projects/DPA_courses/sickSickerRcpp")`
2. ensure the file 'src/define_microsim.cpp' is in the working directory
`"define_microsim.cpp" %in% list.files("src/")` 
3. install the "Rcpp" and "RcppArmadillo" packages
```
install.packages("Rcpp")
install.packages("RcppArmadillo")
```

Run the example code in the R script `example_script.R` or the Rmarkdown vignette `vignette/run_microsim.rmd`.

These files source the functions defined in 'src/define_microsim.cpp' (also in this repository) and then pass some parameters to the `MicroSimV_Cpp` R function before collating some results. This is intended as a barebones example, a full course is under development.

Code developed by Wael Mohammed & Robert Smith at [Dark Peak Analytics](https://darkpeakanalytics.com/).

[Dark Peak Analytics](https://darkpeakanalytics.com/) provide consulting services at the intersection of health economics and data-science. Established by Paul Schneider & Robert Smith in 2019, we have experienced high demand and continue to grow the team. Interested in collaborating on a research, teaching or consulting project with us, contact us via [email](contact@darkpeakanalytics.com).
