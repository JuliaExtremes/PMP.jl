module PMP

using Distributions, Optim, Statistics, SpecialFunctions, MambaLite, ForwardDiff
using DataFrames, Dates, CSV, Extremes, RollingFunctions, LogExpFunctions

import Distributions: @check_args, location, scale, shape, params
import Distributions: cdf, insupport, logpdf, minimum, maximum, quantile, rand 
import Distributions: mean, var, std, modes, mode, skewness, kurtosis, entropy
import Distributions: fit_mle, Beta
import SpecialFunctions: loggamma, logbeta
import LogExpFunctions: logit, logistic, xlogy, xlog1py

import Random

import DataFrames

import Extremes: returnlevel

# distributions
include("distributions/pearsontype1.jl");
include("distributions/pearsontype1b.jl");
include("distributions/pearsontype1c.jl");

# usual PMP 
include("moisture_maximization.jl");
include("other_PMP_methods.jl");

# documentation example
include("data_example.jl");

export
    # distribution types
    PearsonType1,
    PearsonType1b,
    PearsonType1c,

    # methods
    params,      # get the tuple of parameter
    location,
    scale,
    shape,
    minimum,
    maximum,
    insupport,   # predicate, is x in the support of the distribution?
    cdf,         # cumulative distribution function
    logpdf,      # log probability density
    pdf,         # probability density function
    quantile,    # inverse of cdf (defined for p in (0,1))
    rand,
    
    mean, 
    var, 
    std, 
    modes, 
    mode, 
    skewness,
    kurtosis, 
    entropy,

    # distribution fitting
    fit_mme,
    fit_mle,
    getinitialvalues,
    fit_bayes,

    # moisture maxmization
    total_precipitation,
    storm_selection,
    get_max_persisting_dew,
    dewpoint_to_PW,
    PW_storm,
    PW_max,
    PW_return_period,
    storm_maximization,
    PMP_mm,

    # other PMP estimation methods
    PMP_GEV,
    PMP_Hershfield

end # module PMP
