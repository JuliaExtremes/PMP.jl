module PMP

using Distributions, Optim
using DataFrames, Dates, CSV, Extremes

import Distributions: @check_args, location, scale, shape, params
import Distributions: cdf, insupport, logpdf, minimum, maximum, quantile, rand 
import Distributions: mean, var, std, modes, mode, skewness, kurtosis, entropy
import Distributions: fit_mle

import Random

import DataFrames

import Extremes: getcluster, returnlevel

# distributions
include("distributions/pearsontype1.jl");

include("moisture_maximization/mm_observed_data.jl");

include("other_PMP_methods.jl");

export
    # distribution types
    PearsonType1,

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
    entropy

    # distribution fitting
    fit_mme
    #fit_mle
    #fit_bayes

    # moisture maxmization
    max_rain
    total_precipitation
    storm_selection
    get_max_persisting_dew
    dewpoint_to_PW
    PW_max
    PW_return_period
    PMP_mm

    # other PMP estimation methods
    PMP_GEV
    PMP_Hershfield

end # module PMP
