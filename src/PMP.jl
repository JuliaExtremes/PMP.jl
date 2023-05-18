module PMP

using Distributions

import Distributions: @check_args, location, scale, shape, params
import Distributions: cdf, insupport, logpdf, minimum, maximum, quantile, rand

import Random

# distributions
include("functions.jl")
include("distributions/pearsontype1.jl");

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
    rand

end # module PMP
