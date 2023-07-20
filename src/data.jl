"""
    dataset(name::String)::DataFrame

Load the dataset associated with `name`.

Datasets available:
 - `rain`: observed precipitations (in mm) recorded at the Montréal-Trudeau International Airport;
 - `dewpoint`: observed dew point (in °C) recorded at the Montréal-Trudeau International Airport.


# Examples
```julia-repl
julia> PMP.dataset("rain")
```
"""
function dataset(name::String)::DataFrame

    filename = joinpath(dirname(@__FILE__), "..", "data", string(name, ".csv"))
    if isfile(filename)
        # return DataFrame!(CSV.File(filename))
        return CSV.read(filename, DataFrame)
    end
    error("There is no dataset with the name '$name'")

end