
#####
##### Input file setup: has to be done only once!
#####

# Save a file named `filename`, containing a time series of surface heat flux and surface stress
# at times `times` on a grid `grid` 
function generate_forcing_data!(grid, times, filename; ρ_ocean = 1000, cp_ocean = 4000)


    T_forcing_tmp = Field{Center,Center,Center}(grid) 
    T_forcing = FieldTimeSeries{Center, Center, Center}(grid, times; backend = OnDisk(), path = forcing_file, name = "T_forcing")

    for (t, time) in enumerate(T_forcing.times)
        @info "writing down data for timestep $t and time $time"
        set!(T_forcing_tmp, (x, y, z) -> sin(π * x) * time )
        set!(T_forcing,T_forcing_tmp,t) #set!(fts::InMemoryFTS, value, n::Int) = set!(fts[n], value)
    end

end