# MWE of some FTS errors

using Oceananigans
using Oceananigans.Utils: Time 

# First we generate some data and write to disk
grid = RectilinearGrid(size=(2, 2, 2), extent=(1, 1, 1))
times = 0:0.1:3
filename = joinpath(@__DIR__, "MWE_data_file.jld2")
f_tmp = Field{Center,Center,Center}(grid) 
f = FieldTimeSeries{Center, Center, Center}(grid, times; backend = OnDisk(), path = filename, name = "f")

for (it, time) in enumerate(f.times)
    set!(f_tmp, (x, y, z) -> sin(2* π * time/3) )
    set!(f,f_tmp,it)
end

# Now we load the FTS partly in memory
N_in_mem = 5
f_fts = FieldTimeSeries(filename, "f"; backend = InMemory(N_in_mem))
Nt = length(f_fts.times)

# Bug 1: We can't access elements from Nt - N_in_mem + 2 to Nt unless we first access a previous element that loads it in memory.
# This is due to line 268 of field_time_series_indexing.jl, which tries to load into memory times that are out of bounds
function try_accessing_element(FTS,n)
    println("Trying to access element $n")
    try
        println(FTS[n])
    catch e
        println(e)
    end
end

try_accessing_element(f_fts,Nt-N_in_mem +2 ) 
try_accessing_element(f_fts,Nt-N_in_mem +1 ) 
try_accessing_element(f_fts,Nt-N_in_mem +2 ) 

# Bug 2: field_time_series_indexing.jl

# When linearly interpolating a FTS and going out of range of what's in memory, the update backend logic isn't right. 
# Suppose we are interpolating between indices 4 and 5. Currently, only (2,3,4) are in memory (that is, InMemory(2,3)). Then 
# fts[4] is loaded first (line 146, no problem) before fts[5] is loaded, which requires a new backend. This new backend is 
# InMemory(5,3) which is a) not efficient as we'll probably need to go back to InMemory(4,3) for the next timestep and b)
# causes an update of what's in memory mid-computation, which gives an error:

for i in [0,0.39,0.4,0.41,0.42,0.43]
    println("\nTime = $i")
    println(maximum(f_fts[Time(i)]))
    println(f_fts.backend)
end
# See value at time 0.42

# Bug 3: ≈ isn't good enough for values near zero (affects line 251 field_time_series.jl and line 10 set_field_time_series.jl)
time_range = range(0, 1, length=11)
time_arr = [0:0.1:1;]
time_arr[1] = 1e-16 
println(all(time_range .≈ time_arr))
println(isapprox(time_arr,time_range))