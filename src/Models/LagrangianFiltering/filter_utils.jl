using JLD2
using JLD2: Group

function set_times_on_disk!(fields_filename; direction = "backward")
    # We want this to check if the tracers are on center points too
    jldopen(fields_filename,"r+") do file

        if !haskey(file, "direction")
            # It should currently be forward, because direction hasn't been assigned yet
            file["direction"] = "forward"
        end    
        current_direction = file["direction"]
        println("Current direction is $current_direction")
        if current_direction != direction
            Base.delete!(file, "direction")
            file["direction"] = direction
            iterations = parse.(Int, keys(file["timeseries/t"]))
            times = [file["timeseries/t/$iter"] for iter in iterations]
            T = maximum(times)
            Base.delete!(file, "timeseries/t")
            g = Group(file, "timeseries/t")
            for (i,iter) in enumerate(reverse(iterations))
                g["$iter"] = times[i]
            end

            # We should also reorder all elements so they are indexed correctly.
            #Assuming the actual iteration number doesn't matter, just the order in the jld2 group
            for var in keys(file["timeseries"])[1:end-1] # not t
                for iter in reverse(keys(file["timeseries/$var"])[2:end]) # not serialized
                    data = file["timeseries/$var/$iter"]
                    Base.delete!(file, "timeseries/$var/$iter")
                    file["timeseries/$var/$iter"] = data
                end
            end

            println("New direction is $direction")
        else
            iterations = parse.(Int, keys(file["timeseries/t"]))
            times = [file["timeseries/t/$iter"] for iter in iterations]
            T = maximum(times)
            println("No need to update direction")
        end
        return T    
    end
    
end

function create_tracers(;N=1)
    
    gC = ()  # Start with an empty tuple
    gS = ()  # Start with an empty tuple
    for i in 1:N
        new_gC = Symbol("gC", i) 
        new_gS = Symbol("gS", i) 
        gC = (gC..., new_gC)  
        gS = (gS..., new_gS)  
    end
    return (gC...,gS...)
end

function set_forcing_params(;N=1,freq_c=1)
    #TODO add N=0 case 
    N_coeffs = 2^(N-1)
    a = Array{Float64}(undef,N_coeffs)
    b = Array{Float64}(undef,N_coeffs)
    c = Array{Float64}(undef,N_coeffs)
    d = Array{Float64}(undef,N_coeffs)
    for i in 1:N_coeffs
        a[i] = (freq_c/2^N)*sin(pi/(2^(N+1))*(2*i-1))
        b[i] = (freq_c/2^N)*cos(pi/(2^(N+1))*(2*i-1))
        c[i] = freq_c*sin(pi/(2^(N+1))*(2*i-1))
        d[i] = freq_c*cos(pi/(2^(N+1))*(2*i-1))
    end
    forcing_params = (a=a,b=b,c=c,d=d)
    return forcing_params
end