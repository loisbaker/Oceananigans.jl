using JLD2
using JLD2: Group
using Oceananigans
using Oceananigans.Units: Time

function set_data_on_disk!(fields_filename; direction = "backward", T_start = nothing, T_end = nothing)

    # First check that a valid direction has been given
    if (direction !== "backward") && (direction !== "forward")
        error("Invalid direction: $direction")
    end
    
    # Open the file
    jldopen(fields_filename,"r+") do file

        # Checking if this is the first time we've edited the file
        if !haskey(file, "direction")
            # It should currently be forward, because direction hasn't been assigned yet
            file["direction"] = "forward"

            # This is also the first time we've edited the file, so we'll create some data to store original simulation times
            iterations = parse.(Int, keys(file["timeseries/t"]))
            times = [file["timeseries/t/$iter"] for iter in iterations]
            g = Group(file, "timeseries/t_simulation")
            for (i,iter) in enumerate(iterations)
                g["$iter"] = times[i]
            end
        end    

        # Now load in current direction and simulation times
        current_direction = file["direction"]
        println("Current direction is $current_direction")
        iterations = parse.(Int, keys(file["timeseries/t"]))
        t_simulation = [file["timeseries/t_simulation/$iter"] for iter in iterations]
        Nt = length(t_simulation)

        # If already the right direction, we just need to make sure the times correspond to the filter period we want
        if current_direction == direction
            println("No need to reverse order of data")
            # We might need to update times
            if (current_direction == "forward") 

                if isnothing(T_start) # If we're not given a starting time, set it to simulation start time 
                    T_start = t_simulation[1]
                    
                end

                if isnothing(T_end) 
                    T_end = t_simulation[Nt] # If we're not given a end time, set it to simulation end time 
                    
                end

                # Make a new time coordinate that is zero at T_start
                Base.delete!(file, "timeseries/t")
                g = Group(file, "timeseries/t")
                for (i,iter) in enumerate((iterations))
                    g["$iter"] = t_simulation[i] - T_start
                end
                
            else # current direction is backward

                if isnothing(T_start) # If we're not given a starting time, set it to simulation start time 
                    T_start = t_simulation[Nt]
                    
                end

                if isnothing(T_end) 
                    T_end = t_simulation[1] # If we're not given a end time, set it to simulation end time 
                    
                end

                # Make a new time coordinate that is zero at T_end
                Base.delete!(file, "timeseries/t")
                g = Group(file, "timeseries/t")
                for (i,iter) in enumerate(iterations)
                    g["$iter"] = T_end - t_simulation[i]
                end 
            end

        else # current direction is wrong, so we need to reverse everything and also update times
            Base.delete!(file, "direction")
            file["direction"] = direction
            println("Reversing order of data")
            for var in keys(file["timeseries"])
                if (var == "t_simulation") || (var == "t") # Then no serialized entry
                    backward_iterations = reverse(keys(file["timeseries/$var"])) 
                else
                    backward_iterations = reverse(keys(file["timeseries/$var"])[2:end]) # don't include serialized
                end
                # for velocities, we reverse all entries and negate them
                if (var == "u") || (var == "v")
                    for iter in backward_iterations 
                        data = file["timeseries/$var/$iter"] # Load in data
                        Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                        file["timeseries/$var/$iter"] = -data # Write it again, but negative                  
                    end
                    println("Reversed order of and switched sign of $var")

                # For time, we need to reverse the data and also make sure its correctly aligned with start and end points
                elseif var == "t"
                    if current_direction == "forward"
                        if isnothing(T_start) # If we're not given a starting time, set it to simulation start time 
                            T_start = t_simulation[1]
                            
                        end
        
                        if isnothing(T_end) 
                            T_end = t_simulation[Nt] # If we're not given a end time, set it to simulation end time 
                            
                        end

                        for (i, iter) in enumerate(backward_iterations)
                            Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                            file["timeseries/$var/$iter"] = T_end - t_simulation[Nt+1-i]
                        end

                    else # current_direction is backward

                        if isnothing(T_start) # If we're not given a starting time, set it to simulation start time 
                            T_start = t_simulation[Nt]
                            println("G new T_start is $T_start")
                        end
        
                        if isnothing(T_end) 
                            T_end = t_simulation[1] # If we're not given a end time, set it to simulation end time 
                            println("H new T_end is $T_end")
                        end

                        for (i, iter) in enumerate(backward_iterations)
                            Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                            file["timeseries/$var/$iter"] = t_simulation[Nt+1-i] - T_start
                        end
                    end
                    println("Reversed order of and shifted $var")
                # if var is t_simulation, the data doesn't change
                elseif var == "t_simulation"
                    for iter in backward_iterations 
                        data = file["timeseries/$var/$iter"] # Load in data
                        Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                        file["timeseries/$var/$iter"] = data # Write it again 
                    end
                    println("Reversed order of $var")
                # if var is a scalar, we need to make sure that it's on the correct (center) grid   
                else
                    location = file["timeseries/$var/serialized/location"]
                    if location !== (Center, Center, Center)
                        @warn "$var is on grid $location, but should be output on (Center, Center, Center). Continuing anyway."
                    end 
                    for iter in backward_iterations 
                        data = file["timeseries/$var/$iter"] # Load in data
                        Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                        file["timeseries/$var/$iter"] = data # Write it again 
                    end
                    println("Reversed order of $var")
                end
            end
            println("New direction is $direction")
        end
        T_filter = T_end - T_start # Update total filter time
        return T_filter 
    end
      
end


function set_forcing_params(;N=1,freq_c=1) # Try not having arrays as params
    N_coeffs = 2^(N-1)
    forcing_params = NamedTuple()
    for i in 1:N_coeffs
        
        a = (freq_c/2^N)*sin(pi/(2^(N+1))*(2*i-1))
        b = (freq_c/2^N)*cos(pi/(2^(N+1))*(2*i-1))
        c = freq_c*sin(pi/(2^(N+1))*(2*i-1))
        d = freq_c*cos(pi/(2^(N+1))*(2*i-1))

        temp_params = NamedTuple{(Symbol("a$i"), Symbol("b$i"),Symbol("c$i"),Symbol("d$i"))}([a,b,c,d])
        forcing_params = merge(forcing_params,temp_params)
    end
    
    return forcing_params
end

function create_tracers(forcing_params)
    
    N_coeffs = Int(length(forcing_params)/4)
    gC = ()  # Start with an empty tuple
    gS = ()  # Start with an empty tuple
    for i in 1:N_coeffs
        new_gC = Symbol("gC", i) 
        new_gS = Symbol("gS", i) 
        gC = (gC..., new_gC)  
        gS = (gS..., new_gS)  
    end
    return (gC...,gS...)
end

function make_gC_forcing_func(i)
    c_index = (i-1)*4 + 3
    d_index = (i-1)*4 + 4
    # Return a new function 
    return (x,y,t,gC,gS,p) -> -p[c_index]*gC - p[d_index]*gS
end

function make_gS_forcing_func(i)
    c_index = (i-1)*4 + 3
    d_index = (i-1)*4 + 4
    # Return a new function 
    return (x,y,t,gC,gS,p) -> -p[c_index]*gS + p[d_index]*gC  
end

function create_forcing(scalar_fts, forcing_params)

    N_coeffs = Int(length(forcing_params)/4)
    scalar_forcing = scalar_fts

    # Initialize dictionary
    gCdict = Dict()
    gSdict = Dict()

    for i in 1:N_coeffs
        
        gCkey = Symbol("gC$i")   # Dynamically create a Symbol for the key
        gSkey = Symbol("gS$i")   # Dynamically create a Symbol for the key

        # Store in dictionary
        gC_forcing_i = Forcing(make_gC_forcing_func(i), parameters =forcing_params, field_dependencies = (gCkey,gSkey))
        gCdict[gCkey] = (scalar_forcing,gC_forcing_i)

        gS_forcing_i = Forcing(make_gS_forcing_func(i), parameters =forcing_params, field_dependencies = (gCkey,gSkey))
        gSdict[gSkey] = gS_forcing_i
        
    end

    forcing = (; NamedTuple(gSdict)..., NamedTuple(gCdict)...)

    return forcing
end

function sum_tracers(model, forcing_params)
    N_coeffs = Int(length(forcing_params)/4)
    g_total = (forcing_params[1]*model.tracers[1] + forcing_params[2]*model.tracers[N_coeffs+1])
    
    for i in 2:N_coeffs
        a_index = (i-1)*4 + 1
        b_index = (i-1)*4 + 2
        gC_index = i
        gS_index = N_coeffs + i
        g_total = g_total + (forcing_params[a_index]*model.tracers[gC_index] + forcing_params[b_index]*model.tracers[gS_index])
    end
    return g_total 
end

function update_velocities!(sim, fts_velocities)
    model = sim.model
    time = sim.model.clock.time
    u_fts, v_fts = fts_velocities
    set!(model, u=u_fts[Time(time)], v= v_fts[Time(time)])
    return nothing
end
