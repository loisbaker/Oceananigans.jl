# Now try to load back in the data

using Oceananigans
using Oceananigans.Grids
using NCDatasets, Printf, CairoMakie

filename_stem1 = "forward_LF"
filename_stem2 = "forward_LF"
fields_filename_nc1 = joinpath(@__DIR__, filename_stem1 * ".nc")
fields_filename_nc2 = joinpath(@__DIR__, filename_stem2 * ".nc")
ds1 = NCDataset(fields_filename_nc1, "r")
ds2 = NCDataset(fields_filename_nc2, "r")

times = ds1["time"][:]
times_check = ds2["time"][:]
print("$(length(times)) and $(length(times_check))" )
# Visualize the results
# Load required packages to read output and plot.

# # Define the coordinates for plotting.
xF, yF = ds1["xF"], ds1["yF"]
xC, yC = ds1["xC"], ds1["yC"]
#xu, yu = xnodes(grid, Face()), ynodes(grid, Center())

# Read in the `output_writer` for the two-dimensional fields and then create an animation 
# showing both the perturbation vorticities.

fig = Figure(size = (1200, 750))

axis_kwargs = (xlabel = "x", ylabel = "y")
ax_ω  = Axis(fig[2, 1]; title = "Total vorticity, ω", axis_kwargs...)
ax_c = Axis(fig[2, 3]; title = "Tracer, gC", axis_kwargs...)

n = Observable(1)

ω = @lift ds1["ω"][:, :, 1, $n]
hm_ω = heatmap!(ax_ω, xF, yF, ω, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[2, 2], hm_ω)

# g = @lift (ds1["g_total"][:, :, 1, $n] + ds2["g_total"][:, :, 1, $n])
g = @lift ds1["g"][:, :, 1, $n]
hm_c = heatmap!(ax_c, xC, yC, g, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[2, 4], hm_c)

title = @lift @sprintf("t = %.1f", times[$n])
fig[1, 1:4] = Label(fig, title, fontsize=24, tellwidth=false)

current_figure() #hide
fig

# Finally, we record a movie.

frames = 1:length(times)

record(fig, "check_filter.mp4", frames, framerate=5) do i
    n[] = i
end
nothing #hide
close(ds1)
close(ds2 )
