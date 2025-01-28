# Now try to load back in the data

using Oceananigans
using Oceananigans.Grids
using NCDatasets, Printf, CairoMakie

filename_stem1 = "forward_LF"
filename_stem2 = "backward_LF"
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

fig = Figure(size = (1200, 1000))

axis_kwargs = (xlabel = "x", ylabel = "y")
ax_ω1  = Axis(fig[2, 1]; title = "Total vorticity, ω, forward", axis_kwargs...)
ax_c1 = Axis(fig[2, 3]; title = "Filtered tracer, g_total forward", axis_kwargs...)
ax_ω2  = Axis(fig[3, 1]; title = "Total vorticity, ω, backward", axis_kwargs...)
ax_c2 = Axis(fig[3, 3]; title = "Filtered tracer, g_total backward", axis_kwargs...)

n = Observable(1)

ω1 = @lift ds1["ω"][:, :, 1, $n]
hm_ω = heatmap!(ax_ω1, xF, yF, ω1, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[2, 2], hm_ω)

g1 = @lift ds1["g_total"][:, :, 1, $n]
hm_c = heatmap!(ax_c1, xC, yC, g1, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[2, 4], hm_c)

ω2 = @lift ds2["ω"][:, :, 1, $n]
hm_ω = heatmap!(ax_ω2, xF, yF, ω2, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[3, 2], hm_ω)

g2 = @lift ds2["g_total"][:, :, 1, $n]
hm_c = heatmap!(ax_c2, xC, yC, g2, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[3, 4], hm_c)


title = @lift @sprintf("t = %.1f", times[$n])
fig[1, 1:4] = Label(fig, title, fontsize=24, tellwidth=false)

current_figure() #hide
fig

# Finally, we record a movie.

frames = 1:length(times)
movie_filename = joinpath(@__DIR__,  "check_filter.mp4")
record(fig, movie_filename, frames, framerate=5) do i
    n[] = i
end
nothing #hide
close(ds1)
close(ds2 )
