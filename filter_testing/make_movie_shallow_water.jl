# Now try to load back in the data

using Oceananigans
using Oceananigans.Grids
using NCDatasets, Printf, CairoMakie

filename_stem = "SW_vort"
fields_filename = joinpath(@__DIR__, filename_stem * ".nc")
ds = NCDataset(fields_filename, "r")

times = ds["time"][:]
# Visualize the results
# Load required packages to read output and plot.

# # Define the coordinates for plotting.
xF, yF = ds["xF"], ds["yF"]
xC, yC = ds["xC"], ds["yC"]
#xu, yu = xnodes(grid, Face()), ynodes(grid, Center())

# Read in the `output_writer` for the two-dimensional fields and then create an animation 
# showing both the perturbation vorticities.

fig = Figure(size = (550, 500))

axis_kwargs = (xlabel = "x", ylabel = "y")
ax_ω  = Axis(fig[2, 1]; title = "Total relative vorticity, ω", axis_kwargs...)

n = Observable(1)

ω = @lift ds["ω"][:, :, 1, $n]
hm_ω = heatmap!(ax_ω, xF, yF, ω, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[2, 2], hm_ω)

title = @lift @sprintf("t = %.1f", times[$n])
fig[1, 1:2] = Label(fig, title, fontsize=24, tellwidth=false)

current_figure() #hide
fig

# Finally, we record a movie.

frames = 1:length(times)
movie_filename = joinpath(@__DIR__,  "SW_vorticity.mp4")
record(fig, movie_filename, frames, framerate=5) do i
    n[] = i
end
nothing #hide
close(ds)

