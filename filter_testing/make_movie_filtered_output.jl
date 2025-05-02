# Now try to load back in the data

using Oceananigans
using Oceananigans.Grids
using NCDatasets, Printf, CairoMakie

filename_stem = "combined_LF"

fields_filename_nc = joinpath(@__DIR__, filename_stem * ".nc")

ds = NCDataset(fields_filename_nc, "r")


times = ds["time_fwd"][:]

# Visualize the results
# Load required packages to read output and plot.

# # Define the coordinates for plotting.
xF, yF = ds["xF"], ds["yF"]
xC, yC = ds["xC"], ds["yC"]
#xu, yu = xnodes(grid, Face()), ynodes(grid, Center())

# Read in the `output_writer` for the two-dimensional fields and then create an animation 
# showing both the perturbation vorticities.

fig = Figure(size = (1200, 450))

axis_kwargs = (xlabel = "x", ylabel = "y")
ax_1  = Axis(fig[2, 1]; title = "Total vorticity, ω", axis_kwargs...)
ax_2 = Axis(fig[2, 3]; title = "Total filtered vorticity, no regridding", axis_kwargs...)
ax_3  = Axis(fig[2, 5]; title = "Total filtered vorticity, regridded", axis_kwargs...)

n = Observable(1)

# f1 = @lift ds["ω_filtered"][:, :, 1, $n]
f1 = @lift ds["ω"][:, :, $n]
hm_ω = heatmap!(ax_1, xF, yF, f1, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[2, 2], hm_ω)

f2 = @lift ds["ω_filtered"][:, :, $n]
# f2 = @lift ds["ω_filtered_regrid"][:, :, 1, $n]

hm_c = heatmap!(ax_2, xC, yC, f2, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[2, 4], hm_c)

# f3 = @lift ds["u_filtered_regrid"][:, :, 1, $n]
f3 = @lift ds["ω_filtered_regrid"][:, :, $n]

hm_ω = heatmap!(ax_3, xF, yF, f3, colorrange = (-1, 1), colormap = :balance)
Colorbar(fig[2, 6], hm_ω)

# g2 = @lift ds2["g_total"][:, :, 1, $n]
# hm_c = heatmap!(ax_c2, xC, yC, g2, colorrange = (-1, 1), colormap = :balance)
# Colorbar(fig[3, 4], hm_c)


title = @lift @sprintf("t = %.1f", times[$n])
fig[1, 1:6] = Label(fig, title, fontsize=24, tellwidth=false)

current_figure() #hide
fig

# Finally, we record a movie.

frames = 1:length(times)
movie_filename = joinpath(@__DIR__,  "filtered_output.mp4")
record(fig, movie_filename, frames, framerate=5) do i
    n[] = i
end
nothing #hide
close(ds)

