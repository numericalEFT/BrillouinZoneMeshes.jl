#using .PlotlyJS
#import .PlotlyJS: plot
using LinearAlgebra: norm, dot
using StaticArrays
using .WignerSeitz: face_normal, merge_coplanar!
# ---------------------------------------------------------------------------------------- #
# CONSTANTS

# default layout
const DEFAULT_PLOTLY_LAYOUT_3D = Layout(
    showlegend=false,
    scene=attr(
        xaxis=attr(tickvals=[], zeroline=false,
            showgrid=false, showbackground=false,
            title=attr(text=""),
        ),
        yaxis=attr(tickvals=[], zeroline=false,
            showgrid=false, showbackground=false,
            title=attr(text=""),
        ),
        zaxis=attr(tickvals=[], zeroline=false,
            showgrid=false, showbackground=false,
            title=attr(text=""),
        ),
        aspectmode="data",
        camera=attr(up=attr(x=0, z=1, y=0), center=attr(x=0, y=0, z=0)),
        dragmode="turntable",
    ),
    margin=attr(l=0, r=0, b=0, t=0),
    autosize=false,
    #width=200, height=200,
    plot_bgcolor=TRANSPARENT_COL[], paper_bgcolor=TRANSPARENT_COL[],
)

# ---------------------------------------------------------------------------------------- #
# 3D

function plot(c::Cell{3}, layout::Layout=Layout();
    config::PlotConfig=PlotConfig(responsive=true, displaylogo=false))

    layout = merge(DEFAULT_PLOTLY_LAYOUT_3D, layout)

    setting(c) !== CARTESIAN && (c = cartesianize(c))
    scale = maximum(norm, basis(c))

    # BZ
    tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(c))
    for (i, poly) in enumerate(c)
        tbz[i] = PlotlyJS.scatter3d(
            x=push!(getindex.(poly, 1), poly[1][1]),
            y=push!(getindex.(poly, 2), poly[1][2]),
            z=push!(getindex.(poly, 3), poly[1][3]);
            mode="lines", hovertext="Cell", hoverinfo="text+x+y+z",
            line=attr(color=BZ_COL[], width=3)
        )
    end

    # lattice vectors
    tgs = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 6)
    tgtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)
    intersects = axis_intersections(c, basis(c))
    for (i, V) in enumerate(basis(c))
        V′ = V ./ norm(V)
        name = "<b>v</b><sub>$(i)</sub>"
        for j in (1, 2) # in/outside BZ
            start = j == 1 ? 0.0 : intersects[i]
            stop = j == 1 ? intersects[i] : 1.0
            V₀ = start * V .+ (j == 1 ? 0.0 : 0.025) * scale * V′
            V₁ = stop * V .- (j == 1 ? 0.025 : 0.0) * scale * V′
            tgs[i+(j-1)*3] = PlotlyJS.scatter3d(
                x=[V₀[1], V₁[1]], y=[V₀[2], V₁[2]], z=[V₀[3], V₁[3]];
                mode="lines", hovertext=name, hoverinfo="text",
                line=attr(color=ifelse(j == 1, BASIS_LIGHT_COL[], BASIS_COL[]),
                    width=ifelse(j == 1, 5, 6))
            )
        end
        tgtips[i] = PlotlyJS.cone(
            x=[V[1]], u=[V′[1]], y=[V[2]], v=[V′[2]], z=[V[3]], w=[V′[3]],
            sizeref=0.1 * scale, showscale=false, anchor="tail",
            colorscale=[[0, BASIS_COL[]], [1, BASIS_COL[]]],
            hovertext=name, hoverinfo="text+x+y+z")
    end

    # Cartesian axes
    cart_basis = cartesian_axes(Val(3))
    name_axs = ("x", "y", "z")
    taxs = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 6)
    taxtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)
    intersects = axis_intersections(c, cart_basis)
    for (i, V) in enumerate(cart_basis)
        name = "<b>" * string(name_axs[i]) * "</b>"
        V′ = V ./ norm(V)
        for j in (1, 2) # in/outside BZ
            start = j == 1 ? 0.0 : intersects[i]
            stop = j == 1 ? intersects[i] : intersects[i]
            V₀ = start * V .+ (j == 1 ? 0.0 : 0.025) * scale * V′
            V₁ = stop * V .- (j == 1 ? 0.025 : -0.2) * scale * V′
            taxs[i+(j-1)*3] = PlotlyJS.scatter3d(
                x=[V₀[1], V₁[1]], y=[V₀[2], V₁[2]], z=[V₀[3], V₁[3]];
                mode="lines", hovertext=name, hoverinfo="text",
                line=attr(color=ifelse(j == 1, AXIS_LIGHT_COL[], AXIS_COL[]),
                    width=ifelse(j == 1, 5, 6))
            )
            if j == 2
                taxtips[i] = PlotlyJS.cone(
                    x=[V₁[1]], u=[V[1]],
                    y=[V₁[2]], v=[V[2]],
                    z=[V₁[3]], w=[V[3]],
                    sizeref=0.1 * scale, showscale=false, anchor="tail",
                    colorscale=[[0, AXIS_COL[]], [1, AXIS_COL[]]],
                    hovertext=name, hoverinfo="text")
            end
        end
    end

    # combine traces and plot
    ts = vcat(tbz, tgs, tgtips, taxs, taxtips)
    return PlotlyJS.plot(ts, layout; config=config)
end
 
function plot(clist::Array{Cell{3}}, idx_center::Int, layout::Layout=Layout(); ibz=nothing,
    config::PlotConfig=PlotConfig(responsive=true, displaylogo=false))

    layout = merge(DEFAULT_PLOTLY_LAYOUT_3D, layout)

    # scale = maximum(norm, basis(c))
    c_center = clist[idx_center]
    setting(c_center) !== CARTESIAN && (c_center = cartesianize(c_center))
    scale = maximum(norm, basis(c_center))
    tbz_list = []
    tgs = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 6)
    tgtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)
    taxs = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 6)
    taxtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)

    if !isnothing(ibz)
        facecolor = repeat([
    	      "rgb(50, 200, 200)",
    	      "rgb(100, 200, 255)",
    	      "rgb(150, 200, 115)",
    	      "rgb(200, 200, 50)",
    	      "rgb(230, 200, 10)",
    	      "rgb(255, 140, 0)"
        ], inner=[2])
        setting(ibz) !== CARTESIAN && (ibz = cartesianize(ibz))        
        faces = PlotlyJS.mesh3d(x=[ibz.verts[i][1] for i in 1:length(ibz.verts)],
                       y=[ibz.verts[i][2] for i in 1:length(ibz.verts)],
                       z=[ibz.verts[i][3] for i in 1:length(ibz.verts)],
                       i=[ibz.faces[i][1]-1 for i in 1:length(ibz.faces)],
                       j=[ibz.faces[i][2]-1 for i in 1:length(ibz.faces)],
                       k=[ibz.faces[i][3]-1 for i in 1:length(ibz.faces)],
                       facecolor = facecolor, opacity=0.7                               
                       )
        merge_coplanar!(ibz) # Coplanar triangles have to be merged after calling mesh3d.
        # BZ
        tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(ibz)+1)
        for (i, poly) in enumerate(ibz)
            tbz[i] = PlotlyJS.scatter3d(
                x=push!(getindex.(poly, 1), poly[1][1]),
                y=push!(getindex.(poly, 2), poly[1][2]),
                z=push!(getindex.(poly, 3), poly[1][3]);
                mode="lines", hovertext="Cell", hoverinfo="text+x+y+z",
                line=attr(color=BZ_COL[], width=3))
        end
        tbz[length(ibz)+1] = faces
        push!(tbz_list, tbz)
    end
    for (ci, c) in enumerate(clist)
        if (!isnothing(c))
            setting(c) !== CARTESIAN && (c = cartesianize(c))
            # BZ
            tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(c))
            for (i, poly) in enumerate(c)
                tbz[i] = PlotlyJS.scatter3d(
                    x=push!(getindex.(poly, 1), poly[1][1]),
                    y=push!(getindex.(poly, 2), poly[1][2]),
                    z=push!(getindex.(poly, 3), poly[1][3]);
                    mode="lines", hovertext="Cell", hoverinfo="text+x+y+z",
                    line=attr(color=BZ_COL[], width=3)
                )
            end
            push!(tbz_list, tbz)
            if ci == idx_center
                # lattice vectors
                intersects = axis_intersections(c, basis(c))
                for (i, V) in enumerate(basis(c))
                    V′ = V ./ norm(V)
                    name = "<b>v</b><sub>$(i)</sub>"
                    for j in (1, 2) # in/outside BZ
                        start = j == 1 ? 0.0 : intersects[i]
                        stop = j == 1 ? intersects[i] : 1.0
                        V₀ = start * V .+ (j == 1 ? 0.0 : 0.025) * scale * V′
                        V₁ = stop * V .- (j == 1 ? 0.025 : 0.0) * scale * V′
                        tgs[i+(j-1)*3] = PlotlyJS.scatter3d(
                            x=[V₀[1], V₁[1]], y=[V₀[2], V₁[2]], z=[V₀[3], V₁[3]];
                            mode="lines", hovertext=name, hoverinfo="text",
                            line=attr(color=ifelse(j == 1, BASIS_LIGHT_COL[], BASIS_COL[]),
                                width=ifelse(j == 1, 5, 6))
                        )
                    end
                    tgtips[i] = PlotlyJS.cone(
                        x=[V[1]], u=[V′[1]], y=[V[2]], v=[V′[2]], z=[V[3]], w=[V′[3]],
                        sizeref=0.1 * scale, showscale=false, anchor="tail",
                        colorscale=[[0, BASIS_COL[]], [1, BASIS_COL[]]],
                        hovertext=name, hoverinfo="text+x+y+z")
                end

                # Cartesian axes
                cart_basis = cartesian_axes(Val(3))
                name_axs = ("x", "y", "z")
                intersects = axis_intersections(c, cart_basis)
                for (i, V) in enumerate(cart_basis)
                    name = "<b>" * string(name_axs[i]) * "</b>"
                    V′ = V ./ norm(V)
                    for j in (1, 2) # in/outside BZ
                        start = j == 1 ? 0.0 : intersects[i]
                        stop = j == 1 ? intersects[i] : intersects[i]
                        V₀ = start * V .+ (j == 1 ? 0.0 : 0.025) * scale * V′
                        V₁ = stop * V .- (j == 1 ? 0.025 : -0.2) * scale * V′
                        taxs[i+(j-1)*3] = PlotlyJS.scatter3d(
                            x=[V₀[1], V₁[1]], y=[V₀[2], V₁[2]], z=[V₀[3], V₁[3]];
                            mode="lines", hovertext=name, hoverinfo="text",
                            line=attr(color=ifelse(j == 1, AXIS_LIGHT_COL[], AXIS_COL[]),
                                width=ifelse(j == 1, 5, 6))
                        )
                        if j == 2
                            taxtips[i] = PlotlyJS.cone(
                                x=[V₁[1]], u=[V[1]],
                                y=[V₁[2]], v=[V[2]],
                                z=[V₁[3]], w=[V[3]],
                                sizeref=0.1 * scale, showscale=false, anchor="tail",
                                colorscale=[[0, AXIS_COL[]], [1, AXIS_COL[]]],
                                hovertext=name, hoverinfo="text")
                        end
                    end
                end
            end
        end
    end

    # combine traces and plot
    ts = tbz_list[1]
    for i in 2:length(tbz_list)
        ts = vcat(ts, tbz_list[i])
    end
    ts = vcat(ts, tgs, tgtips, taxs, taxtips)
    # println(length(ts))
    return PlotlyJS.plot(ts, layout; config=config)
end



# ---------------------------------------------------------------------------------------- #
# 2D

# default layout
const DEFAULT_PLOTLY_LAYOUT_2D = Layout(
    showlegend=false,
    xaxis=attr(tickvals=[], zeroline=false,
        showgrid=false, showbackground=false,
        title=attr(text=""),
    ),
    yaxis=attr(tickvals=[], zeroline=false,
        showgrid=false, showbackground=false,
        title=attr(text=""),
        scaleanchor="x", scaleratio=1
    ),
    aspectmode="data",
    hovermode="closest",
    margin=attr(l=0, r=0, b=0, t=0),
    autosize=false,
    plot_bgcolor=TRANSPARENT_COL[], paper_bgcolor=TRANSPARENT_COL[],
    annotations=PlotlyBase.PlotlyAttribute[]
)

function plot(c::Cell{2}, layout::Layout=Layout();
    config::PlotConfig=PlotConfig(responsive=true, displaylogo=false))

    layout = merge(DEFAULT_PLOTLY_LAYOUT_2D, layout)

    setting(c) !== CARTESIAN && (c = cartesianize(c))

    scale = maximum(norm, basis(c))
    max_x, max_y = maximum(v -> abs(v[1]), basis(c)), maximum(v -> abs(v[2]), basis(c))
    get!(layout[:xaxis], :range, [-max_x - scale / 15, max_x + scale / 15])
    get!(layout[:yaxis], :range, [-max_y - scale / 15, max_y + scale / 15])

    # Cell boundaries
    tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(c))
    for (i, poly) in enumerate(c)
        tbz[i] = PlotlyJS.scatter(
            x=push!(getindex.(poly, 1), poly[1][1]),
            y=push!(getindex.(poly, 2), poly[1][2]);
            mode="lines", hovertext="Cell", hoverinfo="text+x+y",
            line=attr(color=BZ_COL[], width=3)
        )
    end

    # lattice vectors
    tgs = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 4)
    tgtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 2)
    intersects = axis_intersections(c, basis(c))
    for (i, V) in enumerate(basis(c))
        V′ = V ./ norm(V)
        name = "<b>v</b><sub>$(i)</sub>"
        for j in (1, 2) # in/outside BZ
            start = j == 1 ? 0.0 : intersects[i]
            stop = j == 1 ? intersects[i] : 1.0
            V₀ = start * V .+ (j == 1 ? 0.0 : 0.025) * scale * V′
            V₁ = stop * V .- (j == 1 ? 0.025 : 0.0) * scale * V′
            tgs[i+(j-1)*2] = PlotlyJS.scatter(
                x=[V₀[1], V₁[1]], y=[V₀[2], V₁[2]];
                mode="lines", hovertext=name, hoverinfo=ifelse(j == 1, "text", "text+x+y"),
                line=attr(color=ifelse(j == 1, BASIS_LIGHT_COL[], BASIS_COL[]),
                    width=ifelse(j == 1, 5, 6))
            )
        end
        # arrow heads have to be added as annotations to layout in 2D :/
        haskey(layout, :annotations) || (layout[:annotations] = PlotlyBase.PlotlyAttribute[])
        push!(layout[:annotations],
            attr(x=V[1] + 0.05V′[1] * scale, y=V[2] + 0.05V′[2] * scale,   # awful fidgeting; plotly's
                ax=V[1] - 0.05V′[1] * scale, ay=V[2] - 0.05V′[2] * scale, # arrows are stupid
                xref="ax", yref="ay", axref="x", ayref="y",
                showarrow=true, arrowhead=2, arrowwidth=6, arrowsize=0.5,
                arrowcolor=BASIS_COL[]))
    end

    # Cartesian axes
    cart_basis = cartesian_axes(Val(2))
    name_axs = ("x", "y")
    taxs = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 4)
    taxtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 2)
    intersects = axis_intersections(c, cart_basis)
    for (i, V) in enumerate(cart_basis)
        name = "<b>" * string(name_axs[i]) * "</b>"
        V′ = V ./ norm(V)
        for j in (1, 2) # in/outside BZ
            start = j == 1 ? 0.0 : intersects[i]
            stop = j == 1 ? intersects[i] : intersects[i]
            V₀ = start * V .+ (j == 1 ? 0.0 : 0.025) * scale * V′
            V₁ = stop * V .- (j == 1 ? 0.025 : -0.2) * scale * V′

            taxs[i+2(j-1)] = PlotlyJS.scatter(
                x=[V₀[1], V₁[1]], y=[V₀[2], V₁[2]];
                mode="lines", hovertext=name, hoverinfo="text",
                line=attr(color=ifelse(j == 1, AXIS_LIGHT_COL[], AXIS_COL[]),
                    width=ifelse(j == 1, 5, 6))
            )
            if j == 2
                # arrow heads have to be added as annotations to layout in 2D :/
                haskey(layout, :annotations) || (layout[:annotations] = PlotlyBase.PlotlyAttribute[])
                push!(layout[:annotations],
                    attr(x=V₁[1] + 0.05V′[1] * scale, y=V₁[2] + 0.05V′[2] * scale,   # awful fidgeting; plotly's
                        ax=V₁[1] - 0.05V′[1] * scale, ay=V₁[2] - 0.05V′[2] * scale, # arrows are stupid
                        xref="ax", yref="ay", axref="x", ayref="y",
                        showarrow=true, arrowhead=2, arrowwidth=6, arrowsize=0.5,
                        arrowcolor=AXIS_COL[]))
            end
        end
    end

    # combine traces and plot
    ts = vcat(tbz, tgs, taxs)
    return PlotlyJS.plot(ts, layout; config=config)
end

function plot(clist::Array{Cell{2}}, idx_center::Int, layout::Layout=Layout(); ibz =nothing,
    config::PlotConfig=PlotConfig(responsive=true, displaylogo=false))

    layout = merge(DEFAULT_PLOTLY_LAYOUT_2D, layout)
    tbz_list = []
    tgs = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 4)
    tgtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 2)
    taxs = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 4)
    taxtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 2)

    c_center = clist[idx_center]
    setting(c_center) !== CARTESIAN && (c_center = cartesianize(c_center))
    scale = maximum(norm, basis(c_center))
    max_x, max_y = maximum(v -> abs(v[1]), basis(c_center)), maximum(v -> abs(v[2]), basis(c_center))
    get!(layout[:xaxis], :range, [-max_x - scale / 15, max_x + scale / 15])
    get!(layout[:yaxis], :range, [-max_y - scale / 15, max_y + scale / 15])


    if !isnothing(ibz)
        setting(ibz) !== CARTESIAN && (ibz = cartesianize(ibz))        
        # merge_coplanar!(ibz) # Coplanar triangles have to be merged after calling mesh3d.
        # BZ
        merge_coplanar!(ibz)
        tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(ibz))
        for (i, poly) in enumerate(ibz)
            tbz[i] = PlotlyJS.scatter(
                x=push!(getindex.(poly, 1), poly[1][1]),
                y=push!(getindex.(poly, 2), poly[1][2]);
                mode="lines", hovertext="Cell", hoverinfo="text+x+y",
                line=attr(color=BZ_COL[], width=3 ), fill = "toself", fillcolor = "rgb(50, 200, 200)"
            )
        end
        push!(tbz_list, tbz)
    end

    for (ci, c) in enumerate(clist)
        if (!isnothing(c))
            # Cell boundaries
            setting(c) !== CARTESIAN && (c = cartesianize(c))
            tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(c))
            for (i, poly) in enumerate(c)
                tbz[i] = PlotlyJS.scatter(
                    x=push!(getindex.(poly, 1), poly[1][1]),
                    y=push!(getindex.(poly, 2), poly[1][2]);
                    mode="lines", hovertext="Cell", hoverinfo="text+x+y",
                    line=attr(color=BZ_COL[], width=3)
                )
            end
            push!(tbz_list, tbz)
            # lattice vectors
            if ci == idx_center

                intersects = axis_intersections(c, basis(c))
                for (i, V) in enumerate(basis(c))
                    V′ = V ./ norm(V)
                    name = "<b>v</b><sub>$(i)</sub>"
                    for j in (1, 2) # in/outside BZ
                        start = j == 1 ? 0.0 : intersects[i]
                        stop = j == 1 ? intersects[i] : 1.0
                        V₀ = start * V .+ (j == 1 ? 0.0 : 0.025) * scale * V′
                        V₁ = stop * V .- (j == 1 ? 0.025 : 0.0) * scale * V′
                        tgs[i+(j-1)*2] = PlotlyJS.scatter(
                            x=[V₀[1], V₁[1]], y=[V₀[2], V₁[2]];
                            mode="lines", hovertext=name, hoverinfo=ifelse(j == 1, "text", "text+x+y"),
                            line=attr(color=ifelse(j == 1, BASIS_LIGHT_COL[], BASIS_COL[]),
                                width=ifelse(j == 1, 5, 6))
                        )
                    end
                    # arrow heads have to be added as annotations to layout in 2D :/
                    haskey(layout, :annotations) || (layout[:annotations] = PlotlyBase.PlotlyAttribute[])
                    push!(layout[:annotations],
                        attr(x=V[1] + 0.05V′[1] * scale, y=V[2] + 0.05V′[2] * scale,   # awful fidgeting; plotly's
                            ax=V[1] - 0.05V′[1] * scale, ay=V[2] - 0.05V′[2] * scale, # arrows are stupid
                            xref="ax", yref="ay", axref="x", ayref="y",
                            showarrow=true, arrowhead=2, arrowwidth=6, arrowsize=0.5,
                            arrowcolor=BASIS_COL[]))
                end

                # Cartesian axes
                cart_basis = cartesian_axes(Val(2))
                name_axs = ("x", "y")
                intersects = axis_intersections(c, cart_basis)
                for (i, V) in enumerate(cart_basis)
                    name = "<b>" * string(name_axs[i]) * "</b>"
                    V′ = V ./ norm(V)
                    for j in (1, 2) # in/outside BZ
                        start = j == 1 ? 0.0 : intersects[i]
                        stop = j == 1 ? intersects[i] : intersects[i]
                        V₀ = start * V .+ (j == 1 ? 0.0 : 0.025) * scale * V′
                        V₁ = stop * V .- (j == 1 ? 0.025 : -0.2) * scale * V′

                        taxs[i+2(j-1)] = PlotlyJS.scatter(
                            x=[V₀[1], V₁[1]], y=[V₀[2], V₁[2]];
                            mode="lines", hovertext=name, hoverinfo="text",
                            line=attr(color=ifelse(j == 1, AXIS_LIGHT_COL[], AXIS_COL[]),
                                width=ifelse(j == 1, 5, 6))
                        )
                        if j == 2
                            # arrow heads have to be added as annotations to layout in 2D :/
                            haskey(layout, :annotations) || (layout[:annotations] = PlotlyBase.PlotlyAttribute[])
                            push!(layout[:annotations],
                                attr(x=V₁[1] + 0.05V′[1] * scale, y=V₁[2] + 0.05V′[2] * scale,   # awful fidgeting; plotly's
                                    ax=V₁[1] - 0.05V′[1] * scale, ay=V₁[2] - 0.05V′[2] * scale, # arrows are stupid
                                    xref="ax", yref="ay", axref="x", ayref="y",
                                    showarrow=true, arrowhead=2, arrowwidth=6, arrowsize=0.5,
                                    arrowcolor=AXIS_COL[]))
                        end
                    end
                end
                # combine traces and plot
                #ts = vcat(tbz, tgs, taxs)
            end
        end
    end
    #ts=vcat(tb for tb in tbz_list)
    ts = tbz_list[1]
    for i in 2:length(tbz_list)
        ts = vcat(ts, tbz_list[i])
    end
    ts = vcat(ts, tgs, taxs)
    #ts = tbz_list[1]
    return PlotlyJS.plot(ts, layout; config=config)
end




# ---------------------------------------------------------------------------------------- #
# UTILITIES 

# the "outward" intersections of of lines with direction `Vs` though origo with a cell face
function axis_intersections(c::Cell{3},
    Vs::AbstractVector{<:AbstractVector}=cartesian_axes(Val(3)))

    intersects = MVector{3,Float64}(undef)
    fill!(intersects, Inf)
    # intersection between plane (r-rₒ) ⋅ n = 0 and line r(t) = l₀ + lt:
    #   t = (r₀ - l₀) ⋅ n / (l ⋅ n)     [https://wikipedia.org/wiki/Line–plane_intersection]
    for (i, poly) in enumerate(c)
        r₀ = sum(poly) / length(poly)
        n = face_normal(c, i)
        for (j, V) in enumerate(Vs)
            t = dot(r₀, n) / dot(V, n)
            if t > 0.0 && t < intersects[j]
                intersects[j] = t
            end
        end
    end
    return intersects
end

# the "outward" intersections of of lines with direction `Vs` though origo with a cell segment
function axis_intersections(c::Cell{2},
    Vs::AbstractVector{<:AbstractVector}=cartesian_axes(Val(2)))

    intersects = MVector{2,Float64}(undef)
    fill!(intersects, Inf)
    # intersection between line rₐ(t) a + nₐt and line rᵦ(t) = β + nᵦt can be gotten by 
    # linear algebra: [-nₐ|nᵦ][tₐ,tᵦ] = a-β. Here we, pick rᵦ for the axes `Vs` (β=0) and
    # let rₐ refer to cell segments. We're then only interested in tᵦ:
    fs = faces(c)
    vs = vertices(c)
    for i in eachindex(vs) # number of vertices and segments are equal in 2D
        if length(fs) == 1 # `merge = true`
            idxs = i ≠ length(vs) ? (i, i + 1) : (i, 1)
        else
            idxs = (fs[i][1], fs[i][2])
        end
        a = vs[idxs[1]]
        nₐ = vs[idxs[2]] - a
        for (j, V) in enumerate(Vs)
            nᵦ = V
            tᵦ = (nₐ[2] * a[1] - nₐ[1] * a[2]) / (-nₐ[1] * nᵦ[2] + nₐ[2] * nᵦ[1]) # ~ inv([-nₐ|nᵦ])
            if tᵦ > 0.0 && tᵦ < intersects[j]
                intersects[j] = tᵦ
            end
        end
    end
    return intersects
end

# a set of cartesian basis vectors in D dimensions
cartesian_axes(Dᵛ::Val{D}) where {D} = SVector(ntuple(i -> SVector(ntuple(j -> i == j ? 1.0 : 0.0, Dᵛ)), Dᵛ))
