# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

struct Truss <: FieldProblem end

FEMBase.get_unknown_field_name(::Problem{Truss}) = "displacement"

function FEMBase.assemble_elements!(problem::Problem{Truss}, assembly::Assembly,
                                    elements::Vector{Element{Seg2}}, time::Float64)

    Ke = zeros(6, 6)
    Me = zeros(6, 6)
    fe = zeros(6)

    for element in elements
        fill!(Ke, 0.0)
        X1, X2 = element("geometry", time)
        v_e = X2 - X1
        L = norm(v_e)

        c = v_e ./ L;

        E = element("youngs modulus", time)
        A = element("cross-section area", time)

        # Assemble stiffness matrix

        for i in 1:3, j in 1:3
            Ke[i,j] = c[i]*c[j]
            Ke[i+3,j] = -c[i]*c[j]
            Ke[i,j+3] = -c[i]*c[j]
            Ke[i+3,j+3] = c[i]*c[j]
        end

        Ke .*= E*A/L

        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
        add!(assembly.M, gdofs, gdofs, Me)
        add!(assembly.f, gdofs, fe)
    end

    return nothing

end

function FEMBase.assemble_elements!(problem::Problem{Truss}, assembly::Assembly,
                                    elements::Vector{Element{Poi1}}, time::Float64)

    for element in elements
        gdofs = get_gdofs(problem, element)
        for i = 1:3
            if haskey(element, "point force $i")
                P = element("point force $i", time)
                add!(assembly.f, gdofs[i], P)
            end
            if haskey(element, "fixed displacement $i")
                g = element("fixed displacement $i", time)
                add!(assembly.C1, gdofs[i], gdofs[i], 1.0)
                add!(assembly.C2, gdofs[i], gdofs[i], 1.0)
                add!(assembly.g, gdofs[i], g)
            end
        end
    end

    return nothing

end
