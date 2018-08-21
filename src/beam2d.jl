# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

type Beam2D <: FieldProblem
end

function FEMBase.get_unknown_field_name(::Problem{Beam2D})
    return "displacement"
end

function FEMBase.assemble_elements!{B}(::Problem{Beam2D}, ::Assembly,
                               elements::Vector{Element{B}}, ::Float64)

    for element in elements
        info("Not doing anything useful right now.")
    end

    return nothing
end