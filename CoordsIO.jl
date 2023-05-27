module CoordsIO

export read_xsf

function read_xsf(filename; read_forces=true)
    f = open(filename)
    for i in 1:2
        readline(f)
    end

    supercell = []
    push!(supercell, parse(Float64, split(readline(f))[1]))
    push!(supercell, parse(Float64, split(readline(f))[2]))

    for i in 1:2
        readline(f)
    end

    na = parse(Int, split(readline(f))[1])
    atoms = zeros(Float64, (na, 3))
    forces = zeros(Float64, (na, 3))
    forces_in_file = false

    for k in 1:na
        split_line = split(readline(f))
        x, y, z = split_line[2:4]
        atoms[k, :] = [parse(Float64, x), parse(Float64, y), parse(Float64, z)]
        if length(split_line) > 4
            forces_in_file = true
            fx, fy, fz = split_line[5:end]
            forces[k, :] = [parse(Float64, fx), parse(Float64, fy), parse(Float64, fz)]
        end
    end

    close(f)

    if forces_in_file && read_forces
        return atoms, forces, supercell
    else
        return atoms, supercell
    end
end

end