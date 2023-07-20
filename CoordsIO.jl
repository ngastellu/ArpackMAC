module CoordsIO

export read_xsf, get_frame, get_frame_bash

function read_xsf(filename; read_forces=false, dump=false)
    
    f = open(filename)
    if dump
        for i=1:3
            readline(f)
        end
        
        na = parse(Int, split(readline(f))[1])
        
        for i=1:5
            readline(f)
        end
    
    else
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

    end

    atoms = zeros(Float64, (na, 3))

    if read_forces
        forces = zeros(Float64, (na, 3))
    end
    
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

    if dump
        return atoms
    else
        if forces_in_file && read_forces
            return atoms, forces, supercell
        else
            return atoms, supercell
        end
    end
end


function get_frame(filename, frame_index)

    nb_non_coord_lines::Int = 9
    fo = open(filename)
    
    for i=1:3
        readline(fo)
    end

    Natoms = parse(Int, readline(fo))

    pos = zeros(Float64, (Natoms,3))

    nlines_per_frame = Natoms + nb_non_coord_lines

    seekstart(fo)

    for i=0:frame_index-1
        for j=1:nlines_per_frame
            readline(fo)
        end
    end

    for i=1:nb_non_coord_lines
        l = readline(fo)
        if i <= 2
            println(l)
        end
    end

    for k in 1:Natoms
        split_line = split(readline(fo))
        x, y, z = split_line[2:4]
        pos[k, :] = [parse(Float64, x), parse(Float64, y), parse(Float64, z)]
    end

    close(fo)

    return pos

end

function get_Natoms_dump(filename)
    fo = open(filename)
    nchars = 39 #number of chars to read before reaching line containing the number of atoms (assuming file is a dump file ⟺ line 1 = 'ITEM: TIMESTEP\n')
    seek(fo,nchars)
    Natoms = parse(Int, readline(fo))
    close(fo)
    return Natoms
end

function get_frame_bash(filename, frame_index; frame_step=1)
    
    nb_non_coord_lines::Int = 9
    Natoms = get_Natoms_dump(filename)
    nlines_per_frame = Natoms + nb_non_coord_lines
    nlines_tail = nlines_per_frame * (Int(frame_index/frame_step) + 1)

    pos = zeros(Float64, (Natoms, 3))

    allpos = split(read(pipeline(`head -n $nlines_tail $filename`, `tail -n $nlines_per_frame`), String), '\n')
    stepnb = parse(Int, allpos[2])
    println("Step number: ", stepnb)

    for i=(nb_non_coord_lines+1):nlines_per_frame
        x, y, z = split(allpos[i])[2:4]
        pos[i-nb_non_coord_lines,:] = [parse(Float64, x), parse(Float64, y), parse(Float64, z)]
    end

    return pos
    
end
    
end