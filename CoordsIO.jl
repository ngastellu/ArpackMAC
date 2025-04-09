module CoordsIO

export read_xsf, read_xyz, read_xyz_supercell, get_frame

function parse_coords(f::IOStream, na::Int)
    symbols = Vector{String}(undef, na)
    pos = zeros(Float64,(3,na))
    
    for k=1:na
        split_line = split(strip(readline(f)))
        symbols[k] = split_line[1]
        x, y, z = split_line[2:4]
        pos[:,k] = [parse(Float64, x), parse(Float64, y), parse(Float64, z)]
    end

    return pos, symbols
end

function read_xyz(filename::String)
    f = open(filename)
    na = parse(Int, strip(readline(f))) # first line contains number of atoms
    
    readline(f) # skip second line
    pos, symbols = parse_coords(f,na)
    
    return pos, symbols
end

function read_xyz_supercell(filename::String)
    f = open(filename)
    na = parse(Int, strip(readline(f))) # first line contains number of atoms

    line2 = readline(f) #second line contains cell info
        
    cell_vec_coords = split(split(split(line2, '=')[2], '"')[2]) # Expects ASE formatting
    supercell = parse.(Float64, cell_vec_coords[[1,5,9]]) # assume orthorhombic cell (all lattice vectors mutually orthogonal)

    supercell[ supercell .== 0 ] .= Inf # only apply PBC in directions where lattice dimensions are nonzero

    pos, symbols = parse_coords(f, na)

    return pos, supercell, symbols

end


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

        supercell = Array{Float64}(undef,2)
        supercell[1] = parse(Float64, split(readline(f))[1])
        supercell[2] = parse(Float64, split(readline(f))[2])

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
        if read_forces && length(split_line) > 4
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


function get_Natoms_dump(filename;frame0_index=0)
    fo = open(filename)
    if frame0_index == 0
        ndigits = 1
    else
        ndigits = Int(floor(log10(frame0_index))) + 1
    end
    nchars = 38 + ndigits #number of chars to read before reaching line containing the number of atoms (assuming file is a dump file ⟺ line 1 = 'ITEM: TIMESTEP\n')
    seek(fo,nchars)
    Natoms = parse(Int, readline(fo))
    close(fo)
    return Natoms
end

function get_frame(filename, frame_index; frame_step=1, frame0_index=0)
    
    nb_non_coord_lines::Int = 9
    Natoms = get_Natoms_dump(filename;frame0_index=frame0_index)
    nlines_per_frame = Natoms + nb_non_coord_lines
    nlines_tail = nlines_per_frame * (Int((frame_index - frame0_index)/frame_step) + 1)

    pos = zeros(Float64, (Natoms, 3))
    rawpos = read(pipeline(`head -n $nlines_tail $filename`, `tail -n $Natoms`), String) #this step is the bottleneck... might need to optimise
    allpos = split(rawpos, '\n')

    for i=1:Natoms
        x, y, z = split(allpos[i])[2:4]
        pos[i,:] = [parse(Float64, x), parse(Float64, y), parse(Float64, z)]
    end

    return pos
    
end
    
end
