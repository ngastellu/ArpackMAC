module CoordsAlloc

include("./CoordsIO.jl")
using .CoordsIO, Base.Filesystem

smolfile = expanduser("~/Desktop/simulation_outputs/MAC_MD_lammps/40x40/frame_10000_40K.xsf")
bigfile = expanduser("~/Desktop/simulation_outputs/MAC_MD_lammps/40x40/40K_norotate_0-100000-500.lammpstrj")

frame_step = 500 #check name of bigfile to make sure this matches
frame_ind = 10000 #check name of smolfile to make sure this matches

@time pos1 = read_xsf(smolfile; dump=true)
print('\n')
@time pos2 = get_frame(bigfile, Int(frame_ind/frame_step))
print('\n')
@time pos3 = get_frame_bash(bigfile,frame_ind;frame_step=frame_step)
println("--------------------")
@time pos1 = read_xsf(smolfile; dump=true)
print('\n')
@time pos2 = get_frame(bigfile, Int(frame_ind/frame_step))
print('\n')
@time pos3 = get_frame_bash(bigfile,frame_ind;frame_step=frame_step)

println(pos1 == pos2)
println(pos1 == pos3)

end