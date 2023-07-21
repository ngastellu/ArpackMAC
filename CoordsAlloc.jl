module CoordsAlloc

include("./CoordsIO.jl")
using .CoordsIO, Base.Filesystem

xsf_file = expanduser("~/Desktop/simulation_outputs/MAC_MD_lammps/40x40/frame_10000_40K.xsf")
#smolfile = expanduser("~/Desktop/simulation_outputs/MAC_MD_lammps/10x10/trajectories/initplanar-276_200K.lammpstrj")
smolfile = "tiny_traj.lammpstrj"
bigfile = expanduser("~/Desktop/simulation_outputs/MAC_MD_lammps/40x40/40K_norotate_0-100000-500.lammpstrj")

frame_step = 500 #check name of bigfile to make sure this matches
frame_ind = 10000 #check name of smolfile to make sure this matches

println("Try 1")
@time pos1 = read_xsf(xsf_file; dump=true)
print('\n')
println("Try 2")
@time pos2 = get_frame(smolfile, 0) #run on smolfile 1st to limit memory allocations during compilation

print('\n')
println("Try 3")
@time pos3 = get_frame_bash(smolfile,0) #run on smolfile 1st to limit memory allocations during compilation

print('\n')
println("--------------------")
print('\n')

println("Try 1")
@time pos1 = read_xsf(xsf_file; dump=true)
print('\n')
println("Try 2")
@time pos2 = get_frame(bigfile, Int(frame_ind/frame_step))
print('\n')
println("Try 3")
@time pos3 = get_frame_bash(bigfile,frame_ind;frame_step=frame_step)

println(pos1 == pos2)
println(pos1 == pos3)

end