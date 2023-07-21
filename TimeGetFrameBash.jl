module TimeGetFrameBash

include("./CoordsIO.jl")
using .CoordsIO

smolfile = "tiny_traj.lammpstrj"
bigfile = expanduser("~/Desktop/simulation_outputs/MAC_MD_lammps/40x40/40K_norotate_0-100000-500.lammpstrj")
frame_step = 500 #check name of bigfile to make sure this matches
frame_ind = 10000 #check name of smolfile to make sure this matches

println("Try 1")
_ = get_frame_bash(smolfile,0)

println('\n')
println("---------------")
println('\n')

println("Try 2")
_ = get_frame_bash(bigfile,frame_ind;frame_step=frame_step)


end