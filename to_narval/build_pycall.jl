module BuildPyCall

using PyCall, Pkg

# The purpose of this script is to tell PyCall to use the version
# of Python and of the desired packages used in the currently loaded
# virtualenv, inside of a job which runs on the compute nodes. I had 
# to write this because Narval's default NumPy version is no longer 
# compatible with SciPy version seen by PyCall.
#
# Usage guide:
# 1. Set up the virtualenv inside of the job (i.e. create the venv, activate
# 	it and install the desired Python packages)
#
# 2. Run this script: i.e. "julia build_pycall.jl"
#
# 3. Run the Julia script of interest.
#
# Somewhere in this repo, you should find a submission script that does all of this.

ENV["PYTHON"] = joinpath(ENV["VIRTUAL_ENV"], "bin", "python")
Pkg.build("PyCall")

end
