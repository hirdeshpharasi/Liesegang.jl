#####################           TEST            ################################
#This program is a test code for the Liesegang package, it is a simple model based
#on the Stochastic Rotation Dynamics model.
#loading the package Liesegang.jl
using Liesegang
using Plots #plotting package
################################################################################
#defining the parameters
Lx = 4; Ly = 4 #size of the space
dim = [Lx,Ly]
a = 1.0 #size of the boxes, default = 1
mass = 1.0
np = 60 #number of particles
Tr = 1/3 #reference temperature
τ = 1.73
tmax = 500
angles = [90.0, 130.0, 180.0]
################################################################################
###########                       INITIALIZING                       ###########

parts = [particle(dim, mass) for _ in 1:np] #initializing the particles.
#normalizing the momentum
norm_momentum!(parts)
#now the temperature to the reference Tr
norm_temperature!(parts, Tr)
#initializing the boxes
boxes = [box() for _ in 1:(Lx * Ly)]

################################################################################
#########################    now the simulation...   ###########################
#computing the momentum of each box
anim = @animate for t in 1:tmax
    #first the grid is shifted
    δ = shift_grid!(parts, a, dim)
    #now label the particles in the boxes
    get_box(parts, Lx)
    #the momentums and rotations are computed
    box_vel(parts,boxes)
    #if t < 100
    #    println([parts[i].vel' for i in 1:np])
    #end
    parts_vels!(parts, boxes, angles)
    #shifting back the particles to their original places
    shiftback_grid!(parts, dim, δ)
    #now getting the new positions of the particles
    getpos_pbc!(parts, τ, dim)
    x = [parts[i].pos[1] for i in 1:20]
    y = [parts[i].pos[2] for i in 1:20]
    x1 = [parts[i].pos[1] for i in 31:np]
    y1 = [parts[i].pos[2] for i in 31:np]
    scatter(x,y)
    scatter!(x1,y1)
end
gif(anim, "testshift.gif", fps = 10)
#x = [parts[i].pos[1] for i in 1:np]
#y = [parts[i].pos[2] for i in 1:np]
#writedlm("endpositions.dat", [x y])
