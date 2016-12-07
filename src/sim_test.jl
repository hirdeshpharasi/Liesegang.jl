#####################           TEST            ################################
#This program is a test code for the Liesegang package, it is a simple model based
#on the Stochastic Rotation Dynamics model.
#loading the package Liesegang.jl
using Liesegang
using Plots #plotting package
################################################################################
#defining the parameters
Lx = 10; Ly = 10 #size of the space
dim = [Lx,Ly]
a = 1.0 #size of the boxes, default = 1
mass = 1.0
np = 1000 #number of particles
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
anim = @animate for _ in 1:tmax
    box_vel(parts,boxes)
    #moving the particles
    α = rand(0.0:180.0)*rand([-1,1]) #random angle
    #rotating...
    parts_vels!(parts, boxes, α, a, dim)
    #now getting the new positions of the particles
    getpos_pbc!(parts, τ, dim)
    x = [parts[i].pos[1] for i in 1:np]
    y = [parts[i].pos[2] for i in 1:np]
    scatter(x,y)
end
gif(anim, "test.gif", fps = 10)
