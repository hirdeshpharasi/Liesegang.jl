#####################           TEST            ################################
#This program is a test code for the Liesegang package, it is a simple model based
#on the Stochastic Rotation Dynamics model.
#loading the package Liesegang.jl
using Liesegang
################################################################################
#defining the parameters
Lx = 10; Ly = 10 #size of the space
dim = [Lx,Ly]
mass = 1.0
np = 100 #number of particles
Tr = 1.0 #reference temperature
τ = 0.01
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
box_vel(parts,boxes)
#moving the particles
α = 130.0 #fixing the angle just for the test
#rotating...
parts_vels!(parts, boxes, α)
#now getting the new positions of the particles
getpos_pbc!(parts, τ, dim)
