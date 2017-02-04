#####################           TEST            ################################
#This program is a test code for the Liesegang package, it is a simple model based
#on the Stochastic Rotation Dynamics model.
#Here the idea is to have two species of particles, A and B
#loading the package Liesegang.jl
using Liesegang
using Plots #plotting package
################################################################################
#defining the parameters
Lx = 50
Ly = 5 #size of the space
dim = [Lx,Ly]
a = 1.0 #size of the boxes, default = 1
m = [1.0, 2.0] #masses
np = [1000,1000] #number of particles
ntp = 3 #number of species.
Tr = 1/3 #reference temperature
τ = 1.73 #1.73
kr = 0.8 # probability of reaction or reaction rate.
tmax = 500
angles = [90.0, 90.0]
################################################################################
###########                       INITIALIZING                       ###########

parts = vcat([particle(1,dim, m[1],1) for _ in 1:np[1]], [particle(Lx,dim,m[1],2) for _ in 1:np[2]]) #initializing the particles.
#normalizing the momentum
norm_momentum!(parts)
#now the temperature to the reference Tr
norm_temperature!(parts, Tr)
#initializing the boxes
boxes = [box(ntp,i) for i in 1:(Lx * Ly)]

################################################################################
#########################    now the simulation...   ###########################

#anim = @animate
for t in 1:tmax
    #streaming step
    getpos_slip!(parts, τ, dim)
    #getpos_pbc!(parts,τ,dim)
    #quiver(x, y, quiver = (vx, vy), xlims =(0,Lx), ylims = (0,Ly))
    #first the grid is shifted
    shift_grid!(parts, a, dim)
    #now label the particles in the boxes
    get_box(parts,boxes, Lx)
    for (i, box) in enumerate(boxes) #cycling over the boxes
        parbox = filter(x-> x.indbox == i, parts) #selecting the particles that are in the box.
        #println(i,'\t',length(parbox))
        if isempty(parbox); continue; end  #ignore next steps if the box is empty
        collide_mc(parbox)
        if countnz(box.np) > 1 #this is if there is more than one type of particle in the box
            collide_sc(parbox, ntp)
            #println("antes",'\t',box.np)
            nr = prob_box(box.np, kr) #compute the probabilites of the transitions
            if nr != 0
                nc = reac_box(parbox, nr)
                println(nc)
                push!(parts, nc...) #this adds the new particles c to the array of particles.
            end
        end
    end
    #shifting back the particles to their original places
    shiftback_grid!(parts)
    #x = grap_pos(parts,1)
    #y = grap_pos(parts,2)
    #z = grap_pos(parts,3)
    #vx = [parts[i].vel[1]/3 for i in 1:np] #dividing the vectors by a factor of 3 just for the visualization.
    #vy = [parts[i].vel[2]/3 for i in 1:np]
    #scatter(x[:,1],x[:,2], xlims = (0,Lx), ylims = (0,Ly), size = (Lx*10,Ly*20))
    #scatter!(y[:,1],y[:,2])# xlims = (0,Lx), ylims = (0,Ly))
    #scatter!(z[:,1],z[:,2])
end

gif(anim, "testmulti$(ARGS[1]).gif", fps = 8)
