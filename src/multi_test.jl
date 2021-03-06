#####################           TEST            ################################
#This program is a test code for the Liesegang package, it is a simple model based
#on the Stochastic Rotation Dynamics model.
#Here the idea is to have two species of particles, A and B
#loading the package Liesegang.jl
using Liesegang
using Plots #plotting package
################################################################################
#defining the parameters
Lx = 100
Ly = 10 #size of the space
dim = [Lx,Ly]
a = 1.0 #size of the boxes, default = 1
m = [1.0, 2.0] #masses
np = [5000,5000] #number of particles
Tr = 1/3 #reference temperature
τ = 1.73 #1.73
tmax = 500
angles = [90.0, 90.0]
ntp = 2
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

anim = @animate for t in 1:tmax
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
        if countnz(box.np) > 1
            collide_sc(parbox, ntp)
        end
    end
    #shifting back the particles to their original places
    shiftback_grid!(parts)
    x = [parts[i].pos[1] for i in 1:np[1]]
    y = [parts[i].pos[2] for i in 1:np[1]]
    x1 = [parts[i].pos[1] for i in (np[1]+1):(np[2]+np[1])]
    y1 = [parts[i].pos[2] for i in (np[1]+1):(np[2]+np[1])]
    #vx = [parts[i].vel[1]/3 for i in 1:np] #dividing the vectors by a factor of 3 just for the visualization.
    #vy = [parts[i].vel[2]/3 for i in 1:np]
    scatter(x,y, xlims = (0,Lx), ylims = (0,Ly), size = (Lx*10,Ly*20))
    scatter!(x1,y1)# xlims = (0,Lx), ylims = (0,Ly))
end
gif(anim, "testmulti.gif", fps = 8)
#gif(anim, "testmulti$(ARGS[1]).gif", fps = 8)
