################################################################################
##################### MODULES FOR MULTI PARTICLE COLLISIONS ####################
################################################################################
#This is to get the number of box where the particle is.
function get_box(parts::Array{particle,1}, boxes::Array{box,1}, Lx::Int64)
    for box in boxes
        box.np[:] = 0
    end
    for p in parts
        p.indbox = ceil(p.pgrid[1]) + Lx * (ceil(p.pgrid[2])-1)
        boxes[p.indbox].np[p.tp] += 1
    end
end
################################################################################
#rotation function
function rotate_vec(v::Array{Float64,1}, α::Float64)
    R = [cos(α) sin(α) ; -sin(α) cos(α)]
    return R*v
end
################################################################################
#computing the momentum of the boxes
function box_vel(parts::Array{particle,1},boxes::Array{box,1})
    for (i, box) in enumerate(boxes) #this is for enumerating the boxes
        box.vel = zeros(2)
        tmass = 0 #initializating the total mass of the box
        for p in filter(x-> x.indbox == i, parts) #the loop is made in the particles inside the box i
            #println("inside the filter ", i)
            box.vel += p.mass * p.vel #sums the momentums of all particles
            tmass += p.mass #sums the mass
        end
        if tmass != 0
            box.vel /= tmass #normalize the momentum of the box with the mass of the particles.
        end
    end
end
################################################################################
#This is for the collisions of the particles of equal mass or same particles
function box_velmc(parts::Array{particle,1},boxes::Array{box,1}, m::Array{Float64,1})
    for (i, box) in enumerate(boxes) #this is for enumerating the boxes
        if countnz(box.np) <= 1; continue; end
        for j in m #cycling over the different masses.
            box.vel = zeros(2)
            tmass = 0 #initializating the total mass of the box
            for p in filter(x-> x.indbox == i && x.mass == j , parts) #the loop is made in the particles inside the box i
                #println("inside the filter ", i)
                box.vel += p.mass * p.vel #sums the momentums of all particles
                tmass += p.mass #sums the mass
            end
            if tmass != 0
            box.vel /= tmass #normalize the momentum of the box with the mass of the particles.
            end

        end
    end
end
################################################################################
#collision in the new way.
function collide_mc(parts::Array{particle,1})
    vel = zeros(2)
    tmass = 0 #initializating the total mass of the box
    for p in parts #the loop is made in the particles inside the box
        vel += p.mass * p.vel #sums the momentums of all particles
        tmass += p.mass #sums the mass
    end
    vel /= tmass #normalize the momentum of the box with the mass of the particles.
    for p in parts #loop over all particles
        v = p.vel - vel #extraction of the velocity of the box
        α = 90.0*rand([-1,1]) #the rotation angle
        vn = rotate_vec(v, α) #rotation of the vector
        p.vel = vn + vel #adding the vector and the velocity of the box
    end
end
function collide_sc(parts::Array{particle,1}, m)
    for j in m
        vel = zeros(2)
        tmass = 0
        sp = filter(x-> x.mass == j, parts)
        if isempty(sp); continue; end
        for p in sp
            vel += p.mass * p.vel
            tmass += p.mass
        end
        vel /= tmass
        for p in parts #loop over all particles
            v = p.vel - vel #extraction of the velocity of the box
            α = 90.0*rand([-1,1]) #the rotation angle
            vn = rotate_vec(v, α) #rotation of the vector
            p.vel = vn + vel #adding the vector and the velocity of the box
        end
    end
end
#function of the shifting, actually you shift the positions of the particles.
function shift_grid!(parts::Array{particle,1},a::Float64, dim::Array{Int64,1})
    δx = rand()*rand(-a/2:a/2)
    δy = rand()*rand(-a/2:a/2)
    for p in parts
        p.pgrid[1] = mod(p.pos[1] + δx, dim[1])
        p.pgrid[2] = mod(p.pos[2] + δy, dim[2])
    end
end
#now we need to shift back the particles.
function shiftback_grid!(parts::Array{particle,1})
    for p in parts
        p.pgrid[1] = p.pos[1]
        p.pgrid[2] = p.pos[2]
    end
end
################################################################################
#computing the new velocities of the particles, collisions.
function parts_vels!(parts::Array{particle,1}, boxes::Array{box,1}, angles::Array{Float64,1})
    for p in parts #loop over all particles
        v = p.vel - boxes[p.indbox].vel #extraction of the velocity of the box
        α = rand(angles)*rand([-1,1]) #the rotation angle
        vn = rotate_vec(v, α) #rotation of the vector
        p.vel = vn + boxes[p.indbox].vel #adding the vector and the velocity of the box
    end
end
################################################################################
#getting new positions of the particles with periodic boundary conditions (pbc)
function getpos_pbc!(parts::Array{particle,1}, τ::Float64, dim::Array{Int64,1})
    for p in parts
        p.pos[1] = mod(p.pos[1] + p.vel[1] * τ , dim[1]) #to get the periodic boundary conditions.
        p.pos[2] = mod(p.pos[2] + p.vel[2] * τ , dim[2])
    end
end
#with slip walls
function getpos_slip!(parts::Array{particle,1}, τ::Float64, dim::Array{Int64,1})
    for p in  parts
        #doing the bouncing on the x axis
        if p.pos[1] + p.vel[1] * τ > dim[1]
            dif = p.pos[1] + p.vel[1] * τ - dim[1]
            p.pos[1] = dim[1] - dif #it moves back
            p.vel[1] = -p.vel[1] #changing direction of velocity in x
        elseif p.pos[1] + p.vel[1] * τ < 0
            dif = abs(p.pos[1] + p.vel[1] * τ)
            p.pos[1] = 0 + dif
            p.vel[1] = -p.vel[1] #changin direction of velocity in x
        else
            p.pos[1] = p.pos[1] + p.vel[1] * τ
        end
        #now the bouncing on the y axis
        if p.pos[2] + p.vel[2] * τ > dim[2]
            dif = p.pos[2] + p.vel[2] * τ - dim[2]
            p.pos[2] = dim[2] - dif
            p.vel[2] = -p.vel[2]
        elseif p.pos[2] + p.vel[2] * τ < 0
            dif = abs(p.pos[2] + p.vel[2] * τ)
            p.pos[2] = 0 + dif
            p.vel[2] = -p.vel[2]
        else
            p.pos[2] = p.pos[2] + p.vel[2] * τ
        end
    end
end
################################################################################
#this is the ininitialization part, first the normalization of the total momentum.
function norm_momentum!(parts::Array{particle,1})
    vt = zeros(2) #sum of velocities
    mt = 0 #total mass
    for p in parts #summing all velocities and masses
        vt += p.vel
        mt += p.mass
    end
    vt /= mt #normalizing the sum
    for p in parts #applying the normalization of momentum
        p.vel = p.vel - vt
    end
end
################################################################################
#now the normalization of temperature
function norm_temperature!(parts::Array{particle,1}, Tr::Float64)
    T = 0 #temperature
    for p in parts #this is the calculation of the Temperature
        T += p.mass * norm(p.vel)^2 / (2 * length(parts))
    end
    for p in parts
        p.vel = p.vel * sqrt(Tr / T) #multiplying each velocity for the factor of sqrt(Ttarget/Tactual)
    end
end
################################################################################
#this function compute probabilities of P(n(t+1)|n(t)) and returns how many particles are going to react.
function prob_box(np::Array{Int64}, kr::Float64)
    nr = min(np[1],np[2]) #estimating minimum or maximum of reactions
    p = zeros(nr+1) #array of probabilities
    p[1] = 1 - kr #no-reaction case
    for i=2:nr
        p[i+1] = kr * ( factorial(np[1]) / factorial(np[1]-i) * factorial(np[2]) / factorial(np[2]-i) * factorial(np[3]) / factorial(np[3]+i) )
    end
    s = sum(p)
    if s != 1
        p = p / s
    end
    cs = cumsum(p); a = rand()
    return findfirst(sort[a;cs],a)-1
end
################################################################################
function col_box(box::box,parts::Array{particle,1}, m::Array{Float64,1})
    parbox = filter(x-> x.indbox == box.ind, parts) #selecting the particles that are in the box.
    #println(i,'\t',length(parbox))
    if isempty(parbox) == false #ignore next steps if the box is empty
        collide_mc(parbox)
        if countnz(box.np) > 1
            collide_sc(parbox, m)
        end
    end
end
