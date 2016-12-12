################################################################################
##################### MODULES FOR MULTI PARTICLE COLLISIONS ####################
################################################################################
#This is to get the number of box where the particle is.
function get_box(parts::Array{particle,1}, Lx::Int64)
    for p in parts
        p.indbox = ceil(p.pos[1]) + Lx * (ceil(p.pos[2])-1)
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
#function of the shifting, actually you shift the positions of the particles.
function shift_grid!(parts::Array{particle,1},a::Float64, dim::Array{Int64,1})
    δx = rand()*rand(-a/2:a/2)
    δy = rand()*rand(-a/2:a/2)
    for p in parts
        p.pos[1] = mod(p.pos[1] + δx, dim[1])
        p.pos[2] = mod(p.pos[2] + δy, dim[2])
    end
    return [δx, δy]
end
#now we need to shift back the particles.
function shiftback_grid!(parts::Array{particle,1}, dim::Array{Int64,1}, δ::Array{Float64,1})
    for p in parts
        p.pos[1] = mod(p.pos[1] + δ[1], dim[1])
        p.pos[2] = mod(p.pos[2] + δ[2], dim[2])
    end
end
################################################################################
#computing the new velocities of the particles
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
            p.pos[1] = p.pos[1] - dif
        elseif p.pos[1] + p.vel[1] * τ < 0
            dif = abs(p.pos[1] + p.vel[1] * τ)
            p.pos[1] = dif
        else
            p.pos[1] = p.pos[1] + p.vel[1] * τ
        end
        #now the bouncing on the y axis
        if p.pos[2] + p.vel[2] * τ > dim[2]
            dif = p.pos[2] + p.vel[2] * τ - dim[2]
            p.pos[2] = p.pos[2] - dif
        elseif p.pos[2] + p.vel[2] * τ < 0
            dif = abs(p.pos[2] + p.vel[2] * τ)
            p.pos[2] = dif
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
