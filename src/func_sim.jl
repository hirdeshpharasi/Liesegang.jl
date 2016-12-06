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
function box_vel(parts::Array{particle,1},boxes::Array{box,1}, Ly::Int64)
    for (i, box) in enumerate(boxes) #this is for enumerating the boxes
        tmass = 0 #initializating the total mass of the box
        for p in filter(x-> x.indbox == i, parts) #the loop is made in the particles inside the box i
            box.vel += p.mass * p.vel #sums the momentums of all particles
            tmass += p.mass #sums the mass
        end
        box.vel /= tmass #normalize the momentum of the box with the mass of the particles.
    end
end
################################################################################
#computing the new velocities of the particles
function parts_vels!(parts::Array{particle,1}, boxes::Array{box,1}, α::Float64)
    for p in parts #loop over all particles
        v = p.vel - boxes[p.indbox].vel #extraction of the velocity of the box
        vn = rotate_vec(v, α) #rotation of the vector
        p.vel = vn + boxes[p.indbox].vel #adding the vector and the velocity of the box
    end
end
################################################################################
#getting new positions of the particles.
function getpos_parts!(parts::Array{particle,1}, τ::Float64, Dim::Array{Int64,1})
    for p in parts
        p.pos[1] += mod(p.vel[1] * τ, Dim[1]) #to get the periodic boundary conditions.
        p.pos[2] += mod(p.vel[2] * τ, Dim[2])
    end
end
