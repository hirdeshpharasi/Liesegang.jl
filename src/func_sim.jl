##############################
#This is to get the number of box where the particle is.
function get_box(p::particle, Ly::Int64)
    p.indbox = ceil(p.pos[1]) + Ly * (ceil(p.pos[2])-1)
end
####################################
#rotation function
function rotate_vec(v::Array{Float64,1}, α)
    R = [cos(α) sin(α) ; -sin(α) cos(α)]
    return R*v
end
##################################################
#computing the momentum of the boxes
function box_m(parts::Array{particle,1},boxes::Array{box,1}, Ly::Int64)
    for (i, box) in enumerate(boxes) #this is for enumerating the boxes
        tmass = 0 #initializating the total mass of the box
        for p in filter(x-> x.indbox == i, parts) #the loop is made in the particles inside the box i
            box.vel += p.mass * p.vel #sums the momentums of all particles
            tmass += p.mass #sums the mass
        end
        box.vel /= tmass #normalize the momentum of the box with the mass of the particles.
    end
end
