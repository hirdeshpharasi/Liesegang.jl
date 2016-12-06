#definition of particles
type particle
    pos::Array{Float64,1} #postitions
    vel::Array{Float64,1} #velocities
    m::Int64 #mass
    indbox::Int64 #the box where the particle is
end
#definition of a box
type box
    pos::Array{Float64,1}
    vel::Array{Float64,1}
    #npart::Int64  not really sure if is neccesary to have the number of particles of each box.
end
