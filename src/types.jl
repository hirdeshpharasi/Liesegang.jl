#definition of particles
type particle
    pos::Array{Float64,1} #postitions
    vel::Array{Float64,1} #velocities
    mass::Float64 #mass
    indbox::Int64 #the box where the particle is
    function particle(dim::Array{Int64,1}, m::Float64)
        this = new()
        this.pos = rand(2) .* dim
        this.vel = rand(2)
        this.mass = m
        this.indbox = ceil(this.pos[1]) + dim[1] * (ceil(this.pos[2])-1)
        return this
    end
end
#definition of a box
type box
    vel::Array{Float64,1} # the velocity of the box.
    #npart::Int64  not really sure if is neccesary to have the number of particles of each box.
    function box()
        this = new()
        this.vel = zeros(2)
        return this
    end
end
