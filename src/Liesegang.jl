module Liesegang
    include("types.jl")

    export get_box, rotate_vec, box_vel, parts_vels!

    include("func_sim.jl")

end
