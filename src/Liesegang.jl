module Liesegang
    include("types.jl")

    export particle, box
    export get_box, rotate_vec, box_vel, parts_vels!, getpos_pbc!, norm_momentum!, norm_temperature!, shift_grid!, shiftback_grid!, getpos_slip!

    include("func_sim.jl")

end
