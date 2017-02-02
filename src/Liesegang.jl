module Liesegang
    include("types.jl")

    export particle, box
    export col_box, prob_box, get_box, rotate_vec, box_vel, parts_vels!, getpos_pbc!, norm_momentum!, norm_temperature!, shift_grid!, shiftback_grid!, getpos_slip!, box_velmc, collide_mc, collide_sc

    include("func_sim.jl")

end
