### =============== ### =============== ### =============== ###
### EXAMPLE OF SIMPLE UNBOUNDED MOTION
### MARTIN ZUMAYA HERNANDEZ
### 10 - NOV - 2016
### =============== ### =============== ### =============== ###

using Plots # Plotting Package

### =============== ### =============== ### =============== ###
### DEFINITION OF INITIAL PARAMETERS
### =============== ### =============== ### =============== ###

ρ  = 5.0 # density
N  = 50 # number of particles
L  = sqrt(N / ρ) # size of box
dt = 1.0 # integration step
v0 = 1.0 # particle's speed
T  = 100 # integration time steps

### =============== ### =============== ### =============== ###
### INITIALIZATION OF PARTICLES INITIAL POSITIONS AND VELOCIDITES
### =============== ### =============== ### =============== ###

pos = [ 2*rand()*L - L for i in 1:2N ] # array of random initial particles' postitions
vel = v0 * vcat([ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:N ]...) # array of  particles' velocities

### =============== ### =============== ### =============== ###
### INITIALIZATION OF PARTICLES INITIAL POSITIONS AND VELOCIDITES
### =============== ### =============== ### =============== ###

anim = Animation()

for i in 1:T

    pos += vel*dt # positions update

    x = [pos[i] for i in 1:2:N]
    y = [pos[i+1] for i in 1:2:N]

    pts = vec(P2[(x[i],y[i]) for i=1:length(x)])

    # overlapping graphs
    quiver!(pts, quiver = ([vel[i] for i in 1:2:N], [vel[i+1] for i in 1:2:N]))

    # overwriting graphs
     #quiver(pts, quiver = ([vel[i] for i in 1:2:N], [vel[i+1] for i in 1:2:N]))

    frame(anim)

end

gif(anim, "motion.gif", fps=12)
