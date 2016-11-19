pos = zeros(10)

for i in 1:10
    pos[i] = rand()
end

pos = [ [rand(), rand(), rand()] for i in 1:10 ]
pos = [ [rand(), rand(), rand()] for i in 1:10 ]

for i in 1:N

    pos[i] += vel[i]*dt

end

for i in 1:3:3N

    for j in i:3:3N

    sqrt((pos[i]-pos[j])^2 + (pos[i+1]-pos[j+1])^2 + (pos[i+2]-pos[j+2])^2)
    
    end
end
