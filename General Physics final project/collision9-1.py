from H2O6 import*
from vpython import *
import numpy as np

N = 200
L = ((24.4E-3/(6E23))*N)**(1/3.0)/2 + size 
Length,Height,Width=2*L,2*L,2*L # 2L is the cubic container's original length, width, and height
k, T = 1.38E-23, 3000 # Boltzmann Constant and initial temperature
t, dt = 0, 5E-16
vrms = (2*k*1.5*T/m)**0.5 # the initial root mean square velocity
H2Os = [] # list to store atoms

#initialization
scene = canvas(width=500, height=500, background=vector(0.2,0.2,0), align = 'left')
container = box(length = Length, height = Height, width = Width, opacity=0.2, color = color.yellow )
p_a, v_a= np.zeros((3*N,3)), np.zeros((3*N,3)) # particle position array and particle velocity array, N particles and 3 for x, y, z
axis1 = vector(d, 0, 0)
axis2 = vector(d*math.cos(theta), d*math.sin(theta), 0)    
i=0
while i < 3*N:
    
    a,b,c=2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L
    p_a[i]=[a,b,c]# random() yields a random number between 0 and 1
    p_a[i+1] = p_a[i] + [d,0,0]
    p_a[i+2] = p_a[i] + [d*math.cos(theta), d*math.sin(theta), 0]
    H2O = H2O_molecule(pos=vec(a,b,c)) # generate one H2O molecule
    ra1, ra2, ra3 = pi*random(), pi*random(), pi*random()
    rb1, rb2, rb3 = 2*pi*random(), 2*pi*random(), 2*pi*random()
    v_a[i] = [vrms*sin(ra1)*cos(rb1), vrms*sin(ra1)*sin(rb1), vrms*cos(ra1)] # particle initially same speed but random direction
    v_a[i+1] = [vrms*sin(ra2)*cos(rb2), vrms*sin(ra2)*sin(rb2), vrms*cos(ra2)] # particle initially same speed but random direction
    v_a[i+2] = [vrms*sin(ra3)*cos(rb3), vrms*sin(ra3)*sin(rb3), vrms*cos(ra3)] # particle initially same speed but random direction
    H2Os.append(H2O) # store this molecule into list H2Os
    i+=3
    
def vcollision(a1p, a2p, a1v,a2v): # the function for handling velocity after collisions between two atoms
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime
            
flag=0 # used to record 1000*dt
while True:
    t += dt
    flag+=1
    rate(1000)
    
    p_a += v_a*dt # calculate new positions for all atoms
    i=0
    for i in range(N): # to display atoms at new positions
        #H2Os[i].pos = vector(p_a[i*3, 0], p_a[i*3, 1], p_a[i*3, 2])
        H2Os[i].O.pos = vector(p_a[i*3, 0], p_a[i*3, 1], p_a[i*3, 2])
        H2Os[i].H1.pos = vector(p_a[i*3+1, 0], p_a[i*3+1, 1], p_a[i*3+1, 2])
        H2Os[i].H2.pos = vector(p_a[i*3+2, 0], p_a[i*3+2, 1], p_a[i*3+2, 2])
        H2Os[i].bond1.pos = H2Os[i].O.pos
        H2Os[i].bond2.pos = H2Os[i].O.pos
        H2Os[i].bond3.pos = H2Os[i].H1.pos
        H2Os[i].bond1.axis = H2Os[i].H1.pos - H2Os[i].O.pos
        H2Os[i].bond2.axis = H2Os[i].H2.pos - H2Os[i].O.pos
        H2Os[i].bond3.axis = H2Os[i].H2.pos - H2Os[i].H1.pos
        
        a,b,c=(H2Os[i].bond_force_on_O()*dt / H2Os[i].O.m).x, (H2Os[i].bond_force_on_O()*dt / H2Os[i].O.m).y, (H2Os[i].bond_force_on_O()*dt / H2Os[i].O.m).z
        v_a[3*i]+=[a,b,c]
        #print(H2Os[i].bond_force_on_O(),H2Os[i].bond_force_on_O()*dt / H2Os[i].O.m,a,b,c)
        a,b,c=(H2Os[i].bond_force_on_H1() / H2Os[i].H1.m*dt).x, (H2Os[i].bond_force_on_H1() / H2Os[i].H1.m*dt).y, (H2Os[i].bond_force_on_H1() / H2Os[i].H1.m*dt).z
        v_a[3*i+1]+=[a,b,c]
        #print(a,b,c)
        a,b,c=(H2Os[i].bond_force_on_H2() / H2Os[i].H2.m*dt).x, (H2Os[i].bond_force_on_H2() / H2Os[i].H2.m*dt).y, (H2Os[i].bond_force_on_H2() / H2Os[i].H2.m*dt).z
        v_a[3*i+2]+=[a,b,c]
        #print(a,b,c)
    
### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:,np.newaxis] # array for vector from one atom to another atom for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) # distance array between atoms for all pairs of atoms
    hit = np.less_equal(rmag,2*size)
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() # change hit to a list
    
    for ij in hitlist: # i,j encoded as i*Natoms+j
        i, j = divmod(ij,3*N) # atom pair, i-th and j-th atoms, hit each other
        hitlist.remove(j*3*N+i) # remove j,i pair from list to avoid handling the collision twice
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0 : # only handling collision if two atoms are approaching each other
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j]) # handle collision

#find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(3*N):
            if abs(p_a[i][0]) >= Length/2 - size and p_a[i][0]*v_a[i][0] > 0 :
                v_a[i][0] = - v_a[i][0]
            
            if abs(p_a[i][1]) >= Height/2 - size and p_a[i][1]*v_a[i][1] > 0 :
                v_a[i][1] = - v_a[i][1]
                
            if abs(p_a[i][2]) >= Width/2 - size and p_a[i][2]*v_a[i][2] > 0 :
                v_a[i][2] = - v_a[i][2]

    if flag==100:
        flag=0
        T=0
        total_v=v_a**2
        for i in range(3*N): # total kinetic energy = 3/2*N*k*T
            T+=(2.0/3.0)*0.5*m*sum(total_v[i])/3/N/k
        for i in range(N):
            T+=H2Os[i].Total_U()
        print(T)       