from vpython import *
from vpython import *
import numpy as np
from histogram import *

N = 200
e=-0.0001
m, size = 4E-3/6E23, 31E-12*10 # He atoms are 10 times bigger for easiear collision but not too big for accuracy
L = ((24.4E-3/(6E23))*N)**(1/3.0)/2 + size 
Length,Height,Width=2*L,2*L,2*L # 2L is the cubic container's original length, width, and height
gamma=5/3
k, T_initial = 1.38E-23, 298.0 # Boltzmann Constant and initial temperature
t, dt = 0, 3E-13
v_W=L/20000.0/dt
vrms = (2*k*1.5*T_initial/m)**0.5 # the initial root mean square velocity
atoms = [] # list to store atoms

# histogram setting
deltav = 50. # slotwidth for v histogram
vdist = graph(x=800, y=0, ymax = N*deltav/1000.,width=500, height=300, xtitle='v', ytitle='dN', align = 'left')
theory_low_T = gcurve(color=color.cyan) # for plot of the curve for the atom speed distribution

dv = 10.
for v in arange(0.,4201.+dv,dv): # theoretical speed distribution (Maxwell-Boltzmann distribution)
    theory_low_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*T_initial))**1.5)*exp((-0.5*m*v**2)/(k*T_initial))*(v**2)*dv))
observation = ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.red) # for the simulation speed distribution
observation2 = ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.blue)
diagram=graph(width=450,align='right')
pv_diagram=gcurve(graph=diagram,color=color.blue)
#initialization
scene = canvas(width=500, height=500, background=vector(0.2,0.2,0), align = 'left')
container = box(length = 2*Length, height = Height, width = Width, opacity=0.2, color = color.yellow )
rectangle1 = box(length = Length/100, height = Height, width = Width, opacity=0.6, color = color.yellow, pos = vec(-Length/2, 0, 0) )
rectangle2 = box(length = Length/100, height = Height, width = Width, opacity=0.6, color = color.yellow, pos = vec(Length/2, 0, 0) )
p_a, v_a = np.zeros((N,3)), np.zeros((N,3)) # particle position array and particle velocity array, N particles and 3 for x, y, z

for i in range(N):
    p_a[i] = [2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L] # particle is initially random positioned in container
    if i== N-1: # the last atom is with yellow color and leaves a trail (retain 50 points)
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=color.yellow, make_trail = True, retain = 50)
    else: # other atoms are with random color and leaves no trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=vector(random(), random(), random()))
    ra = pi*random()
    rb = 2*pi*random()
    v_a[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)] # particle initially same speed but random direction
    atoms.append(atom)

def vcollision(a1p, a2p, a1v,a2v): # the function for handling velocity after collisions between two atoms
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime

T=0
flag=0 # used to record 1000*dt
p=0 # the total momentum impacted on the walls
stage=0
while True:
    #container.length=Length
    t += dt
    flag+=1
    rate(1000)
    rectangle1.pos.x=-Length/2
    rectangle2.pos.x=Length/2
    p_a += v_a*dt # calculate new positions for all atoms
    
    if stage==0 :
        Length+=2*v_W*dt
        e=0
    if Length>=3*L and stage==0 : 
        stage=1
        T_initial=T
    if stage==1 : 
        Length+=2*v_W*dt
         # adiabatic compression
    if Length>=4*L and stage==1 : 
        stage=2 # stop the walls from moving
        v_W*=(-1)
        
    if stage==2 :
        Length+=2*v_W*dt
        e=0 
    if T>=298.0 and stage==2 :
        stage=3 
        T_initial=298.0
    if stage==3 : 
        Length+=2*v_W*dt
        
    if stage==3 and Length <=2*L:
        stage=0
        v_W*=(-1)
         
        
    for i in range(N): atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]) # to display atoms at new positions
    observation.plot(data = np.sqrt(np.sum(np.square(v_a),-1))) # plot the speed distribution of atoms
     
### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:,np.newaxis] # array for vector from one atom to another atom for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) # distance array between atoms for all pairs of atoms
    hit = np.less_equal(rmag,2*size)-np.identity(N) # if smaller than 2*size meaning these two atoms might hit each other
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() # change hit to a list
    
    for ij in hitlist: # i,j encoded as i*Natoms+j
        i, j = divmod(ij,N) # atom pair, i-th and j-th atoms, hit each other
        hitlist.remove(j*N+i) # remove j,i pair from list to avoid handling the collision twice
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0 : # only handling collision if two atoms are approaching each other
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j]) # handle collision

#find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
        if abs(p_a[i][0]) >= Length/2 - size and p_a[i][0]*v_a[i][0] > 0 :
            if p_a[i][0]>0: # collisions between atoms and the right wall
                p+=2*m*(abs(v_a[i][0])+v_W)
                v_a[i][0]=(1+e)*(-v_a[i][0]+2*v_W)
            elif p_a[i][0]<0: # collisions between atoms and the left wall
                p+=2*m*(abs(v_a[i][0])+v_W)
                v_a[i][0]=(1+e)*(-v_a[i][0]-2*v_W)
        if abs(p_a[i][1]) >= Height/2 - size and p_a[i][1]*v_a[i][1] > 0 :
            v_a[i][1] = - (1+e)*v_a[i][1]
            p+=2*m*abs(v_a[i][1])
        if abs(p_a[i][2]) >= Width/2 - size and p_a[i][2]*v_a[i][2] > 0 :
            v_a[i][2] = - (1+e)*v_a[i][2]
            p+=2*m*abs(v_a[i][2])
            
    if flag==100: # after 1000*dt
        flag=0 # initialize
        T=0 # temperature
        area=(Length*Width+Length*Height+Height*Width)*2 # area of the six walls
        P=p/1000/dt/area # pressure = total momentum divided by time and area
        V=Length*Height*Width #volume
        for i in range(N): # total kinetic energy = 3/2*N*k*T
            T+=(2.0/3.0)*0.5*m*(v_a[i,0]**2+v_a[i,1]**2+v_a[i,2]**2)/N/k
        e=0.0001*(T_initial-T)
        print('T:',T,end=' ') # print current macrostates
        print('p:',P,end=' ')
        print('V:',V,end=' ')
        print('P*V=',P*V,end=' ') # ideal gas law : pV=NkT
        print('N*k*T=',N*k*T,end=' ')
        print('P*(V**gamma)=',P*(V**gamma)) # in adiabatic process, p*(V**gamma)=constant
        print('stage=', stage, '\n')
        pv_diagram.plot(pos=(V,P))
        p=0 # initialize