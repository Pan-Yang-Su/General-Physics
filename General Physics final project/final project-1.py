from vpython import *
from vpython import *
import numpy as np
from histogram import *


###some constants
#######################################################################
N = 200
N2 = 200
m, size = 4E-3/6E23, 31E-12*10 # He atoms are 10 times bigger for easiear collision but not too big for accuracy
L = ((24.4E-3/(6E23))*N)**(1/3.0)/2 + size 
Length, Height, Width, pipelength, piperadius = 2*L, 2*L, 2*L, 8*L, 0.3*L
pos = vector(0, 0, 0)#scene center
gamma=5/3
t, dt = 0, 3E-13
v_W=L/20000.0/dt
k, T, T_initial = 1.38E-23, 298.0, 298.0 # Boltzmann Constant and initial temperature
e=0.01
vrms = (2*k*1.5*T/m)**0.5 # the initial root mean square velocity

scene = canvas(width = 500, height = 500, center = vec(0,0,0), background = vec(1,1,1), forward = vec(0,0,-1))

#spring = box(pos=vec(pipelength + 7*L, pipelength + 1/2*Height,0),width=L,length=L,height=L,color=color.yellow)
def keyinput(evt): #keyboard callback function
    global pos, angle
    move = {'left': vector(-L/10, 0, 0), 'right': vector(L/10, 0, 0),
        'up': vector(0, L/10, 0),
        'down': vector(0, L/10, 0), 'i' : vector(0, 0, L/10),
        'o': vector(0, 0, L/10)}
    
    s = evt.key
    if s in move : pos = pos + move[s]
scene.bind('keydown', keyinput) # setting for the binding function

###steam turbine
################################################################
atoms = [] # list to store atoms 
container1 = box(length = Length, height = Height, width = Width, opacity = 0.2, color = color.blue )
container2 = box(pos = vec(pipelength+Length,0,0), length = Length, height = Height, width = Width, opacity = 0.2, color = color.blue )
pipe1 = box( pos = vec(0,Height/2+pipelength/2,0),height =pipelength, length = Length/2,width = Width/2, opacity = 0.2, color = color.green)
pipe2 = box( pos = vec(pipelength/2+Length/2,pipelength+Height/2+Height/4,0),height =Height/2, length =pipelength+Length+Length/2 ,width = Width/2 , opacity = 0.2, color = color.green)
pipe3 = box( pos = vec(Length/2+pipelength/2,0,0),height =Height/2, length =pipelength ,width = Width/2, opacity = 0.2, color = color.green)
pipe4 = box( pos = vec(Length/2+pipelength+Length/2,Height/2+pipelength/2,0),height =pipelength, length = Length/2,width = Width/2, opacity = 0.2, color = color.green)
desk=box(pos=vec(container1.length/2+pipe3.length/2,-container1.height/2-pipe3.height/2,0),height=7*pipe3.height/8,length=2*pipe3.length,width=2*container1.width,texture=textures.stucco)

turbo = cylinder(pos = vec(pipelength + 6*Length/8 , pipelength + Height/2, -Length/2), axis = vec(0, 0, Length), radius = 1/2*Length, thickness = 1,emissive=True, color = color.green, opacity = 0.2)
leave_1 = box(pos = vec(pipelength + 6/8*Length, pipelength + 1/2*Height, 0), length = Length, height = 11/12*Length, width = 0.01*Length, axis = vec(0, 0, Length), m = 0.01, color = color.black)
leave_2 = box(pos = vec(pipelength + 6/8*Length, pipelength + 1/2*Height, 0), length = Length, height = 11/12*Length, width = 0.01*Length, axis = vec(0, 0, Length), m = 0.01, color = color.black)
leave_3 = box(pos = vec(pipelength + 6/8*Length, pipelength + 1/2*Height, 0), length = Length, height = 11/12*Length, width = 0.01*Length, axis = vec(0, 0, Length), m = 0.01, color = color.black)
leave_4 = box(pos = vec(pipelength + 6/8*Length, pipelength + 1/2*Height, 0), length = Length, height = 11/12*Length, width = 0.01*Length, axis = vec(0, 0, Length), m = 0.01, color = color.black)


leave_2.rotate(angle = pi*1/4, axis = vec(0, 0, Length))
leave_3.rotate(angle = pi*2/4, axis = vec(0, 0, Length))
leave_4.rotate(angle = pi*3/4, axis = vec(0, 0, Length))

#bar1 = box(pos=vec(pipelength + Length, pipelength + Height/2, -Length/2), length = Length/4, height = 0.1*Length, width = 0.01*Length, axis = vec(0, Length, 0), color = color.black)
#bar2 = box(pos=vec(pipelength + Length, pipelength + 1/4*Height, -Length/2), length = 2.4*Length, height = 0.1*Length, width = 0.01*Length, axis = vec(-Length, 0, 0), color = color.black)

cylinder1=cylinder(pos = vec(pipelength + 9/8*Length , pipelength + Height/2, 0), axis = vec(2*Length, 0, 0), radius = 1/10*Length, thickness = 1,emissive=True,texture=textures.stucco )

rotational_inertia = 4*1/12*leave_1.m*leave_1.height**2
R = turbo.radius*11/12
M = 4*leave_1.m
###carnot refrigerator
######################################################################### 
atoms2 = []
container3 = box(pos=vec(pipelength + 4*Length, pipelength + 1/2*Height, 0), length = 3*Length, height = Height, width = Width, opacity=0.2, color = color.yellow )
rectangle1 = box(length = Length/100, height = Height, width = Width, opacity=0.6, color = color.yellow, pos = vec(pipelength + 7*L, pipelength + 1/2*Height, 0) )
rectangle2 = box(length = Length/100, height = Height, width = Width, opacity=0.6, color = color.yellow, pos = vec(pipelength + 9*L, pipelength + 1/2*Height, 0) )
p_a2, v_a2 = np.zeros((N2,3)), np.zeros((N2,3)) # particle position array and particle velocity array, N particles and 3 for x, y, z

for i in range(N2):
    p_a2[i] = [L*random() +8*L + pipelength, 2 * L*random() - L + pipelength+Height/2, 2 * L*random() - L] # particle is initially random positioned in container
    if i== N2-1: # the last atom is with yellow color and leaves a trail (retain 50 points)
        atom = sphere(pos=vector(p_a2[i, 0], p_a2[i, 1], p_a2[i, 2]), radius = size, color=color.yellow, make_trail = True, retain = 50, trail_radius = size)
    else: # other atoms are with random color and leaves no trail
        atom = sphere(pos=vector(p_a2[i, 0], p_a2[i, 1], p_a2[i, 2]), radius = size, color=vector(random(), random(), random()))
    ra = pi*random()
    rb = 2*pi*random()
    v_a2[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)] # particle initially same speed but random direction
    atoms2.append(atom)

###histogram setting
####################################################################### 
deltav = 50. # slotwidth for v histogram
vdist = graph(x=800, y=0, ymax = N*deltav/1000.,width=500, height=300, xtitle='v', ytitle='dN', align = 'left')
theory_low_T = gcurve(color=color.cyan) # for plot of the curve for the atom speed distribution
diagram=graph(width=450,align='right')
pv_diagram=gcurve(graph=diagram,color=color.blue)

dv = 10.
for v in arange(0.,4201.+dv,dv): # theoretical speed distribution (Maxwell-Boltzmann distribution)
    theory_low_T.plot(pos=(v,(deltav/dv)*N2*4.*pi*((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))
observation = ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.red) # for the simulation speed distribution
observation2 = ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.blue)

#initialization
p_a, v_a = np.zeros((N,3)), np.zeros((N,3)) # particle position array and particle velocity array, N particles and 3 for x, y, z
for i in range(N):
    p_a[i] = [2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L] # particle is initially random positioned in container
    if i== N-1: # the last atom is with yellow color and leaves a trail (retain 50 points)
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=color.yellow, make_trail = True, retain = 50, trail_radius = size, flag=0)
    else: # other atoms are with random color and leaves no trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=vector(random(), random(), random()), flag=0)
    ra = pi*random()
    rb = 2*pi*random()
    v_a[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)] # particle initially same speed but random direction
    atoms.append(atom)

def vcollision(a1p, a2p, a1v,a2v): # the function for handling velocity after collisions between two atoms
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime
e2=0
omega=0
dV=0
flag=0 # used to record 1000*dt
p=0 # the total momentum impacted on the walls
stage=0
T=0
P=100000000000000000


#spring.axis = vec(pipelength+Length, pipelength+1/2*Height, 0) - rectangle1.pos
while True:
    #spring.axis = vec(pipelength+Length, pipelength+1/2*Height, 0)- rectangle1.pos
    scene.center = pos
    cylinder1.axis.x = rectangle1.pos.x - (pipelength + 2.25*L)
    t += dt
    flag+=1
    rate(10000)
    #print("angular_momentum = ", angular_momentum)
    p_a += v_a*dt # calculate new positions for all atoms
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
    rotational_energy = 0    
    for i in range(N):
        if (atoms[i].pos.x-pipelength-6/8*Length)**2 + (atoms[i].pos.y-pipelength-Height/2)**2 <= R**2 and v_a[i][0] >= 0 and atoms[i].pos.y >= pipelength+Height/2:
            v_a[i][0] = v_a[i][0]/2
            rotational_energy += abs(1/2*m*(v_a[i][0]/2)**2)
        elif (atoms[i].pos.x-pipelength-6/8*Length)**2 + (atoms[i].pos.y-pipelength-Height/2)**2 <= R**2 and v_a[i][1] <= 0 and atoms[i].pos.x >= pipelength+7/8*Length:
            v_a[i][1] = v_a[i][1]/2
            rotational_energy += abs(1/2*m*(v_a[i][0]/2)**2)
            
    omega += sqrt(2*rotational_energy/rotational_inertia)
    leave_1.rotate(angle = -omega * dt*10**15 , axis = vec(0, 0, Length))
    leave_2.rotate(angle = -omega * dt*10**15 , axis = vec(0, 0, Length))
    leave_3.rotate(angle = -omega * dt*10**15 , axis = vec(0, 0, Length))
    leave_4.rotate(angle = -omega * dt*10**15 , axis = vec(0, 0, Length))
    #print(omega)
    if omega >= 0.01:#1/2*I*omega**2 = PdV
        dV=1/2*rotational_inertia*(omega**2-0.01**2)/P
        omega=0.01
    
    
    
    
    
#find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
    
        ##container1 and pipe1/container1(y)
        if (p_a[i][1]) >= (container1.height)/2 - size and abs(p_a[i][0])<(pipe1.length)/2 and abs(p_a[i][2])<(pipe1.width)/2 and atoms[i].flag==0:
                atoms[i].flag=1#container1 to pipe1
        if abs(p_a[i][1]) >= (container1.height)/2 - size and p_a[i][1]* v_a[i][1] > 0 and atoms[i].flag==0:v_a[i][1] = - (1+e)*v_a[i][1] #in container1
     
        if (p_a[i][1]) <= (container1.height)/2 + size and abs(p_a[i][0])<(pipe1.length)/2 and abs(p_a[i][2])<(pipe1.width)/2 and atoms[i].flag==1:
            atoms[i].flag=0#pipe1 to container1 
    
       ##pipe1 and pipe2/pipe2(y)/pipe2 and pipe4
        if (p_a[i][1]) >= (container1.height)/2+pipe1.height- size and atoms[i].flag==1:
            atoms[i].flag=2#pipe1 to pipe2
        if (p_a[i][1]) >= (container1.height)/2+(pipe2.height)+pipe1.height- size  and p_a[i][1]* v_a[i][1] > 0 and atoms[i].flag==2:
            v_a[i][1] = -v_a[i][1]#in pipe2 with upper wall
        if (p_a[i][1]) <= (container1.height)/2+pipe1.height+ size and abs(p_a[i][0])<(pipe1.length)/2 and abs(p_a[i][2])<(pipe1.width)/2 and atoms[i].flag==2:
            atoms[i].flag=1#pipe2 to pipe1
    
        if (p_a[i][1]) <= (container2.height)/2+pipe4.height+ size and (p_a[i][0])>(container1.length/2+pipe3.length+container2.length/4) and abs(p_a[i][2])<(pipe4.width)/2 and atoms[i].flag==2:
            atoms[i].flag=5#pipe2 to pipe4
        if (p_a[i][1]) <= (container1.height)/2+pipe1.height+ size and (p_a[i][0])>=(pipe1.length)/2 and (p_a[i][0])<=(container1.length/2+pipe3.length+container2.length/4) and atoms[i].flag==2:
            v_a[i][1] = -v_a[i][1]#in pipe2 with down wall
        ### pipe3(y)
        if abs(p_a[i][1]) >= (pipe3.height)/2-size and p_a[i][1]* v_a[i][1] > 0 and atoms[i].flag==3:
            v_a[i][1] = -v_a[i][1]#in pipe3
        ### container2 and pipe4/ container2 (y)
            
        if (p_a[i][1]) >= (container2.height)/2-size and (p_a[i][0])>(container1.length/2+pipe3.length+container2.length/4) and (p_a[i][0])<(container1.length/2+pipe3.length+3*container2.length/4) and abs(p_a[i][2])<(pipe4.width)/2 and atoms[i].flag==4:
            atoms[i].flag=5#container2 to pipe4
        if abs(p_a[i][1]) >= (container2.height)/2-size and p_a[i][1]* v_a[i][1] > 0 and atoms[i].flag==4:
            v_a[i][1] = -(1-e)*v_a[i][1]#in container2
            
        ### pipe4(y)
        if (p_a[i][1]) >= (container2.height)/2+pipe4.height- size and atoms[i].flag==5:
            atoms[i].flag=2#pipe4 to pipe2
        if (p_a[i][1]) <= (container2.height)/2+size and atoms[i].flag==5:
            atoms[i].flag=4#pipe4 to container2
        ### container1 and pipe3/container1(x)
        if (p_a[i][0]) >= (container1.length)/2 - size and abs(p_a[i][1])<(pipe3.height)/2 and abs(p_a[i][2])<(pipe3.width)/2 and atoms[i].flag==0:
            atoms[i].flag=3#container1 to pipe3
             
        if abs(p_a[i][0]) >= (container1.length)/2 - size and p_a[i][0]*v_a[i][0] > 0 and atoms[i].flag==0:v_a[i][0] = - (1+e)*v_a[i][0] #in container1
            
        if (p_a[i][0]) <= (container1.length)/2 + size and abs(p_a[i][1])<(pipe3.height)/2 and abs(p_a[i][2])<(pipe3.width)/2 and atoms[i].flag==3:atoms[i].flag=0#pipe3 to container1 
            
        ###pipe3 and container2/container2(x)
        if (p_a[i][0]) >= (container1.length)/2+pipe3.length-size and atoms[i].flag==3:
            atoms[i].flag=4#pipe3 to container2
    
        if (p_a[i][0]) <= (container1.length)/2+pipe3.length + size and abs(p_a[i][1])<(pipe3.height)/2 and abs(p_a[i][2])<(pipe3.width)/2 and atoms[i].flag==4:
            atoms[i].flag=3#container2 to pipe3
    
        if (p_a[i][0]) >= (container1.length)/2+(container2.length)+pipe3.length-size  and p_a[i][0]*v_a[i][0] > 0 and atoms[i].flag==4: v_a[i][0] = -(1-e)*v_a[i][0] #in container2 left side
            
        if (p_a[i][0]) <= (container1.length)/2+pipe3.length+size and p_a[i][0]*v_a[i][0] < 0 and atoms[i].flag==4 :
            v_a[i][0] = -(1-e)*v_a[i][0] #in container2 right side      
        ###pipe1 and pipe2 and pipe4 (x)
        if abs(p_a[i][0]) >= (pipe1.length)/2 - size and p_a[i][0]*v_a[i][0] > 0 and atoms[i].flag==1:
            v_a[i][0] = -v_a[i][0]#in pipe1
        if (p_a[i][0]) <= -(pipe1.length)/2 + size and p_a[i][0]*v_a[i][0] > 0 and atoms[i].flag==2:
            v_a[i][0] = -v_a[i][0]#in pipe2 with left side
        if (p_a[i][0]) >= pipe2.length-(pipe1.length)/2 - size and p_a[i][0]*v_a[i][0] > 0 and atoms[i].flag==2:
            v_a[i][0] = -v_a[i][0]#in pipe2 with right side
            
        if (p_a[i][0]) <=(container1.length/2+pipe3.length+container2.length/4+size) and p_a[i][0]*v_a[i][0] < 0 and atoms[i].flag==5:
            v_a[i][0] = -v_a[i][0]#in pipe4 with left side
        if (p_a[i][0]) >=(container1.length/2+pipe3.length+3*container2.length/4-size) and p_a[i][0]*v_a[i][0] > 0 and atoms[i].flag==5:
            v_a[i][0] = -v_a[i][0]#in pipe4 with right side
            
        ###pipe1 and pipe2 and pipe3 and pipe4(z)
        if abs(p_a[i][2]) >= (pipe1.width)/2 - size and p_a[i][2]* v_a[i][2] > 0 and atoms[i].flag==1:
            v_a[i][2] = -v_a[i][2] #in pipe1
        if abs(p_a[i][2]) >= (pipe2.width)/2 - size and p_a[i][2]* v_a[i][2] > 0 and atoms[i].flag==2:
            v_a[i][2] = -v_a[i][2] #in pipe2
        if abs(p_a[i][2]) >= (pipe3.width)/2 - size and p_a[i][2]* v_a[i][2] > 0 and atoms[i].flag==3:
            v_a[i][2] = -v_a[i][2] #in pipe3
        if abs(p_a[i][2]) >=(pipe4.width)/2 - size and p_a[i][2]* v_a[i][2] > 0 and atoms[i].flag==5:
            v_a[i][2] = -v_a[i][2] #in pipe4
        ###container1 and container2(z)
        if abs(p_a[i][2]) >= (container1.width)/2 - size and p_a[i][2]* v_a[i][2] > 0 and atoms[i].flag==0:
            v_a[i][2] = -v_a[i][2] #in container1   
        if abs(p_a[i][2]) >= (container2.width)/2 - size and p_a[i][2]* v_a[i][2] > 0 and atoms[i].flag==4:
            v_a[i][2] = -(1-e)*v_a[i][2] #in container2  
            
###carnot refrigerator            
##############################################################################    
    rectangle1.pos.x=-Length/2+pipelength+8*L
    rectangle2.pos.x=Length/2+pipelength+8*L
    p_a2 += v_a2*dt # calculate new positions for all atoms
    
    v_W=dV/Height/Width/dt#L/20000.0/dt
   
    if stage==0 :
        Length+=v_W*dt
        e2=0
    if Length>=3*L and stage==0 : 
        stage=1
        T_initial=T
    if stage==1 : 
        Length+=v_W*dt
         # adiabatic compression
    if Length>=4*L and stage==1 : 
        stage=2 # stop the walls from moving
        v_W*=(-1)
    if stage==2 :
        v_W*=(-1)
        Length+=v_W*dt
        e2=0 
    if T>=298.0 and stage==2 :
        stage=3 # stop the walls from moving
        T_initial=298.0
    if stage==3 : 
        v_W*=(-1)
        Length+=v_W*dt
        
    if Length<=2*L and stage==3 : 
        stage=0
    
    
    for i in range(N2): atoms2[i].pos = vector(p_a2[i, 0], p_a2[i, 1], p_a2[i, 2]) # to display atoms at new positions
    observation.plot(data = np.sqrt(np.sum(np.square(v_a2),-1))) # plot the speed distribution of atoms
     
### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a2-p_a2[:,np.newaxis] # array for vector from one atom to another atom for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) # distance array between atoms for all pairs of atoms
    hit = np.less_equal(rmag,2*size)-np.identity(N2) # if smaller than 2*size meaning these two atoms might hit each other
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() # change hit to a list
    
    for ij in hitlist: # i,j encoded as i*Natoms+j
        i, j = divmod(ij,N2) # atom pair, i-th and j-th atoms, hit each other
        hitlist.remove(j*N2+i) # remove j,i pair from list to avoid handling the collision twice
        if sum((p_a2[i]-p_a2[j])*(v_a2[i]-v_a2[j])) < 0 : # only handling collision if two atoms are approaching each other
            v_a2[i], v_a2[j] = vcollision(p_a2[i], p_a2[j], v_a2[i], v_a2[j]) # handle collision

#find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N2):
        if abs(p_a2[i][0]-8*L-pipelength) >= Length/2 - size and (p_a2[i][0]-8*L-pipelength)*v_a2[i][0] > 0 :
            if (p_a2[i][0]-8*L-pipelength)>0: # collisions between atoms and the right wall
                p+=2*m*(abs(v_a2[i][0])+v_W/2)
                v_a2[i][0]=(1+e2)*(-v_a2[i][0]+v_W)
            elif (p_a2[i][0]-8*L-pipelength)<0: # collisions between atoms and the left wall
                p+=2*m*(abs(v_a2[i][0])+v_W/2)
                v_a2[i][0]=(1+e2)*(-v_a2[i][0]-v_W)
        if abs(p_a2[i][1]-(pipelength+Height/2)) >= Height/2 - size and (p_a2[i][1]-(pipelength+Height/2))*v_a2[i][1] > 0 :
            v_a2[i][1] = - (1+e2)*v_a2[i][1]
            p+=2*m*abs(v_a2[i][1])
        if abs(p_a2[i][2]) >= Width/2 - size and p_a2[i][2]*v_a2[i][2] > 0 :
            v_a2[i][2] = - (1+e2)*v_a2[i][2]
            p+=2*m*abs(v_a2[i][2])
            
    if flag==50: # after 50*dt
        flag=0 # initialize
        T=0 # temperature
        area=(Length*Width+Length*Height+Height*Width)*2 # area of the six walls
        P=p/50/dt/area # pressure = total momentum divided by time and area
        V=Length*Height*Width #volume
        for i in range(N2): # total kinetic energy = 3/2*N*k*T
            T+=(2.0/3.0)*0.5*m*(v_a2[i,0]**2+v_a2[i,1]**2+v_a2[i,2]**2)/N2/k
        e2=0.0001*(T_initial-T)
        print(Length/L,stage,e2,T,dV,v_W)
        print('T:',T,end=' ') # print current macrostates
        print('p:',P,end=' ')
        print('V:',V,end=' ')
        print('P*V=',P*V,end=' ') # ideal gas law : pV=NkT
        print('N*k*T=',N2*k*T,end=' ')
        print('P*(V**gamma)=',P*(V**gamma),'\n') # in adiabatic process
        pv_diagram.plot(pos=(V,P))
        p=0 # initialize
