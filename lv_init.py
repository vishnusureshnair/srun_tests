import os
import numpy as np



path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
print(os.getcwd)

#Import Mission String

MS = np.genfromtxt('input_s1',names=True,dtype=float,delimiter='\t')


#--- Atmospheric data---
atmdat = np.genfromtxt('std_atm_data_1000.txt',names=True,dtype=float,delimiter='\t')
# (Altitude(km), Temperature(K), Pressure(Pa), Density(kg/m3))
alt_data = [atmdat[x][0] for x in range(len(atmdat))]
T_data = [atmdat[x][1] for x in range(len(atmdat))]
P_data = [atmdat[x][2] for x in range(len(atmdat))]
rho_data = [atmdat[x][3] for x in range(len(atmdat))]

delt = 0.01 # timestep, s
delt23 = 0.1 # timestep, s


# Constants
ftm=0.3048                  #Feet to meter conversion factor
gamma_air = 1.4
R_air = 287

#-------------------------------------Outputs-----------------------------------

pl = ['t','Q','h','am','vrm','fpa','r','aoa','steer[:,1]','event','mach','pr']

#-------------------------------Orbit parameters--------------------------------
r0 = [x*10**3 for x in [400,400]]               # [m] Orbit h
inclination = np.radians(90)                     # [deg] Orbit inclination

#------------------------------Vehicle parameters-------------------------------
mpl_req = 500                           # [kg] Payload mass
m_PLF = 150
m_PLA = 45
m_EB = 40
m_pl = mpl_req+m_PLA+m_EB+m_PLF
pload =[m_pl-m_PLF-m_PLA-m_EB,m_PLF,m_PLA,m_EB]
mis = [0, 0, 0]
Sp_imp = np.array([270, 310, 310])          # [s] Isp of stages
strf = np.array([0.15, 0.15, 0.15])  # [] Structural factors
m_lo = np.array([[], 50, 15])               # left out mass in kg
# segn = np.array([[], 20, 20])               # no. of pitchrate segments
# Thr = np.array([[], 20*1000*9.8055, 10*1000*9.8055])
tw_ratio = 2
dia = 2                                # [m] diameter of LV
area = np.pi*(dia/2)**2                    # [m2] CS area of LV
roll = 0
yaw = 0
pitch = 0
#-------------------------Vehicle Inputs---------------------------
#prate_in = ([[],np.radians(np.zeros(segn[1])),np.radians(np.zeros(segn[2]))])
Thr = np.array([[], 20*1000*9.8055, 4.2*1000*9.8055])
Ns = 3
#Thrust_in = [[],np.ones(segn[1])*20*1000*9.8055,np.ones(segn[2])*4.2*1000*9.8055]
#----------------------------Kickrate Maneuvre---------------------------------
#krate = radians(-1.5) #deg/s
kstart = 7.5 # s
kend = 20 # s
#prate_in = ([[],np.radians(np.random.rand(20)),np.radians(np.random.rand(20))])
#------------------------------Gravity Turn-------------------------------------
gtstart = 20 # s
gtend = 125.5 # s
#------------------------------Launch Pad parameters----------------------------
# Kulasekharapatnam
lat_L = np.radians(8.3989)                  # [rad], latitude
long_L = np.radians(78.0524)                # [rad], longitude
h0 = 8                             # [km] LP h above sea level (KSP)
# Sriharikotta
#lat_L = np.radians(13.7259)                 # [rad], latitude
#long_L = np.radians(80.2266)                # [rad], longitude
#h0 = 26.40                             # [km] LP h above sea level (SHAR)

#----------------------------- Earth Parameters --------------------------------
g0 = 9.80665 # acceleration due to gravity at surface of earth, m/s**2
Re = 2.0925741*ftm*(10**7) # Equitorial radius of earth in m
Rp = 2.0855590*ftm*(10**7) #Polar radius of earth in m
k_e = (Re/Rp)**2
Om_p = [0,0,7.2921159*(10**-5)] #Angular Velocity of earth in rad/sec
mu = 3.986004418*10**14                 # [m3/sec2] gravity parameter
j2 = 1.0823*10**-3                     # [] gravity parameter
f = (Re-Rp)/Re                         # [] Oblateness/Flattening
#-------------------------------Derived quantities------------------------------

apogee = max(r0)                        # [m] Apogee h
perigee = min(r0)                     # [m] Perigee h
k_e = (Re/Rp)**2                        
PhiL = np.arctan((1-f)**2*np.tan(lat_L));
Rs = Re/np.sqrt(1+(k_e-1)*(np.sin(PhiL))**2);
Ri = Rs+h0;    
ri = [Ri*x for x in [np.cos(lat_L)*np.cos(long_L),np.cos(lat_L)*np.sin(long_L),np.sin(lat_L)]];
vi = np.cross(Om_p,ri);
alti = np.linalg.norm(ri) - Rs
IGL = np.array( [[-np.sin(lat_L)*np.cos(long_L), -np.sin(lat_L)*np.sin(long_L), np.cos(lat_L)],
            [-np.sin(long_L), np.cos(long_L), 0],
            [-np.cos(lat_L)*np.cos(long_L), -np.cos(lat_L)*np.sin(long_L), -np.sin(lat_L)]])
#AzL = np.arcsin(np.cos(inclination)/np.cos(lat_L)); # As per computation

AzL = np.radians(179.9)
# Azl = radians(135) # Launch azimuth for SHAR polar launches

#----------------------------Trajectory Determination---------------------------

# An estimate of the velocity required, as the actual velocity required 
# depends on the final trajectory.

v_p = 0                                # [m/s] propulsive loss
margin_v = 3                       # [m/s] velocity margin for safety
v_sy = 0                               # [m/s] steering, yaw maneuver loss 
a = (apogee+perigee)/2+Re              # [m] semi-major axis of orbit
v_o = np.sqrt(mu*(2/(Re+perigee)-1/a))   # [m/s] orbital velocity
v_rot = np.sqrt(np.dot(vi[0], vi[0]))     # [m/s] Earth's rotational velocity
v_rgain = v_o - np.sqrt(np.square(v_rot+v_o*np.sin(AzL))+np.square(v_o*np.cos(AzL)))

# as per Loftus and Teixeira empirical relation for gravity loss
v_g = 81.006*np.square(tw_ratio)-667.62*(tw_ratio)+1505.4

# as per Loftus and Teixeira empirical relation for drag loss
v_d = -32.692*np.square(tw_ratio)+258.86*tw_ratio-226.57

# ve_req = v_o+v_g+v_d+v_p+v_sy+v_rgain+margin_v
ve_req = 10.09*1e3
# print('Estimated velocity required : %f km/s\n',ve_req/1000)

# print('\n\nReading Mach no. vs Coeff of Drag data...\n')
#f=open("data.txt", "r")
#if f.mode == 'r':
#    data = f.read()
data = np.genfromtxt('data.txt',
                    names=True,
                    dtype=None,
                    delimiter='  ')
mach_data = [data[x][0] for x in range(len(data))]
cd_data = [data[x][1] for x in range(len(data))]
#mach_data = append([0],mach_data)
#cd_data = append([0],cd_data)
# print('Aerodynamic Data assimilated!\n\n')