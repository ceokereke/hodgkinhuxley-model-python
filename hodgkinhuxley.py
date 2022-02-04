import numpy as np
import array as a
import matplotlib.pyplot as plt
import math 


#===simulation time===
simulationTime =100 # in milliseconds
deltaT=0.01 #space out simulation time

t= np.arange(0,simulationTime+deltaT,deltaT) #spaced out simulation time

#===specify the external current I===
#Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)
currentLevels = 50

#Set externally applied current across time
#Here, first 500 timesteps are at current of 50, next 1500 timesteps at
#current of zero (resets resting potential of neuron), and the rest of
#timesteps are at constant current

I =[currentLevels]*len(t)
I[500:1999]=[0]*1500

#Comment out the above line and uncomment the line below for constant current, and observe effects on voltage timecourse
#I = [currentLevels] *len(t)

#===constant parameters===
#All of these can be found in Table 3

gbar_K=36
gbar_Na=120 
g_L=.3
E_K = -12 
E_Na=115
E_L=10.6
C=1


#===set the initial states===
V=[0] # Baseline Volatge


 #===set the initial states===
	#Equation 12
alp_n = 0.01 * ( (10-V[0]) / (math.exp((10-V[0])/10)-1))
        #Equation 13
bet_n = 0.125 * math.exp(-V[0]/80) 
	# Equation 20
alp_m = 0.1*( (25-V[0]) / (math.exp((25-V[0])/10)-1) )

	#Equation 21
bet_m = 4 * math.exp(-V[0]/18)
	#Equation 23
alp_h = 0.07* math.exp(-V[0]/20)
	#Equation 24
bet_h = 1/(math.exp((30-V[0])/10)+1) 
	
n = [alp_n/(alp_n+bet_n)] # Equation 9
m = [alp_m/(alp_m+bet_m)] # Equation 18
h = [alp_h/(alp_h+bet_h)] # Equation 18

alpha_n = []

beta_n = []

alpha_m = []

beta_m =  []

alpha_h = []

beta_h = []


for i in range(len(t)-1):

     #---calculate the coefficients---%
    #Equations here are same as above, just calculating at each time step
        alp_n = 0.01 * ( (10-V[i]) / (math.exp((10-V[i])/10)-1) )
        alpha_n.append(alp_n)
        
        bet_n = 0.125*math.exp(-V[i]/80)
        beta_n.append(bet_n)
        
        alp_m= 0.1*( (25-V[i]) / (math.exp((25-V[i])/10)-1) )
        alpha_m.append(alp_m)
        
        bet_m = 4*math.exp(-V[i]/18)
        beta_m.append(bet_m)
        
        alp_h = 0.07*math.exp(-V[i]/20);
        alpha_h.append(alp_h)
       
        bet_h = 1/(math.exp((30-V[i])/10)+1);
        beta_h.append(bet_h)
   
    #%---calculate the currents---%
        I_Na = (m[i]**3) * gbar_Na * h[i] * (V[i]-E_Na) #; %Equations 3 and 14
        I_K = (n[i]**4) * gbar_K * (V[i]-E_K)#; %Equations 4 and 6
        I_L = g_L *(V[i]-E_L) #; %Equation 5
        I_ion = I[i] - I_K - I_Na - I_L
   
   
    #%---calculate the derivatives using Euler first order approximation---%
        Vv = V[i] + deltaT*I_ion/C
        
        V.append(Vv)
        
        nn = n[i] + deltaT*(alpha_n[i] *(1-n[i]) - beta_n[i] * n[i])#; %Equation 7
        n.append(nn)
        
        mm = m[i] + deltaT*(alpha_m[i] *(1-m[i]) - beta_m[i] * m[i])#; %Equation 15
        m.append(mm)
        
        hh = h[i] + deltaT*(alpha_h[i] *(1-h[i]) - beta_h[i] * h[i])#; %Equation 16	
        h.append(hh)

#Set resting potential to -70mv

v = np.array(V)
Vo = v-70
print (Vo)
#%===plot Voltage===%
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 5))


ax[0].plot(t, Vo, label='Voltage')
ax[0].legend(loc="upper right")
ax[0].set_xlabel('time (ms)')
ax[0].set_ylabel('Voltage (mv)')
ax[0].set_title("Voltage over Time in Simulated Neuron")

na = np.array(n)

ma = np.array(m)

ha = np.array(h)

gKn=float(gbar_K)*na**4

gNam= float(gbar_Na)*(ma**3)*ha

ax[1].plot(t,gKn, label ='Conductance for Potassium');
ax[1].plot(t,gNam, label ='Conductance for Sodium');
ax[1].legend(loc="upper right")
ax[1].set_xlabel('time (ms)')
ax[1].set_ylabel('Conductance')
ax[1].set_title("Conductance for Potassium and Sodium Ions in Simulated Neuron")
#p2 = plot(t,gbar_Na*(m.^3).*h,'r','LineWidth',2);
#legend([p1, p2], 'Conductance for Potassium', 'Conductance for Sodium')
#ylabel('Conductance')
#xlabel('time (ms)')
#title('Conductance for Potassium and Sodium Ions in Simulated Neuron')
plt.show()








