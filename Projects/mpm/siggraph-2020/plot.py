import math
import sys
import matplotlib.pyplot as plt
import os
import numpy as np
import random

#############################################################
## Plot time && time step && velocity
#############################################################
lim="3.402823466385289e+38"
path=sys.argv[1]
if path[-1] == r'/':
    path=path[:-1]
velocity_string=os.popen('cat '+path+r'/common/log.txt | grep "max velo" | sed "s/.*velocity: \([^<]*\)<.*/\1/g" ').read().split('\n')
velocity_string.pop()
if velocity_string[-1]=="Binary file (standard input) matches":
    velocity_string.pop()
velocity=list(map(float, velocity_string))
velocity=list(filter(lambda a: a != 0.0, velocity))
velocity=[(math.log10(x)) for x in velocity]


# predicted_dt_string=os.popen('cat '+path+r'/common/log.txt | grep "dx/sound" | sed "s/.*dx\/soundspeed: \([^<]*\)<.*/\1/g" ').read().split('\n')
# predicted_dt_string.pop()
# predicted_dt_string=list(filter(lambda a: a != lim, predicted_dt_string))
# if len(predicted_dt_string)>0 and predicted_dt_string[-1]=="Binary file (standard input) matches":
#     predicted_dt_string.pop()
# predicted_dt=list(map(float, predicted_dt_string))
# predicted_dt=[x * .9 for x in predicted_dt]
# predicted_dt=[math.log10(x) for x in predicted_dt ]
# if len(predicted_dt)==0:
#     print("no sound speed")

# actual_dt_string=os.popen('cat '+path+r'/common/log.txt | grep "actual dt" | sed "s/.*dt: \([^<]*\)<.*/\1/g" ').read().split('\n')
# actual_dt_string.pop()
# actual_dt_string=list(filter(lambda a: a != lim, actual_dt_string))
# if actual_dt_string[-1]=="Binary file (standard input) matches":
#     actual_dt_string.pop()
# actual_dt=list(map(float, actual_dt_string))
# actual_dt=[math.log10(x) for x in actual_dt ]
# fig.subplots_adjust(hspace=0.5)
# ax1.plot([x for x in range(len(predicted_dt))],predicted_dt)
# ax1.plot([x for x in range(len(actual_dt))],actual_dt)
# ax1.set_title(path+": dt, with blue the predicted and orange the actual")
plt.plot([x for x in range(len(velocity))],velocity)
title=''
for tit in sys.argv[2:]:
    title=title+' '+tit
plt.title(title)
plt.ylabel("Velocity Magnitude")
plt.xlabel("Time Step")
plt.show()
# ax2.set_title(path+": velocity")


#############################################################
## Plot sound speed wrt scale_E
#############################################################
# scale_e=[]
# sound_speed=[]
# sph_res="2k"
# sim_res=32
# for scale in [500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2300,2500,2700,3000,4000,6000]:
#     scale_e.append(scale)
#     ss=os.popen(r'cat Test_3d_13_'+sph_res+'_'+str(scale)+'_'+str(sim_res)+r'/common/log.txt | grep "max sound" | tail -1 | sed "s/.*max sound speed: \([^<]*\)<.*/\1/g"')
#     sound_speed.append(float(ss.read()[:-1]))
# plt.scatter(scale_e, sound_speed, c='red', label='2k_32')
# co=np.polyfit(np.sqrt(scale_e), sound_speed, 1)
# print(co)
# plt.plot(scale_e, np.poly1d(co)(np.sqrt(scale_e)))
# plt.legend()
# plt.show()
# exit(0)



# demand=sys.argv[1]
# random.seed()
# if demand=='sph':
# # ############################################################
# # # Normalize stable dt  ==== sph res is variable
# # ############################################################
#     sim_res=sys.argv[2]
#     scale_e=sys.argv[3]
#     fig=plt.figure()
#     ax=plt.subplot(2,1,1)
#     ax1=plt.subplot(2,1,2)
#     sph_resv=[]
#     stable_dt=[]
#     for sph_res in ['192','390','726','1k','2k','5k','9k']:
#         sph_resv.append(int(sph_res[0]))
#         ss=os.popen(r'head Test_3d_13_'+sph_res+'_'+str(scale_e)+'_'+str(sim_res)+r'/common/log.txt | grep command | sed "s/.*min_dt \([^-]*\) -max_dt.*/\1/g"')
#         stable_dt.append(float(ss.read()))
#     ax.scatter(sph_resv,stable_dt,c=[[random.random(),random.random(),random.random()]])
#     sph_resv=[192,390,726,1000,2000,5000,9000]
#     ax1.plot(sph_resv,np.asarray(stable_dt)/np.asarray(np.log(sph_resv)))
#     print(sph_resv)
#     print(stable_dt)
#     plt.xlabel("Sph Res")
#     plt.ylabel("Stable dt")
#     plt.title("Stable dt wrt sph res")
#     plt.legend()
#     plt.show()
# elif demand=='sim':
# ############################################################
# # Normalize stable dt  ==== sim res is variable
# ############################################################
#     sph_res=sys.argv[2]
#     scale_e=sys.argv[3]
#     fig=plt.figure()
#     ax=plt.subplot(2,1,1)
#     ax1=plt.subplot(2,1,2)
#     sim_resv=[]
#     stable_dt=[]
#     for sim_res in [32,48,64,96,128]:
#         sim_resv.append(sim_res)
#         ss=os.popen(r'head Test_3d_13_'+sph_res+'_'+str(scale_e)+'_'+str(sim_res)+r'/common/log.txt | grep command | sed "s/.*min_dt \([^-]*\) -max_dt.*/\1/g"')
#         stable_dt.append(float(ss.read()))
#     ax.scatter(sim_resv,stable_dt,c=[[random.random(),random.random(),random.random()]],label=str(''))
#     co=np.polyfit(sim_resv,stable_dt,2)
#     print(co)
#     ax1.plot(sim_resv,np.asarray(stable_dt)/np.asarray(np.log(sim_resv)))

#     plt.xlabel("Sim Res")
#     plt.ylabel("Stable dt")
#     plt.title("Stable dt wrt sim res")
#     plt.legend()
#     plt.show()
# else:
# ############################################################
# # Normalize stable dt  ==== scale_e is variable
# ############################################################
#     # sound speed = Asqrt(scale_e)+b
#     # dt = 1/ (Asqrt(scale_e)+b)
#     sph_res=sys.argv[2]
#     sim_res=sys.argv[3]
#     fig=plt.figure()
#     ax=plt.subplot(2,1,1)
#     ax1=plt.subplot(2,1,2)
#     scale_ev=[]
#     stable_dt=[]
#     for scale_e in [500,1000,2000]:
#         scale_ev.append(scale_e)
#         ss=os.popen(r'head Test_3d_13_'+sph_res+'_'+str(scale_e)+'_'+str(sim_res)+r'/common/log.txt | grep command | sed "s/.*min_dt \([^-]*\) -max_dt.*/\1/g"')
#         stable_dt.append(float(ss.read()))
#     ax.scatter(scale_ev,stable_dt,c=[[random.random(),random.random(),random.random()]],label=str(''))
#     # co=np.polyfit(np.sqrt(scale_ev), 1/np.asarray(stable_dt), 1)
#     # print(co)
#     ax1.plot(scale_ev,np.asarray(stable_dt)/np.asarray(np.power(scale_ev,-0.5)))
#     # plt.plot(scale_ev,stable_dt*np.poly1d(co)(np.sqrt(scale_ev))/np.poly1d(co)(np.sqrt(500)))
#     plt.xlabel("Scale_E")
#     plt.ylabel("Stable dt")
#     plt.title("Stable dt wrt Scale_E")
#     plt.legend()
#     plt.show()