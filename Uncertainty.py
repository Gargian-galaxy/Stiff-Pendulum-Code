
import numpy as np
import pandas as pd

#measured constants
Mr=0.27187
Mp=0.74483
Lr=0.786
h=0.532

Ts=pd.read_csv('C:/Users/prana/OneDrive - Imperial College London/Labs/Stiff Pendulum/Stiff_pendulum_results.csv',usecols=['Ts'])
Ts=Ts.dropna()
Ts=Ts.drop(index=11)
Tr=pd.read_csv('C:/Users/prana/OneDrive - Imperial College London/Labs/Stiff Pendulum/Stiff_pendulum_results.csv',usecols=['Tr'])
Tr=Tr.drop(index=6)
Tr=Tr.dropna()
Tp=pd.read_csv('C:/Users/prana/OneDrive - Imperial College London/Labs/Stiff Pendulum/Stiff_pendulum_results.csv',usecols=['Tp'])
Tp=Tp.drop(index=4)
Tp=Tp.dropna()

Ts=Ts.to_numpy()
Ts=Ts.reshape(-1)
Tr=Tr.to_numpy()
Tr=Tr.reshape(-1)
Tp=Tp.to_numpy()
Tp=Tp.reshape(-1)

Tp=Tp/5
Tp_avg=np.mean(Tp)
Tr=Tr/5
Tr_avg=np.mean(Tr)
Ts=Ts/10
Ts_avg=np.mean(Ts)


#calculate g
Ir=(1/12)*(Mr*((Lr)*(Lr)))
Ip=(Ir*(Tp_avg/Tr_avg)*(Tp_avg/Tr_avg))

I=Ip+Mp*h*h
pi=np.pi
g=(4*pi*pi*I)/(Mp*h*Ts_avg*Ts_avg)
print(g)

#uncertainty is masses is 0.000005 and in lenght it is 0.0005
#uncertainties- sigma_x and error_x are both the error on mean
sigma_Lsquared=(Lr**2)*np.sqrt(2*((0.0005/Lr)**2))
sigma_Ir=np.sqrt((((1/12)*(Lr**2))**2)*(0.000005)**2+(sigma_Lsquared**2)*((1/12)*2*Lr*Mr)**2)           

sd_Ts_temp=0
for i in range(0,len(Ts)):
    sd_Ts_temp=(Ts[i]-Ts_avg)**2
sd_Ts=(1/(len(Ts)-1)*sd_Ts_temp)
error_Ts=sd_Ts/np.sqrt(len(Ts))

sd_Tp_temp=0
for i in range(0,len(Tp)):
    sd_Tp_temp=(Tp[i]-Tp_avg)**2
sd_Tp=(1/(len(Tp)-1)*sd_Tp_temp)
error_Tp=sd_Tp/np.sqrt(len(Tp))

sd_Tr_temp=0
for i in range(0,len(Tr)):
    sd_Tr_temp=(Tr[i]-Tr_avg)**2
sd_Tr=(1/(len(Tr)-1)*sd_Tr_temp)
error_Tr=sd_Tr/np.sqrt(len(Tr))

sigma_Ip=Ip*np.sqrt((error_Tp/Tp_avg)**2+(error_Tr/Tr_avg)**2+(sigma_Ir/Ir)**2)

sigma_Mphsquare=((Mp*h**2)*np.sqrt((0.000005/Mp)**2 + 2*((0.0005/h)**2)))

sigmaI=np.sqrt((sigma_Ip**2)+(sigma_Mphsquare**2))


error_g=np.sqrt(((((4*pi*pi)/(Mp*h*Ts_avg*Ts_avg))**2)*((sigmaI)**2))+((error_Ts)**2)*(((4*pi*pi*I)/(-2*Mp*h*Ts_avg*Ts_avg*Ts_avg))**2))
print(error_g)

print("g= %.2e +/- %.0e"
      %(g,error_g))