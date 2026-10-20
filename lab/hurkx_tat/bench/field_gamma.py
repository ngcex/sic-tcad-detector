import numpy as np, devsim
from etna.api.device import DeviceConfig, device_kwargs
from etna.core.dark_current import create_dark_current_device, dark_current_at_bias
hbar=1.0546e-34; m0=9.109e-31; q=1.602e-19; kB=8.617e-5
for T in (150.,300.,450.):
    FG = np.sqrt(24*0.25*m0*(kB*T*q)**3)/(q*hbar)/100
    print("T",T,"F_Gamma V/cm %.3e"%FG)
cfg=DeviceConfig()
info=create_dark_current_device(N_t=cfg.N_t,E_t=cfg.E_t,**device_kwargs(cfg))
for V in (-10.,-60.,-100.):
    r=dark_current_at_bias(info,V,area=cfg.area_cm2)
    E=np.array(devsim.get_node_model_values(device=info["device_name"],region=info["region_name"],name="E_field_node"))
    x=np.array(devsim.get_node_model_values(device=info["device_name"],region=info["region_name"],name="x"))
    FG=np.sqrt(24*0.25*m0*(kB*300*q)**3)/(q*hbar)/100
    f=E/FG; G=2*np.sqrt(3*np.pi)*f*np.exp(f**2)
    print(V, "Emax %.3e"%E.max(), "Gmax %.3f"%G.max(), "I %.4e"%r["I_total"], "ISRH %.3e ITAT %.3e"%(r["I_SRH"],r["I_TAT"]))
    m = E>3.5e4
    print("   mean Gamma over E>Eref nodes %.3f"%G[m].mean(), "x-range of gen", x[m].min()*1e4, x[m].max()*1e4)
