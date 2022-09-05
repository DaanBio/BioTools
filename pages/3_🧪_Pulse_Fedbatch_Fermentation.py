import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import math

col1, col2, col3, =st.columns(3)


tf=col1.number_input('The time that the reactor is running (hours)',min_value=1,max_value=240, value=48, step=1)
t_act_S=col1.number_input('The interval at which substrate is added to the reactor (hours)',min_value=0.,max_value=240., value=2., step=0.1)
Fin_S=col1.number_input('The amount of volume added to the reactor when pulse-feeding (ml)',min_value=0.,max_value=100000., value=1., step=0.1)
K_Pace=col3.number_input('The inhibition constant of growth inhibition molecule 1 (kg/m^3)',min_value=0.,max_value=400., value=10., step=0.1)
K_Peth=col3.number_input('The inhibition constant of growth inhibition molecule 2 (kg/m^3)',min_value=0.,max_value=400., value=1., step=0.1)
Ypco2=col3.number_input('The Yield of CO2 on the substrate (kg_CO2/kg_S)',min_value=0.,max_value=10., value=1., step=0.01)
C_X0 = col1.number_input('Biomass concentration at t=0 (g/L)',min_value=0.,max_value=100., value=0.1, step=0.01)  # gX/L or 1e-4 mass fraction
C_S0 = col2.number_input('Substrate concentration at t=0 (g/L)',min_value=0., max_value=530.,value=0., step=0.1)  # gS/L or 0.1 mass fraction
Y = col3.number_input('Biomass yield over substrate(kgX/kgS)',min_value=0., max_value=1.,value=0.4, step=0.01)   # kgX/kgS
mu_max_hour = col1.number_input('Maximum growth rate (1/h)',min_value=0.,value=0.1, step=0.01)   # 1/h
k_d_hour = col1.number_input('Death rate (1/h)',min_value=0., value=0.01, step=0.001) # 1/h
K_S = col2.number_input('Substrate specificity constant Ks (kgS/m3)',min_value=0., value=3., step=0.01)   # kgS/m3
C_S_in = col2.number_input('Substrate concentration in the feed (g/L)',min_value=0., max_value=530.,value=250., step=0.1)# kg/m3 (pure glucose in)
V_0 = col2.number_input('Volume of reactor at t=0 (ml)',min_value=0.1,value=200., step=0.1)   # ml
F_0 = 0     # ml/s
activation = 0
set_point = 0.25     # kg/m3
Yp = col3.number_input('Yield of product 1 over substrate (kgP1/kgS)',min_value=0., max_value=1.,value=0.2, step=0.01) #kgP/kgS
Ypace =col3.number_input('Yield of product 2 over substrate (kgP2/kgS)',min_value=0., max_value=1.,value=0.1, step=0.01) #kgP/kgS
C_P0 = col2.number_input('Concentration of product 1 at t=0 (g/L)',min_value=0.,value=0., step=0.1) #gP/L
C_Pace0 = col2.number_input('Concentration of product 2 at t=0 (g/L)',min_value=0.,value=0., step=0.1)  #gP/L

t0 = 0
dt = 1
tf = tf * 3600
molarmass_acid = 60.05
Ka_acid = 1.8 * 10 ** -5
Kb_ace = 5.6 * 10 ** -10
C_X = C_X0
C_S = C_S0
C_P = C_P0
C_Pace = C_Pace0
V_co2 = 0
C_co2 = 0
t = t0
V = V_0 / 10 ** 6
F = F_0 / 10 ** 6
F_S = F_0 / 10 ** 6
mu_max = mu_max_hour / 3600
k_d = k_d_hour / 3600
activation = 0
check = -1
warning = 0
density_co2 = 1.808  # kg/m3
KH = 0.037 * 44  # g/l*atmosphere

C_X_v, C_S_v, t_v, V_v, F_v, F_S_v, C_P_v, C_Pace_v, V_co2_v, C_co2_v = [C_X], [C_S], [t], [V], [F], [F_S], [C_P], [
    C_co2], [C_Pace], [V_co2]
mu_v, rX_v, rS_v, rXg_v, rXd_v = [], [], [], [], []
while t < tf:
    check = check + 1
    # Check if there is activation of the pump1
    if t % (t_act_S * 3600) == 0:
        activation = 1
        F_S = Fin_S / 10 ** 6  # kg/m3/s
    else:
        activation = 0
        F_S = 0
    F = F_S

    # Compute rates
    mu = mu_max * C_S / (C_S + K_S) * (K_Pace / (C_Pace + K_Pace)) * (K_Peth / (C_P + K_Peth))
    rXg = mu * C_X
    rXd = k_d * C_X
    rX = rXg - rXd
    rS = -1 / Y * rXg
    rP = -Yp * rS
    rPace = -Ypace * rS
    rPco2 = -Ypco2 * rS

    # Compute derivatives
    dXdt = rX - F / V * C_X
    dSdt = F / V * (C_S_in - C_S) + rS
    dVdt = F
    dPdt = rP - F / V * C_P
    dPacedt = rPace - F / V * C_Pace
    dPco2dt = rPco2

    # Update results
    t += dt
    C_X += dXdt * dt
    C_S += dSdt * dt
    V += dVdt * dt
    C_P += dPdt * dt
    C_Pace += dPacedt * dt
    C_co2 += dPco2dt * dt
    if C_co2 < KH:
        V_co2 += 0
    else:
        V_co2 += dPco2dt / density_co2

    # C_NAOH+=dNAOHdt*dt
    # C_ph=C_Pace/molarmass_acid
    # OH-
    # C_pH_NAOH=C_NAOH/molarmass_NAOH
    # H3O+
    # sqrt1=math.sqrt(Ka_acid*C_ph)
    # sqrt2=math.sqrt(Kb_ace*C_ph)

    # if C_ph==0 and C_NAOH==0:
    #    pH=7
    # elif C_ph==0 and C_NAOH!=0:
    #    pH=14-math.log10(C_pH_NAOH)
    # elif C_ph!=0 and C_NAOH==0:
    #    pH=-1*math.log10(sqrt1)
    # elif C_pH_NAOH>C_ph:
    #    pH=14-1*math.log10(C_pH_NAOH-C_ph)
    # elif C_pH_NAOH==C_ph:
    #    pH=14-math.log10(sqrt)
    # else:
    #    pH=-1*math.log10(Ka_acid)+math.log10(C_pH_NAOH/sqrt1)

    # if C_pH_NAOH>C_ph:
    #    warning=warning+1
    # else:
    #    pass

    # Store results
    C_X_v.append(C_X)
    C_S_v.append(C_S)
    C_P_v.append(C_P)
    C_Pace_v.append(C_Pace)
    V_co2_v.append(V_co2)
    C_co2_v.append(C_co2)
    # C_NAOH_v.append(C_NAOH)
    # pH_v.append(pH)
    t_v.append(t)
    V_v.append(V * 10 ** 6)
    F_v.append(F * 10 ** 6)
    F_S_v.append(F_S * 10 ** 6)
    # F_NAOH_v.append(F_NAOH*10**6)
    mu_v.append(mu)
    rS_v.append(rS)
    rX_v.append(rX)
    rXg_v.append(rXg)
    rXd_v.append(rXd)

# mu_v=mu_v[::3600]

t_plot = [t / 24 / 3600 for t in t_v]
mu_plot = [mu * 3600 for mu in mu_v]
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
ax1.plot(t_plot, C_S_v, 'r-', t_plot, C_X_v, t_plot, C_P_v, t_plot, C_Pace_v, 'g-', linewidth=2)
ax1.legend(['Concentration glucose', 'Concentration biomass', 'Concentration ethanol', 'Concentration acetate'])
ax1.set(xlabel='Time (days)', ylabel='Concentration (kg/$m^3$)',
        title='Concentration of different compounds in N.I.N.O.')

# ax2.plot(t_plot, pH_v,linewidth=2)
# ax2.plot([0,tf/(3600*24)],[6,6], color='r')
# ax2.plot([0,tf/(3600*24)],[8,8], color='r')
# ax2.set(xlabel='Time (days)', ylabel='pH', title='pH in the reactor over time'+title_ph)
# ax2.legend(['pH of reactor','pH threshold'])

ax2.plot(t_plot, V_co2_v, linewidth=2)
ax2.set(xlabel='Time (days)', ylabel='Total off gas volume (ml)', title='Total offgass volume')

ax3.plot(t_plot[:-1], mu_plot, 'y-', linewidth=2)
ax3.set(xlabel='Time (days)', ylabel='$\mu$ (1/h)', title='Growth rate over time')

ax4.set_title('Volume and flowrate of N.I.N.O')
ax4.scatter(t_plot, F_v, linewidth=2, label='Flowrate')
ax4.set(xlabel='Time (days)', ylabel='Flowrate (ml/s)')
ax4.legend(loc=6)
ax5 = ax4.twinx()
ax5.plot(t_plot[1:], V_v[1:], 'g-', linewidth=2, label='Volume')
ax5.set_ylabel('Volume (mL)')
ax5.ticklabel_format(axis="y")
ax5.legend(loc=5)
# lns = lns1+lns2
# labs = [l.get_label() for l in lns]
# ax4.legend(lns,labs,loc=0)

st.pyplot(fig)
