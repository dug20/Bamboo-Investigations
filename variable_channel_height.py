import bamboo as bam
import bamboo.cooling as cool
import bamboo.materials

import numpy as np
import matplotlib.pyplot as plt
import pypropep as ppp
import bamboo.plot
import thermo
import time

'''Engine dimensions'''
Ac = np.pi*0.1**2               #Chamber cross-sectional area (m^2)
L_star = 1.5                    #L_star = Volume_c/Area_t
wall_thickness = 2e-3

'''Chamber conditions'''
pc = 15e5               #Chamber pressure (Pa)
p_tank = 50e5           #Tank / inlet coolant stagnation pressure (Pa) - used for cooling jacket
mdot = 5.4489           #Mass flow rate (kg/s)
p_amb = 1.01325e5       #Ambient pressure (Pa). 1.01325e5 is sea level atmospheric.
OF_ratio = 3.5          #Oxidiser/fuel mass ratio

'''We want to investigate adding water to the isopropyl alcohol'''
water_mass_fraction = 0.10  #Fraction of the fuel that is water, by mass

'''Coolant jacket'''
wall_material = bam.materials.CopperC700
mdot_coolant = mdot/(OF_ratio + 1) 
inlet_T = 298.15                    #Coolant inlet temperature

'''Get combustion properties from pypropep'''
#Initialise and get propellants
ppp.init()
e = ppp.Equilibrium()
ipa = ppp.PROPELLANTS['ISOPROPYL ALCOHOL']
water = ppp.PROPELLANTS['WATER']
n2o = ppp.PROPELLANTS['NITROUS OXIDE']

#Add propellants by mass fractions (note the mass fractions can add up to more than 1)
e.add_propellants_by_mass([(ipa, 1-water_mass_fraction), 
                        (water, water_mass_fraction), 
                        (n2o, OF_ratio)])

#Adiabatic combustion using chamber pressure                      
e.set_state(P = pc/1e5, type='HP')                      

gamma = e.properties.Isex   #I don't know why they use 'Isex' for gamma. 
cp = 1000*e.properties.Cp   #Cp is given in kJ/kg/K, we want J/kg/K
Tc = e.properties.T

'''Choose the models we want to use for transport properties of the coolant and exhaust gas'''
thermo_coolant = thermo.mixture.Mixture(['isopropanol', 'water'], ws = [1 - water_mass_fraction, water_mass_fraction])
thermo_gas = thermo.mixture.Mixture(['N2', 'H2O', 'CO2'], zs = [e.composition['N2'], e.composition['H2O'], e.composition['CO2']])   

gas_transport = cool.TransportProperties(model = "thermo", thermo_object = thermo_gas, force_phase = 'g')
coolant_transport = cool.TransportProperties(model = "thermo", thermo_object = thermo_coolant, force_phase = 'l')

'''Create the engine object'''
perfect_gas = bam.PerfectGas(gamma = gamma, cp = cp)    #Gas for frozen flow
chamber_conditions = bam.ChamberConditions(pc, Tc, mdot)
nozzle = bam.Nozzle.from_engine_components(perfect_gas, chamber_conditions, p_amb, type = "rao", length_fraction = 0.8)
white_dwarf = bam.Engine(perfect_gas, chamber_conditions, nozzle)
chamber_length = L_star*nozzle.At/Ac

'''Add the cooling system to the engine'''
white_dwarf.add_geometry(chamber_length, Ac, wall_thickness)
white_dwarf.add_exhaust_transport(gas_transport)
white_dwarf.add_ablative(bam.materials.Graphite, bam.materials.CopperC700, regression_rate = 0.0033e-3, xs = [-100, 100], ablative_thickness = None)

def get_data(channel_height):
    white_dwarf.add_cooling_jacket(wall_material, inlet_T, p_tank, coolant_transport, mdot_coolant, configuration = "vertical", channel_height = channel_height)

    '''Run the heating analysis'''
    #print(f"Sea level thrust = {white_dwarf.thrust(1e5)/1000} kN")
    #print(f"Sea level Isp = {white_dwarf.isp(1e5)} s")

    return white_dwarf.steady_heating_analysis(number_of_points = 250, to_json = False)

heights = np.linspace(0.5e-3, 5e-3, 20)
isps = np.zeros(len(heights))
max_inner_wall_Ts = np.zeros(len(heights))
pressure_drop = np.zeros(len(heights))

for i in range(len(heights)):
    data = get_data(heights[i])
    max_inner_wall_Ts[i] = np.amax(data["T_wall_inner"])
    pressure_drop[i] = data["p_coolant"][0] - data["p_coolant"][-1] 
    print(f"i = {i}, max T = {max_inner_wall_Ts[i]} K, pressure drop = {pressure_drop[i]/1e5} bar")

fig, T_axs = plt.subplots()
p_axs = T_axs.twinx() 

T_axs.plot(heights/1e-3, max_inner_wall_Ts - 273.15, label = "Temperature (°C)", color = 'orange')
p_axs.plot(heights/1e-3, pressure_drop/1e5, label = "Pressure drop (bar)", color = 'green')

T_axs.grid()
T_axs.set_xlabel("Cooling channel height (mm)")
T_axs.set_ylabel("Maximum inner liner temperature (°C)")
p_axs.set_ylabel("Pressure drop (bar)")
T_axs.legend(loc='lower right')
p_axs.legend(loc='upper left')
plt.title("Vertical channels, copper walls, graphite insert")
plt.show()
