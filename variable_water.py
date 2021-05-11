import bamboo as bam
import bamboo.cooling as cool
import bamboo.materials

import numpy as np
import matplotlib.pyplot as plt
import pypropep as ppp
import bamboo.plot
import thermo
import time

def white_dwarf_cooling_data(water_mass_fraction):
    '''Engine dimensions'''
    Ac = np.pi*0.1**2               #Chamber cross-sectional area (m^2)
    L_star = 1.5                    #L_star = Volume_c/Area_t
    wall_thickness = 2e-3

    '''Chamber conditions'''
    pc = 15e5               #Chamber pressure (Pa)
    p_tank = 20e5           #Tank / inlet coolant stagnation pressure (Pa) - used for cooling jacket
    mdot = 5.4489           #Mass flow rate (kg/s)
    p_amb = 1.01325e5       #Ambient pressure (Pa). 1.01325e5 is sea level atmospheric.
    OF_ratio = 3.5          #Oxidiser/fuel mass ratio

    '''Coolant jacket'''
    wall_material = bam.materials.CopperC700
    mdot_coolant = mdot/(OF_ratio + 1) 
    inlet_T = 298.15                    #Coolant inlet temperature

    '''Get combustion properties from pypropep'''
    #Initialise and get propellants
    ppp.init()
    e =  ppp.Equilibrium()
    p =  ppp.ShiftingPerformance()

    ipa = ppp.PROPELLANTS['ISOPROPYL ALCOHOL']
    water = ppp.PROPELLANTS['WATER']
    n2o = ppp.PROPELLANTS['NITROUS OXIDE']

    #Add propellants by mass fractions (note the mass fractions can add up to more than 1)
    e.add_propellants_by_mass([(ipa, 1), 
                            (water, water_mass_fraction), 
                            (n2o, OF_ratio)])

    p.add_propellants_by_mass([(ipa, 1), 
                            (water, water_mass_fraction), 
                            (n2o, OF_ratio)])

    #Adiabatic combustion
    e.set_state(P = pc/p_amb, type = 'HP')  

    #Set chamber pressure and exit pressure in atmospheres                           
    p.set_state(P = pc/p_amb, Pe = 1)           

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

    #Spiral channels
    #white_dwarf.add_cooling_jacket(wall_material, inlet_T, p_tank, coolant_transport, mdot_coolant, 
    #                               configuration = "spiral", channel_shape = "semi-circle", channel_width = 0.020)

    #Or vertical channels
    white_dwarf.add_cooling_jacket(wall_material, inlet_T, p_tank, coolant_transport, mdot_coolant, configuration = "vertical", channel_height = 0.001)

    '''Run the heating analysis'''
    cooling_data = white_dwarf.steady_heating_analysis(number_of_points = 250, to_json = False)

    return cooling_data, white_dwarf.isp(p_amb)


water_mass_fractions = np.linspace(0, 0.8, 20)
isps = np.zeros(len(water_mass_fractions))
max_inner_wall_Ts = np.zeros(len(water_mass_fractions))

for i in range(len(water_mass_fractions)):
    cooling_data, isps[i] = white_dwarf_cooling_data(water_mass_fractions[i])
    max_inner_wall_Ts[i] = np.amax(cooling_data["T_wall_inner"])
    print(f"i = {i}, max T = {max_inner_wall_Ts[i]} K, isp = {isps[i]} s")

fig, T_axs = plt.subplots()
isp_axs = T_axs.twinx() 

T_axs.plot(water_mass_fractions, max_inner_wall_Ts - 273.15, label = "Temperature (°C)", color = 'orange')
isp_axs.plot(water_mass_fractions, isps, label = "Isp (s)", color = 'green')

for i in range(len(isps)):
    print(f"{water_mass_fractions[i]} {isps[i]}")

T_axs.grid()
T_axs.set_xlabel("Water 'proportional' added")
T_axs.set_ylabel("Max. liner temperature (°C)")
isp_axs.set_ylabel("Specific impulse (from bamboo) (s)")
T_axs.legend(loc='lower right')
isp_axs.legend(loc='upper left')
plt.title("Isp and maximum liner temperature against water mass fraction")
plt.show()
