import rocketcea
import os
import numpy as np
import scipy
from scipy.optimize import minimize
import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib
print(matplotlib.matplotlib_fname())
from rocketcea.cea_obj import CEA_Obj
import pandas as pd
import math
import datetime
import rocketpy
from rocketpy import Environment, Rocket, Flight, Function
from rocketpy.motors import LiquidMotor



#Fluid Properties (SI units)
rho_LOX = 1140.0
rho_ETH = 798.0

#System Test Data Results
CdA_inj_LOX = 0.00001134 #faked for testing
CdA_inj_ETH = 0.00001078 #faked for testing

CdA_feed_LOX = 0.0000305333 #faked for testing
CdA_feed_ETH = 0.0000244267 #faked for testing

#Hydraulic Resistance Terms
R_ox_inj = 1/(2*(CdA_inj_LOX**2)) #dP=Rhyd*mdot^2/rho
R_eth_inj = 1/(2*(CdA_inj_ETH**2)) #dP=Rhyd*mdot^2/rho

R_ox_feed = 1/(2*(CdA_feed_LOX**2)) #dP=Rhyd*mdot^2/rho
R_eth_feed = 1/(2*(CdA_feed_ETH**2)) #dP=Rhyd*mdot^2/rho

R_ox = R_ox_inj + R_ox_feed #Equivalent Hydraulic System Resistance
R_eth = R_eth_inj + R_ox_feed #Equivalent Hydraulic System Resistance

#Tank Properties
gamma_tanks = 1.41 #1.41=GN2, 1.67=GHe
V_oxtank = 6.926 #L
V_ethtank = 7.57 #L

V_oxinit = 3.95 #OPTMIMIZE THIS
V_ethinit = 3.75 #OPTIMIZE THIS

V_oxgas = V_oxtank-V_oxinit
V_ethgas = V_ethtank - V_ethinit

#Initial Tank Pressures
P_tank_ox_psi = 550.0 #psia
P_oxtank = P_tank_ox_psi*6895 #Pa

P_tank_eth_psi = 475.0 #psia
P_ethtank = P_tank_eth_psi*6895 #Pa

#define cstar efficiency: completeion of energy release. See RPE Pg64
Efficiency = 0.925
chamber = CEA_Obj(propName="", oxName="LOX", fuelName="C2H5OH") #initializs CEA object

#define Throat Diameter, Area
Dt = 27.54/1000 #m
At = (Dt**2)/4*math.pi

def Calculate_Residual(Pc, P_oxtank, P_ethtank):
    #solve mass flow from fluid resistances
    mdot_ox = ((1/R_ox)*rho_LOX*(P_oxtank - (Pc*6895)))**(1/2)
    mdot_eth = ((1/R_eth)*rho_ETH*(P_ethtank - (Pc*6895)))**(1/2)
    mdot_fluid = mdot_ox + mdot_eth
    #get OF ratio
    OF_ratio = mdot_ox/mdot_eth
    #solve mass flow from CEA
    Cstar_fps = chamber.get_Cstar(Pc = Pc[0], MR = OF_ratio[0]) #see RPE pg64
    Cstar = float(Cstar_fps*0.3048) #m/s
    if Cstar == 0:
        print("CSTAR ZERO")
        print(f"PC = {Pc[0]}")
        mdot_CEA_res = Pc*6895*At/(abs(Cstar)*Efficiency)
    else:
        mdot_CEA_res = (Pc*6895)*At/(Cstar*Efficiency) #kg/s

    #compare residual, append to array
    error = mdot_CEA_res-mdot_fluid

    global OF_ratio_glob
    OF_ratio_glob = float(OF_ratio[0])

    global mdot_total_glob
    mdot_total_glob = float(mdot_fluid[0])

    if OF_ratio > 1.8: #tank pressure drop eqs breaking
        error = 500

    return abs(error)

def GradientDescent(guess, P_oxtank, P_ethtank):
    # Use scipy optimize minimize with residual function to find Chamber Pressure
    result = minimize(
        Calculate_Residual,
        guess,
        args = (P_oxtank, P_ethtank),
        bounds = [(100, min([P_oxtank/6895, P_ethtank/6895]))],
    )
    P_chamber = result.x[0]
    # t = type(P_chamber)
    #print(f"Pchamb TYPE {t}")
    OF = OF_ratio_glob
    isp = chamber.estimate_Ambient_Isp(Pc=P_chamber,MR=OF,eps=4.35)[0]
    thrust = 9.8*isp*(mdot_total_glob)/1000 #kN
    print(f"Thrust {thrust}")
    massflow_total = mdot_total_glob
    print(f"massflow {massflow_total}")
    #print(f"PC {P_chamber/6895} and MR {OF_ratio} at {i*dt}")

    return P_chamber, thrust, OF, massflow_total

def thrust(V_oxgas, V_ethgas, P_oxtank, P_ethtank):
  iterations = 200
  time = np.linspace(0, 20, iterations) #200 pts from 0 to 15 seconds
  dt = float(time[1]-time[0])
  print(f"TIMESTEP {dt}")
  OF_array = []
  Thrust_array = []
  P_chamber_array = []
  mdtot_array = []
  fin = 0
  for i in range(len(time)): #perform this for every timestep in the profile

      if i == 0:
          Pc_guess = 400
      else:
          Pc_guess = P_chamber_last-5

      P_chamber, Thrust, OF, md_tot = GradientDescent(Pc_guess, P_oxtank, P_ethtank)
      md_ox = md_tot/(
          1+1/OF)
      md_eth = md_tot-md_ox

      #print(f"Timestep {dt}")
      masslost_ox = md_ox*dt
      masslost_eth = md_eth*dt
      #print(f"masslostox {masslost_ox}")

      #print(f"Voxgas(L) {V_oxgas}")
      V_oxgas_next = V_oxgas + (masslost_ox/(rho_LOX*0.001))
      V_ethgas_next = V_ethgas + (masslost_eth/(rho_ETH*0.001))
      #print(f"Voxgasnext(L) {V_oxgas_next}")

      P_oxtank = P_oxtank*((V_oxgas/V_oxgas_next)**gamma_tanks)
      P_ethtank = P_ethtank*((V_ethgas/V_ethgas_next)**gamma_tanks)
      # print(f"Oxtank = {P_oxtank/6895}[psi] ... Ethtank = {P_ethtank/6895}[psi] at {i*dt}")

      V_oxgas = V_oxgas_next
      V_ethgas = V_oxgas_next

      mdtot_array.append(md_tot)
      OF_array.append(OF)
      Thrust_array.append(Thrust)
      P_chamber_array.append(P_chamber)
      assert len(Thrust_array) == len(P_chamber_array)
      P_chamber_last = P_chamber

      if P_chamber>(0.90*P_oxtank/6895) or P_chamber>(0.90*P_oxtank/6895):
          OxDrop = (P_oxtank/6895)/P_chamber
          ETHDrop = (P_ethtank/6895)/P_chamber
          print(f"Flow Stability Violated with {OxDrop}% LOXratio and {ETHDrop}% ETHratio")
          break

      if V_oxgas>=(V_oxtank-V_oxtank/250) or V_ethgas>=(V_ethtank-V_ethtank/250):
          oxrem = V_oxtank-V_oxgas
          ethrem = V_ethtank-V_ethgas
          print(f"Burn finished with {oxrem}L LOX and {ethrem}L ETH at Time {i*dt}s")
          break
  
  file_path = "../LE2/LE2 Simulation and Analysis/Flight sim/data/motors/LE2.eng"  # Specify the desired file path

  if not os.path.exists(file_path):
      f = open(file_path, "x")
      f.close()

  # Open the file for writing
  with open(file_path, "w") as f:
      f.write("; ALULA - LE2 \n")
      f.write("; 8/1/2023 ver. \n")
      f.write("; created by UCB STAR \n")
      f.write("LE2 98 732 0 6.325 8.98822 ALULA\n")

      for i in range(len(Thrust_array)):
          f.write(f"{str(time[i])} {str(Thrust_array[i] * 1000)}\n")  # Thrust_array values in N
  
  env = Environment(
    latitude=32.9901,
    longitude=-106.9751,
    elevation=1400.556
    )

  tomorrow = datetime.date.today() + datetime.timedelta(days=1)

  env.set_date(
    (tomorrow.year, tomorrow.month, tomorrow.day, 12), timezone="America/Denver"
    ) # Tomorrow's date in year, month, day, hour UTC format

  env.set_atmospheric_model(type='Forecast', file='GFS')  
  LiquidMotor.propellant_initial_mass = 4.845+3.3915

  LE2 = LiquidMotor(
                thrust_source="../LE2/LE2 Simulation and Analysis/Flight sim/data/motors/LE2.eng",
                dry_mass=12.685,
                center_of_dry_mass=1.107,
                dry_inertia=(7.332,7.333,0.0318,-0.00153,0.0219,0.0284),
                nozzle_radius=0.0515/2,
                #burn_time,
                nozzle_position=0,
                reshape_thrust_curve=False,
                interpolation_method="linear",
                coordinate_system_orientation="nozzle_to_combustion_chamber",
    )

  LE2.center_of_propellant_mass = 1 ##UPDATE

  ALULA = Rocket(
        #self,
        radius=0.0785,
        mass=36.038,     # total dry
        inertia=(40.32089762,40.32707615,0.16968255,-0.00501791,-0.33114339,0.44168709),
        power_off_drag="../LE2/LE2 Simulation and Analysis/Flight sim/Alula_Cd_PowerOn.csv", ##UPDATE
        power_on_drag="../LE2/LE2 Simulation and Analysis/Flight sim/Alula_Cd_PowerOff.csv", ##UPDATE
        center_of_mass_without_motor=1.996,
        coordinate_system_orientation="tail_to_nose",
)

#LE2.center_of_dry_mass_position=3
#LE2.center_of_propellant_mass=1 ##UPDATE
#ALULA.total_mass=36.038+4.845+3.3915

  ALULA.add_motor(LE2, position=0) # origin = nozzle outlet

  railButtons = ALULA.set_rail_buttons(
    upper_button_position=0.18,
    lower_button_position=-1.4246,
    angular_position=60,
)

  NoseCone = ALULA.add_nose(length=0.762, kind="vonKarman", position=3.99)

  FinSet = ALULA.add_trapezoidal_fins(
    n=3,
    root_chord =0.305,
    tip_chord=0.102,
    span=0.152,
    position=0.4,
    sweep_angle=33.7
  )


  main = ALULA.add_parachute(
    name="main", #120''d
    cd_s=7.07050353*2.2,
    trigger=457.2,  # ejection altitude: 1500ft
    sampling_rate=105,
    lag=1.5,
    noise=(0, 8.3, 0.5),
  )

  drogue = ALULA.add_parachute(
    name="drogue", #60''d
    cd_s=1.767638271*2.2,
    trigger="apogee",  # ejection at apogee
    sampling_rate=105,
    lag=1.5,
    noise=(0, 8.3, 0.5),
  )
  test_flight = Flight(
  rocket=ALULA, environment=env, rail_length=18.288, inclination=85, heading=0
  )

  return test_flight.apogee


res = minimize(lambda x: (thrust(x[0],x[1],x[2],x[3]) - 10000)**2,np.array([3.95, 3.75,475*6895,500*6895]))
print(res.x)
#LE2.total_mass=8.236+12.685


