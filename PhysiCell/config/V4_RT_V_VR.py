import numpy as np 
from matplotlib import pyplot as plt
from scipy.integrate import odeint 

class Model: 
  def __init__(self): 
    # model basics 
    self.ns = {} 
    self.time_series = {} 
    self.time_pt = {} 
    self.time_ctr = 0 
    #compartments 
    self.compartments = { "cell": 1.0}
    # parameters 
    self.parameters = { "TIME_CONV": 60.0,
              "ACE2_0": 10000.0,
              "r_rep_half": 200.0,
              "max_infected_half_max": 250.0,
              "apoptosis_hill": 1.0,
              "virus_fraction_released_at_death": 0.0,
              "r_bind": 0.6,
              "r_endo": 0.6,
              "r_release": 0.6,
              "r_recycle": 0.6,
              "r_uncoat": 0.6,
              "r_prep": 0.6,
              "r_synth": 0.6,
              "r_rep_max": 180.0,
              "r_rna_deg": 0.6,
              "r_assemble": 0.6,
              "r_exo": 0.6,
              "max_infected_apoptosis_rate": 0.12}
    # boundary species 
    self.bspecies= {} 
    # species 
    self.species= { "S1": 0,
              "S2": 0,
              "S3": 0,
              "S4": 0,
              "S5": 0,
              "S6": 0,
              "S7": 0,
              "S8": 0,
              "S9": 0,
              "S10": 0,
              "S11": 0,
              "S12": 0}
    # initial species array 
    self.init_species = np.array([10.0, 10000.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    # species name map 
    self.species_name_map = { "S1": "nV()",
              "S2": "R_EU()",
              "S3": "Cell()",
              "S4": "E()",
              "S5": "R_EB()",
              "S6": "R_IB()",
              "S7": "R_IU()",
              "S8": "V()",
              "S9": "U()",
              "S10": "R()",
              "S11": "P()",
              "S12": "A()"}
    # assignment rules 
    self.arules= { "nV": "0 + S1",
              "R_EB": "0 + S5",
              "V": "0 + S8",
              "U": "0 + S9",
              "R": "0 + S10",
              "P": "0 + S11",
              "A": "0 + S12",
              "E": "0 + S4",
              "Cell": "0 + S3",
              "infected_apoptosis_rate": "max_infected_apoptosis_rate * (pow(A, apoptosis_hill) / (pow(A, apoptosis_hill) + pow(max_infected_half_max, apoptosis_hill)))",
              "rna_rep_rate": "r_rep_max / (R + r_rep_half)",
              "virus_export_rate": "Cell * r_exo",
              "vout_day": "Cell * r_exo * A * 24",
              "dead_cell_release": "virus_fraction_released_at_death * infected_apoptosis_rate"}
    # reactions 
    self.rules= { "R1": "r_bind * S1 * S2",
              "R2": "infected_apoptosis_rate * S3",
              "R3": "r_endo * S5",
              "R4": "r_release * S6",
              "R5": "r_recycle * S7",
              "R6": "r_uncoat * S8",
              "R7": "r_prep * S9",
              "R8": "r_synth * S10",
              "R9": "rna_rep_rate * S10",
              "R10": "r_rna_deg * S10",
              "R11": "r_assemble * S11",
              "R12": "virus_export_rate * S12",
              "R13": "dead_cell_release * S12"}
    # variables to calc derivatives for 
    self.variables = { "S1": "-R1",
              "S2": "-R1+R5",
              "S3": "-R2",
              "S4": "+R12+R13",
              "S5": "+R1-R3",
              "S6": "+R3-R4",
              "S7": "+R4-R5",
              "S8": "+R4-R6",
              "S9": "+R6-R7",
              "S10": "+R7-R8+R8-R9+R9+R9-R10",
              "S11": "+R8-R11",
              "S12": "+R11-R12-R13"}

  def update_namespace(self): 
    ''' 
    this function starts up a new namespace 
    which will be used to evaluate formulas
    ''' 
    self.ns = {} 
    self.ns.update(self.compartments) 
    self.ns.update(self.parameters) 
    self.ns.update(self.bspecies) 
    self.ns.update(self.species) 

  def calc_assignment_rules(self): 
    ''' 
    this function evaluates assignment rule 
    formulas and keeps the namespace updated
    ''' 
    arule_vals = {} 
    for arule in self.arules: 
      formula = self.arules[arule] 
      arule_vals[arule] = eval(formula, self.ns)
      self.ns[arule] = arule_vals[arule]
    return arule_vals 

  def calc_reaction_rates(self): 
    ''' 
    this function evaluates kinetics to calculate 
    reaction rates using parameter and evaluated 
    assignment rule values 
    ''' 
    reaction_rates = {} 
    for rxn in self.rules: 
      formula = self.rules[rxn] 
      reaction_rates[rxn] = eval(formula, self.ns)
    self.ns.update(reaction_rates) 
    return reaction_rates 

  def calc_deriv(self, spec_vals, time=None): 
    ''' 
    this function calls all the helper functions 
    to start up a new name space, evaluate assignment rules 
    calculate reaction rates and finally calculate 
    the derivatives of the variables. A time point  
    dictionary will be created for each call to this 
    function and saved under self.time_series. 
    ''' 
    # reset time point dict 
    self.time_pt = {} 
    # set current species values 
    for i,spec in enumerate(self.species): 
      self.species[spec] = spec_vals[i] 
    # get new namespace and update 
    self.update_namespace() 
    # update time point 
    self.time_pt["compartments"] = self.compartments 
    self.time_pt["parameters"] = self.parameters 
    self.time_pt["bspecies"] = self.bspecies 
    self.time_pt["species"] = self.species 
    # calculate assignment rules and save 
    ar = self.calc_assignment_rules() 
    self.time_pt["arules"] = ar 
    # calculate reaction rates and save 
    rr = self.calc_reaction_rates() 
    self.time_pt["rrates"] = rr 
    # calculate derivatives and save 
    deriv = [] 
    var_dict = {} 
    for i,var in enumerate(self.variables): 
      dxdt = eval(self.variables[var], self.ns)
      deriv.append(dxdt)
      var_dict["{0}_deriv".format(var)] = dxdt 
    self.time_pt["deriv"] = var_dict 
    # save time point dictionary 
    if time is not None: 
      self.time_series[time] = self.time_pt 
    else: 
      self.time_series[self.time_ctr] = self.time_pt 
      self.time_ctr += 1 
    # return derivatives 
    return np.array(deriv) 

  def __call__(self, vals, time): 
    ''' 
    this allows us to call the object like a function
    ''' 
    return self.calc_deriv(vals, time) 

if __name__ == "__main__": 
  import argparse 
  parser = argparse.ArgumentParser(description="A python model object generated from SBML") 
  parser.add_argument("-t0", dest="t0", type=float, default=0.0, 
          help="Initial time point for integration") 
  parser.add_argument("-tEnd", dest="tEnd", type=float, default=1.0, 
          help="Final time point for integration") 
  parser.add_argument("-nt", dest="nt", type=int, default=100, 
          help="Number of time points for integration") 
  parser.add_argument("-o", "--output", dest="output", type=str, default="test.png", 
          help="Output plot file") 
  args = parser.parse_args() 
  # get time points 
  time = np.linspace(args.t0, args.tEnd, args.nt)

  # initialize model 
  m = Model()
  # integrate the model 
  result = odeint(m, m.init_species, time)

  # do some plotting 
  fig = plt.figure()
  ax = plt.subplot(111)
  plt.plot(time,result[:,0], label=m.species_name_map["S1"], lw=1.5)
#   plt.plot(time,result[:,1], label=m.species_name_map["S2"], lw=1.5)
  plt.plot(time,result[:,2], label=m.species_name_map["S3"], lw=1.5)
  plt.plot(time,result[:,3], label=m.species_name_map["S4"], lw=1.5)
  plt.plot(time,result[:,4], label=m.species_name_map["S5"], lw=1.5)
  plt.plot(time,result[:,5], label=m.species_name_map["S6"], lw=1.5)
  plt.plot(time,result[:,6], label=m.species_name_map["S7"], lw=1.5)
  plt.plot(time,result[:,7], label=m.species_name_map["S8"], lw=1.5)
  plt.plot(time,result[:,8], label=m.species_name_map["S9"], lw=1.5)
  plt.plot(time,result[:,9], label=m.species_name_map["S10"], lw=1.5)
  plt.plot(time,result[:,10], label=m.species_name_map["S11"], lw=1.5)
  plt.plot(time,result[:,11], label=m.species_name_map["S12"], lw=1.5)
  box = ax.get_position()
  ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
  plt.xlabel("time")
  plt.ylabel("concentration")
  plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
  # plt.show()
  plt.savefig(args.output)
