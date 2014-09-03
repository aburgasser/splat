"""
.. note::
         sUsing a suite of evolutionary models, this code translates 
         between the following brown dwarf parameters: mass, age, 
         temperature, radius, surface gravity, and luminosity. We allow 
         the user to choose a set of evolutionary model 
         (Baraffe, Burrows, or Saumon) and two parameters, then output
         the rest of the interpolated parameters. 
"""

# Standard library imports.
from sys import exit 
from urllib2 import urlopen

# Related third party imports.
from scipy.interpolate import interp1d 
from astropy import units as u
from astropy.io import ascii
from math import exp

###############################################################################
###############################################################################
class ReadModel(object):
    """
    :Description: This class reads in evolutionary models that are defined 
                  in the methods below, and their data is acquired `here
                  <http://pono.ucsd.edu/~adam/splat/EvolutionaryModels/>`_.
                  Units are the following: masses are in M/Msun, 
                  luminosities in log L/Lsun, radius in R/Rsun, surface 
                  gravities in log g (cm/s^2), temperatures in kelvins, 
                  and ages in Gyr. Each evolutionary model gives 
                  snapshots of Brown Dwarfs with different masses and 
                  other physical properties as a function of time.

    :param model: 
        - **baraffe:** 
          Isochrones from Baraffe et. la models (2003), described in the 
	  following paper: "Evolutionary models for cool brown dwarfs and 
	  extrasolar giant planets. The case of HD 20945": `Here 
	  <http://arxiv.org/abs/astro-ph/0302293>`_. Original model's `URL 
          <https://perso.ens-lyon.fr/isabelle.baraffe/COND03_models>`_.
          Ages (in Gyr) used for interpolation were the following: 

	     >>> from bdevopar import *
             >>> print ReadModel('baraffe')['age']
             [0.001, 0.005, 0.01, 0.05, 0.1, 0.12, 0.5, 1.0, 5.0, 10.0]
             >>> print ReadModel('BaRaFfE')['age']
             [0.001, 0.005, 0.01, 0.05, 0.1, 0.12, 0.5, 1.0, 5.0, 10.0]

	 .. note:: Capilatizing some letters or none won't cause any problems.
	           However, incorrect spelling will raise a NameError, so
		   please make sure you spell correctly.

        - **burrow:**

        - **saumon:**
          Isochrones from Saumon & Marley models (2008), described in 
          `paper <http://adsabs.harvard.edu/abs/2008ApJ...689.1327S>`_. 
          Original models' URL and the README used to differentiate
          between metallicites are found `here.
          <https://laws.lanl.gov/x7/dsaumon/BD_evolution/>`_
          Brown dwarfs ages used here are as follows:

	     >>> from bdevopar import *
             >>> print ReadModel('saumon')['age']
             [0.003, 0.004, 0.006, 0.008, 0.01, 0.015, 0.02, 0.03, 0.04, 
	     0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 
	     2.0, 3.0, 4.0, 6.0, 8.0, 10.0]

           :metallicity:
             - **nc_solar:** The atmosphere model is cloudless with [M/H]=0
             - **nc+0.3:** The atmosphere model is cloudless with [M/H]=+0.3
             - **nc-0.3:** The atmosphere model is cloudless with [M/H]=-0.3
             - **f2_solar:** Atmosphere model is cloudless Fsed=2) with [M/H]=0
             - **hybrid_solar** Atmosphere cloudless(Fsed=2 to nc) with [M/H]=0

              >>> from bdevopar import *
	      >>> saumon = ReadModel('saumon',metallicity='nc_solar')
	      >>> saumon = ReadModel('saumon',metallicity='f2_solar')
	      >>> saumon = ReadModel('saumon',metallicity='hybri_solar')

    :Return: A dictionary where each key maps to a 3 dimensional matrix. The 
             dimensions are as follows: (Age x Mass x OtherParam). This
	     method is used in the subclass Parameters, which is 
	     further discussed in the section belowed.
    """
    def __new__(cls, *model, **z): 
        Emodels_URL = 'http://pono.ucsd.edu/~adam/splat/EvolutionaryModels/'

        assert len(model) <= 1, "Only one argument (model) allowed."
        assert len(z.keys()) <= 1, "Only one keyword (metallicity) allowed."

        try: model = model[0].lower()
  	except TypeError: raise TypeError("Model can't be a number.")
        except IndexError: model = 'baraffe'
	finally: print "You are using " + model + "'s model."

        assert model in ['baraffe','burrow','saumon'], """
	    Incorrect model or bad spelling. Please choose from the 
	    following models: 'baraffe', 'burrow', or 'saumon'.\n"""

        ##################### BARAFFE OR BURROW MODEL #########################
        if model == 'baraffe' or model == 'burrow':
            ages = ['0.001', '0.005', '0.010', '0.050', '0.100',
                    '0.120', '0.500', '1.000', '5.000', '10.000']

            if model == 'baraffe': Emodels = Emodels_URL + 'Baraffe/cond_'
            else: Emodels = Emodels_URL + 'Burrow/b97_'
            num_files = len(ages); v = 0
        ########################### SAUMON MODEL ##############################
        else:
            try: metallicity = z[z.keys()[0]]
  	    except TypeError: raise TypeError("Metallicity can't be a number.")
            except IndexError: metallicity = 'Hybrid_solar'
	    finally: print "You are using the metallicity " + metallicity + "."

            if metallicity == 'Hybrid_solar': Z = 'hybrid_solar_age_'
            elif metallicity == 'No_Clouds_+0.3': Z = 'nc+0.3_age_'
            elif metallicity == 'No_Clouds_-0.3': Z = 'nc-0.3_age_'
            elif metallicity == 'No_Clouds_Solar': Z = 'no_solar_age_'
            elif metallicity == 'f2_solar': Z = 'f2_solar_age_'
            else: raise NameError("The inputted metallicity '" + metallicity +\
	     "' does not exist. Make sure spelling \n           is correct."+\
	     " Please choose from the following: 'Hybrid_solar',\n           "\
	     +"'No_Clouds_+0.3', 'No_Clouds_-0.3', 'No_Clouds_Solar',\n "+\
	     "          or 'f2_solar'.\n")

            ages = ['0.003','0.004','0.006','0.008','0.010','0.015','0.020',
                    '0.030','0.040','0.060','0.080','0.100','0.150','0.200',
                    '0.300','0.400','0.600','0.800','1.000','1.500','2.000',
                    '3.000','4.000','6.000','8.000','10.000']

            Emodels = Emodels_URL + 'Saumon/' + metallicity + '/' + Z
	    num_files = len(ages); v = 1
        #######################################################################

        n_tables = range(num_files)
  
        masses = [[] for i in n_tables]
        temperatures = [[] for i in n_tables]
        luminosities = [[] for i in n_tables]
        gravities = [[] for i in n_tables]
        radiuses = [[] for i in n_tables]

        for age in n_tables:
	    try:data =ascii.read(urlopen(Emodels+ages[age]).read(),comment=';')
	    except: assert False, "The servers are offline or yours is."
            for line,value in enumerate(data):
                masses[age].append(float(value[0+v]))
                temperatures[age].append(float(value[1+v]))
                luminosities[age].append(float(value[2+v]))
                gravities[age].append(float(value[3+v]))
                radiuses[age].append(float(value[4+v]))

        return {'masses':masses, 'gravities':gravities, 'radiuses':radiuses,
                'luminosities':luminosities, 'temperatures':temperatures, 
		'ages':[float(i) for i in ages]}

############################### End of class ReadModel ########################
###############################################################################
class Parameters(ReadModel):
    """
    :Description: Checks that the input parameters are valid (e.g user 
                  inputs must be capable of being interpolated). 
                  Interpolation between evolutionary models is done here 
                  as well, and we used the method interp1d from astropy.

    :param *args: If you're using Saumon's models, then u must specified the
                  metallicity since there's 5 total: hybrid_solar, f2_solar, 
                  nc+0.3, nc-0.3, and nc_solar.

    :param *kwargs: You must input any of the two following parameters: age, 
                    mass, temperature, luminosity, gravity, or radius.

    :returns: Returns the interpolated parameters.

     - **Examples:**
      
       >>> model = ReadModel('saumon', z='hybrid_solar')
       >>> saumon = Parameters(model, mass=0.4, age=0.01)
       >>> model = ReadModel('saumon', z='nc+0.3')
       >>> saumon = Parameters(a=0.01, mss=0.4)
       >>> model = ReadModel('saumon', z='nc-0.3')
       >>> saumon = Parameters(mas=0.4, RiUS=0.01)
       >>> model = ReadModel('saumon', z='nc_solar')
       >>> saumon = Parameters(M=0.4, AgE=0.01)
       >>> model = ReadModel('saumon', z='f2_solar')
       >>> saumon = Parameters(adf=0.01, MASSSS=0.4)

       >>> model = ReadModel('baraffe', z='f2_solar')
       >>> burrow = Parameters(model, L=0.05, Gradfty=0.4)
       >>> model = ReadModel('burrow', z='f2_solar')
       >>> baraffe = Parameters(model, temp=2000, age=0.4)

    .. note:: Keyword spelling doesn't matter as long as the first 
              letter is right. Also, ordering of the keywords is not 
              important. However, metallicites must be spelled exactly 
              as mentioned aboved--there are no execptions.
    """
    def __new__(cls, *args, **kwargs):
      ################# Check User Inputted Keywords ##########################
      keywords = kwargs.keys()
      numberKeys = len(keywords)
      params = {'temperature':0, 'luminosity':0, 
                'age':0, 'gravity':0, 'mass':0, 'radius':0}

      try: model = args[0]
      except IndexError: model = ReadModel('baraffe')

      for i in range(numberKeys):
         if keywords[i][0].upper().startswith('T'): 
            if params['temperature'] == 0: 
               params['temperature'] = float(kwargs[keywords[i]])
            else: raise NameError('Only 1 temperature allowed.\n')
         elif keywords[i][0].upper().startswith('A'): 
            if params['age'] == 0:
               params['age'] = float(kwargs[keywords[i]])
            else: raise NameError('Only 1 age allowed.\n')
         elif keywords[i][0].upper().startswith('G'): 
            if params['gravity'] == 0: 
               params['gravity'] = float(kwargs[keywords[i]])
            else: raise NameError('Only 1 gravity allowed.\n')
         elif keywords[i][0].upper().startswith('R'): 
            if params['radius'] == 0: 
               params['radius'] = float(kwargs[keywords[i]])
            else: raise NameError('Only 1 radius allowed.\n')
         elif keywords[i][0].upper().startswith('M'): 
            if params['mass'] == 0: 
               params['mass'] = float(kwargs[keywords[i]])
            else: raise NameError('Only 1 mass allowed.\n')
         elif keywords[i][0].upper().startswith('L'): 
            if params['luminosity'] == 0: 
               params['luminosity'] = float(kwargs[keys[i]])
            else: raise NameError('Only 1 luminosity allowed.\n')
         else: raise NameError("Keyword '"+keywords[i]+"' is nonexistent.\n")

      if params['temperature'] != 0: 
         assert min(min(model['temperatures'][:])) <= \
              params['temperature'] <= max(max(model['temperatures'][:])),\
              "Temperature is out of model's range"
      if params['luminosity'] != 0: 
         assert min(min(model['luminosities'][:])) <= \
              params['luminosity'] <= max(max(model['luminosities'][:])), \
              "Luminosity is out of model's range"
      if params['gravity'] != 0: 
         assert min(min(model['gravities'][:])) <= params['gravity'] <= \
            max(max(model['gravities'][:])), "Gravity is out of model's range."
      if params['radius'] != 0: 
         assert min(min(model['radiuses'][:])) <= params['radius'] <= \
              max(max(model['radiuses'][:])), "Radius is out of model's range."
      if params['mass'] != 0: 
         assert min(min(model['masses'][:])) <= params['mass'] <= \
              max(max(model['masses'][:])), "Mass is out of model's range."
      if params['age'] != 0: 
         assert min(model['ages'][:]) <= params['age'] <= \
              max(model['ages'][:]), "Age is out of model's range."

      #########################################################################
      Ag, Ma, Te, Le, Ge, Re = [],[],[],[],[],[]
      input_type = 'mass_age'
      valid_ages = []
      n_tables = range(len(model['ages']))

      ############### WITH TWO KNOWN PARAMETERS, SPIT OUT AGE #################
      if (params['mass'] == False) and (params['age'] == False):

          input_type = 'two_params'; P = []; Lno, Gno, Rno = True, True, True
          for i in [0,1]: # Because there's two input parameters
              if (params['luminosity'] != False) and Lno:  
                  P.append(['luminosity', params['luminosity']]); Lno=0
              elif (params['gravity'] != False) and Gno:   
                  P.append(['gravity', params['gravity']]); Gno=0
              elif (params['radius'] != False) and Rno:      
                  P.append(['radius', params['radius']]); Rno=0
              else: P.append(['temperature', params['temperature']])

          for i in n_tables:
              if min(model[P[0][0]][i]) <= P[0][1] <= max(model[P[0][0]][i]) \
               and min(model[P[1][0]][i]) <= P[1][1] <= max(model[P[1][0]][i]):
                  Ag.append(model['age'][i]) 
                  valid_ages.append(i)
                  f = interp1d(model[P[0][0]][i], model['mass'][i])
                  Ma = f(P[0][1])
                  f = interp1d(model['mass'][i], model[P[1][0]][i])
                  Ge.append(f(Ma))

          f = interp1d(Ge, Ag)
          try:
             params['age'] = f(P[1][1])
          except ValueError:
             assert False, P[1][0] + " is above the interpolation range."

          del Ge[:], Ag[:]; Ma = []
 
      ################ WITH KNOWN AGE OR MASS AND ANOTHER PARAMETER ###########
      ################### SPIT OUT MASS OR AGE, RESPECTIVELY ##################
      if ((params['age'] != False) and (params['mass'] == False)) \
         or ((params['mass'] != False) and (params['age'] == False)):

          if input_type != 'two_params': 
              input_type = 'one_param'; P = []
              if params['temperature'] != False: 
                  P.append(['temperature', params['temperature']])
              elif params['luminosity'] != False:  
                  P.append(['luminosity', params['luminosity']])
              elif params['gravity'] != False:   
                  P.append(['gravity', params['gravity']])
              elif params['radius'] != False:      
                  P.append(['radius', params['radius']])
              P.append(['', 0])

          if input_type == 'two_params': 
              n_t = valid_ages
              valid_ages = []
          else: n_t = model['n_tables']

          for i in n_t:
             if min(model[P[0][0]][i]) <= P[0][1] <= max(model[P[0][0]][i]):
                 valid_ages.append(i)
                 Ag.append(model['age'][i]) 
                 f = interp1d(model[P[0][0]][i], model['mass'][i])
                 Ma.append(f(P[0][1]))
          
          if params['mass'] == False:
              f = interp1d(Ag, Ma)
              params['mass'] = round(f(params['age']), 5)
          else:
              f = interp1d(Ma, Ag)
              params['age'] = round(f(params['mass']), 5)
          del Ag[:]

      ###################### WITH KNOWN MASS AND AGE ##########################
      if (params['mass'] != False) and (params['age'] != False):

          if input_type == 'mass_age': valid_ages = n_tables

          for i in valid_ages:
              if min(model['masses'][i]) <= params['mass'] \
                                                  <= max(model['masses'][i]):
                  Ag.append(model['ages'][i])
                  f =interp1d(model['masses'][i],model['temperatures'][i])
                  Te.append(f(params['mass']))
                  f = interp1d(model['masses'][i],model['luminosities'][i])
                  Le.append(f(params['mass']))
                  f = interp1d(model['masses'][i],model['gravities'][i])
                  Ge.append(f(params['mass']))
                  f = interp1d(model['masses'][i],model['radiuses'][i])
                  Re.append(f(params['mass']))
      
      f = interp1d(Ag, Te) 
      params['temperature'] = round(f(params['age']), 8)
      f = interp1d(Ag, Le)
      params['luminosity'] = round(f(params['age']), 8)
      f = interp1d(Ag, Ge) 
      params['gravity'] = round(f(params['age']), 8)
      f = interp1d(Ag, Re)
      params['radius'] = round(f(params['age']), 8)
      
      if input_type == 'one_param': params[P[0][0]] = P[0][1]
      elif input_type == 'two_params': 
          params[P[0][0]] = P[0][1]
          params[P[1][0]] = P[1][1]

#      print "\n------------------------------------------"
#      print "| BD's age (Gyr):           | " + str(params['age'])
#      print "------------------------------------------"
#      print "| BD's mass (Msun):         | " + str(params['mass'])
#      print "------------------------------------------"
#      print "| Effective temperature (K) | " + str(params['temperature'])
#      print "------------------------------------------"
#      print "| log L_bd/L_sun:           | " + str(params['luminosity'])
#      print "------------------------------------------"
#      print "| log surface gravity:      | " + str(params['gravity'])
#      print "------------------------------------------"
#      print "| BD's radius (Rsun):       | " + str(params['radius'])
#      print "------------------------------------------\n"

      return {'temperature':params['temperature']*u.K, 
              'mass':params['mass']*u.solMass, 'age':params['age']*u.Gyr, 
              'luminosity':exp(params['luminosity'])*u.solLum,
	      'gravity':exp(params['gravity'])*u.centimeter/u.second**2, 
	      'radius':params['radius']*u.solRad}
########################## End of the class: bdevopar #########################
###############################################################################

class ParamsList(Parameters, ReadModel):
    """
    :Description: Allows the user to input a list of parameters.
    :Returns: 
       >>> model = ReadModel('baraffe')
       >>> m = [0.04,0.06,0.07]
       >>> a = [0.001,0.003,0.004]
       >>> params = ParamsList(model,masses=m,ages=a)
       >>> print params
    """
    def __new__(cls,*model,**kwargs):
	 keywords = kwargs.keys()
         assert len(kwargs[keywords[0]]) == len(kwargs[keywords[1]]), """
	     Number of elements in both input lists must be equal."""
	 assert len(keywords) == 2, "Only two keywords (lists) allowed."
	 
	 try: model = model[0]
	 except IndexError: model = 'baraffe'

         if type(model) is not dict: model = ReadModel(model)

         params = {'temperatures':[],'ages':[],'gravities':[],
	           'radiuses':[],'masses':[],'luminosities':[]}
         T, A, G, R, M, L = False,False,False,False,False,False

         for i in range(2):
             if keywords[i].upper().startswith('T'): 
                 T = True
	         temperature = kwargs[keywords[i]]
             elif keywords[i].upper().startswith('A'): 
	         A = True
	         age = kwargs[keywords[i]]
             elif keywords[i].upper().startswith('G'): 
	         G = True
	         gravity = kwargs[keywords[i]]
             elif keywords[i].upper().startswith('R'): 
	         R = True
	         radius = kwargs[keywords[i]]
             elif keywords[i].upper().startswith('M'): 
	         M = True
	         mass = kwargs[keywords[i]]
             elif keywords[i].upper().startswith('L'): 
	         L = True 
	         luminosity = kwargs[keywords[i]]
         
	 numberValues = len(kwargs[keywords[0]])
         p = [[] for i in range(numberValues)]

	 if A == True and M == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,m=mass[i],a=age[i])
	 elif A == True and L == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,l=luminosity[i],a=age[i])
	 elif A == True and T == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,t=temperature[i],a=age[i])
	 elif A == True and G == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,g=gravity[i],a=age[i])
	 elif A == True and R == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,r=radius[i],a=age[i])
	 elif M == True and R == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,r=radius[i],m=mass[i])
	 elif M == True and L == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,l=luminosity[i],m=mass[i])
	 elif M == True and T == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,t=temperature[i],m=mass[i])
	 elif M == True and G == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,m=mass[i],g=gravity[i])
	 elif R == True and T == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,r=radius[i],t=temperature[i])
	 elif R == True and G == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,r=radius[i],g=gravity[i])
	 elif R == True and L == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,r=radius[i],l=luminosity[i])
	 elif G == True and T == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,g=gravity[i],t=temperature[i])
	 elif G == True and L == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,g=gravity[i],l=luminosity[i])
	 elif L == True and T == True:
  	    for i in range(numberValues):
	        p[i] = Parameters(model,r=luminosity[i],t=temperature[i])

	 for i in range(numberValues):
	     params['temperatures'].append(p[i]['temperature'])
	     params['ages'].append(p[i]['age'])
             params['gravities'].append(p[i]['gravity'])
	     params['radiuses'].append(p[i]['radius'])
	     params['masses'].append(p[i]['mass'])
	     params['luminosities'].append(p[i]['luminosity'])

         return params
