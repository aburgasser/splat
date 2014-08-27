# FIX:
#    make sure that all the values between models are consistent in unit types.
#
#    if parameter (not mass/age) is the min/max, then there's no interplation
#    because we need two files. THus, if input param value is min/max, then just
#    read out of model. 
#       ex:   if self.params['Temp'] == min(min(model['Temp'])):

"""
.. note::
              Using a suite of evolutionary models, this code translates 
              between the following brown dwarf parameters: stellar mass, age, 
              temperature, radius, surface gravity, and luminosity. We allow 
              the user to choose a set of evolutionary model 
              (Baraffe, Burrows, and Saumon) and two parameters, then give back
              the rest of the interpolated parameters. 
"""

from sys import exit 
from warnings import filterwarnings

from astropy.io import ascii
from scipy.interpolate import interp1d 
from urllib2 import urlopen

Emodelss = 'http://pono.ucsd.edu/~adam/splat/EvolutionaryModels/'
filterwarnings("ignore")

################################################################################
################################################################################
class EvolutionaryModels(object):
    """
    :Description: ``This class reads in evolutionary models, that are defined in
                    the methods below, and their data was acquired from the URL:
                    http://pono.ucsd.edu/~adam/splat/EvolutionaryModels/.
                    Units are the following: masses are in M/Msun, luminosities
                    in log L/Lsun, radius in R/Rsun, surface gravities in log 
                    g (cm/s^2), temperatures in kelvins, and ages in Gyr.
                    Each evolutionary model gives snapshots of Brown Dwarfs 
                    with different masses and other physical properties as a 
                    function of time.``
    """
    ############################## BARAFFE MODEL ###############################
    def baraffe_model(self):
        """
        :Summary: ``Isochrones from Baraffe et. la models (2003), described in
                    the following paper: "Evolutionary models for cool brown
                    dwarfs and extrasolar giant planets. The case of HD 20945":
                    arxiv.org./abs/astro-ph/0302293. Original model's URL: 
                    perso.ens-lyon.fr/isabelle.baraffe/COND03_models.
                    Ages (in Gyr) used for interpolation were the following:
                    0.001, 0.005, 0.01, 0.05, 0.1, 0.12, 0.5, 1, 5, and 10.``
       
        :Parameters: ``None.``

        :Return: ``A dictionary where each key maps to a 3 dimensional matrix.
                   The dimensions are as follows: (Age x Mass x OtherParam).
                   This method is used in the subclass Bdevopar, which is 
                   discussed in the section belowed.``
                 
        """
        bd_ages = ['0.001', '0.005', '0.010', '0.050', '0.100',
                   '0.120', '0.500', '1.000', '5.000', '10.000']
        n_tables = range(len(bd_ages))

        # Declare list where each element is a list containing values for a file
        masses = [[] for i in n_tables]
        temperatures = [[] for i in n_tables]
        luminosities = [[] for i in n_tables]
        gravities = [[] for i in n_tables]
        radiuses = [[] for i in n_tables]
       
        global Emodelss
        Emodels = Emodelss + 'Baraffe/cond_'
        for age in n_tables:
            data = ascii.read(urlopen(Emodels+bd_ages[age]).read(),comment=';')
            for line in range(len(data)):
                masses[age].append(data[line][0])
                temperatures[age].append(data[line][1])
                luminosities[age].append(data[line][2])
                gravities[age].append(data[line][3])
                radiuses[age].append(data[line][4])

        # Convert string list 'bd_ages' into a list of float numbers.
        bd_ages = [float(i) for i in bd_ages] 

        return {'n_tables':n_tables, 'Mass':masses, 'Temp':temperatures,
                'Age':bd_ages, 'Radi':radiuses, 'Lumi':luminosities,
                'Grav':gravities}

    ############################## BURROW MODEL ###############################
    def burrow_model(self): 
        """
        :Summary: ``Isochrones from Burrow.``
        """
        bd_ages = ['0.001', '0.005', '0.010', '0.050', '0.100',
                   '0.120', '0.500', '1.000', '5.000', '10.000']
        n_tables = range(len(bd_ages))

        # Declare list where each element is a list containing values for a file
        masses = [[] for i in n_tables]
        temperatures = [[] for i in n_tables]
        luminosities = [[] for i in n_tables]
        gravities = [[] for i in n_tables]
        radiuses = [[] for i in n_tables]

        global Emodelss
        Emodels = Emodelss + 'Burrow/b97_'
        for age in n_tables:
            data = ascii.read(urlopen(Emodels+bd_ages[age]).read(),comment=';')
            for line in range(len(data)):
                masses[age].append(data[line][0])
                temperatures[age].append(data[line][1])
                luminosities[age].append(data[line][2])
                gravities[age].append(data[line][3])
                radiuses[age].append(data[line][4])

        # Convert string list 'bd_ages' into a list of float numbers.
        bd_ages = [float(i) for i in bd_ages] 

        return {'n_tables':n_tables, 'Mass':masses, 'Temp':temperatures,
                'Age':bd_ages, 'Radi':radiuses, 'Lumi':luminosities,
                'Grav':gravities}

    ############################## SAUMON MODELS #############################
    def saumon_model(self, metallicity):
        """
        :Summary: ``Isochrones from Saumon & Marley models (2008), described 
                    in this paper: adsabs.harvard.edu/abs/2008ApJ...689.1327S. 
                    Original models' URL and the README used to differentiate
                    between metallicites are found here: 
                    https://laws.lanl.gov/x7/dsaumon/BD_evolution/
                    Brown dwarfs ages used here are as follows: 0.003, 0.004,
                    0.006, 0.008, 0.010, 0.015, 0.020, 0.030, 0.040, 0.060,
                    0.080, 0.100, 0.150, 0.200, 0.300, 0.400, 0.600, 0.800,
                    1.000, 1.500, 2.000, 3.000, 4.000, 6.000, 8.000, 10.000``

        :param nc_solar: ``The atmosphere model is cloudless with [M/H]=0``
        :param nc+0.3: ``The atmosphere model is cloudless with [M/H]=+0.3``
        :param nc-0.3: ``The atmosphere model is cloudless with [M/H]=-0.3``
        :param f2_solar: ``Atmosphere model is cloudless (Fsed=2) with [M/H]=0``
        :param hybrid_solar: ``Atmosphere cloudless(Fsed=2 to nc) with [M/H]=0``
 
        :returns: ``A dictionary with similar format as baraffe's, but with 
                    saumon's model values.``
        """
        bd_ages = ['0.003','0.004','0.006','0.008','0.010','0.015','0.020',
                   '0.030','0.040','0.060','0.080','0.100','0.150','0.200',
                   '0.300','0.400','0.600','0.800','1.000','1.500','2.000',
                   '3.000','4.000','6.000','8.000','10.000']
        n_tables = range(len(bd_ages))

        if metallicity == 'hybrid_solar': 
            Z = 'hybrid_solar_age_'
        elif metallicity == 'nc+0.3':
            Z = 'nc+0.3_age_'
        elif metallicity == 'nc-0.3':
            Z = 'nc-0.3_age_'
        elif metallicity == 'nc_solar':
            Z = 'no_solar_age_'
        elif metallicity == 'f2_solar':
            Z = 'f2_solar_age_'
        else: print('This metallicity does not exist.'); NameError()

        masses = [[] for i in n_tables]
        temperatures = [[] for i in n_tables]
        luminosities = [[] for i in n_tables]
        gravities = [[] for i in n_tables]
        radiuses = [[] for i in n_tables]

        global Emodelss
        type(metallicity)
        Emodels = Emodelss + 'Saumon/' + metallicity + '/' + Z
        print Emodels
        for age in n_tables:
            data = ascii.read(urlopen(Emodels+bd_ages[age]).read(),comment=';')
            for line in range(len(data)):
                masses[age].append(data[line][1])
                temperatures[age].append(data[line][3])
                luminosities[age].append(data[line][2])
                gravities[age].append(data[line][4])
                radiuses[age].append(data[line][5])

        # Convert string list 'bd_ages' into a list of float numbers.
        bd_ages = [float(i) for i in bd_ages] 

        return {'n_tables':n_tables, 'Mass':masses, 'Temp':temperatures,
                'Age':bd_ages, 'Radi':radiuses, 'Lumi':luminosities,
                'Grav':gravities}

###################### End of class EvolutionaryModels #########################
################################################################################
class Bdevopar(EvolutionaryModels):
   """
   :Description: ``Checks that the input parameters are valid (e.g user inputs
                   must be capable of being interpolated). Interpolation between
                   evolutionary models is done here as well, and we used the 
                   method interp1d from astropy.``

   :param *args:  ``If you're using Saumon's models, then you must specified the
                    metallicity since there's 5 total:  hybrid_solar, f2_solar, 
                    nc+0.3, nc-0.3, and nc_solar.``

   :param *kwargs: ``You must input any of the two following parameters: age, 
                     mass, temperature, luminosity, gravity, or radius.`` 

   :returns:       ``Returns an instance of the class object bdevopar. In
                     other words, you will receive interpolated parameters.``

   - **Examples:**
      >>> saumon = Bdevopar('hybrid_solar', E='saumon', mass=0.4, age=0.01)
      >>> saumon = Bdevopar('nc+0.3', a=0.01, mss=0.4, E='saumon')
      >>> saumon = Bdevopar('nc-0.3', E='saumon', mas=0.4, RiUS=0.01)
      >>> saumon = Bdevopar('nc_solar', M=0.4, esdf='saumon', AgE=0.01)
      >>> saumon = Bdevopar('f2_solar', e='saumon', adf=0.01, MASSSS=0.4)

      >>> burrow = Bdevopar(E='burrow', L=0.05, Gradfty=0.4)
      >>> baraffe = Bdevopar(E='baraffe', temp=2000, age=0.4)

   .. note:: ``Keyword spelling doesn't matter as long as the first letter is
               right. Also, ordering of the keywords is not important. However,
               metallicites must be spelled exactly as mentioned aboved--there 
               are no execptions.``  
   """
   def __init__(self, *args, **kwargs):
      ################# Check User Inputted Keywords ###########################
      self.params = {'Emodel':0, 'Temp':0, 'Age':0, 'Grav':0, 'Mass':0, 
                     'Lumi':0, 'Radi':0}
      keys = kwargs.keys() # Make a list of all the user-inputted keywords.
      if len(keys) != 3: raise NameError('Exactly three keywords required.\n')
      if 'e' not in [keys[i][0].lower() for i in [0,1,2]]: 
         raise NameError('Need an evolutionary model: baraffe,burrow,saumon\n')
      for i in range(3): # 3 keywords required. 
         if list(keys[i])[0].upper() == 'T': # Check Temperature keyword. 
            if type(kwargs[keys[i]])==str:kwargs[keys[i]]=float(kwargs[keys[i]])
            if self.params['Temp'] == 0: self.params['Temp'] = kwargs[keys[i]]
            else: NameError('Only 1 temperature allowed.\n')
         elif list(keys[i])[0].upper() == 'A': # Check Age keyword.
            if type(kwargs[keys[i]])==str:kwargs[keys[i]]=float(kwargs[keys[i]])
            if self.params['Age'] == 0: self.params['Age'] = kwargs[keys[i]]
            else: raise NameError('Only 1 age allowed.\n')
         elif list(keys[i])[0].upper() == 'G': # Check Gravity keyword.
            if type(kwargs[keys[i]])==str:kwargs[keys[i]]=float(kwargs[keys[i]])
            if self.params['Grav'] == 0: self.params['Grav'] = kwargs[keys[i]]
            else: raise NameError('Only 1 gravity allowed.\n')
         elif list(keys[i])[0].upper() == 'R': # Check Radius keyword.
            if type(kwargs[keys[i]])==str:kwargs[keys[i]]=float(kwargs[keys[i]])
            if self.params['Radi'] == 0: self.params['Radi'] = kwargs[keys[i]]
            else: raise NameError('Only 1 radius allowed.\n')
         elif list(keys[i])[0].upper() == 'M': # Check Mass keyword.
            if type(kwargs[keys[i]])==str:kwargs[keys[i]]=float(kwargs[keys[i]])
            if self.params['Mass'] == 0: self.params['Mass'] = kwargs[keys[i]]
            else: raise NameError('Only 1 mass allowed.\n')
         elif list(keys[i])[0].upper() == 'L': # Check Luminosity keyword.
            if type(kwargs[keys[i]])==str:kwargs[keys[i]]=float(kwargs[keys[i]])
            if self.params['Lumi'] == 0: self.params['Lumi'] = kwargs[keys[i]]
            else: raise NameError('Only 1 luminosity allowed.\n')
         elif list(keys[i])[0].upper() == 'E': # Check Emodel keyword.
            if type(kwargs[keys[i]]) != str: 
               raise ValueError('Evolutionary model must be a string.\n')
            if self.params['Emodel'] == 0: 
               self.params['Emodel'] = kwargs[keys[i]]
               if self.params['Emodel'].lower() == 'baraffe': 
                  model = super(Bdevopar, self).baraffe_model()
               elif self.params['Emodel'].lower() == 'burrow': 
                  model = super(Bdevopar, self).burrow_model()
               elif self.params['Emodel'].lower()  == 'saumon': 
                  model = super(Bdevopar, self).saumon_model(args[0])
               else: raise NameError('This model does not exist.\n')
            else: raise NameError('Only 1 evolutionary model allowed.')
         else: raise NameError("Keyword '"+ keys[i] + "' is nonexistent.\n")

      if self.params['Temp'] != 0: 
         assert min(min(model['Temp'][:])) <= self.params['Temp'] <= \
                max(max(model['Temp'][:])),"Temperature is out of model's range"
      if self.params['Lumi'] != 0: 
         assert min(min(model['Lumi'][:])) <= self.params['Lumi'] <= \
                max(max(model['Lumi'][:])), "Luminosity is out of model's range"
      if self.params['Grav'] != 0: 
         assert min(min(model['Grav'][:])) <= self.params['Grav'] <= \
                max(max(model['Grav'][:])), "Gravity is out of model's range."
      if self.params['Radi'] != 0: 
         assert min(min(model['Radi'][:])) <= self.params['Radi'] <= \
                max(max(model['Radi'][:])), "Radius is out of model's range."
      if self.params['Mass'] != 0: 
         assert min(min(model['Mass'][:])) <= self.params['Mass'] <= \
                max(max(model['Mass'][:])), "Mass is out of model's range."
      if self.params['Age'] != 0: 
         assert min(model['Age'][:]) <= self.params['Age'] <= \
                max(model['Age'][:]), "Age is out of model's range."

      if self.params['Emodel'] == 'baraffe' and \
                    (self.params['Age'] != 0 and self.params['Mass'] != 0):
         if 0.001 <= self.params['Age'] <= 1.:
            if 0.0005 <= self.params['Mass'] <= 0.1: pass
            else: assert False, 'mass must be in between 0.0005 and 0.1'
         elif 1. <= self.params['Age'] <= 5.:
            if 0.002 <= self.params <= 0.1: pass
            else: assert False, 'mass must be in between 0.002 and 0.1'
         elif 5 <= self.params['Age'] <= 10:
            if 0.003 <= self.params['Mass'] <= 0.1: pass
            else: assert False, 'mass must be in between 0.003 and 0.1'
         else: assert False, 'age must be in between 0.001 and 10'

      #########################################################################
      Ag, Ma, Te, Le, Ge, Re = [],[],[],[],[],[]
      input_type = 'mass_age'
      valid_ages = []
      n_tables = model['n_tables']

      ############### WITH TWO KNOWN PARAMETERS, SPIT OUT AGE ##################
      if (self.params['Mass'] == False) and (self.params['Age'] == False):

          input_type = 'two_params'; P = []; Lno, Gno, Rno = True, True, True
          for i in [0,1]: # Because there's two input parameters
              if (self.params['Lumi'] != False) and Lno:  
                  P.append(['Lumi', self.params['Lumi']]); Lno=0
              elif (self.params['Grav'] != False) and Gno:   
                  P.append(['Grav', self.params['Grav']]); Gno=0
              elif (self.params['Radi'] != False) and Rno:      
                  P.append(['Radi', self.params['Radi']]); Rno=0
              else: P.append(['Temp', self.params['Temp']])

          for i in n_tables:
              if min(model[P[0][0]][i]) <= P[0][1] <= max(model[P[0][0]][i]) \
                and min(model[P[1][0]][i]) <= P[1][1] <= max(model[P[1][0]][i]):
                  Ag.append(model['Age'][i]) 
                  valid_ages.append(i)
                  f = interp1d(model[P[0][0]][i], model['Mass'][i])
                  Ma = f(P[0][1])
                  f = interp1d(model['Mass'][i], model[P[1][0]][i])
                  Ge.append(f(Ma))

          f = interp1d(Ge, Ag)
          try:
             self.params['Age'] = f(P[1][1])
          except ValueError:
             assert False, P[1][0] + " is above the interpolation range."

          del Ge[:], Ag[:]; Ma = []
 
      ################ WITH KNOWN AGE OR MASS AND ANOTHER PARAMETER ###########
      ################### SPIT OUT MASS OR AGE, RESPECTIVELY ##################
      if ((self.params['Age'] != False) and (self.params['Mass'] == False)) \
          or ((self.params['Mass'] != False) and (self.params['Age'] == False)):

          if input_type != 'two_params': 
              input_type = 'one_param'; P = []
              if self.params['Temp'] != False: 
                  P.append(['Temp', self.params['Temp']])
              elif self.params['Lumi'] != False:  
                  P.append(['Lumi', self.params['Lumi']])
              elif self.params['Grav'] != False:   
                  P.append(['Grav', self.params['Grav']])
              elif self.params['Radi'] != False:      
                  P.append(['Radi', self.params['Radi']])
              P.append(['', 0])

          if input_type == 'two_params': 
              n_t = valid_ages
              valid_ages = []
          else: n_t = model['n_tables']

          for i in n_t:
             if min(model[P[0][0]][i]) <= P[0][1] <= max(model[P[0][0]][i]):
                 valid_ages.append(i)
                 Ag.append(model['Age'][i]) 
                 f = interp1d(model[P[0][0]][i], model['Mass'][i])
                 Ma.append(f(P[0][1]))
          
          if self.params['Mass'] == False:
              f = interp1d(Ag, Ma)
              self.params['Mass'] = round(f(self.params['Age']), 5)
          else:
              f = interp1d(Ma, Ag)
              self.params['Age'] = round(f(self.params['Mass']), 5)
          del Ag[:]

      ###################### WITH KNOWN MASS AND AGE ##########################
      if (self.params['Mass'] != False) and (self.params['Age'] != False):

          if input_type == 'mass_age': valid_ages = model['n_tables']

          for i in valid_ages:
              if min(model['Mass'][i]) <= self.params['Mass'] \
                                      <= max(model['Mass'][i]):
                  Ag.append(model['Age'][i])
                  f =interp1d(model['Mass'][i],model['Temp'][i])
                  Te.append(f(self.params['Mass']))
                  f = interp1d(model['Mass'][i],model['Lumi'][i])
                  Le.append(float(f(self.params['Mass'])))
                  f = interp1d(model['Mass'][i],model['Grav'][i])
                  Ge.append(f(self.params['Mass']))
                  f = interp1d(model['Mass'][i],model['Radi'][i])
                  Re.append(f(self.params['Mass']))
      
      f = interp1d(Ag, Te) 
      self.params['Temp'] = round(f(self.params['Age']), 8)
      f = interp1d(Ag, Le)
      self.params['Lumi'] = round(f(self.params['Age']), 8)
      f = interp1d(Ag, Ge) 
      self.params['Grav'] = round(f(self.params['Age']), 8)
      f = interp1d(Ag, Re)
      self.params['Radi'] = round(f(self.params['Age']), 8)
      
      if input_type == 'one_param': self.params[P[0][0]] = P[0][1]
      elif input_type == 'two_params': 
          self.params[P[0][0]] = P[0][1]
          self.params[P[1][0]] = P[1][1]

   ############################################################################
   def all_params(self):
       """
       :Returns: ``This method returns a dictionary with all the parameters.
                   You may also retrieve only a parameter. See example below.``
   
       - **Example:**
          >>> baraffe = Bdevopar(E='baraffe', a=0.4, ma= 0.05)
          >>> params = baraffe.all_params()
          >>> print(params)
          {'Emodel': 'baraffe', 'Temp': 2023.94736842, 
          'Lumi': -3.76789474, 'Age': 0.4, 'Grav': 5.06442105, 
          'Mass': 0.05, 'Radi': 0.10994737}
          >>> print params['Mass'], params['Lumi']
          0.05 -3.76789474
       """
       return self.params

   ############################################################################
   def __str__(self):
       print "\n------------------------------------------"
       print "| BD's age (Gyr):            | " + str(self.params['Age'])
       print "------------------------------------------"
       print "| BD's mass (Msun):          | " + str(self.params['Mass'])
       print "------------------------------------------"
       print "| Effective temperature (K): | " + str(self.params['Temp'])
       print "------------------------------------------"
       print "| log L_bd/L_sun:            | " + str(self.params['Lumi'])
       print "------------------------------------------"
       print "| log surface gravity:       | " + str(self.params['Grav'])
       print "------------------------------------------"
       print "| BD's radius (Rsun):        | " + str(self.params['Radi'])
       print "------------------------------------------"
       return ''

########################## End of the class: bdevopar #########################
###############################################################################
#from random import uniform
#
#class Distribution(Bdevopar, EvolutionaryModels): 
#    def __init__(self): pass
#    def power_law(self, xmin=0.03, xmax=0.1, alpha=-1.3, nstars=10)
#        output = []
#        for i in range(nstars):
#           x = ((xmax**(alpha+1) - xmin**(alpha+1))**uniform(0,1) + \
#                xmin**(alpha+1))**(1./(alpha+1))
#           output.append(x)
#        return output

