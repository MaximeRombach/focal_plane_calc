#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import norm
from numpy.polynomial import Polynomial
from astropy.table import Table, vstack
import json
import logging
from shapely import affinity, MultiPolygon, MultiPoint, GeometryCollection
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, box
from shapely.ops import unary_union
from shapely.validation import make_valid
from scipy import optimize
from scipy.interpolate import interp1d
from shapely.plotting import plot_polygon, plot_points
from types import SimpleNamespace  # Python 3.3+ only.
import os
from datetime import datetime
import geopandas as gpd
import pandas as pd
import math
from shapely.geometry import mapping
import ezdxf
import ezdxf.addons.geo
import updown_tri as tri
import pandas as pd

logging.basicConfig(level=logging.INFO)
tan30 = np.tan(np.deg2rad(30))
sqr3 = np.sqrt(3)

## Parameters ##
# NOTE: FocalSurf class
""" Focal plane parameters """
class FocalSurf():
     def __init__(self, project: str) -> None:

          self.project = project
          self.focal_surf_param = self.set_surface_parameters()

          self.curvature_R = self.focal_surf_param['curvature_R'] #curvature radius of focal plane
          self.vigR = self.focal_surf_param['vigR'] # vignetting radius
          self.BFS = self.focal_surf_param['BFS']
          self.f_number = self.focal_surf_param['f-number']
          
          self.asph_formula = self.focal_surf_param['asph_formula'] # Checks the presence or not of aspherical coefficients
          if self.asph_formula: # Store spherical focal plane parameters if known
               ## TODO: use that to unwrap focal surf param --> replace everything call in the code
               # self.asph_coeff = SimpleNamespace(**self.focal_surf_param)
               self.a2 = self.focal_surf_param['a2']
               self.a4 = self.focal_surf_param['a4']
               self.a6 = self.focal_surf_param['a6']
               self.a8 = self.focal_surf_param['a8']
               self.k = self.focal_surf_param['k']
               self.c = 1/self.curvature_R
          else:
               self.read_focal_plane_data()

          if 'WST' in self.project:
               self.donutR = self.focal_surf_param['donutDiam']/2

          self.surfaces_polygon = {} # dictionnary to store the polygons of the focal plane surfaces (main vigR, GFA, donut hole, etc)
          self.surf_name = self.focal_surf_param['name']
          self.FoV = self.focal_surf_param['FoV']
          self.focus_tolerance_width = self.focal_surf_param['focus_tolerance_width'] # [mm] width of the focus tolerance band around the focal surface

          # Transfer functions
          self.R2Z, self.R2CRD, self.R2NORM, self.R2S, self.R2NUT, self.S2R = self.transfer_functions()

     def set_surface_parameters(self):

          if self.project == 'MUST':
               ## Last update of MUST focal plane parameters (update from: 2024-01-29)
               focal_surf_param = {
                                   'name': 'MUST',
                                   'curvature_R': -11918, # [mm], curvature radius
                                   'vigD': 1184.7, # [mm], vignetting diameter
                                   'vigR': 1184.7/2, # [mm], vignetting radius
                                   'asph_formula': False,
                                   'k': 0,
                                   'a2': 0, 
                                   'a4': -3.913848e-11, 
                                   'a6': 5.905507e-17,
                                   'a8': 0,
                                   'f-number': 3.72059,
                                   'BFS': 10477.594, # [mm], radius of BFS # Value calculated with Joe's cal_BFS function with MUST data from 2023-09-06
                                   'FoV': None, # [deg]
                                   'focus_tolerance_width': 0.02 # [mm]
                                   }
               ## Previous values for MUST focal plane parameters (update from: 2023-11-08)
               # focal_surf_param = {
               #                     'name': 'MUST',
               #                     'curvature_R': -11365, # [mm], curvature radius
               #                     'vigD': 1068, # [mm], vigneting diameter: last updated 2023-11-08 (previous value: 589.27 * 2)
               #                     'vigR': 1068/2, # [mm], vignetting radius
               #                     'asph_formula': True,
               #                     'k': 0,
               #                     'a1': 0, 
               #                     'a2': -6e-12, 
               #                     'a3': 0,
               #                     'f-number': 3.699,
               #                     'BFS': 10992.7, # [mm], radius of BFS # Value calculated with Joe's cal_BFS function with MUST data from 2023-09-06
               #                     'FoV': None, # [deg]
               #                     'focus_tolerance_width': 0.02 # [mm]
               #                     }
          elif self.project == 'MegaMapper':
               focal_surf_param = {
                                   'name': 'MegaMapper',
                                   'curvature_R': -11067, 
                                   'vigR': 613.2713,
                                   # 'vigR': 440,
                                   'asph_formula': False,
                                   'BFS': 11045.6, # [mm], radius of BFS
                                   'f-number': 3.57,
                                   'FoV': None, # [deg]
                                   'focus_tolerance_width': 0.1 # [mm]
                                   }

          elif self.project == 'DESI':
               focal_surf_param = {
                                   'name': 'DESI',
                                   'curvature_R': -11067, 
                                   'vigR': 406.,
                                   'asph_formula': False,
                                   'BFS': 11067, # [mm], radius of BFS,
                                   'f-number': 3.699,
                                   'FoV': None, # [deg]
                                   'focus_tolerance_width': 0.1 # [mm]
                                   }
               
          elif self.project == 'WST1':
               focal_surf_param = {
                                   'name': r'$\bf{WST - Focal Plane \varnothing: 1.2 m - Center \varnothing: 0.17m}$',
                                   'curvature_R': -11067, 
                                   'vigR': 600, # [mm], vignetting radius
                                   'donutDiam': 20, # [arcmin], diameter of donut hole
                                   'asph_formula': False,
                                   'BFS': 11067, # [mm], radius of BFS,
                                   'f-number': 3.699, # assumption based on previous telescope designs
                                   'FoV': 1.8, # [deg]
                                   'focus_tolerance_width': None # [mm]
                                   }
               
          elif self.project == 'WST2':
               focal_surf_param = {
                              'name': r'$\bf{WST - Focal Plane \varnothing: 1.4 m - Center \varnothing: 0.2m}$',
                              'curvature_R': -11067, 
                              'vigR': 700,# [mm], vignetting radius
                              'donutDiam': 20, # [arcmin], diameter of donut hole
                              'asph_formula': False,
                              'BFS': 11067, # [mm], radius of BFS,
                              'f-number': 3.7, # assumption based on previous telescope designs
                              'FoV': 1.8, # [deg]
                              'focus_tolerance_width': None # [mm]
                              }

          elif self.project == 'WST3':
               focal_surf_param = {
                              'name': r'$\bf{WST - Focal Plane \varnothing: 1.3 m - Center \varnothing: 0.24m (20 arcmin)}$',
                              'curvature_R': -11067, 
                              'vigR': 650,# [mm], vignetting radius
                              'donutDiam': 20, # [arcmin], diameter of donut hole
                              'asph_formula': False,
                              'BFS': 11067, # [mm], radius of BFS,
                              'f-number': 3.7, # assumption based on previous telescope designs
                              'FoV': 1.8, # [deg]
                              'focus_tolerance_width': None # [mm]
                              }
               
          elif self.project == 'Spec-s5-old':
               focal_surf_param = {
                    'name': r'$\bf{Spec-s5}$',
                    'curvature_R': -13000, # [mm], main radius of curvature (roughly 13m from Claire's info 20.10.2023)
                    'vigR': 4.084421E+02,# [mm], vignetting radius
                    'asph_formula': False,
                    'BFS': 1.277364E+04, # [mm], radius of BFS,
                    'f-number': 3.64, # assumption based on previous telescope designs
                    'FoV': None, # [deg],
                    'focus_tolerance_width': None # [mm]
                    }
               
          elif self.project == 'Spec-s5':
               #2024.02.16 update: New focal surface parameters for Spec-s5 (Pat Jelinsky new design - SPHERICAL focal surface youpi)
               focal_surf_param = {
                    'name': r'$\bf{Spec-s5}$',
                    'curvature_R': None, # [mm], main radius of curvature (roughly 13m from Claire's info 20.10.2023)
                    'vigR': 409.4, # [mm], vignetting radius
                    'asph_formula': False,
                    'BFS': 12657, # [mm], radius of BFS,
                    'f-number': 3.62002, # assumption based on previous telescope designs
                    'FoV': None, # [deg],
                    'focus_tolerance_width': None # [mm]
                    }
          else: 
               logging.error(f'Setting focal surface parameters: {self.project} project is not defined')
               raise Exception(f'Error in setting focal surface parameters: {self.project} project is not defined')

          return focal_surf_param
     
     def read_focal_plane_data(self):

               filename = f"./Data_focal_planes/{self.project}.csv" # optics data from Zemax
               comment_character = "#"  # The character that indicates a commented line
               # Read CSV file and ignore commented lines
               if self.project == 'MUST':
                    sep = ',' # keeping distinction as delimiter might change from time to time
               elif self.project == 'MegaMapper':
                    sep = ','
               else :
                    sep = ','
               optics_data = pd.read_csv(filename, comment=comment_character, sep=sep)

               # Set headers using the first non-comment line
               with open(filename, 'r') as file:
                    for line in file:
                         if not line.startswith(comment_character):
                              headers = line.strip().split(sep)
                              break
               optics_data.columns = headers
               
               print(f"{self.project} focal plane data successfully read from {filename}")
               print(optics_data.head())
               self.optics_data = optics_data
     
     def transfer_functions(self):
          # Those transfer functions and there implementation later are inspired from the work of Joseph Silber (LBNL) in the generate_raft_layout.py code in https://github.com/joesilber/raft-design 

          """ 
          Input: 
               - optics_data: pandas dataframe containing the focal plane data from zmax csv file
               - analytical: [bool] flag to use analytical solution of R2Z instead of data 

          Output:
               - R2Z: functions to convert radius to height on focal surface
               - R2CRD: function to get chief ray deviation interms of radial position on focal surface
          
          ATTENTION: Needs csv file to contain columns named: R, Z, CRD """
          r = np.linspace(0,self.vigR,2000) # [mm] radial positions on focal plane
          if not self.asph_formula:
               logging.info('Transfer functions made from sampled data in csv')
               R = self.optics_data['R']
               # DEBUG check: check for duplicate in R that causes interp1d to fail, can't have similar x values
               if R.duplicated().any():
                    R_dup = np.sum(R.duplicated())
                    raise Exception(f'Found duplicate in R: {R_dup} rows of csv file. Check csv file for duplicate entries and remove them.')
               
               # R2Z maps the radial position on the focal plane to the height of the focal plane
               if 'Z' in self.optics_data.keys():
                    Z = self.optics_data['Z']
                    R2Z = interp1d(R,Z,kind='cubic', fill_value = "extrapolate") #leave 'cubic' interpolation for normal vectors calculations
               else:
                    R2Z = lambda r: -self.BFS + np.sqrt(self.BFS**2 - r**2) # Spherical focal surface
                    logging.warning('No Z data available in samples - assuming SHERICAL surface')

               # R2CRD maps the radial position on the focal plane to the chief ray deviation (CRD) on the focal plane
               if 'CRD' in self.optics_data.keys():
                    CRD = self.optics_data['CRD']
                    R2CRD = interp1d(R,CRD,kind='cubic', fill_value = "extrapolate")
                    
                    logging.info('CRD data available in samples - using it for CRD transfer function')
               else:
                    CRD = np.zeros_like(R)
                    R2CRD = interp1d(R,CRD,kind='cubic', fill_value = "extrapolate")
                    logging.warning('No CRD data available - CRD = 0 everywhere on focal plane')
               crd = R2CRD(r)
          else:
               logging.info('Transfer functions: using aspherical formula cofficients, no CRD data available')
               R2Z = self.analytical_R2Z
               R2CRD = None
               crd = np.zeros_like(r)

          z = R2Z(r) # Calculate focal plane curve from csv data
          dr = np.diff(r)
          dz = np.diff(z)
          dzdr = dz/dr # Get derivative of focal plane curve at the sample points

          # Assert in slope data is available in csv file
          if 'Slope' in self.optics_data.keys():
               norm = self.optics_data['Slope'] # Calculate normal angles at the sample points
               R2NORM = interp1d(R, norm,kind='cubic', fill_value = "extrapolate")
               logging.info('Slope data available in samples - using it for R2NORM transfer function')
          else:
               norm = np.degrees(np.arctan(dzdr))# Calculate normal angles at the sample points
               R2NORM = interp1d(r[:-1], norm,kind='cubic', fill_value = "extrapolate")
               NORM2R = interp1d(norm, r[:-1])
               logging.info('No slope data available in samples - deriving normal angles from focal plane curve')

          nut = -(R2NORM(r) + crd) # Joe: -(R2NORM(r) + crd)
          R2NUT = interp1d(r, nut, kind='cubic', fill_value = "extrapolate")
          NUT2R = interp1d(nut, r)

          # Parametric path length along the focal surface
          ds = (1 + dzdr**2)**0.5 * dr
          s = np.cumsum(ds)
          s = np.insert(s, 0, 0.)  # first value of path length is 0
          Z2R = interp1d(z, r)
          R2S = interp1d(r, s)
          S2R = interp1d(s, r, kind='cubic', fill_value = "extrapolate")
     
          return R2Z, R2CRD, R2NORM, R2S, R2NUT, S2R

     def analytical_R2Z(self, r: np.array):
          """ 
          Returns the aspherical focal plane curve from the aspherical coefficients

          Input:
                    - r: [np.array] radial positions on focal plane
          Output:   
                    - z: [np.array] height positions on focal plane
          """
          if self.focal_surf_param['asph_formula']:
               main_term = self.c * np.power(r, 2) / (1 + np.sqrt(1 - (1 + self.k) * self.c**2 * np.power(r, 2)))
               secondary_terms = self.a2 * np.power(r, 2) + self.a4 * np.power(r, 4) + self.a6 * np.power(r, 6) + self.a8 * np.power(r, 8)
               return  main_term + secondary_terms
          else:
               return print('Aspherical coefficients not defined for this project \n Taking sampled data from csv file instead')
          
     
     def calc_BFS(self, r, z):

          # Courtesy of Joseph Silber from generate_raft_layout.py code in https://github.com/joesilber/raft-design 
          # best-fit sphere
          calc_sphR = lambda z_ctr: (r**2 + (z - z_ctr)**2)**0.5
          def calc_sphR_err(z_ctr):
               sphR_test = calc_sphR(z_ctr) 
               errors = sphR_test - np.mean(sphR_test) 
               scalar_error = np.sum(np.power(errors, 2)) # metric for optimization 
               return scalar_error
          typical_fov = 3.7  # deg
          z_guess = np.sign(np.mean(z)) * np.max(r) / np.radians(typical_fov/2)
          result = optimize.least_squares(fun=calc_sphR_err, x0=z_guess)
          z_ctr = float(result.x) # signed BFS radius
          sphR = abs(z_ctr) # unsigned BFS radius
          # is_convex = np.sign(z_ctr) < 1  # convention where +z is toward the fiber tips
          # sphR_sign = -1 if is_convex else +1
          logging.info(f"Best-fit sphere radius: {sphR:.3f} mm")

          return sphR
     
     def make_vigR_polygon(self, trimming_angle = 360, n_vigR = 500, changed_vigR = None):
          
          if changed_vigR is not None:
               vigR = changed_vigR # [mm] allows to change the vignetting radius for specific cases
          else:
               vigR = self.vigR # [mm] vignetting radius, nominal case declared for each project

          vigR_lim_x = vigR * np.cos(np.deg2rad(np.linspace(0,trimming_angle,n_vigR)))
          vigR_lim_y = vigR * np.sin(np.deg2rad(np.linspace(0,trimming_angle,n_vigR)))
          if trimming_angle == 360:
               end_point = [vigR_lim_x[0], vigR_lim_y[0]]
          else:
               end_point = [0, 0]
          vigR_lim_x = np.insert(vigR_lim_x, 0, end_point[0])
          vigR_lim_y = np.insert(vigR_lim_y, 0, end_point[1])
          pizza = Polygon(to_polygon_format(vigR_lim_x, vigR_lim_y))

          if self.project == 'WST' or self.project == 'WST1' or self.project == 'WST2' or self.project == 'WST3':
               # WST has a donut hole in the center of the focal plane for big fiber bundle
               donut = Polygon(to_polygon_format(vigR_lim_x*self.arcmin2mm(self.donutR)/vigR, vigR_lim_y*self.arcmin2mm(self.donutR)/vigR))
               pizza = pizza.difference(donut)
               self.surfaces_polygon['donut'] = donut
               
          self.surfaces_polygon['pizza'] = pizza

          return pizza

     def plot_vigR_poly(self, pizza, label = None, ax = None):
          plot_polygon(pizza, ax = ax, add_points = False, edgecolor = 'black', linestyle = '--', facecolor= 'None', label = label)

     def arcmin2mm(self, arcmin):
          return arcmin * 2 * self.vigR / (self.FoV * 60)
     
     def mm2arcmin(self, mm):
          return mm * (self.FoV * 60) / (2 * self.vigR)
     
     def mm2arcsec(self, mm):
          return mm * (self.FoV * 3600) / (2 * self.vigR)
     
     def arcsec2mm(self, arcsec):
          return arcsec * 2 * self.vigR / (self.FoV * 3600)
     
     def arcmin2arcsec(self, arcmin):
          return arcmin / 60
     
     def arcsec2arcmin(self, arcsec):
          return arcsec * 60

# NOTE : Saving results class
class SavingResults:
     """
     Class for saving the various results produced with focal_plane_coverage.py
     Initialized with "saving_df" dictionnary containing booleans for:

     - Saving final plots (save_plots) of the focal plane
     - Saving the layout of the frame in dxf format (save_dxf) for further SolidWorks results
     - Saving csv files produced along the way (e.g. xy positions of all individual robots)

     It also creates the "Results/" directory, if not already existing, to store eveything in one place
     """

     def __init__(self, saving_df, project_name = None) -> None:

          self.project_name = project_name
          self.results_dir_path = self.path_to_results_dir()
          if 'save_plots' in saving_df.keys():
               self.save_plots = saving_df['save_plots']

               
          if 'save_dxf' in saving_df.keys():
               self.save_dxf = saving_df['save_dxf']


          if 'save_csv' in saving_df.keys():
               self.save_csv = saving_df['save_csv']


          if 'save_txt' in saving_df.keys():     
               self.save_txt = saving_df['save_txt']   

     def path_to_results_dir(self):

          script_dir = os.path.dirname(__file__)
          # results_dir_path = os.path.join(script_dir, 'Results_examples/')
          results_dir_path = os.path.join(script_dir, 'Results/')

          # Creates a subfolder corresponding to the project name in Results/
          if self.project_name is not None: 
               results_dir_path = os.path.join(results_dir_path, self.project_name + '/')

          # Checks if the directory exists, if not creates it
          if not os.path.isdir(results_dir_path):
               os.makedirs(results_dir_path)
          
          return results_dir_path
     
     def save_dxf_to_dir(self, geometries: dict, suffix_name):

          if not self.save_dxf:
               return
          
          now = datetime.now()
          name_frame = now.strftime("%Y-%m-%d--%H-%M-%S_") + suffix_name

          doc = ezdxf.new()
          msp = doc.modelspace()
          for key, geom in geometries.items():

               geoproxy = ezdxf.addons.geo.GeoProxy.parse(mapping(geom))

               # Use LWPOLYLINE instead of hatch.
               for entity in geoproxy.to_dxf_entities(polygon=2):
                    msp.add_entity(entity)
                    entity.set_dxf_attrib('layer', f"{key}")

          doc.saveas(self.results_dir_path + f"{name_frame}.dxf")

     def save_figures_to_dir(self, image_name: str, file_format: str = 'png', save_eps: bool = False, dpi: int = 400):

          if not self.save_plots:
               return

          now = datetime.now()
          today_filename = now.strftime("%Y-%m-%d-%H-%M-%S_") + image_name + f'.{file_format}'
          plt.savefig(self.results_dir_path + today_filename, bbox_inches = 'tight', format=file_format, dpi = dpi)
          if save_eps:
               plt.savefig(self.results_dir_path + today_filename, bbox_inches = 'tight', format='eps')

          logging.info(f'{image_name}.png successfully saved in in {self.results_dir_path}')

     def save_grid_to_txt(self, grid, filename, direct_SW = False):

          if not self.save_txt:
               return
          now = datetime.now()
          # grid MUST be a (N,4) numpy array
          #TODO: modify that function to accept any number of columns and take headers as input for column names
          x = grid[:,0]
          y = grid[:,1]
          z = grid[:,2]
          if not direct_SW:
               up_tri = grid[:,3]
          with open(self.results_dir_path + now.strftime("%Y-%m-%d-%H-%M-%S_") + f'{filename}.txt', 'w') as file:
               if direct_SW: # removes orientation column to directly read cloud point in solidworks
                    file.write("x[mm] y[mm] z[mm]\n")
                    for (dx,dy,dz) in zip(x,y,z):
                         file.write(f"{dx:.3f} {dy:.3f} {dz:.3f}\n")
               else:
                    file.write("x[mm] y[mm] z[mm] upward_tri [bool]\n")
                    for (dx,dy,dz,up) in zip(x,y,z, up_tri):
                              file.write(f"{dx:.3f} {dy:.3f} {dz:.3f} {int(up)}\n")
          
          logging.info(f'{filename}.txt succesfully saved in {self.results_dir_path}')

     def save_grid_to_txt2(self, grid: pd, filename: str, **kwargs):

          columns = kwargs.get('columns', None)
          index = kwargs.get('index', False)
          if not self.save_txt:
               return
          now = datetime.now()
          text_representation = grid.to_string(header= True, columns=columns, index=index)
          # Write string to file
          with open(self.results_dir_path + now.strftime("%Y-%m-%d-%H-%M-%S_") + f'{filename}.txt', 'w') as file:
               file.write(text_representation)

          logging.info(f'{filename}.txt succesfully saved in {self.results_dir_path}')

""" Module parameters """ 

# NOTE: Individual Module class
class Module(SavingResults):

     def __init__(self, nb_robots, saving_df, BFS: float, is_wall: bool = True, width_increase = 0, chanfer_length = 7.5):

          """Robot parameters""" 

          self.nb_robots = nb_robots
          self.la = 1.8 # alpha arm length [mm] /!\ lb > la /!\
          self.lb = 1.8 # beta arm length [mm]
          self.pitch = 6.2 # pitch [mm]
          self.alpha = np.linspace(-180,180,180) # [deg] rotational range of alpha arm
          self.beta = np.linspace(-180,180,180) # [deg] rotational range of beta arm (holding the fiber)
          self.beta2fibre = 1 # [mm] distance fiber center to edge of beta arm
          
          
          @property
          def robots_module_table(self):
               """ Table of robots centroids within the module
               Output: [nb_robots x 3] array with columns: x, y, z """
               return self._robots_module_table


          @property
          def robots_module_table(self):
               """ Table of robots centroids within the module
               Output: [nb_robots x 3] array with columns: x, y, z """
               return self._robots_module_table

          """Module parameters"""

          self.chanfer_length = chanfer_length # [mm] Length of chanfer on module vertices
          self.width_increase = width_increase # [mm] Increase module side length w/o changing nb of pos
          
          self.x_first = 5.9 + self.width_increase/2 # [mm] Horizontal pos of first robot (bottom left)
          self.y_first = 3.41 + self.width_increase/2 * tan30 # [mm] Vertical pos of first robot (bottom left)
          
          self.x_inc = 3.1 # [mm] Horizontal increment at each row
          self.y_inc = 5.369 # [mm] Vertical increment at each row
          self.test_pitch = norm(np.array([self.x_inc,self.y_inc]))

          self.safety_distance = 0.2 # [mm] physical distance kept between shields and edge of beta arm
          self.offset_from_module_edges = self.safety_distance + self.beta2fibre
          self.start_offset_x = 6.2 # [mm]
          self.start_offset_y = 3.41 # [mm]

          self.is_wall = is_wall # flag for protective shields or not on modules

          self.module_length = 590 # [mm] last dimension updae for length of module unit
          self.BFS = BFS # [mm] radius of BFS

          @property
          def module_centroid(self):
               return self._module.centroid
          
          @property
          def base_z(self):
               return self.BFS

          # Add 1 row of positioners from one case to another
          if self.nb_robots == 42:

               self.module_width = 62.4 + self.width_increase # [mm] triangle side length
               self.nb_rows = 8 # number of rows of positioners
               self.key = 'n42'
          
          elif self.nb_robots == 52:

               self.module_width = 67.6 + self.width_increase # [mm] triangle side length
               self.nb_rows = 9 # number of rows of positioners
               self.key = 'n52'

          elif self.nb_robots == 63:
               self.module_width = 73.8 + self.width_increase # [mm] triangle side length
               self.width_hole_in_frame = self.module_width + self.width_increase # [mm] ONGOING TEST: width of the module hole in the frame, represents the clearance between the module and the frame --> influences coverage
               self.nb_rows = 10 # number of rows of positioners
               self.key = 'n63'

          elif self.nb_robots == 75:

               self.module_width = 80 + self.width_increase # [mm] triangle side length
               self.nb_rows = 11 # number of rows of positioners
               self.key = 'n75'

          elif self.nb_robots == 88:
               self.module_width = 86.2 + self.width_increase # [mm] triangle side length
               self.nb_rows = 12 # number of rows of positioners
               self.key = 'n88'

          elif self.nb_robots == 102:
               
               self.module_width = 92.4 + self.width_increase # [mm] triangle side length
               self.nb_rows = 13 # number of rows of positioners
               self.key = 'n102'

          else:
               raise Exception('Error: only 42, 52, 63, 75, 88, 102 robots per module supported')
          
          self.module_collection, self.multi_wks_list, self.coverages = self.create_module()
          self.module = self.module_collection.geoms[0]
          self.module_w_beta_and_safety_dist = self.module_collection.geoms[1]
          self.effective_wks = self.module_collection.geoms[2]
          self.triang_meshgrid = self.module_collection.geoms[3]

          super().__init__(saving_df)
          
     def equilateral_vertices(self, triangle_width, xoff = 0, yoff = 0):
          """
          Create the raw triangular base of a module w/ the correct width

          Output: sets of x,y coords of the triangle's vertices
          """
          x_vertices = np.ones(4)
          offset = [xoff, yoff]
          x_vertices[0] = offset[0]
          x_vertices[1] = x_vertices[0] + triangle_width
          x_vertices[2] = x_vertices[0] + triangle_width * np.cos(np.deg2rad(60))
          x_vertices[3] = x_vertices[0]


          y_vertices = np.ones(4)
          y_vertices[0] = offset[1]
          y_vertices[1] = y_vertices[0]
          y_vertices[2] = y_vertices[0] + triangle_width * np.sin(np.deg2rad(60))
          y_vertices[3] = y_vertices[0]

          return x_vertices,y_vertices

     def chanfered_base(self):
          """
          1) Creation of small triangles that will serve to cut the raw triangle for chanfering
          2) Cuts triangular base with small triangles to create the chanfers
          
          Output: sets of x,y coords of the chanfered module's vertices
          """
          
          x_vertices,y_vertices = self.equilateral_vertices(self.module_width)
          module_triangle = Polygon(to_polygon_format(x_vertices,y_vertices))

          x_chanfer,y_chanfer = self.equilateral_vertices(self.chanfer_length)

          # First chanfer created from the origin
          chanfer_triangle = Polygon(to_polygon_format(x_chanfer-0.0001,y_chanfer-0.0001)) # tiny diff for geometry to intersect in one point rather tha infinite points (there is a better solution)

          chanfers = [chanfer_triangle]
          angles = [120, 240]
          module_centroid = module_triangle.centroid
          
          for angle in angles:
               # Move the 2 remaining triangles to their corresponding location at the triangle's edges
               rot_chanfer_triangle = rotate_and_translate(chanfer_triangle, angle, 0, 0, origin = module_centroid)
               chanfers.append(rot_chanfer_triangle)
               chanfered_base = module_triangle.difference(rot_chanfer_triangle)
               

          multi_chanfers = MultiPolygon(chanfers)
          # Apply cut on triangular base
          chanfered_base = module_triangle.difference(multi_chanfers)

          # TODO: finish proper implementation of chanfered base
          # chanfered_vertices = np.zeros((7,3))
          # chanfered_vertices[0] = [-self.chanfer_length/2, sqr3 * (self.module_width/3 - self.chanfer_length/2), 0]
          # chanfered_vertices[:] = chanfered_vertices[0] #start and end point must be the same for polygon creation
          # chanfered_vertices[1] = [-chanfered_vertices[0][0], chanfered_vertices[0][1], 0]

          # for idx, ang in enumerate([240, 120]):
          #      idx = idx + 1
          #      chanfered_vertices[2*idx] = rotation3D(chanfered_vertices[0].T, 0, np.deg2rad(ang))
          #      chanfered_vertices[2*idx + 1] = rotation3D(chanfered_vertices[1].T, 0, np.deg2rad(ang))

          # print(chanfered_vertices)
          # figure = plt.figure(figsize=(8,8))
          # chanfered_base2 = Polygon(to_polygon_format(chanfered_vertices[:,0], chanfered_vertices[:,1]))
          # plot_polygon(chanfered_base2)
          # plt.show()
          
          return chanfered_base
     
     def remove_positioner(self, xx, yy, list_to_remove):
          """ Input:
               - xx, yy: grid of center of positioners
               - list_to_remove: list of positioners position to remove

               Output:
               - xx_sliced, yy_sliced: grid of center of positioners minus the ones removed
          """
          xx_sliced = np.delete(xx,list_to_remove)
          yy_sliced = np.delete(yy,list_to_remove)

          return xx_sliced, yy_sliced
     
     def create_module(self):
          
          module = self.chanfered_base() # Create chamfered module shape as shapely Polygon
          # coords_module_x, coords_module_y = module.exterior.coords.xy
          reference_centroid = module.centroid

          # Their norm corresponds to the pitch defined in parameters
          xx1 = np.ones(self.nb_rows + 1)
          yy1 = np.ones(self.nb_rows + 1)
          zz1 = self.BFS * np.ones(self.nb_rows + 1)

          for idx in range(self.nb_rows): # Create the grid of positioner's center points
          
               # Place the first robot of each row
               start_offset_x = self.x_inc * idx + self.x_first
               start_offset_y = self.y_inc * idx + self.y_first

               # Generate each row of robot: nb robots/row decreases as going upward in triangle
               xx_new = np.linspace(start_offset_x, (self.nb_rows - idx ) * self.pitch + start_offset_x, self.nb_rows - idx + 1)
               yy_new = start_offset_y * np.ones(len(xx_new))

               if idx == 0:
                    xx1 = xx_new
                    yy1 = yy_new
               else: 
                    xx1 = np.hstack((xx1, xx_new))
                    yy1 = np.hstack((yy1, yy_new))


               if idx == self.nb_rows - 1:
                    xx_new = np.array([self.x_inc * (idx + 1)])
                    yy_new = self.y_inc * (idx + 1) * np.ones(len(xx_new))
                    xx1 = np.hstack((xx1, xx_new))
                    yy1 = np.hstack((yy1, yy_new))

          list_to_remove = [0, self.nb_rows, -1] # remove positioners at the edges of the triangle
          # list_to_remove = [] # remove no positioner
          xx1, yy1 = self.remove_positioner(xx1, yy1, list_to_remove)
          nbots = len(xx1)

          triang_meshgrid = MultiPoint(to_polygon_format(xx1, yy1)) # convert the meshgrid into shapely standard for later manipulation

          """ 1)b) Define coverage for 1 modules """

          c1 = np.cos(np.deg2rad(self.alpha))
          s1 = np.sin(np.deg2rad(self.alpha))

          c2 = np.cos(np.deg2rad(self.beta))
          s2 = np.sin(np.deg2rad(self.beta))

          xa, ya, xab, yab = (self.lb-self.la)*c1, (self.lb-self.la)*s1, (self.lb+self.la)*c2, (self.la+self.lb)*s2

          wks_list = []

          for idx, (dx, dy) in enumerate(zip(xx1, yy1)):

               xa1, ya1, xab1, yab1 = xa + dx, ya + dy, xab + dx, yab + dy

               coords_int = to_polygon_format(xa1, ya1)
               coords_ext = to_polygon_format(xab1, yab1)

               interior = coords_int[::-1]
               poly_c1 = Polygon(coords_ext, [interior])
               wks_list.append(poly_c1)

          multi_wks = MultiPolygon(wks_list)
          total_positioners_workspace = unary_union(wks_list)
          module_w_beta_and_safety_dist = module.buffer(-self.offset_from_module_edges)

          if self.is_wall:
               effective_wks = module_w_beta_and_safety_dist.intersection(total_positioners_workspace)
          else:
               effective_wks = total_positioners_workspace


          module_collection = GeometryCollection([module, module_w_beta_and_safety_dist, effective_wks, triang_meshgrid])
          module_collection = affinity.translate(module_collection, xoff = -reference_centroid.x, yoff = -reference_centroid.y)
          multi_wks_list = affinity.translate(multi_wks, xoff = -reference_centroid.x, yoff = -reference_centroid.y)
          coverage_with_walls = round(effective_wks.area/module.area * 100,1)
          coverage_no_walls = round(total_positioners_workspace.area/module.area * 100,1)
          coverages = [coverage_with_walls, coverage_no_walls]

          logging.info(f'Module object created with width {self.module_width} mm & {nbots} robots')
          return module_collection, multi_wks_list, coverages
     
     def plot_raw_module(self):

          fig = plt.figure(figsize=(8,8))
          figtitle = f"Module coverage raw - {self.nb_robots} robots per module"
          plt.title(figtitle)
          plot_polygon(self.module, facecolor='None', edgecolor='black', add_points=False)
          plot_polygon(self.module_w_beta_and_safety_dist, facecolor='None', edgecolor='red', linestyle = '--', add_points=False
                       , label = "Safety dist = {} mm".format(self.offset_from_module_edges))
          plot_points(self.triang_meshgrid, marker='.', color='k', label = "{} robots".format(self.nb_robots))
          
          self.save_figures_to_dir(self.save_plots, figtitle)

# NOTE: Intermediate Triangle class
class IntermediateTriangle: 

     def __init__(self, module_collection, module_width, intermediate_frame_thick):

          self.module_collection = module_collection
          self.dist_inter = 2*module_width*np.sqrt(3)/6 + intermediate_frame_thick # distance between each neighbor from center module
          self.angles = np.array([-30, 90, 210])
          self.upward_tri = [False,True,True,True]

          self.x_grid_inter = np.cos(np.deg2rad(self.angles))*self.dist_inter
          self.x_grid_inter = np.insert(self.x_grid_inter, 0, 0)
          self.y_grid_inter = np.sin(np.deg2rad(self.angles))*self.dist_inter
          self.y_grid_inter = np.insert(self.y_grid_inter, 0, 0)

     def create_intermediate_triangle(self):
          
          boundaries = []
          boundary_xx = []
          boundary_yy = []
          intermediate_collection = []
          inter_coverage = []
          covered_area = 0
          inter_boundaries_df = {'name':[],'geometry':[], 'color': []}
          inter_modules_df = {'geometry':[], 'centroids': [], 'upward_tri':[]}
          inter_coverage_df = {'geometry':[]}
          inter_robots_df = {'geometry':[]}
          

          for idx, (up_tri, dx, dy) in enumerate(zip(self.upward_tri, self.x_grid_inter, self.y_grid_inter)):

               if not up_tri:
                    angle = 180
               else:
                    angle = 0

               transformed_all = rotate_and_translate(self.module_collection, angle, dx, dy, origin = "centroid") # place module at the correct pos and orientation
               boundaries.append(transformed_all.geoms[0]) # log the exterior boundary points of the transformed module
               inter_modules_df['geometry'].append(transformed_all.geoms[0])
               inter_modules_df['centroids'].append([dx, dy])
               inter_modules_df['upward_tri'].append(up_tri)
               inter_coverage.append(transformed_all.geoms[2])
               inter_coverage_df['geometry'].append(transformed_all.geoms[2])
               intermediate_collection.append(transformed_all)
               inter_robots_df['geometry'].append(transformed_all.geoms[3])
               
               xx,yy = transformed_all.geoms[0].exterior.coords.xy
               boundary_xx.append(xx.tolist())
               boundary_yy.append(yy.tolist())

          
          inter_modules_df['centroids'] = np.array(inter_modules_df['centroids'])
          inter_robots_df['geometry'] = GeometryCollection(inter_robots_df['geometry'])
          modules_polygon_intermediate = MultiPolygon(boundaries)
          # Convex hull of boundaries
          bounding_polygon_intermediate = modules_polygon_intermediate.convex_hull
          inter_boundaries_df['name'].append('Inter_convex_hull')
          inter_boundaries_df['geometry'].append(modules_polygon_intermediate.convex_hull)
          inter_boundaries_df['color'].append('green')
          # Concave hull of boundaries
          int_cent = np.array(bounding_polygon_intermediate.centroid.xy).reshape((1,2))
          boundary_xx = flatten_list(boundary_xx)
          boundary_yy = flatten_list(boundary_yy)
          pol_points = sort_points_for_polygon_format(boundary_xx,boundary_yy, int_cent)
          concave_hull = Polygon(pol_points)
          inter_boundaries_df['name'].append('Inter_concave_hull')
          inter_boundaries_df['geometry'].append(concave_hull)
          inter_boundaries_df['color'].append('orange')
          coverage_polygon_intermediate = MultiPolygon(inter_coverage) # save coverage as whole to speedup global calculation
          # plot_polygon(coverage_polygon_intermediate, add_points=False)
          covered_area_inter = coverage_polygon_intermediate.area
          intermediate_collection.append(bounding_polygon_intermediate)
          intermediate_collection = GeometryCollection(intermediate_collection)
          intermediate_collection_speed = GeometryCollection([bounding_polygon_intermediate, modules_polygon_intermediate, coverage_polygon_intermediate, concave_hull, inter_robots_df['geometry']])
          available_intermediate_area = round(bounding_polygon_intermediate.area,1)
          intermediate_coverage = round(covered_area_inter/available_intermediate_area*100,1)
          
          inter_df = {'inter_boundaries': inter_boundaries_df, 'inter_modules': inter_modules_df, 'inter_coverage': inter_coverage_df, 'inter_robots': inter_robots_df, 'intermediate_coverage': intermediate_coverage}

          return intermediate_collection, intermediate_collection_speed, intermediate_coverage, inter_df

# NOTE: GFA class
class GFA(SavingResults):
     def __init__(self, length: float, width: float, nb_gfa: int, saving_df, vigR: float, trimming_angle = 360, trimming_geometry = None) -> None:
          self.length = length
          self.width = width
          self.vigR = vigR
          self.nb_gfa = nb_gfa
          self.trimming_angle = trimming_angle
          self.trimming_geometry = trimming_geometry
          self.gdf_gfa = self.make_GFA_array()

          super().__init__(saving_df)

     def make_GFA(self):
          # Dummy shape, to be replaced with true footprint
          minx = -self.length/2
          miny = -self.width/2
          maxx = self.length/2
          maxy = self.width/2
          gfa_footprint = Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy), (minx, miny)])
          return gfa_footprint
     
     def make_GFA_array(self):

          gfa_df = {'gfa_index':[], 'center': [], 'orientation': [], 'geometry':[], 'color':[], 'label':[]}
          # gfa_pos_on_vigR = make_vigR_polygon(n_vigR = self.nb_gfa + 1).exterior.coords.xy
          Dangle = 360/self.nb_gfa
          angles = np.linspace(Dangle,Dangle * self.nb_gfa, self.nb_gfa)
          gfa_pos_on_vigR_x = self.vigR * np.cos(np.deg2rad(angles))
          gfa_pos_on_vigR_y = self.vigR * np.sin(np.deg2rad(angles))

          for i in range(self.nb_gfa):
               x = gfa_pos_on_vigR_x[i]
               y = gfa_pos_on_vigR_y[i]
               theta = angles[i]
               gfa = self.make_GFA()
               placed_gfa = rotate_and_translate(gfa, angles[i], x, y)

               if self.trimming_angle != 360 and not placed_gfa.intersects(self.trimming_geometry):
                    continue

               gfa_df['gfa_index'].append(i)
               gfa_df['center'].append([x,y])
               gfa_df['orientation'].append(angles[i])
               gfa_df['geometry'].append(placed_gfa)
               gfa_df['label'].append('GFA')
               gfa_df['color'].append('brown')
          
          gdf_gfa = gpd.GeoDataFrame(gfa_df)

          return gdf_gfa
     
     def closest_gfa(self, pos):
          """
          Finds the closest gfa to a given point
          
          Input: pos = [x,y] # [mm, mm]
          Ouput: index of closest gfa
          """

          faraway = np.ones(self.nb_gfa)
          for idx, gfa_center in enumerate(self.gdf_gfa['center']):
               faraway[idx] = norm2d(pos, gfa_center)
          closest_gfa = np.where(faraway == np.amin(faraway))

          return closest_gfa[0][0]
     
     def GFA_to_csv(self):

          if not self.save_csv:
               return

          now = datetime.now()
          today_filename = now.strftime("%Y-%m-%d-%H-%M-%S_") + "GFA.csv"
          self.gdf_gfa.to_csv(self.results_dir_path() + today_filename)

     # NOTE: Grid generation class
class Grid(FocalSurf):

     def __init__(self, project_surface: str, module_width: float, inter_frame_thick: float, global_frame_thick: float, trimming_angle: float = 360, centered_on_triangle: bool  = False) -> None:
          
          
          self.centered_on_triangle = centered_on_triangle
          self.module_width = module_width
          self.inter_frame_thick = inter_frame_thick
          self.global_frame_thick = global_frame_thick
          self.inter_frame_width = 2*self.module_width + 2*self.inter_frame_thick*np.cos(np.deg2rad(30)) + 2*self.global_frame_thick*np.cos(np.deg2rad(30))
          
          super().__init__(project_surface)
          self.pizza = FocalSurf.make_vigR_polygon(self, trimming_angle = trimming_angle)

          self.grid_flat_init = {'geometry': []}
          self.grid_BFS = {'front':{}, 'back':{}}
          self.fiducials = {'x': [], 'y': [], 'z': [], 'geometry': []}
          self.x_grid, self.y_grid, self.z_grid = self.create_flat_grid()
          
          
     def create_flat_grid(self):
                    # Make global grid out of triangular grid method credited in updown_tri.py
          # Its logic is also explained in updown_tri.py
          n = 7
          center_coords = []
          a_max = n
          b_max = n
          c_max = n

          x_grid = []
          y_grid = []
          flip_global = [] # stores the orientation of each triangle (either up or downward)

          if self.centered_on_triangle:
               origin = 'triangle'
               valid = [0,1]
          else:
               origin = 'vertex'
               valid = [1,2]
          
          vigR_tresh = 150

          for a in np.arange(-a_max,a_max):
               for b in np.arange(-b_max,b_max):
                    for c in np.arange(-c_max,c_max):
                         sum_abc = a + b + c
                         # if valid == 1 or valid == 2: 
                         if valid.count(sum_abc): # check if sum abc corresponds to a valid triangle depending on the centering case
                              x,y = tri.tri_center(a,b,c,self.inter_frame_width)
                              #TODO: make smaller grid of fiducials to fill up more spaces
                              fiducial = tri.tri_corners(a,b,c, self.inter_frame_width)
                              # Intermediate triangles placement
                              if Point(x,y).within(self.pizza.buffer(vigR_tresh)): # allow centroid of inter modules to go out of vigR for further filling purpose
                                   center_coords.append((a,b,c))
                                   x_grid.append(x)
                                   y_grid.append(y)
                                   self.grid_flat_init['geometry'].append(Point(x,y))
                                   if tri.points_up(a,b,c, origin = origin): # flip the modules depending on there position on the grid
                                        flip_global.append(0)
                                   else: 
                                        flip_global.append(1)
                              # Fiducials placement
                              for fid in fiducial:
                                   if Point(fid[0],fid[1]).within(self.pizza.buffer(10)): # limit fiducials outside vigR to inter modules with center inside vigR
                                   
                                        self.fiducials['x'].append(fid[0])
                                        self.fiducials['y'].append(fid[1])
                                        self.fiducials['z'].append(0)
                                        self.fiducials['geometry'].append(Point(fid))
                                   # self.fiducials['xyz'] = np.vstack((self.fiducials['xyz'], np.asarray(fiducial)))
          
          x_grid = np.array(x_grid)
          y_grid = np.array(y_grid)
          z_grid = -np.sqrt(self.BFS**2 - (self.vigR)**2)*np.ones(len(x_grid))
          # z_grid = 0*np.ones(len(x_grid))

          #TODO: make use of the built-in strucutre there instead of returning just arrays
          # might be more clean and easier to manipulate

          self.grid_flat_init['x'] = x_grid
          self.grid_flat_init['y'] = y_grid
          self.grid_flat_init['z'] = z_grid
          self.grid_flat_init['tri_orientation'] = np.array(flip_global)
          # self.grid_flat_init['geometry'] = MultiPoint(to_polygon_format(x_grid, y_grid, z_grid))

          return x_grid, y_grid, z_grid
     
     def trim_grid(self, final_grid: dict , trimming_angle = 360):

          """
          Input:

          - final_grid: [dict] contains final grids of both modules and intermediate triangles center points
               - final_grid['indiv']['xyz']: (N,3) [numpy array] contains x,y,z coordinates of each module center point
               - final_grid['indiv']['upward_tri']: (N,1) [numpy array] contains 0 or 1 depending on the orientation of the module (0 = up, 1 = down)
          - trim_angle: [float] angle of the pizza slice to trim the grid with (default = 360 --> no trimming, full grid taken)


          Output:

          - grid_points: (N,3) [numpy array] contains x,y,z coordinates of each module center point within the trimming angle
          - module_up: (N,1) [numpy array] contains 0 or 1 depending on the orientation of the module (0 = up, 1 = down) within the trimming angle
          """

          if trimming_angle == 360:
               grid_points = np.asarray(final_grid['indiv']['xyz'])
               module_up = np.asarray(final_grid['indiv']['upward_tri'])
          else:
               extra_vigR = 20 # [mm] extra radius to include modules whose centroid are outside vigR but fulfills the out_allowance condition
               pizza_slice = FocalSurf.make_vigR_polygon(self, trimming_angle = trimming_angle, changed_vigR = self.vigR + extra_vigR)
               grid_points = []
               module_up = []
               for idx, point in enumerate(final_grid['indiv']['geometry']):
                    if point.within(pizza_slice):
                         grid_points.append(final_grid['indiv']['xyz'][idx])
                         module_up.append(final_grid['indiv']['upward_tri'][idx])
               grid_points = np.asarray(grid_points)
               module_up = np.asarray(module_up)

          return grid_points, module_up

     def project_grid_on_sphere(self, grid_points: np.array, sphere_radius: float, module_length: float, module_up: np.array):

          """
               Input:

               - grid_points: (N,3) [numpy array]
               - sphere_radius: [float] radius of the sphere onto which the grid is projected
               - module_length: [float] length of the module (used to project the grid on the back of the sphere)
               - module_up: (N,1) [numpy array] contains 0 or 1 depending on the orientation of the module (0 = up, 1 = down)

               Doing:
               1) Takes flat grid and projectes it on a sphere of given radius (use case: focal surface BFS radius)
               2) Saves result in "projection" dictionnary in 3 different formats developped below for further manipulations

               Output:

               - front_proj: (N,4) [numpy array] projected grid on the front of the sphere
               - back_proj: (N,4) [numpy array] projected grid on the back of the sphere (where the back of the modules is)
               - proj: (2*N,4) [numpy array] contains both front and back projections
          """
          projection = {}
          ## Normalize 3D flat grid so that every point lie on the unit sphere
          norm_points = norm(grid_points, axis=1)
          normalized = grid_points/norm_points[:, np.newaxis]
          ## Scale the unit sphere to the desired sphere with radius R
          front_proj = normalized * sphere_radius
          front_proj[:,2] = front_proj[:,2] + sphere_radius # brings back z coordinates to centered around 0 instead of BFS
          back_proj = normalized * (sphere_radius + module_length)
          back_proj[:,2] = back_proj[:,2] + sphere_radius # brings back z coordinates to centered around 0 instead of BFS

          projection['front'] = np.hstack((front_proj, module_up.reshape(len(module_up),1)))
          projection['back'] = np.hstack((back_proj, module_up.reshape(len(module_up),1)))
          projection['proj'] = np.vstack((projection['front'], projection['back']))

          for key in ['front', 'back']:
               self.grid_BFS[key]['x'] = projection[key][:,0]
               self.grid_BFS[key]['y'] = projection[key][:,1]
               self.grid_BFS[key]['z'] = projection[key][:,2]

          self.grid_BFS['geometry'] = MultiPoint(to_polygon_format(projection['front'][:,0], projection['front'][:,1], projection['front'][:,2]))

          return projection
     
     def orientation_vector(self, phi: float, theta: float):

          """
          Input:

          - phi: [float] azimuthal angle in spherical coordinates
          - theta: [float] polar angle in spherical coordinates

          Output:

          - orientation_vector: [3x1 numpy array] contains the x,y,z coordinates of the orientation vector
          """

          x = - np.sin(theta) * np.cos(phi)
          y = - np.sin(theta) * np.sin(phi)
          z = np.cos(theta)

          return np.array([x,y,z]).T

     def plot_3D_grid(self,ax, x, y ,z,  color = 'red', label=None):

          ax.scatter(x, y, z , label=label, color=color)
          ax.set_box_aspect((5,5,1))
          plt.legend()
          ax.set_xlabel('X')
          ax.set_ylabel('Y')
          ax.set_zlabel('Z')

     def plot_2D_grid(self):

          fig = plt.figure('2D grid', figsize=(8,8))
          ax = fig.add_subplot()
          ax.scatter(self.grid_flat_init['x'], self.grid_flat_init['y'], label=f'Projected', color='red')
          # ax.scatter(self.grid_df['x_grid_flat'], self.grid_df['y_grid_flat'], label=f'Flat', color='blue')
          ax.set_box_aspect(1)
          plt.legend()
          ax.set_xlabel('X')
          ax.set_ylabel('Y')

def to_polygon_format(x,y, z = None):
     """ Input:
          - x,y: 2 sets of coordinates for polygon creation

          Output:
          - coords: list of tuples for each polygon vertices (just for convenience) """
     coords = []
     if z is not None:
          for (i,j,k) in zip(x,y,z):
               coords.append((i,j,k))
     else:
          for (i,j) in zip(x,y):
               coords.append((i,j))
               
     return coords

def rotate_and_translate(geom, angle, dx, dy, dz = None, origin = 'centroid'):

     rotated_geom = affinity.rotate(geom, angle, origin=origin)
     transformed_geom = affinity.translate(rotated_geom, dx, dy, dz)

     return transformed_geom

def rotation3D(vector, phi, theta, rotation_axis='z'):
     """
     Rotate a 3D vector around a given axis by an angle theta.

     Parameters:
     - vector: numpy array of shape (3,) representing the 3D vector to be rotated
     - phi: azimuthal angle in spherical coordinates
     - theta: polar angle in spherical coordinates
     - rotation_axis: axis of rotation, can be 'x', 'y', or 'z' (default is 'z')

     Returns:
     - rotated_vector: numpy array of shape (3,) representing the rotated 3D vector
     """

     # Define the rotation matrix based on the rotation axis
     if rotation_axis == 'x':
          rotation_matrix = np.array([[1, 0, 0],
                                             [0, np.cos(theta), -np.sin(theta)],
                                             [0, np.sin(theta), np.cos(theta)]])
     elif rotation_axis == 'y':
          rotation_matrix = np.array([[np.cos(theta), 0, np.sin(theta)],
                                             [0, 1, 0],
                                             [-np.sin(theta), 0, np.cos(theta)]])
     elif rotation_axis == 'z':
          rotation_matrix = np.array([[np.cos(theta), -np.sin(theta), 0],
                                             [np.sin(theta), np.cos(theta), 0],
                                             [0, 0, 1]])
     else:
          raise ValueError("Invalid rotation axis. Please choose 'x', 'y', or 'z'.")

     # Perform the rotation
     rotated_vector = np.dot(rotation_matrix, vector)

     return rotated_vector



def norm2d(p,q):
     return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)


def plot_module(module_collection, label_coverage, label_robots, ignore_points):
     for jdx, geo in enumerate(module_collection.geoms):
          # plot_polygon(geometry[jdx], add_points=False, facecolor='None' , edgecolor='black')
          if (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 0:
               plot_polygon(module_collection.geoms[jdx], add_points=False, facecolor='None' , edgecolor='black')
          elif (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 1:
               plot_polygon(module_collection.geoms[jdx], add_points=False, facecolor='None', linestyle = '--')
          elif (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 2:
               plot_polygon(module_collection.geoms[jdx], add_points=False, alpha=0.2, edgecolor='black', label = label_coverage)
          elif (isinstance (module_collection.geoms[jdx], MultiPoint) and not ignore_points):
               plot_points(module_collection.geoms[jdx], marker='.', color='k', label = label_robots)

def plot_intermediate(intermediate_collection, nb_robots, ignore_points, intermediate_coverage=None, available_intermediate_area=None, draw_legend = False):
     for idx, mod_collection in enumerate(intermediate_collection.geoms):
          if idx == 0 and draw_legend:
               label_coverage = 'Coverage: {} %'.format(intermediate_coverage)
               # label_coverage = 'Coverage: {} %'.format(area_to_cover)
               label_robots = "{} robots".format(nb_robots)
          else: 
               label_coverage = None
               label_robots = None
          if (isinstance (mod_collection, Polygon)):
               plot_polygon(mod_collection, add_points=False, facecolor='None', linestyle = '-.', color = 'green')

               # plot_polygon(mod_collection, add_points=False, facecolor='None', linestyle = '-.', color = 'green', label = 'Available area: {} mm$^2$'.format(available_intermediate_area))
               continue
          plot_module(mod_collection, label_coverage, label_robots, ignore_points)
               
               
def plot_module(module_collection, label_coverage, label_robots, ignore_points):
     for jdx, geo in enumerate(module_collection.geoms):
          # plot_polygon(geometry[jdx], add_points=False, facecolor='None' , edgecolor='black')
          if (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 0:
               plot_polygon(module_collection.geoms[jdx], add_points=False, facecolor='None' , edgecolor='black')
          elif (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 1:
               plot_polygon(module_collection.geoms[jdx], add_points=False, facecolor='None', linestyle = '--')
          elif (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 2:
               plot_polygon(module_collection.geoms[jdx], add_points=False, alpha=0.2, edgecolor='black', label = label_coverage)
          elif (isinstance (module_collection.geoms[jdx], MultiPoint) and not ignore_points):
               plot_points(module_collection.geoms[jdx], marker='.', color='k', label = label_robots)

def plot_intermediate(intermediate_collection, nb_robots, ignore_points, intermediate_coverage=None, available_intermediate_area=None, draw_legend = False):
     for idx, mod_collection in enumerate(intermediate_collection.geoms):
          if idx == 0 and draw_legend:
               label_coverage = 'Coverage: {} %'.format(intermediate_coverage)
               # label_coverage = 'Coverage: {} %'.format(area_to_cover)
               label_robots = "{} robots".format(nb_robots)
          else: 
               label_coverage = None
               label_robots = None
          if (isinstance (mod_collection, Polygon)):
               plot_polygon(mod_collection, add_points=False, facecolor='None', linestyle = '-.', color = 'green')

               # plot_polygon(mod_collection, add_points=False, facecolor='None', linestyle = '-.', color = 'green', label = 'Available area: {} mm$^2$'.format(available_intermediate_area))
               continue
          plot_module(mod_collection, label_coverage, label_robots, ignore_points)
               
def final_title(project: str, vigR:float, nb_robots: int, total_modules: int, total_robots: int, inter_frame_thick, global_frame_thick, allow_small_out: bool, out_allowance: float, disp_robots_info: bool = True):

     project_info = r"$\bf{Project:}$" + project + f" - vigR: {vigR:0.1f} mm"
     
     if allow_small_out:
          small_out_info = f"Out allowance: {out_allowance * 100:0.1f} %"
     else:
          small_out_info = ''

     if disp_robots_info:
          robots_info = f"\n Total # modules: {total_modules} - Total # robots: {total_robots}"
          modules_info = f" {nb_robots} robots per module"
     else:
          robots_info = ''
          modules_info = ''
     if inter_frame_thick != global_frame_thick:
          figtitle = f"{project_info}\n Semi frameless - {modules_info} \n Inner gap: {inter_frame_thick} mm - Global gap: {global_frame_thick} mm {robots_info} \n {small_out_info} \n"
     elif inter_frame_thick == global_frame_thick and global_frame_thick == 0:
          figtitle = f"{project_info}\n Frameless - {modules_info} {robots_info} \n {small_out_info} \n"
     else:
          figtitle = f"{project_info}\n Framed - {modules_info} \n Gap: {inter_frame_thick} mm {robots_info} \n {small_out_info} \n"

     return figtitle

def flatten_list(l:list):
     return [item for sublist in l for item in sublist]

def makeKey(key_int: int):
     return f'n{key_int}'

def sort_points_for_polygon_format(x: list, y: list, centroid):
     """
     Sort points (np array) counterclockwise for later polygon creation
     """
     points = []
     for (xi,yi) in zip(x, y):
          points.append([xi,yi])
     points = np.array(points)
     vec = points - centroid # get vector connecting centroid of points to them
     angles = np.arctan2(vec[:,1],vec[:,0]) # calculate angle of each vector
     d = np.hstack((vec, angles.reshape(len(angles),1))) # put everything in one array
     d = d[d[:, 2].argsort()] # sort the points by ascending order array
     return to_polygon_format(d[:,0] + centroid[0,0], d[:,1]+ centroid[0,1])