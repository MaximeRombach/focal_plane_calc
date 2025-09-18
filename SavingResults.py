from datetime import datetime
import os
import logging
import matplotlib.pyplot as plt
import pandas as pd
import ezdxf
import ezdxf.addons.geo
from ezdxf import units
from shapely.geometry import mapping



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
     
     def save_dxf_to_dir(self, geometries: list, suffix_name):

          if not self.save_dxf:
               return
          
          now = datetime.now()
          name_frame = now.strftime("%Y-%m-%d--%H-%M-%S_") + suffix_name

          doc = ezdxf.new()
          doc.units = units.MM
          msp = doc.modelspace()
          for index, geom in enumerate(geometries):

               geoproxy = ezdxf.addons.geo.GeoProxy.parse(mapping(geom))

               # Use LWPOLYLINE instead of hatch.
               for entity in geoproxy.to_dxf_entities(polygon=2):
                    msp.add_entity(entity)
                    entity.set_dxf_attrib('layer', f"{index}")

          doc.saveas(self.results_dir_path + f"{name_frame}.dxf")
          logging.info(f'{name_frame}.dxf successfully saved in {self.results_dir_path}')

     def save_figures_to_dir(self, image_name: str, file_format: str = 'png', save_eps: bool = False, dpi: int = 400):

          if not self.save_plots:
               return

          now = datetime.now()
          today_filename = now.strftime("%Y-%m-%d-%H-%M-%S_") + image_name + f'.{file_format}'
          plt.savefig(self.results_dir_path + today_filename, bbox_inches = 'tight', format=file_format, dpi = dpi)
          if save_eps:
               plt.savefig(self.results_dir_path + today_filename + ".eps", bbox_inches = 'tight', format='eps')
               logging.info(f'{image_name}.eps successfully saved in in {self.results_dir_path}')

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
          #TODO: more optimal than the previous function, needs to replace all calls to the previous function
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

     def save_grid_to_csv(self, grid: pd, filename: str, results_string: str = None):
          if not self.save_txt:
               return
          now = datetime.now()
          now_str = now.strftime("%Y-%m-%d-%H-%M_")
          file_path = self.results_dir_path + now_str + self.project_name + '_' + filename + '.csv'
          if results_string is not None:
               with open(file_path, 'w') as file:
                    file.write(results_string)
          grid.to_csv(file_path, sep = ",", decimal = ".", index = False)
          
          logging.info(f'{filename}.csv succesfully saved in {self.results_dir_path}')