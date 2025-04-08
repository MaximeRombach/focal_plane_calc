import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-4s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
import pandas as pd
import numpy as np
from shapely.geometry import Polygon


class FocalSurf:
    """
    Defines the focal surface properties of the telescope
    The convention is the following:
    - coordinate system origin in the center of the focal surface
    - z axis points upward
    - r: cartesian norm in 3D
    - theta: declination between r vector and polar axis (here z). Starts at 0 when parallel to z pointing upwards and goes up to 180°
    - phi: (azimuthal angle) angular position about z axis from 0 to 360°

    """

    def __init__(self, project_name: str, **kwargs): 

        self.project_name = project_name
        self.curvature_R = kwargs.get('curvature_R', None) # [mm] radius of curvature of the focal surface
        self.vigD = kwargs.get('vigD', 0) # [mm] diameter of the focal surface
        self.vigR = self.vigD / 2 # [mm] vignetting radius, nominal case declared for each project
        self.asph_formula = kwargs.get('asph_formula', False) # True if the focal surface is aspheric, False if spherical
        self.k = kwargs.get('k', None) # aspheric coefficient
        self.a2 = kwargs.get('a2', None) # [mm] spherical aberration coefficient
        self.a4 = kwargs.get('a4', None) # [mm] spherical aberration coefficient
        self.a6 = kwargs.get('a6', None) # [mm] spherical aberration coefficient
        self.a8 = kwargs.get('a8', None) # [mm] spherical aberration coefficient
        self.fnumber = kwargs.get('fnumber', None) # average f-number of the telescope
        self.BFS = kwargs.get('BFS', None) # [mm] Best Fit Sphere radius
        self.FoV = kwargs.get('FoV', None) # [deg] Field of View of the telescope
        self.focus_tolerance_width = kwargs.get('focus_tolerance_width', None) # [mm] tolerance on the focal surface
        self.plate_scale = kwargs.get('focus_tolerance_width', None) # [um/arcsec] corresponding distance on the focal surface mapped to angular position on sky
        self.donut_diam = kwargs.get('donut_diam', None) # [um/arcsec] corresponding distance on the focal surface mapped to angular position on sky


    
    @property    
    def optics_data(self):

        filename = f"./Data_focal_planes/{self.project_name}.csv" # optics data from Zemax
        comment_character = "#"  # The character that indicates a commented line
        # Read CSV file and ignore commented lines
        if self.project_name == 'MUST':
            sep = ';' # keeping distinction as delimiter might change from time to time
        elif self.project_name == 'MegaMapper':
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

        logging.info(f"{self.project_name} focal plane data successfully read from {filename}")

        """ Check if Chief Ray Deviation data is available (optional)"""
        if 'CRD' in optics_data.columns:
            self.is_CRD = True
            logging.info(f"CRD data found")
        else:
            self.is_CRD = False
            logging.warning(f"CRD data NOT found")    

        return optics_data
    
    @property
    def vignetting_disk(self, trimming_angle: int = 360, changed_vigR: float = None):
        """
        Returns the vignetting disk as a shapely Polygon object
        
        Input: 
        - trimming_angle [deg]: angle to trim the vignetting disk at
        - changed_vigR [mm]: allows to change the vignetting radius for specific cases
        
        Output:
        - pizza: shapely Polygon object representing the vignetting disk
        """
        if changed_vigR is not None:
            vigR = changed_vigR # [mm] allows to change the vignetting radius for specific cases
        else:
            vigR = self.vigR # [mm] vignetting radius, nominal case declared for each project

        vigR_lim_x = vigR * np.cos(np.deg2rad(np.linspace(0,trimming_angle,500)))
        vigR_lim_y = vigR * np.sin(np.deg2rad(np.linspace(0,trimming_angle,500)))

        pizza = Polygon([[x,y] for x,y in zip(vigR_lim_x, vigR_lim_y)])

        if 'WST' in self.project_name:

            # WST project has a hole in the center of MOS for the IFU mode
            donutR = self.donut_diam / 2 # [mm] donut radius
            donut_lim_x = donutR * np.cos(np.deg2rad(np.linspace(0,360,500)))
            donut_lim_y = donutR * np.sin(np.deg2rad(np.linspace(0,360,500)))
            donut = Polygon([[x,y] for x,y in zip(donut_lim_x, donut_lim_y)])
            
            pizza = pizza.difference(donut) # remove the donut hole from the vignetting disk

        return pizza
    
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