import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-4s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
import pandas as pd
import numpy as np
from shapely.geometry import Polygon
import json
from scipy.interpolate import interp1d
from shapely.plotting import plot_polygon
import matplotlib.pyplot as plt


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
    def optics_data(self, alt_file_name: str = None):

        if alt_file_name is not None:
            fileName = alt_file_name
        else:
            fileName = self.project_name
        filePath = f"./Data_focal_planes/{fileName}.csv" # optics data from Zemax
        comment_character = "#"  # The character that indicates a commented line
        # Read CSV file and ignore commented lines
        if self.project_name == 'MUST':
            sep = ';' # keeping distinction as delimiter might change from time to time
        elif self.project_name == 'MegaMapper':
            sep = ','
        else :
            sep = ','
        optics_data = pd.read_csv(filePath, comment=comment_character, sep=sep)

        # Set headers using the first non-comment line
        with open(filePath, 'r') as file:
            for line in file:
                    if not line.startswith(comment_character):
                        headers = line.strip().split(sep)
                        break
        optics_data.columns = headers

        logging.info(f"{self.project_name} focal plane data successfully read from {filePath}")

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
            pizza = pizza.difference(self.donut_hole)  # remove the central region for IFU mode in case of WST layouts
        
        return pizza
    
    @property
    def donut_hole(self):
            
        # WST project has a hole in the center of MOS for the IFU mode
        if  'WST' in self.project_name:
            donutR = self.arcmin2mm(self.donut_diam) / 2
            donut_lim_x = donutR * np.cos(np.deg2rad(np.linspace(0,360,500)))
            donut_lim_y = donutR * np.sin(np.deg2rad(np.linspace(0,360,500)))
            donut = Polygon([[x,y] for x,y in zip(donut_lim_x, donut_lim_y)])

        else:
            donut = Polygon()

        return donut
    
    def transfer_functions(self):
        # Those transfer functions and there implementation later are inspired from the work of Joseph Silber (LBNL) in the generate_raft_layout.py code in https://github.com/joesilber/raft-design 

        """ 
        Input: 
            - optics_data: pandas dataframe containing the focal plane data from zmax csv file
            - analytical: [bool] flag to use analytical solution of R2Z instead of data 

        Output:
            - R2Z: functions to convert radius to height on focal surface
            - R2CRD: function to get chief ray deviation interms of radial position on focal surface
            - R2NORM: maps the radial position to the angle of the normal vector of the focal surface
            - R2NUT: maps the radial position to the actual nutation angle which is: R2NORM + CRD
        
        ATTENTION: Needs csv file to contain columns named: R, Z, CRD """
        r = np.linspace(0,self.vigR,2000) # [mm] radial positions on focal plane

        """ R2Z """
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
                R2Z = lambda r: self.BFS - np.sqrt(self.BFS**2 - r**2) # Spherical focal surface
                logging.warning('No Z data available in samples - assuming SPHERICAL surface')

        else:
            logging.info('Transfer functions: using aspherical formula cofficients, no CRD data available')
            R2Z = self.analytical_R2Z

        """ R2CRD """
        if self.is_CRD:
            print(self.optics_data.keys())
            CRD = self.optics_data['CRD']
            R = self.optics_data['R']
            R2CRD = interp1d(R,CRD,kind='cubic', fill_value = "extrapolate")
        else:
            CRD = np.zeros_like(R)
            R2CRD = interp1d(R,CRD,kind='cubic', fill_value = "extrapolate")
        crd = R2CRD(r)

        """ Calculate derivative of focal surface """
        z = R2Z(r)
        dr = np.diff(r)
        dz = np.diff(z)
        dzdr = dz/dr

        """ R2NORM """
        norm = np.degrees(np.arctan(dzdr))# Calculate normal angles at the sample points
        R2NORM = interp1d(r[:-1], norm,kind='cubic', fill_value = "extrapolate")

        """ R2NORM """
        nut = R2NORM(r) + crd # Joe: -(R2NORM(r) + crd)
        R2NUT = interp1d(r, nut, kind='cubic', fill_value = "extrapolate")

        # Parametric path length along the focal surface
        ds = (1 + dzdr**2)**0.5 * dr
        s = np.cumsum(ds)
        s = np.insert(s, 0, 0.)  # first value of path length is 0
        S2R = interp1d(s, r, kind='cubic', fill_value = "extrapolate")

        return R2Z, R2CRD, R2NORM, R2NUT, S2R
    
    def trimming_polygon(self, geometry = 'hex', trim_diff_to_vigR: float = 0):
        """
        Returns the trimming polygon as a shapely Polygon object

        Input:
        - geometry: 'hex' or 'square'
        - trim_diff_to_vigR [mm]: distance between the trimming polygon and the vignetting disk

        Output:
        - trim_polygon: shapely Polygon object representing the trimming polygon
        """

        if geometry == 'square':
            points = 5 #one more point to close the polygon
        elif geometry == 'hex':
            points = 7 #one more point to close the polygon
        elif geometry == 'circle':
            points = 20
        else:
            raise ValueError("geometry must be 'hex', 'square' or 'circle'")

        trimR = self.vigR - trim_diff_to_vigR
        trimR_lim_x = trimR * np.cos(np.deg2rad(np.linspace(0,360,points)))  # [mm] trimming radius
        trimR_lim_y = trimR * np.sin(np.deg2rad(np.linspace(0,360,points)))  # [mm] trimming radius
        trim_polygon = Polygon([[x,y] for x,y in zip(trimR_lim_x, trimR_lim_y)])
        
        if 'WST' in self.project_name:
            # WST project has a hole in the center of MOS for the IFU mode
            donutR = self.arcmin2mm(self.donut_diam) / 2
        
        return trim_polygon
    
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
    
    def mm2arcsec(self, mm):
        # Converts 
        return mm * (self.FoV * 60 * 60) / self.vigD
    
    def mm2deg(self, mm):
        return mm * self.FoV / self.vigD
    
    def mm2arcmin(self, mm):
        return mm * (self.FoV * 60) / self.vigD
    
    def arcmin2mm(self, arcmin):
        return arcmin * self.vigD / (self.FoV * 60)

if __name__ == "__main__":
    project = 'MUST'  # Example project name, change as needed
    project_parameters = json.load(open('projects.json', 'r'))
    surf = FocalSurf(project, **project_parameters[project])
    R2Z, R2CRD, R2NORM, R2NUT, S2R = surf.transfer_functions()

    print(optical_data := surf.optics_data)  # Should print the optics data as a pandas DataFrame

    plot_polygon(surf.vignetting_disk) # Should print the vignetting disk as a shapely Polygon object
    plot_polygon(surf.trimming_polygon(geometry='hex')) # Should print the trimming polygon as a shapely Polygon object
    
    figure, ax = plt.subplots()
    r = np.linspace(0,surf.vigR,500)
    z = R2Z(r)
    plt.plot(r, z)
    plt.grid()
    
    plt.show()

    