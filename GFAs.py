from shapely.geometry import Polygon, MultiPolygon
from shapely import affinity
import numpy as np
import datetime
import math
import geopandas as gpd
from SavingResults import SavingResults


class GFA():
    def __init__(self, length: float, width: float, nb_gfa: int, vigR: float, trimming_angle = 360, trimming_geometry = None, angle_offset = None) -> MultiPolygon:
        self.length = length
        self.width = width
        self.vigR = vigR
        self.nb_gfa = nb_gfa
        self.trimming_angle = trimming_angle
        self.trimming_geometry = trimming_geometry
        self.angle_offset = angle_offset
        self.gdf_gfa = self.make_GFA_array()

        # super().__init__(saving_df)

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
        angles = np.linspace(Dangle,Dangle * self.nb_gfa, self.nb_gfa) + self.angle_offset
        gfa_pos_on_vigR_x = self.vigR * np.cos(np.deg2rad(angles))
        gfa_pos_on_vigR_y = self.vigR * np.sin(np.deg2rad(angles))

        for i in range(self.nb_gfa):
            x = gfa_pos_on_vigR_x[i]
            y = gfa_pos_on_vigR_y[i]
            theta = angles[i]
            gfa = self.make_GFA()
            placed_gfa = self.rotate_and_translate(gfa, angles[i], x, y)

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
            faraway[idx] = self.norm2d(pos, gfa_center)
        closest_gfa = np.where(faraway == np.amin(faraway))

        return closest_gfa[0][0]
     
    # def GFA_to_csv(self):

    #     if not self.save_csv:
    #         return

    #     now = datetime.now()
    #     today_filename = now.strftime("%Y-%m-%d-%H-%M-%S_") + "GFA.csv"
    #     self.gdf_gfa.to_csv(self.results_dir_path() + today_filename)

    def rotate_and_translate(self, geom, angle, dx, dy, dz = None, origin = 'centroid'):

            rotated_geom = affinity.rotate(geom, angle, origin=origin)
            transformed_geom = affinity.translate(rotated_geom, dx, dy, dz)

            return transformed_geom
    
    def norm2d(p,q):
     return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)