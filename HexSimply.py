#############################################################################
# Filename: HexSimply.py                                                    #
# Author: Karolina Mamczarz                                                 #
# Institution: AGH University of Science and Technology in Cracow, Poland   #
# Last update: 2017-07-24                                                   #
# Version: 1.0.0                                                            #
# Description: #                                                            #
# Class: HexSimply                                                          #
# Methods: #                                                                #
# Parameters: #                                                             #
# Result: #                                                                 #
#############################################################################

import arcpy
import os
from math import sqrt, trunc, sin, cos, pi, atan2
from sys import exit


class HexSimply(object):

    def __init__(self, original, method, l, s, simplified):
        self.original = original
        self.simplified = simplified
        self.l = l / 1000
        self.s = s
        self.method = method
        self.workspace = os.path.dirname(self.original)  # POTEM DO WYRZUCENIA
        self.choose_method()

    def sort(self):
        #ewentualnie skrypt na dissolve NIEEEEEEEEEEE
        return

    def read_geom(self, shape):
        with arcpy.da.SearchCursor(shape, ["SHAPE@"]) as cursor:
            self.shape_coords = []
            partnum = 0
            for row in cursor:
                for part in row[0]:
                    pntnum = 0
                    for pnt in part:
                        self.shape_coords.append([partnum, pntnum, pnt.X, pnt.Y])
                        pntnum += 1
                    partnum += 1
        return self.shape_coords

    @staticmethod
    def read_geom_attr_rectangle(rectangle):
        with arcpy.da.SearchCursor(rectangle, ["SHAPE@"]) as cursor:
            rect_area_poly_coords = []
            for row in cursor:
                partnum = 0
                for part in row[0]:
                    pntnum = 0
                    for pnt in part:
                        rect_area_poly_coords.append([partnum, pntnum, pnt.X, pnt.Y])
                        pntnum += 1
                    partnum += 1
        return rect_area_poly_coords

    def tessera_width(self):
        self.r = 5*self.l*self.s
        return self.r

    def largest_diagonal_half(self):
        self.d = round(self.r/sqrt(3), 4)
        return self.d

    """
    The algorithm is known as "Ray Casting Method"
    Implemented to the code based on the code from website:
    http://geospatialpython.com/2011/01/point-in-polygon.html
    """
    @staticmethod
    def ray_casting_method(hex_coords, iter_polyline_coords):
        n = len(hex_coords)
        inside = False
        p1x, p1y = hex_coords[0]
        for k in range(n+1):
            p2x, p2y = hex_coords[k % n]
            if iter_polyline_coords[3] > min(p1y, p2y):
                if iter_polyline_coords[3] <= max(p1y, p2y):
                    if iter_polyline_coords[2] <= max(p1x, p2x):
                        if p1y != p2y:
                            xints = (iter_polyline_coords[3]-p1y) * (p2x-p1x) / (p2y-p1y) + p1x
                        if p1x == p2x or iter_polyline_coords[2] <= xints:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    def vertex_clustering(self, points_list):
        points_list.sort(key=lambda id_point: id_point[1])
        self.cluster = []
        prev = 0
        for point in points_list:
            if point[0] - prev != 0:
                prev = point[0]
                self.cluster.append([point])
            elif point[0] == 0:
                idx = points_list.index(point)
                if idx == 0:
                    prev = point[0]
                    self.cluster.append([point])
                else:
                    self.cluster[-1].append(point)
            else:
                self.cluster[-1].append(point)
        return self.cluster

    def spatial_mean(self, cluster_list, points_list):
        self.mean_xy_list = []
        self.mean_xy_list.append(arcpy.Point(points_list[0][2], points_list[0][3]))
        for one_cluster in cluster_list:
            for xy in one_cluster:
                n = one_cluster.index(xy) + 1
            sums = [sum(i) for i in zip(*one_cluster)]
            mean_x = sums[2]/n
            mean_y = sums[3]/n
            self.mean_xy_list.append(arcpy.Point(mean_x, mean_y))
        self.mean_xy_list.append(arcpy.Point(points_list[-1][2], points_list[-1][3]))
        return self.mean_xy_list

    def eliminate_self_crossing(self):
        return

    def statistics(self):
        return

    def set_path(self):
        path, file = os.path.split(self.simplified)
        filename, file_ext = os.path.splitext(file)
        arcpy.CreateFeatureclass_management(path, filename, "POLYLINE", "", "", "", self.original)
        container, container_ext = os.path.splitext(path)
        if container_ext == ".gdb" or container_ext == ".mdb" or container_ext == ".sde":
            file = filename
        self.full_pathname = path + "\\" + file
        return self.full_pathname

    def create_new_feature(self, new_polyline_coords):
        with arcpy.da.InsertCursor(self.full_pathname, ["SHAPE@"]) as cursor:
            cursor.insertRow([arcpy.Polyline(arcpy.Array(new_polyline_coords))])
        return

    """
    Direction of tessellation consistent horizontally to minimal area bounding box starting from upper-left corner
    of this bounding box
    """
    def bounding_box(self):
        self.set_path()
        polyline_coords = self.read_geom(self.original)
        az = [30, 90, 150, 210, 270, 330]
        max_min = []
        for coord in zip(*polyline_coords):
            max_min.append(max(coord))
            max_min.append(min(coord))
        x_ur = max_min[4]  # x of upper-right corner
        y_ur = max_min[6]  # y of upper-right corner
        x_dl = max_min[5]  # x of down-left corner
        y_dl = max_min[7]  # y of down-left corner
        x_ul = x_dl        # x of upper-left corner
        y_ul = y_ur        # y of upper-left corner
        a = sqrt((x_ul-x_ur)**2 + (y_ul-y_ur)**2)
        b = sqrt((x_ul-x_dl)**2 + (y_ul-y_dl)**2)
        vertical_cover = trunc((b-(self.tessera_width()/2)) /self.tessera_width()) + 1
        horizontal_cover = trunc((a-0.5*self.largest_diagonal_half()) / (1.5*self.largest_diagonal_half())) + 2
        points_in_hex_coords = []
        id_hex = 0
        for i in range(horizontal_cover + 1):
            for j in range(vertical_cover + 1):
                hex_coords_temp = []
                for i_az in az:
                    if (i % 2 == 0):
                        x = round((x_ul + self.largest_diagonal_half()*sin((i_az*pi)/180)) + 1.5*self.largest_diagonal_half()*i, 4)
                        y = round((y_ul + self.largest_diagonal_half()*cos((i_az*pi)/180)) - self.tessera_width()*j, 4)
                        hex_coords_temp.append([x, y])
                    else:
                        x = round((x_ul + self.largest_diagonal_half()*sin((i_az*pi)/180)) + 1.5*self.largest_diagonal_half()*i, 4)
                        y = round((y_ul + self.largest_diagonal_half()*cos((i_az*pi)/180)) - self.tessera_width()*j - self.tessera_width()/2, 4)
                        hex_coords_temp.append([x, y])
                for point_coords in polyline_coords:
                    # Using Ray Casting Method
                    if self.ray_casting_method(hex_coords_temp, point_coords) == True:
                        points_in_hex_coords.append([id_hex, point_coords[1], point_coords[2], point_coords[3]])
                id_hex += 1
        # Using Vertex Clustering
        cluster = self.vertex_clustering(points_in_hex_coords)
        # Using Spatial Mean
        mean_xy = self.spatial_mean(cluster, polyline_coords)
        self.create_new_feature(mean_xy)
        return

    def minimal_rectangle_general(self, method):
        self.set_path()
        polyline_coords = self.read_geom(self.original)
        az = [30, 90, 150, 210, 270, 330]
        temp_rect_area = "in_memory\\rect_area"
        arcpy.MinimumBoundingGeometry_management(self.original, temp_rect_area, method, "ALL")
        data = self.read_geom(temp_rect_area)
        dx01 = data[1][2]-data[0][2]
        dy01 = data[1][3]-data[0][3]
        dx12 = data[2][2]-data[1][2]
        dy12 = data[2][3]-data[1][3]
        a = sqrt(dx01**2 + dy01**2)
        b = sqrt(dx12**2 + dy12**2)
        if a > b:
            if (dx01 > 0 and dy01 > 0) or (dx01 < 0 and dy01 > 0):
                orient = atan2(dy01, dx01)
            else:
                orient = atan2(dy01, dx01) + 2*pi
            vertical_cover = trunc((b-(self.tessera_width()/2)) /self.tessera_width()) + 1
            horizontal_cover = trunc((a-0.5*self.largest_diagonal_half()) / (1.5*self.largest_diagonal_half())) + 2
            points_in_hex_coords = []
            id_hex = 0
            for i in range(horizontal_cover + 1):
                for j in range(vertical_cover + 1):
                    hex_coords_temp = []
                    for i_az in az:
                        if (i % 2 == 0):
                            xpp = round((data[0][2] + self.largest_diagonal_half()*sin((i_az*pi)/180)) + 1.5*self.largest_diagonal_half()*i, 4)
                            ypp = round((data[0][3] + self.largest_diagonal_half()*cos((i_az*pi)/180)) - self.tessera_width()*j, 4)
                            x_prim = xpp - data[0][2]
                            y_prim = ypp - data[0][3]
                            x_bis = x_prim*cos(orient) - y_prim*sin(orient)
                            y_bis = x_prim*sin(orient) + y_prim*cos(orient)
                            x = x_bis + data[0][2]
                            y = y_bis + data[0][3]
                            hex_coords_temp.append([x, y])
                        else:
                            xpp = round((data[0][2] + self.largest_diagonal_half()*sin((i_az*pi)/180)) + 1.5*self.largest_diagonal_half()*i, 4)
                            ypp = round((data[0][3] + self.largest_diagonal_half()*cos((i_az*pi)/180)) - self.tessera_width()*j - self.tessera_width()/2, 4)
                            x_prim = xpp - data[0][2]
                            y_prim = ypp - data[0][3]
                            x_bis = x_prim*cos(orient) - y_prim*sin(orient)
                            y_bis = x_prim*sin(orient) + y_prim*cos(orient)
                            x = x_bis + data[0][2]
                            y = y_bis + data[0][3]
                            hex_coords_temp.append([x, y])
                    for point_coords in polyline_coords:
                        # Using Ray Casting Method
                        if self.ray_casting_method(hex_coords_temp, point_coords) == True:
                            points_in_hex_coords.append([id_hex, point_coords[1], point_coords[2], point_coords[3]])
                    id_hex += 1
        else:
            if (dx12 > 0 and dy12 > 0) or (dx12 < 0 and dy12 > 0):
                orient = atan2(dy12, dy12)
            else:
                orient = atan2(dy12, dy12) + 2*pi
            vertical_cover = trunc((b-(self.tessera_width()/2)) /self.tessera_width()) + 1
            horizontal_cover = trunc((a-0.5*self.largest_diagonal_half()) / (1.5*self.largest_diagonal_half())) + 2
            points_in_hex_coords = []
            id_hex = 0
            for i in range(horizontal_cover + 1):
                for j in range(vertical_cover + 1):
                    hex_coords_temp = []
                    for i_az in az:
                        if (i % 2 == 0):
                            xpp = round((data[1][2] + self.largest_diagonal_half()*sin((i_az*pi)/180)) + 1.5*self.largest_diagonal_half()*i, 4)
                            ypp = round((data[1][3] + self.largest_diagonal_half()*cos((i_az*pi)/180)) - self.tessera_width()*j, 4)
                            x_prim = xpp - data[1][2]
                            y_prim = ypp - data[1][3]
                            x_bis = x_prim*cos(orient) - y_prim*sin(orient)
                            y_bis = x_prim*sin(orient) + y_prim*cos(orient)
                            x = x_bis + data[1][2]
                            y = y_bis + data[1][3]
                            hex_coords_temp.append([x, y])
                        else:
                            xpp = round((data[1][2] + self.largest_diagonal_half()*sin((i_az*pi)/180)) + 1.5*self.largest_diagonal_half()*i, 4)
                            ypp = round((data[1][3] + self.largest_diagonal_half()*cos((i_az*pi)/180)) - self.tessera_width()*j - self.tessera_width()/2, 4)
                            x_prim = xpp - data[1][2]
                            y_prim = ypp - data[1][3]
                            x_bis = x_prim*cos(orient) - y_prim*sin(orient)
                            y_bis = x_prim*sin(orient) + y_prim*cos(orient)
                            x = x_bis + data[1][2]
                            y = y_bis + data[1][3]
                            hex_coords_temp.append([x, y])
                    for point_coords in polyline_coords:
                        # Using Ray Casting Method
                        if self.ray_casting_method(hex_coords_temp, point_coords) == True:
                            points_in_hex_coords.append([id_hex, point_coords[1], point_coords[2], point_coords[3]])
                    id_hex += 1
        # Using Vertex Clustering
        cluster = self.vertex_clustering(points_in_hex_coords)
        # Using Spatial Mean
        mean_xy = self.spatial_mean(cluster, polyline_coords)
        self.create_new_feature(mean_xy)
        return

    """
    Direction of tessellation consistent to direction of minimal rectangle area polygon starting from the corner
    of the rectangle
    """
    def minimal_rectangle_area(self):
        self.minimal_rectangle_general("RECTANGLE_BY_AREA")
        return

    """
    Direction of tessellation consistent to direction of minimal rectangle width polygon starting from the corner
    of the rectangle
    """
    def minimal_rectangle_width(self):
        self.minimal_rectangle_general("RECTANGLE_BY_WIDTH")
        return

    """
    Direction of tesselation consistent with the direction perpendicular to the direction of bisector of an angle,
    which vertex is simultaneously a global maximum of the original polyline
    """
    def global_maximum(self):
        return

    def choose_method(self):
        if self.method == 'FROM BOUNDING BOX':
            self.bounding_box()
        elif self.method == 'FROM MINIMAL RECTANGLE AREA':
            self.minimal_rectangle_area()
        elif self.method == 'FROM MINIMAL RECTANGLE WIDTH':
            self.minimal_rectangle_width()
        elif self.method == 'FROM GLOBAL MAXIMUM':
            self.global_maximum()


if __name__ == '__main__':
    polyline = HexSimply(arcpy.GetParameterAsText(0), arcpy.GetParameterAsText(1), arcpy.GetParameter(2),
                         arcpy.GetParameter(3), arcpy.GetParameterAsText(4))
    exit(polyline)
