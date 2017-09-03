#######################################################################
# Filename: HexSimply.py                                              #
# Author: Karolina Mamczarz                                           #
# Institution: AGH University of Science and Technology in Cracow,    #
#               Poland                                                #
# Faculty: Mining Surveying and Environmental Engineering             #
# Department: Integrated Geodesy and Cartography                      #
# Last update: 2017-09-02                                             #
# Version: 1.0.0                                                      #
# Description: Implementation of line simplification methods based on #
#              hexagonal tessellation generated on bounding box and   #
#              oriented rectangles according to:                      #
#               - minimal area of rectangle                           #
#               - minimum width of input data                         #
#               - the furthest point from line which joins first and  #
#                 last point of input data                            #
# Class: HexSimply                                                    #
# Methods: __init__, tessera_width, largest_diagonal_half, set_path,  #
#          create_new_feature, bounding_box, oriented_rectangle,      #
#          minimal_rectangle_area, minimal_rectangle_area_or_width,   #
#          furthest_point, choose_method                              #
# Result: Object - Simplified polyline feature class with             #
#         no self-crossing                                            #
#######################################################################

import os
from math import trunc, sin, cos
from sys import exit
from HexTools import *


class HexSimply(object):

    def __init__(self, original, method, l, s, simplified):
        self.original = original
        self.simplified = simplified
        self.l = l / 1000
        self.s = s
        self.method = method
        self.choose_method()

    def tessera_width(self):
        return 5*self.l*self.s

    def largest_diagonal_half(self):
        return round(self.tessera_width()/sqrt(3), 4)
        
    def set_path(self):
        path, my_file = os.path.split(self.simplified)
        filename, file_ext = os.path.splitext(my_file)
        arcpy.CreateFeatureclass_management(path, filename, "POLYLINE", "", "",
                                            "", self.original)
        container, container_ext = os.path.splitext(path)
        if container_ext == ".gdb" or container_ext == ".mdb" \
                or container_ext == ".sde":
            my_file = filename
        return path + "\\" + my_file

    def create_new_feature(self, new_polyline_coords):
        with arcpy.da.InsertCursor(self.set_path(), ["SHAPE@"]) as cursor:
            arc_point_list = []
            for point in new_polyline_coords:
                arc_point = arcpy.Point(point[0], point[1])
                arc_point_list.append(arc_point)
            cursor.insertRow([arcpy.Polyline(arcpy.Array(arc_point_list))])
        return

    """
    Direction of tessellation consistent horizontally to minimal area
    bounding box starting from upper-left corner of this bounding box.
    """
    def bounding_box(self):
        self.set_path()
        polyline_coords = HexTools.read_geom(self.original)
        az = [30, 90, 150, 210, 270, 330]
        max_min = []
        for coord in zip(*polyline_coords):
            max_min.append(max(coord))
            max_min.append(min(coord))
        x_lr = max_min[4]
        y_lr = max_min[7]
        x_ll = max_min[5]
        y_ll = max_min[7]
        x_ul = max_min[5]
        y_ul = max_min[6]
        a = HexTools.calculate_distance(x_lr, x_ll, y_lr, y_ll)
        b = HexTools.calculate_distance(x_ll, x_ul, y_ll, y_ul)
        vertical_cover = trunc((b-(self.tessera_width()/2))
                               / self.tessera_width()) + 1
        horizontal_cover = trunc((a-0.5*self.largest_diagonal_half())
                                 / (1.5*self.largest_diagonal_half())) + 2
        points_in_hex_coords = []
        id_hex = 0
        for i in range(horizontal_cover + 1):
            for j in range(vertical_cover + 1):
                hex_coords_temp = []
                for i_az in az:
                    if i % 2 == 0:
                        x = round((x_ul + self.largest_diagonal_half()
                                  * sin((i_az*pi)/180))
                                  + 1.5*self.largest_diagonal_half()*i, 4)
                        y = round((y_ul + self.largest_diagonal_half()
                                  * cos((i_az*pi)/180))
                                  - self.tessera_width()*j, 4)
                        hex_coords_temp.append([x, y])
                    else:
                        x = round((x_ul + self.largest_diagonal_half()
                                  * sin((i_az*pi)/180))
                                  + 1.5*self.largest_diagonal_half()*i, 4)
                        y = round((y_ul + self.largest_diagonal_half()
                                  * cos((i_az*pi)/180))
                                  - self.tessera_width()*j
                                  - self.tessera_width()/2, 4)
                        hex_coords_temp.append([x, y])
                for point_coords in polyline_coords:
                    # Using Ray Casting Method
                    if HexTools.ray_casting_method(hex_coords_temp,
                                                   point_coords) is True:
                        points_in_hex_coords.append([id_hex, point_coords[1],
                                                     point_coords[2],
                                                     point_coords[3]])
                id_hex += 1
        # Using Vertex Clustering
        cluster = HexTools.vertex_clustering(points_in_hex_coords)
        # Using Spatial Mean
        mean_xy = HexTools.spatial_mean(cluster, polyline_coords)
        # Detecting self-crossing and creating simplified polyline
        self.create_new_feature(HexTools.eliminate_self_crossing(mean_xy))
        return

    def oriented_rectangle(self, a, b, x0, x1, x2, y0, y1, y2,
                           polyline_coords):
        az = [30, 90, 150, 210, 270, 330]
        if a > b:
            orient = HexTools.azimuth(x1-x0, y1-y0)
            oriented_vertical_cover = trunc((b-(self.tessera_width()/2))
                                   / self.tessera_width()) + 1
            oriented_horizontal_cover = \
                trunc((a-0.5*self.largest_diagonal_half())
                      / (1.5*self.largest_diagonal_half())) + 2
            points_in_hex_coords = []
            id_hex = 0
            for i in range(oriented_horizontal_cover + 1):
                for j in range(oriented_vertical_cover + 1):
                    hex_coords_temp = []
                    for i_az in az:
                        if i % 2 == 0:
                            xpp = round((x0 + self.largest_diagonal_half()
                                        * sin((i_az*pi)/180))
                                        + 1.5*self.largest_diagonal_half()*i,
                                        4)
                            ypp = round((y0 + self.largest_diagonal_half()
                                        * cos((i_az*pi)/180))
                                        - self.tessera_width()*j, 4)
                            x_prim = xpp - x0
                            y_prim = ypp - y0
                            x_bis = x_prim*cos(orient) - y_prim*sin(orient)
                            y_bis = x_prim*sin(orient) + y_prim*cos(orient)
                            x = x_bis + x0
                            y = y_bis + y0
                            hex_coords_temp.append([x, y])
                        else:
                            xpp = round((x0 + self.largest_diagonal_half()
                                        * sin((i_az*pi)/180))
                                        + 1.5*self.largest_diagonal_half()*i,
                                        4)
                            ypp = round((y0 + self.largest_diagonal_half()
                                        * cos((i_az*pi)/180))
                                        - self.tessera_width()*j
                                        - self.tessera_width()/2, 4)
                            x_prim = xpp - x0
                            y_prim = ypp - y0
                            x_bis = x_prim*cos(orient) - y_prim*sin(orient)
                            y_bis = x_prim*sin(orient) + y_prim*cos(orient)
                            x = x_bis + x0
                            y = y_bis + y0
                            hex_coords_temp.append([x, y])
                    for point_coords in polyline_coords:
                        # Using Ray Casting Method
                        if HexTools.ray_casting_method(hex_coords_temp,
                                                       point_coords) is True:
                            points_in_hex_coords.append([id_hex,
                                                         point_coords[1],
                                                         point_coords[2],
                                                         point_coords[3]])
                    id_hex += 1
        else:
            orient = HexTools.azimuth(x2-x1, y2-y1)
            oriented_vertical_cover = trunc((a-(self.tessera_width()/2))
                                   / self.tessera_width()) + 1
            oriented_horizontal_cover = \
                trunc((b-0.5*self.largest_diagonal_half())
                      / (1.5*self.largest_diagonal_half())) + 2
            points_in_hex_coords = []
            id_hex = 0
            for i in range(oriented_horizontal_cover + 1):
                for j in range(oriented_vertical_cover + 1):
                    hex_coords_temp = []
                    for i_az in az:
                        if i % 2 == 0:
                            xpp = round((x1 + self.largest_diagonal_half()
                                        * sin((i_az*pi)/180))
                                        + 1.5*self.largest_diagonal_half()*i,
                                        4)
                            ypp = round((y1 + self.largest_diagonal_half()
                                        * cos((i_az*pi)/180))
                                        - self.tessera_width()*j, 4)
                            x_prim = xpp - x1
                            y_prim = ypp - y1
                            x_bis = x_prim*cos(orient) - y_prim*sin(orient)
                            y_bis = x_prim*sin(orient) + y_prim*cos(orient)
                            x = x_bis + x1
                            y = y_bis + y1
                            hex_coords_temp.append([x, y])
                        else:
                            xpp = round((x1 + self.largest_diagonal_half()
                                        * sin((i_az*pi)/180))
                                        + 1.5*self.largest_diagonal_half()*i,
                                        4)
                            ypp = round((y1 + self.largest_diagonal_half()
                                        * cos((i_az*pi)/180))
                                        - self.tessera_width()*j
                                        - self.tessera_width()/2, 4)
                            x_prim = xpp - x1
                            y_prim = ypp - y1
                            x_bis = x_prim*cos(orient) - y_prim*sin(orient)
                            y_bis = x_prim*sin(orient) + y_prim*cos(orient)
                            x = x_bis + x1
                            y = y_bis + y1
                            hex_coords_temp.append([x, y])
                    for point_coords in polyline_coords:
                        # Using Ray Casting Method
                        if HexTools.ray_casting_method(hex_coords_temp,
                                                       point_coords) is True:
                            points_in_hex_coords.append([id_hex,
                                                         point_coords[1],
                                                         point_coords[2],
                                                         point_coords[3]])
                    id_hex += 1
        # Using Vertex Clustering
        cluster = HexTools.vertex_clustering(points_in_hex_coords)
        # Using Spatial Mean
        mean_xy = HexTools.spatial_mean(cluster, polyline_coords)
        # Detecting self-crossing and creating simplified polyline
        self.create_new_feature(HexTools.eliminate_self_crossing(mean_xy))
        return

    """
    Direction of tessellation consistent to direction of minimal area
    rectangle polygon starting from the corner of the rectangle
    """
    def minimal_rectangle_area(self):
        self.set_path()
        polyline_coords = HexTools.read_geom(self.original)
        temp_rect_area = "in_memory\\rect_area"
        arcpy.MinimumBoundingGeometry_management(self.original, temp_rect_area,
                                                 "RECTANGLE_BY_AREA", "ALL")
        data = HexTools.read_geom(temp_rect_area)
        x0, x1, x2, y0, y1, y2 = data[0][2], data[1][2], data[2][2], \
                                 data[0][3], data[1][3], data[2][3]
        a = HexTools.calculate_distance(x0, x1, y0, y1)
        b = HexTools.calculate_distance(x1, x2, y1, y2)
        self.oriented_rectangle(a, b, x0, x1, x2, y0, y1, y2, polyline_coords)
        return

    """
    Direction of tessellation consistent to direction of minimal width
    rectangle starting from the corner of the rectangle. It is not a
    pure method for this algorithm, because it is divided with
    conditions for minimal area rectangle.
    """
    def minimal_rectangle_area_or_width(self):
        self.set_path()
        polyline_coords = HexTools.read_geom(self.original)
        temp_rect_width = "in_memory\\rect_width"
        arcpy.MinimumBoundingGeometry_management(self.original,
                                                 temp_rect_width,
                                                 "RECTANGLE_BY_WIDTH", "ALL")
        data = HexTools.read_geom(temp_rect_width)
        x0, x1, x2, y0, y1, y2 = data[0][2], data[1][2], data[2][2], \
                                 data[0][3], data[1][3], data[2][3]
        a = HexTools.calculate_distance(x0, x1, y0, y1)
        b = HexTools.calculate_distance(x1, x2, y1, y2)
        first_last_distance = HexTools.calculate_distance(
            polyline_coords[0][2], polyline_coords[-1][2],
            polyline_coords[0][3], polyline_coords[-1][3])
        if a < b:
            if first_last_distance < a:
                self.oriented_rectangle(a, b, x0, x1, x2, y0, y1, y2,
                                        polyline_coords)
            else:
                self.minimal_rectangle_area()
        else:
            if first_last_distance < b:
                self.oriented_rectangle(a, b, x0, x1, x2, y0, y1, y2,
                                        polyline_coords)
            else:
                self.minimal_rectangle_area()
        return

    """
    Direction of tessellation consistent with the directions of two
    bisectors of vertically opposite angles, which vertex is
    simultaneously a furthest point of the original polyline to the
    line which joins first and last point. Angles are created with
    lines:
        - first point - furthest vertex,
        - furthest vertex - last point.
    """
    def furthest_point(self):
        self.set_path()
        polyline_coords = HexTools.read_geom(self.original)
        x_first = polyline_coords[0][2]
        y_first = polyline_coords[0][3]
        x_last = polyline_coords[-1][2]
        y_last = polyline_coords[-1][3]
        a_m, b_m = HexTools.coefficients_linear_function(x_first, y_first,
                                                         x_last, y_last)
        d_major_line, x_major_line, y_major_line = \
            HexTools.point_to_line_distance(x_first, y_first, polyline_coords,
                                            a_m, -1, b_m)
        a_first, b_first = \
            HexTools.coefficients_linear_function(x_first, y_first,
                                                  x_major_line, y_major_line)
        a_last, b_last = \
            HexTools.coefficients_linear_function(x_major_line, y_major_line,
                                                  x_last, y_last)
        a_b1, b_b1, a_b2, b_b2 = \
            HexTools.coefficients_general_equation(a_first, a_last, -1, -1,
                                                   b_first, b_last)
        result_b1 = \
            HexTools.point_to_line_distance_with_sides(x_first, y_first,
                                                       polyline_coords,
                                                       a_b1, -1, b_b1)
        result_b2 = \
            HexTools.point_to_line_distance_with_sides(x_first, y_first,
                                                       polyline_coords,
                                                       a_b2, -1, b_b2)
        lines_coefficients = [[a_b1, -1, result_b1["y_point_one_side"]
                               - a_b1*result_b1["x_point_one_side"]],
                              [a_b1, -1, result_b1["y_point_other_side"]
                               - a_b1*result_b1["x_point_other_side"]],
                              [a_b2, -1, result_b2["y_point_one_side"]
                               - a_b2*result_b2["x_point_one_side"]],
                              [a_b2, -1, result_b2["y_point_other_side"]
                               - a_b2*result_b2["x_point_other_side"]]]
        xy_inter = HexTools.intersection_vertices(lines_coefficients)
        xy_for_side_temp = []
        for xy1 in xy_inter:
            for xy2 in xy_inter:
                dxy_o = HexTools.calculate_distance(xy1[0], xy2[0], xy1[1],
                                                    xy2[1])
                if dxy_o != 0.0:
                    xy_for_side_temp.append([dxy_o, xy1[0], xy1[1], xy2[0],
                                             xy2[1]])
            break
        del_id = max(xy_for_side_temp, key=lambda item: item[0])
        xy_for_side = []
        for element in xy_for_side_temp:
            if element == del_id:
                xy_for_side_temp.remove(element)
            else:
                xy_for_side.append([0, 0, element[1], element[2]])
                xy_for_side.append([0, 0, element[3], element[4]])
        xy_for_side.pop(0)
        dx1 = xy_for_side[1][2]-xy_for_side[0][2]
        dy1 = xy_for_side[1][3]-xy_for_side[0][3]
        dx2 = xy_for_side[2][2]-xy_for_side[1][2]
        dy2 = xy_for_side[2][3]-xy_for_side[1][3]
        angle = HexTools.azimuth(dx1, dy1) - HexTools.azimuth(dx2, dy2)
        if angle > pi/2:
            xy_for_side[2], xy_for_side[0] = xy_for_side[0], xy_for_side[2]
        x0, x1, x2, y0, y1, y2 = xy_for_side[0][2], xy_for_side[1][2], \
                                 xy_for_side[2][2], xy_for_side[0][3], \
                                 xy_for_side[1][3], xy_for_side[2][3]
        a = HexTools.calculate_distance(x0, x1, y0, y1)
        b = HexTools.calculate_distance(x1, x2, y1, y2)
        self.oriented_rectangle(a, b, x0, x1, x2, y0, y1, y2, polyline_coords)
        return

    def choose_method(self):
        if self.method == 'FROM BOUNDING BOX':
            self.bounding_box()
        elif self.method == \
                'FROM MINIMAL RECTANGLE WIDTH OR MINIMAL RECTANGLE AREA':
            self.minimal_rectangle_area_or_width()
        elif self.method == 'FROM THE FURTHEST POINT OF POLYLINE':
            self.furthest_point()
        return

if __name__ == '__main__':
    polyline = HexSimply(arcpy.GetParameterAsText(0),
                         arcpy.GetParameterAsText(1), arcpy.GetParameter(2),
                         arcpy.GetParameter(3), arcpy.GetParameterAsText(4))
    exit(polyline)
