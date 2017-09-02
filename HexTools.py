#######################################################################
# Filename: HexTools.py                                               #
# Author: Karolina Mamczarz                                           #
# Institution: AGH University of Science and Technology in Cracow,    #
#               Poland                                                #
# Faculty: Mining Surveying and Environmental Engineering             #
# Department: Integrated Geodesy and Cartography                      #
# Last update: 2017-09-02                                             #
# Version: 1.0.0                                                      #
# Description: Tools used by HexSimply script, encloses mathematical  #
#              functions or algorithms and methods for data management#
# Class: HexTools                                                     #
# Methods: read_geom, read_geom_attr_rectangle, ray_casting_method,   #
#          vertex_clustering, spatial_mean, point_to_line_distance,   #
#          point_to_line_distance_with_sides,                         #
#          coefficients_linear_function,                              #
#          coefficients_general_equation, intersection_vertices,      #
#          calculate_distance, azimuth, eliminate_self_crossing       #
# Result: Data passed with lists, variables, dictionaries, boolean    #
#######################################################################

import arcpy
from math import sqrt, pi, atan2, fabs


class HexTools:

    @staticmethod
    def read_geom(shape):
        with arcpy.da.SearchCursor(shape, ["SHAPE@"]) as cursor:
            shape_coords = []
            partnum = 0
            for row in cursor:
                for part in row[0]:
                    pntnum = 0
                    for pnt in part:
                        shape_coords.append([partnum, pntnum, pnt.X, pnt.Y])
                        pntnum += 1
                    partnum += 1
        return shape_coords

    @staticmethod
    def read_geom_attr_rectangle(rectangle):
        with arcpy.da.SearchCursor(rectangle, ["SHAPE@"]) as cursor:
            rect_area_poly_coords = []
            for row in cursor:
                partnum = 0
                for part in row[0]:
                    pntnum = 0
                    for pnt in part:
                        rect_area_poly_coords.append([partnum, pntnum, pnt.X,
                                                      pnt.Y])
                        pntnum += 1
                    partnum += 1
        return rect_area_poly_coords

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
                            xints = (iter_polyline_coords[3]-p1y) * (p2x-p1x) \
                                    / (p2y-p1y) + p1x
                        if p1x == p2x or iter_polyline_coords[2] <= xints:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    @staticmethod
    def vertex_clustering(points_list):
        points_list.sort(key=lambda id_point: id_point[1])
        cluster = []
        prev = 0
        for point in points_list:
            if point[0] - prev != 0:
                prev = point[0]
                cluster.append([point])
            elif point[0] == 0:
                idx = points_list.index(point)
                if idx == 0:
                    prev = point[0]
                    cluster.append([point])
                else:
                    cluster[-1].append(point)
            else:
                cluster[-1].append(point)
        return cluster

    @staticmethod
    def spatial_mean(cluster_list, points_list):
        mean_xy_list = []
        mean_xy_list_clean = []
        mean_xy_list.append([points_list[0][2], points_list[0][3]])
        for one_cluster in cluster_list:
            for xy in one_cluster:
                n = one_cluster.index(xy) + 1
            sums = [sum(i) for i in zip(*one_cluster)]
            mean_x = sums[2]/n
            mean_y = sums[3]/n
            mean_xy_list.append([mean_x, mean_y])
        mean_xy_list.append([points_list[-1][2], points_list[-1][3]])
        for duplicate in mean_xy_list:
            if duplicate not in mean_xy_list_clean:
                mean_xy_list_clean.append(duplicate)
        return mean_xy_list_clean

    @staticmethod
    def point_to_line_distance(initial_x, initial_y, set_of_coords, a, b, c):
        d = 0.0
        x_point = initial_x
        y_point = initial_y
        for point in set_of_coords:
            distance = fabs(a*point[2]+b*point[3]+c) / sqrt(a**2 + b**2)
            if distance > d:
                d = distance
                x_point = point[2]
                y_point = point[3]
        return d, x_point, y_point

    @staticmethod
    def point_to_line_distance_with_sides(initial_x, initial_y,
                                          set_of_coords, a, b, c):
        result = {"d_one_side": 0.0,
                  "x_point_one_side": initial_x,
                  "y_point_one_side": initial_y,
                  "d_other_side": 0.0,
                  "x_point_other_side": initial_x,
                  "y_point_other_side": initial_y}
        distance_exception = None
        x_exception = None
        y_exception = None
        for point in set_of_coords:
            equation = a*point[2]+b*point[3]+c
            if equation > 0:
                distance = fabs(equation) / sqrt(a**2 + b**2)
                if distance > result.get("d_one_side"):
                    result.update({"d_one_side": distance,
                                   "x_point_one_side": point[2],
                                   "y_point_one_side": point[3]})
            elif equation < 0:
                distance = fabs(equation) / sqrt(a**2 + b**2)
                if distance > result.get("d_other_side"):
                    result.update({"d_other_side": distance,
                                   "x_point_other_side": point[2],
                                   "y_point_other_side": point[3]})
            elif equation == 0:
                distance_exception = fabs(equation) / sqrt(a**2 + b**2)
                x_exception = point[2]
                y_exception = point[3]
        for k, v in dict(result).items():
            if k == "d_one_side" and v == 0.0:
                result.update({"d_one_side": distance_exception,
                               "x_point_one_side": x_exception,
                               "y_point_one_side": y_exception})
            if k == "d_other_side" and v == 0.0:
                result.update({"d_other_side": distance_exception,
                               "x_point_other_side": x_exception,
                               "y_point_other_side": y_exception})
        return result

    @staticmethod
    def coefficients_linear_function(x1, y1, x2, y2):
        a = (y2-y1) / (x2-x1)
        b = y1 - a*x1
        return a, b

    @staticmethod
    def coefficients_general_equation(a1, a2, b1, b2, c1, c2):
        k = sqrt((a2**2 + b2**2) / (a1**2 + b1**2))
        a_bisector_one = -(k*a1-a2) / (k*b1-b2)
        b_bisector_one = -((k*c1-c2) / (k*b1-b2))
        a_bisector_two = -(k*a1+a2) / (k*b1+b2)
        b_bisector_two = -((k*c1+c2) / (k*b1+b2))
        return a_bisector_one, b_bisector_one, a_bisector_two, b_bisector_two

    @staticmethod
    def intersection_vertices(coefficient_array):
        xy_temp = []
        for i in coefficient_array:
            for j in coefficient_array:
                w = i[0]*j[1]-j[0]*i[1]
                if w != 0:
                    wx = (-i[2])*j[1] - (-j[2])*i[1]
                    wy = i[0]*(-j[2]) - j[0]*(-i[2])
                    x = wx/w
                    y = wy/w
                    xy_temp.append([x, y])
        xy = []
        for pair in xy_temp:
            if pair not in xy:
                xy.append(pair)
        return xy

    @staticmethod
    def calculate_distance(x1, x2, y1, y2):
        d = sqrt((x2-x1)**2 + (y2-y1)**2)
        return d

    @staticmethod
    def azimuth(dx, dy):
        if (dx > 0 and dy > 0) or dy > 0 > dx:
            azimuth = atan2(dy, dx)
        else:
            azimuth = atan2(dy, dx) + 2*pi
        return azimuth

    @staticmethod
    def eliminate_self_crossing(maybe_tangled_line):
        points_to_revert = []
        points = len(maybe_tangled_line)
        i_point = 0
        while i_point < points - 3:
            j_point = i_point + 2
            while j_point < points - 1:
                dx_ab = maybe_tangled_line[i_point + 1][0] \
                        - maybe_tangled_line[i_point][0]
                dx_ac = maybe_tangled_line[j_point][0] \
                        - maybe_tangled_line[i_point][0]
                dx_cd = maybe_tangled_line[j_point + 1][0] \
                        - maybe_tangled_line[j_point][0]
                dy_ab = maybe_tangled_line[i_point + 1][1] \
                        - maybe_tangled_line[i_point][1]
                dy_ac = maybe_tangled_line[j_point][1] \
                        - maybe_tangled_line[i_point][1]
                dy_cd = maybe_tangled_line[j_point + 1][1] \
                        - maybe_tangled_line[j_point][1]
                k = (dx_ac * dy_cd - dx_cd * dy_ac) \
                    / (dx_ab * dy_cd - dx_cd * dy_ab)
                xp = maybe_tangled_line[i_point][0] + k*dx_ab
                yp = maybe_tangled_line[i_point][1] + k*dy_ab
                if maybe_tangled_line[i_point + 1][0] \
                        - maybe_tangled_line[i_point][0] > 0:
                    xp_range_first_seg = maybe_tangled_line[i_point][0] \
                                         < xp \
                                         < maybe_tangled_line[i_point + 1][0]
                else:
                    xp_range_first_seg = maybe_tangled_line[i_point][0] \
                                         > xp \
                                         > maybe_tangled_line[i_point + 1][0]
                if maybe_tangled_line[j_point + 1][0] \
                        - maybe_tangled_line[j_point][0] > 0:
                    xp_range_second_seg = maybe_tangled_line[j_point][0] \
                                          < xp \
                                          < maybe_tangled_line[j_point + 1][0]
                else:
                    xp_range_second_seg = maybe_tangled_line[j_point][0] \
                                          > xp \
                                          > maybe_tangled_line[j_point + 1][0]
                if maybe_tangled_line[i_point + 1][1] \
                        - maybe_tangled_line[i_point][1] > 0:
                    yp_range_first_seg = maybe_tangled_line[i_point][1] \
                                         < yp \
                                         < maybe_tangled_line[i_point + 1][1]
                else:
                    yp_range_first_seg = maybe_tangled_line[i_point][1] \
                                         > yp \
                                         > maybe_tangled_line[i_point + 1][1]
                if maybe_tangled_line[j_point + 1][1] \
                        - maybe_tangled_line[j_point][1] > 0:
                    yp_range_second_seg = maybe_tangled_line[j_point][1] \
                                          < yp \
                                          < maybe_tangled_line[j_point + 1][1]
                else:
                    yp_range_second_seg = maybe_tangled_line[j_point][1] \
                                          > yp \
                                          > maybe_tangled_line[j_point + 1][1]
                xp_in_range = xp_range_first_seg and xp_range_second_seg
                yp_in_range = yp_range_first_seg and yp_range_second_seg
                if xp_in_range and yp_in_range:
                    points_to_revert.append([i_point + 1, j_point])
                j_point += 1
            i_point += 1
        for reverse_pair in points_to_revert:
            maybe_tangled_line[reverse_pair[0]], \
            maybe_tangled_line[reverse_pair[1]] \
                = maybe_tangled_line[reverse_pair[1]], \
                  maybe_tangled_line[reverse_pair[0]]
        return maybe_tangled_line
