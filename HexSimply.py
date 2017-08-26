#############################################################################
# Filename: HexSimply.py                                                    #
# Author: Karolina Mamczarz                                                 #
# Institution: AGH University of Science and Technology in Cracow, Poland   #
# Last update: 2017-08-24                                                   #
# Version: 1.0.0                                                            #
# Description: #                                                            #
# Class: HexSimply                                                          #
# Methods: #                                                                #
# Parameters: #                                                             #
# Result: #                                                                 #
#############################################################################

import arcpy
import os
from math import sqrt, trunc, sin, cos, pi, atan2, fabs
from sys import exit


class HexSimply(object):

    def __init__(self, original, method, l, s, simplified):
        self.original = original
        self.simplified = simplified
        self.l = l / 1000
        self.s = s
        self.method = method
        self.choose_method()

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
    def point_to_line_distance_with_sides(initial_x, initial_y, set_of_coords, a, b, c):
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
                    result.update({"d_one_side": distance, "x_point_one_side": point[2], "y_point_one_side": point[3]})
            elif equation < 0:
                distance = fabs(equation) / sqrt(a**2 + b**2)
                if distance > result.get("d_other_side"):
                    result.update({"d_other_side": distance, "x_point_other_side": point[2], "y_point_other_side": point[3]})
            elif equation == 0:
                distance_exception = fabs(equation) / sqrt(a**2 + b**2)
                x_exception = point[2]
                y_exception = point[3]
        for k, v in dict(result).items():
            if k == "d_one_side" and v == 0.0:
                result.update({"d_one_side": distance_exception, "x_point_one_side": x_exception, "y_point_one_side": y_exception})
            if k == "d_other_side" and v == 0.0:
                result.update({"d_other_side": distance_exception, "x_point_other_side": x_exception, "y_point_other_side": y_exception})
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
                dx_ab = maybe_tangled_line[i_point + 1][0] - maybe_tangled_line[i_point][0]
                dx_ac = maybe_tangled_line[j_point][0] - maybe_tangled_line[i_point][0]
                dx_cd = maybe_tangled_line[j_point + 1][0] - maybe_tangled_line[j_point][0]
                dy_ab = maybe_tangled_line[i_point + 1][1] - maybe_tangled_line[i_point][1]
                dy_ac = maybe_tangled_line[j_point][1] - maybe_tangled_line[i_point][1]
                dy_cd = maybe_tangled_line[j_point + 1][1] - maybe_tangled_line[j_point][1]
                k = (dx_ac * dy_cd - dx_cd * dy_ac) / (dx_ab * dy_cd - dx_cd * dy_ab)
                xp = maybe_tangled_line[i_point][0] + k*dx_ab
                yp = maybe_tangled_line[i_point][1] + k*dy_ab
                if maybe_tangled_line[i_point + 1][0] - maybe_tangled_line[i_point][0] > 0:
                    xp_range_first_seg = maybe_tangled_line[i_point][0] < xp < maybe_tangled_line[i_point + 1][0]
                else:
                    xp_range_first_seg = maybe_tangled_line[i_point][0] > xp > maybe_tangled_line[i_point + 1][0]
                if maybe_tangled_line[j_point + 1][0] - maybe_tangled_line[j_point][0] > 0:
                    xp_range_second_seg = maybe_tangled_line[j_point][0] < xp < maybe_tangled_line[j_point + 1][0]
                else:
                    xp_range_second_seg = maybe_tangled_line[j_point][0] > xp > maybe_tangled_line[j_point + 1][0]
                if maybe_tangled_line[i_point + 1][1] - maybe_tangled_line[i_point][1] > 0:
                    yp_range_first_seg = maybe_tangled_line[i_point][1] < yp < maybe_tangled_line[i_point + 1][1]
                else:
                    yp_range_first_seg = maybe_tangled_line[i_point][1] > yp > maybe_tangled_line[i_point + 1][1]
                if maybe_tangled_line[j_point + 1][1] - maybe_tangled_line[j_point][1] > 0:
                    yp_range_second_seg = maybe_tangled_line[j_point][1] < yp < maybe_tangled_line[j_point + 1][1]
                else:
                    yp_range_second_seg = maybe_tangled_line[j_point][1] > yp > maybe_tangled_line[j_point + 1][1]
                xp_in_range = xp_range_first_seg and xp_range_second_seg
                yp_in_range = yp_range_first_seg and yp_range_second_seg
                if xp_in_range and yp_in_range:
                    points_to_revert.append([i_point + 1, j_point])
                j_point += 1
            i_point += 1
        for reverse_pair in points_to_revert:
            maybe_tangled_line[reverse_pair[0]], maybe_tangled_line[reverse_pair[1]] = maybe_tangled_line[reverse_pair[1]], maybe_tangled_line[reverse_pair[0]]
        return maybe_tangled_line

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
            arc_point_list = []
            for point in new_polyline_coords:
                arc_point = arcpy.Point(point[0], point[1])
                arc_point_list.append(arc_point)
            cursor.insertRow([arcpy.Polyline(arcpy.Array(arc_point_list))])
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
        a = self.calculate_distance(x_ur, x_ul, y_ur, y_ul)
        b = self.calculate_distance(x_dl, x_ul, y_dl, y_ul)
        vertical_cover = trunc((b-(self.tessera_width()/2)) / self.tessera_width()) + 1
        horizontal_cover = trunc((a-0.5*self.largest_diagonal_half()) / (1.5*self.largest_diagonal_half())) + 2
        points_in_hex_coords = []
        id_hex = 0
        for i in range(horizontal_cover + 1):
            for j in range(vertical_cover + 1):
                hex_coords_temp = []
                for i_az in az:
                    if i % 2 == 0:
                        x = round((x_ul + self.largest_diagonal_half()*sin((i_az*pi)/180)) + 1.5*self.largest_diagonal_half()*i, 4)
                        y = round((y_ul + self.largest_diagonal_half()*cos((i_az*pi)/180)) - self.tessera_width()*j, 4)
                        hex_coords_temp.append([x, y])
                    else:
                        x = round((x_ul + self.largest_diagonal_half()*sin((i_az*pi)/180)) + 1.5*self.largest_diagonal_half()*i, 4)
                        y = round((y_ul + self.largest_diagonal_half()*cos((i_az*pi)/180)) - self.tessera_width()*j - self.tessera_width()/2, 4)
                        hex_coords_temp.append([x, y])
                for point_coords in polyline_coords:
                    # Using Ray Casting Method
                    if self.ray_casting_method(hex_coords_temp, point_coords) is True:
                        points_in_hex_coords.append([id_hex, point_coords[1], point_coords[2], point_coords[3]])
                id_hex += 1
        # Using Vertex Clustering
        cluster = self.vertex_clustering(points_in_hex_coords)
        # Using Spatial Mean
        mean_xy = self.spatial_mean(cluster, polyline_coords)
        # Detecting self-crossing and creating simplified polyline
        self.create_new_feature(self.eliminate_self_crossing(mean_xy))
        return

    def rectangle_general(self, data, polyline_coords):
        az = [30, 90, 150, 210, 270, 330]
        dx01 = data[1][2]-data[0][2]
        dy01 = data[1][3]-data[0][3]
        dx12 = data[2][2]-data[1][2]
        dy12 = data[2][3]-data[1][3]
        a = self.calculate_distance(data[0][2], data[1][2], data[0][3], data[1][3])
        b = self.calculate_distance(data[1][2], data[2][2], data[1][3], data[2][3])
        if a > b:
            orient = self.azimuth(dx01, dy01)
            vertical_cover = trunc((b-(self.tessera_width()/2)) / self.tessera_width()) + 1
            horizontal_cover = trunc((a-0.5*self.largest_diagonal_half()) / (1.5*self.largest_diagonal_half())) + 2
            points_in_hex_coords = []
            id_hex = 0
            for i in range(horizontal_cover + 1):
                for j in range(vertical_cover + 1):
                    hex_coords_temp = []
                    for i_az in az:
                        if i % 2 == 0:
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
                        if self.ray_casting_method(hex_coords_temp, point_coords) is True:
                            points_in_hex_coords.append([id_hex, point_coords[1], point_coords[2], point_coords[3]])
                    id_hex += 1
        else:
            orient = self.azimuth(dx12, dy12)
            vertical_cover = trunc((a-(self.tessera_width()/2)) / self.tessera_width()) + 1
            horizontal_cover = trunc((b-0.5*self.largest_diagonal_half()) / (1.5*self.largest_diagonal_half())) + 2
            points_in_hex_coords = []
            id_hex = 0
            for i in range(horizontal_cover + 1):
                for j in range(vertical_cover + 1):
                    hex_coords_temp = []
                    for i_az in az:
                        if i % 2 == 0:
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
                        if self.ray_casting_method(hex_coords_temp, point_coords) is True:
                            points_in_hex_coords.append([id_hex, point_coords[1], point_coords[2], point_coords[3]])
                    id_hex += 1
        # Using Vertex Clustering
        cluster = self.vertex_clustering(points_in_hex_coords)
        # Using Spatial Mean
        mean_xy = self.spatial_mean(cluster, polyline_coords)
        # Detecting self-crossing and creating simplified polyline
        self.create_new_feature(self.eliminate_self_crossing(mean_xy))
        return

    """
    Direction of tessellation consistent to direction of minimal rectangle area polygon starting from the corner
    of the rectangle
    """
    def minimal_rectangle_area(self):
        self.set_path()
        polyline_coords = self.read_geom(self.original)
        temp_rect_area = "in_memory\\rect_area"
        arcpy.MinimumBoundingGeometry_management(self.original, temp_rect_area, "RECTANGLE_BY_AREA", "ALL")
        data = self.read_geom(temp_rect_area)
        self.rectangle_general(data, polyline_coords)
        return

    """
    Direction of tessellation consistent to direction of minimal rectangle width polygon starting from the corner
    of the rectangle
    """
    def minimal_rectangle_width(self):
        self.set_path()
        polyline_coords = self.read_geom(self.original)
        temp_rect_area = "in_memory\\rect_area"
        arcpy.MinimumBoundingGeometry_management(self.original, temp_rect_area, "RECTANGLE_BY_WIDTH", "ALL")
        data = self.read_geom(temp_rect_area)
        self.rectangle_general(data, polyline_coords)
        return

    """
    Direction of tessellation consistent with the directions of two bisectors of vertically opposite angles,
    which vertex is simultaneously a furthest point of the original polyline to the line which joins first and last
    point. Angles are created with lines:
        - first point - furthest vertex,
        - furthest vertex - last point.
    """
    def furthest_point(self):
        self.set_path()
        polyline_coords = self.read_geom(self.original)
        x_first = polyline_coords[0][2]
        y_first = polyline_coords[0][3]
        x_last = polyline_coords[-1][2]
        y_last = polyline_coords[-1][3]
        a_m, b_m = self.coefficients_linear_function(x_first, y_first, x_last, y_last)
        d_major_line, x_major_line, y_major_line = self.point_to_line_distance(x_first, y_first, polyline_coords, a_m, -1, b_m)
        a_first, b_first = self.coefficients_linear_function(x_first, y_first, x_major_line, y_major_line)
        a_last, b_last = self.coefficients_linear_function(x_major_line, y_major_line, x_last, y_last)
        a_b1, b_b1, a_b2, b_b2 = self.coefficients_general_equation(a_first, a_last, -1, -1, b_first, b_last)
        result_b1 = self.point_to_line_distance_with_sides(x_first, y_first, polyline_coords, a_b1, -1, b_b1)
        result_b2 = self.point_to_line_distance_with_sides(x_first, y_first, polyline_coords, a_b2, -1, b_b2)
        lines_coefficients = [[a_b1, -1, result_b1["y_point_one_side"] - a_b1*result_b1["x_point_one_side"]],
                              [a_b1, -1, result_b1["y_point_other_side"] - a_b1*result_b1["x_point_other_side"]],
                              [a_b2, -1, result_b2["y_point_one_side"] - a_b2*result_b2["x_point_one_side"]],
                              [a_b2, -1, result_b2["y_point_other_side"] - a_b2*result_b2["x_point_other_side"]]]
        xy_inter = self.intersection_vertices(lines_coefficients)
        xy_for_side_temp = []
        for xy1 in xy_inter:
            for xy2 in xy_inter:
                dxy_o = self.calculate_distance(xy1[0], xy2[0], xy1[1], xy2[1])
                if dxy_o != 0.0:
                    xy_for_side_temp.append([dxy_o, xy1[0], xy1[1], xy2[0], xy2[1]])
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
        angle = self.azimuth(dx1, dy1) - self.azimuth(dx2, dy2)
        if angle > pi/2:
            xy_for_side[2], xy_for_side[0] = xy_for_side[0], xy_for_side[2]
        self.rectangle_general(xy_for_side, polyline_coords)
        return

    def choose_method(self):
        if self.method == 'FROM BOUNDING BOX':
            self.bounding_box()
        elif self.method == 'FROM MINIMAL RECTANGLE AREA':
            self.minimal_rectangle_area()
        elif self.method == 'FROM MINIMAL RECTANGLE WIDTH':
            self.minimal_rectangle_width()
        elif self.method == 'FROM THE FURTHEST POINT OF POLYLINE':
            self.furthest_point()


if __name__ == '__main__':
    polyline = HexSimply(arcpy.GetParameterAsText(0), arcpy.GetParameterAsText(1), arcpy.GetParameter(2),
                         arcpy.GetParameter(3), arcpy.GetParameterAsText(4))
    exit(polyline)
