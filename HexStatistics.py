#######################################################################
# Filename: HexStatistics.py                                          #
# Author: Karolina Mamczarz                                           #
# Institution: AGH University of Science and Technology in Cracow,    #
#               Poland                                                #
# Faculty: Mining Surveying and Environmental Engineering             #
# Department: Integrated Geodesy and Cartography                      #
# Last update: 2017-09-12                                             #
# Version: 1.0.0                                                      #
# Description: Statistics used to test quality of line simplification #
# Class: HexStatistics                                                #
# Methods:  __init__, comparison_grid_1_sqkm, point_count,            #
#           std_point_per_area, surface_in_between,                   #
#           surface_in_between_per_area, hausdorff_distance           #
# Result: Report in a text file                                       #
#######################################################################

from sys import exit
import arcpy
from math import sqrt


class HexStatistics(object):

    def __init__(self, original, simplified, report_file_path):
        self.original = original
        self.simplified = simplified
        self.report_file_path = report_file_path
        self.report = open(self.report_file_path + '.txt', 'w')
        self.grid = self.comparison_grid_1_sqkm()
        self.calculate_statistics()

    def comparison_grid_1_sqkm(self):
        temp_merge = "in_memory\\merge"
        arcpy.Merge_management([self.original, self.simplified], temp_merge)
        temp_grid = "in_memory\\grid"
        #temp_grid = r"C:\Student\grid.shp"
        arcpy.GridIndexFeatures_cartography(
            temp_grid, temp_merge, "INTERSECTFEATURE", "NO_USEPAGEUNIT", "#",
            "1 Kilometers", "1 Kilometers")
        arcpy.Delete_management(temp_merge)
        arcpy.AddMessage("Comparison grid (1 sq km) prepared.")
        return temp_grid

    def point_count(self):           # Dodac obliczenie na długosc linii
        self.report.write("\n---Point count per length---\n")
        with arcpy.da.SearchCursor(self.original, ["SHAPE@"]) as cursor:
            original_pntnum = 0
            for row in cursor:
                for part in row[0]:
                    for pnt in part:
                        original_pntnum += 1
        self.report.write("Number of points before simplification: {0}\n"
                          .format(original_pntnum))
        with arcpy.da.SearchCursor(self.simplified, ["SHAPE@"]) as cursor:
            simplified_pntnum = 0
            for row in cursor:
                for part in row[0]:
                    for pnt in part:
                        simplified_pntnum += 1
        self.report.write("Number of points after simplification: {0}\n"
                          .format(simplified_pntnum))
        percentage = (float(simplified_pntnum) / float(original_pntnum))*100
        self.report.write("Percentage change in number of coordinates: "
                          + "{:0.2f}%\n"
                          .format(percentage))
        arcpy.AddMessage("Point count statistics - done.")
        return

    def std_point_per_area(self):           # Dodac obliczenie na długosc linii
        self.report.write("\n---Standard deviation of point count per "
                          + "length and area---\n")
        with arcpy.da.SearchCursor(self.grid, ["PageNumber"]) as cursor:
            count_cells = 1
            for row in cursor:
                count_cells += 1
        self.report.write("Number of cells in grid: {0}\n"
                          .format(count_cells))
        temp_points_original = "in_memory\\points_original"
        temp_points_simplified = "in_memory\\points_simplified"
        arcpy.FeatureVerticesToPoints_management(self.original,
                                                 temp_points_original)
        arcpy.FeatureVerticesToPoints_management(self.simplified,
                                                 temp_points_simplified)
        temp_intersect_original = "in_memory\\intersect_original"
        temp_intersect_simplified = "in_memory\\intersect_simplified"
        arcpy.Intersect_analysis([temp_points_original, self.grid],
                                 temp_intersect_original, "ALL", "#", "POINT")
        arcpy.Intersect_analysis([temp_points_simplified, self.grid],
                                 temp_intersect_simplified, "ALL", "#",
                                 "POINT")
        arcpy.Delete_management(temp_points_original)
        arcpy.Delete_management(temp_points_simplified)
        temp_stats_original = "in_memory\\stats_original"
        temp_stats_simplified = "in_memory\\stats_simplified"
        arcpy.Statistics_analysis(temp_intersect_original, temp_stats_original,
                                  "Id SUM", "PageNumber")
        arcpy.Statistics_analysis(temp_intersect_simplified,
                                  temp_stats_simplified,
                                  "Id SUM", "PageNumber")
        arcpy.Delete_management(temp_intersect_original)
        arcpy.Delete_management(temp_intersect_simplified)
        with arcpy.da.SearchCursor(temp_stats_original,
                                   ["FREQUENCY"]) as cursor:
            points_original_in_cells = []
            for row_fr1 in cursor:
                points_original_in_cells.append(row_fr1[0])
        mean_original = sum(points_original_in_cells) \
                        / len(points_original_in_cells)
        points_original_in_cells[:] = [(x - mean_original)**2
                                       for x in points_original_in_cells]
        std_original = sqrt(float(sum(points_original_in_cells)) \
                       / float(len(points_original_in_cells) - 1))
        self.report.write("Number of cells in grid intersected by original "
                          + "points: {0}\n"
                          .format(len(points_original_in_cells)))
        self.report.write("Standard deviation of original points: {:0.2f} m\n"
                          .format(std_original))
        with arcpy.da.SearchCursor(temp_stats_simplified,
                                   ["FREQUENCY"]) as cursor:
            points_simplified_in_cells = []
            for row_fr2 in cursor:
                points_simplified_in_cells.append(row_fr2[0])
        mean_simplified = sum(points_simplified_in_cells) \
                        / len(points_simplified_in_cells)
        points_simplified_in_cells[:] = [(x - mean_simplified)**2
                                       for x in points_simplified_in_cells]
        std_simplified = sqrt(float(sum(points_simplified_in_cells)) \
                       / float(len(points_simplified_in_cells) - 1))
        self.report.write("Number of cells in grid intersected by simplified "
                          + "points: {0}\n"
                          .format(len(points_simplified_in_cells)))
        self.report.write("Standard deviation of simplified points: "
                          + "{:0.2f} m\n".format(std_simplified))
        percentage = (std_simplified / std_original)*100
        self.report.write("Percentage change in the standard deviation of the "
                          + "number of coordinates: {:0.2f}%\n"
                          .format(percentage))
        arcpy.Delete_management(temp_stats_original)
        arcpy.Delete_management(temp_stats_simplified)
        arcpy.AddMessage("Point count per area statistics - done.")
        return

    def surface_in_between(self):
        temp_merge = "in_memory\\merge"
        arcpy.Merge_management([self.original, self.simplified], temp_merge)
        temp_poly = "in_memory\\poly"
        arcpy.FeatureToPolygon_management(temp_merge, temp_poly)
        arcpy.Delete_management(temp_merge)
        arcpy.AddGeometryAttributes_management(temp_poly, "AREA", "#",
                                               "SQUARE_METERS")
        with arcpy.da.SearchCursor(temp_poly, ["POLY_AREA"]) as cursor:
            summarize = 1
            parts = 1
            for row in cursor:
                summarize += row[0]
                parts += 1
        mean = summarize / parts
        arcpy.Delete_management(temp_poly)
        self.report.write("\n---Differential surface---\n")
        self.report.write("Sum area: {0} sq m\n".format(summarize))
        self.report.write("Mean area: {0} sq m\n".format(mean))
        arcpy.AddMessage("Differential surface statistics - done.")
        return

    def surface_in_between_per_area(self):   # Dodac obliczenie na długosc linii i powierzchni
        temp_merge = "in_memory\\merge"
        arcpy.Merge_management([self.original, self.simplified], temp_merge)
        temp_poly = "in_memory\\poly"
        arcpy.FeatureToPolygon_management(temp_merge, temp_poly)
        arcpy.Delete_management(temp_merge)
        temp_intersect_poly = "in_memory\\intersect_poly"
        arcpy.Intersect_analysis([temp_poly, self.grid],
                                 temp_intersect_poly, "ALL", "#", "INPUT")
        arcpy.AddGeometryAttributes_management(temp_intersect_poly, "AREA",
                                               "#", "SQUARE_METERS")
        temp_stats_poly = "in_memory\\stats_simplified"
        arcpy.Statistics_analysis(temp_intersect_poly, temp_stats_poly,
                                  "POLY_AREA SUM", "PageNumber")

        temp_poly_points = "in_memory\\poly_points"
        arcpy.FeatureToPoint_management(temp_intersect_poly, temp_poly_points,
                                        "INSIDE")
        with arcpy.da.SearchCursor(temp_poly_points, ["SHAPE@", "POLY_AREA",
                                                      "PageNumber"]) as cursor:
            poly_points_array = []
            poly_point_attr_array = []
            for point in cursor:
                poly_points_array.append(point[0])
                poly_point_attr_array.append([point[1], point[2]])
        with arcpy.da.SearchCursor(self.simplified, ["SHAPE@"]) as cursor:
            side_array = []
            for polyline in cursor:
                for i_p in poly_points_array:
                    line = polyline[0]
                    check_point = line.queryPointAndDistance(i_p)
                    side_array.append(check_point[3])
        sum_right = 0
        sum_left = 0
        for elem in side_array:
            if elem is True:
                sum_right += poly_point_attr_array[elem][0]
            else:
                sum_left += poly_point_attr_array[elem][0]
        self.report.write("\n---Differential surface per length and area---\n")
        self.report.write("Sum of positive areal displacement: {0} sq m\n"
                          .format(sum_right))
        self.report.write("Sum of negative areal displacement: {0} sq m\n"
                          .format(sum_left))
        arcpy.AddMessage("Differential surface per area statistics - done.")
        return

    def hausdorff_distance(self):
        arcpy.AddMessage("Hausdorff distance statistics - done.")
        return

    def calculate_statistics(self):
        self.report.write("---------Simplification statistics---------\n")
        arcpy.AddMessage("Started calculating statistics...")
        self.point_count()
        self.std_point_per_area()
        self.surface_in_between()
        self.surface_in_between_per_area()
        self.hausdorff_distance()
        self.report.close()



if __name__ == '__main__':
    calculating_statistics = HexStatistics(arcpy.GetParameterAsText(0),
                                           arcpy.GetParameterAsText(1),
                                           arcpy.GetParameterAsText(2))
    exit(calculating_statistics)
