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
        arcpy.GridIndexFeatures_cartography(
            temp_grid, temp_merge, "INTERSECTFEATURE", "NO_USEPAGEUNIT", "#",
            "1 Kilometers", "1 Kilometers")
        arcpy.Delete_management(temp_merge)
        arcpy.AddMessage("Comparison grid (1 sq km) prepared.")
        return temp_grid

    def point_count(self):
        self.report.write("\n---Point count---\n")
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
        arcpy.AddMessage("Point count statistics - done.")
        return

    def std_point_per_area(self):
        self.report.write("\n---Standard deviation of point count per "
                          + "area---\n")
        with arcpy.da.SearchCursor(self.grid, ["PageNumber"]) as cursor:
            count_cells = 1
            for row in cursor:
                count_cells +=1
        self.report.write("Number of cells in grid: {0}\n"
                          .format(count_cells))

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
        self.report.write("Sum area [square m]: {0}\n".format(summarize))
        self.report.write("Mean area [square m]: {0}\n".format(mean))
        arcpy.AddMessage("Differential surface statistics - done.")
        return

    def surface_in_between_per_area(self):
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
