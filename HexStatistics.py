#######################################################################
# Filename: HexStatistics.py                                          #
# Author: Karolina Mamczarz                                           #
# Institution: AGH University of Science and Technology in Cracow,    #
#               Poland                                                #
# Faculty: Mining Surveying and Environmental Engineering             #
# Department: Integrated Geodesy and Cartography                      #
# Last update: 2017-09-11                                             #
# Version: 1.0.0                                                      #
# Description:                         #
# Class: HexStatistics                                                #
# Methods:                            #
# Result:                                      #
#######################################################################

from sys import exit
import arcpy


class HexStatistics(object):

    def __init__(self, original, simplified, report_file_path):
        self.original = original
        self.simplified = simplified
        self.report_file_path = report_file_path
        self.report = open(self.report_file_path + '.txt', 'w')
        self.report.write("---------Simplification statistics---------\n\n")

        # Methods run on init
        self.grid = self.comparison_grid_1_sqkm()
        self.point_count()
        #self.point_per_square_1_km()
        #self.surface_in_between()

    def comparison_grid_1_sqkm(self):
        temp_merge = "in_memory\\merge"
        arcpy.Merge_management([self.original, self.simplified], temp_merge)
        temp_grid = "in_memory\\grid"
        arcpy.GridIndexFeatures_cartography(
            temp_grid, temp_merge, "INTERSECTFEATURE", "NO_USEPAGEUNIT", "#",
            "1 Kilometers", "1 Kilometers")
        arcpy.Delete_management(temp_merge)
        return temp_grid

    def point_count(self):
        self.report.write("\n---Point count---\n")
        with arcpy.da.SearchCursor(self.original, ["SHAPE@"]) as cursor:
            original_pntnum = 0
            for row in cursor:
                for part in row[0]:
                    for pnt in part:
                        original_pntnum += 1
        self.report.write("Points before simplification: {0}\n"
                          .format(original_pntnum))
        with arcpy.da.SearchCursor(self.simplified, ["SHAPE@"]) as cursor:
            simplified_pntnum = 0
            for row in cursor:
                for part in row[0]:
                    for pnt in part:
                        simplified_pntnum += 1
        self.report.write("Points after simplification: {0}\n"
                          .format(simplified_pntnum))
        return

    def point_per_square_1_km(self):
        return

    def surface_in_between(self):
        temp_merge = "in_memory\\merge"
        arcpy.Merge_management([self.original, self.simplified], temp_merge)
        temp_poly = "in_memory\\poly"
        arcpy.FeatureToPolygon_management(temp_merge, temp_poly)
        arcpy.Delete_management(temp_merge)
        arcpy.AddGeometryAttributes_management(temp_poly, "AREA", "#",
                                               "METERS")
        with arcpy.da.SearchCursor(temp_poly, ["POLY_AREA"]) as cursor:
            sum = 0
            for row in cursor:
                parts = 1
                for part in row[0]:
                    sum += part
                    parts += 1
        mean = sum / parts
        arcpy.Delete_management(temp_poly)
        self.report.write("\n---Surface in between---\n")
        self.report.write("Sum area: {0}\n".format(sum))
        self.report.write("Mean area: {0}\n".format(mean))

if __name__ == '__main__':
    calculate_statistics = HexStatistics(arcpy.GetParameterAsText(0),
                                         arcpy.GetParameterAsText(1),
                                         arcpy.GetParameterAsText(2))
    exit(calculate_statistics)
