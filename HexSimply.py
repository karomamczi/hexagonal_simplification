# -*- coding: utf-8 -*-

############################
# Filename: HexSimply.py
# Author: eng. Karolina Mamczarz
# Institution: AGH University of Science and Technology in Cracow, Poland
# Last update: 2017-03-08
# Version: 1.0.0
# Description: #
# Class: HexSimply
# Methods: #
# Parameters: #
# Result: #
############################



import arcpy
import os
from math import sqrt, trunc, sin, cos, pi

class HexSimply(object):

    def __init__(self,original,method,l,s,simplified):
        self.original = original
        self.simplified = simplified
        self.l = l / 1000
        self.s = s
        self.method = method
        #self.createNewFeature()
        self.workspace = os.path.dirname(self.original) #POTEM DO WYRZUCENIA
        self.chooseMethod()

    def createNewFeature(self):
        path, file = os.path.split(self.simplified)
        filename, fileext = os.path.splitext(file)
        arcpy.CreateFeatureclass_management(path, filename, "POLYLINE","","","",self.original)
        container, containerext = os.path.splitext(path)
        if containerext == ".gdb" or containerext == ".mdb" or containerext == ".sde":
            file = filename
        self.fullpathname = path + "\\" + file
        return self.fullpathname

    def sort(self):
        #ewentualnie skrypt na dissolve NIEEEEEEEEEEE
        return

    def readGeom(self):
        cursor = arcpy.da.SearchCursor(self.original,["SHAPE@"])
        self.polylinecoords = []
        partnum = 0
        for row in cursor:
            for part in row[0]:
                pntnum = 0
                for pnt in part:
                    self.polylinecoords.append([partnum,pntnum,pnt.X,pnt.Y])
                    pntnum += 1
                partnum += 1
        return self.polylinecoords

    def tesseraWidth(self):
        self.r = 5*self.l*self.s
        return self.r

    def largestDiagonalHalf(self):
        self.d = round(self.r/sqrt(3),4)
        return self.d

    # The algorithm is known as "Ray Casting Method"
    # Implemented to the code based on the code from website:
    # http://geospatialpython.com/2011/01/point-in-polygon.html
    def rayCastingMethod(self):
        self.RCM = []
        for hex in self.hexcoords:
            for coords in self.polylinecoords:
                n = len(hex)
                inside = False
                p1x,p1y = hex[0]
                for i in range(n+1):
                    p2x,p2y = hex[i % n]
                    if coords[3] > min(p1y,p2y):
                        if coords[3] <= max(p1y,p2y):
                            if coords[2] <= max(p1x,p2x):
                                if p1y != p2y:
                                     xints = (coords[3]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                                if p1x == p2x or coords[2] <= xints:
                                    inside = not inside
                    p1x,p1y = p2x,p2y
                if inside == True:
                    self.RCM.append([self.hexcoords.index(hex),self.polylinecoords.index(coords),coords[2],coords[3]])
        return self.RCM

    def vertexClustering(self):
        return

    def spatialMean(self):
        return

    def eliminateSelfCrossing(self):
        return

    def statistics(self):
        return

    # direction of tesselation consistent horizontally to minimal area bounding box statring from upper-left corner
    # of this bounding box
    def boundingBox(self):
        Az = [30, 90, 150, 210, 270, 330]
        maxmin= []
        for coord in zip(*self.readGeom()):
            maxmin.append(max(coord))
            maxmin.append(min(coord))
        Xur = maxmin[4] # X of upper-right corner
        Yur = maxmin[6] # Y of upper-right corner
        Xdl = maxmin[5] # X of down-left corner
        Ydl = maxmin[7] # Y of down-left corner
        Xul = Xdl       # X of upper-left corner
        Yul = Yur       # Y of upper-left corner
        a = sqrt((Xul-Xur)**2 + (Yul-Yur)**2)
        b = sqrt((Xul-Xdl)**2 + (Yul-Ydl)**2)
        verticalCover = trunc((b-(self.tesseraWidth()/2))/self.tesseraWidth()) + 1
        horizontalCover = trunc((a-0.5*self.largestDiagonalHalf())/(1.5*self.largestDiagonalHalf())) + 2
        self.hexcoords = []
        idxhex = 0
        for i in range(horizontalCover+1):
            for j in range(verticalCover+1):
                hexcoords_temp = []
                for az in Az:
                    if (i % 2 == 0):
                        X = round((Xul + self.largestDiagonalHalf()*sin((az*pi)/180)) + 1.5*self.largestDiagonalHalf()*i,4)
                        Y = round((Yul + self.largestDiagonalHalf()*cos((az*pi)/180)) - self.tesseraWidth()*j,4)
                        hexcoords_temp.append([X,Y])
                    else:
                        X = round((Xul + self.largestDiagonalHalf()*sin((az*pi)/180)) + 1.5*self.largestDiagonalHalf()*i,4)
                        Y = round((Yul + self.largestDiagonalHalf()*cos((az*pi)/180)) - self.tesseraWidth()*j - self.tesseraWidth()/2,4)
                        hexcoords_temp.append([X,Y])
                for coords in self.polylinecoords:
                    n = len(hexcoords_temp)
                    inside = False
                    p1x,p1y = hexcoords_temp[0]
                    for i in range(n+1):
                        p2x,p2y = hexcoords_temp[i % n]
                        if coords[3] > min(p1y,p2y):
                            if coords[3] <= max(p1y,p2y):
                                if coords[2] <= max(p1x,p2x):
                                    if p1y != p2y:
                                         xints = (coords[3]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                                    if p1x == p2x or coords[2] <= xints:
                                        inside = not inside
                        p1x,p1y = p2x,p2y
                    if inside == True:
                        self.hexcoords.append([idxhex,self.polylinecoords.index(coords),coords[2],coords[3]])
                idxhex +=1
        return arcpy.AddMessage(self.hexcoords)

    # direction of tesselation consistent with the direction of the longest section of the original polyline
    # starting from the first point of this original polyline
    def directionLongestSection(self):
        return

    # direction of tesselation consistent with the direction perpendicular to the direction of bisector of an angle,
    # which vertex is simultaneously a global maximum of the original polyline
    def globalMaximum(self):
        return

    def chooseMethod(self):
        if self.method == 'FROM BOUNDING BOX':
            self.boundingBox()
        elif self.method == 'FROM LONGEST SECTION':
            self.directionLongestSection()
        elif self.method == 'FROM GLOBAL MAXIMUM':
            self.globalMaximum()


if __name__ == '__main__':
    polyline = HexSimply(arcpy.GetParameterAsText(0),arcpy.GetParameterAsText(1),arcpy.GetParameter(2),
                         arcpy.GetParameter(3),arcpy.GetParameterAsText(4))
