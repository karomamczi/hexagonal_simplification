# Hexagonal Simplification ArcToolbox

## Description
Implementation of line simplification methods based on hexagonal tessellation generated on bounding box ([inspired by Raposo's algorithm](https://github.com/paulojraposo/HexQuant)) and oriented rectangles (author's original modification). Implementation of automated tests on efficiency and quality of processed simplification.

## Tools
### Hexagonal Simplification
Tool calculates simplified lines using hexagonal tessellation. You can choose,which method to use in order to orientate the tessellation in Cartesian coordinate system:
* bounding box (parallel to axes),
* minimum input data width rectangle or minimum area rectangle (oriented),
* furthest point of input data from line, which joins first and last point of the data(oriented).
### Statistics
A set of tools useful for testing the quality of simplification.Statistics cover:
* length difference between original and simplified line,
* percentage change in number of coordinates,
* percentage change instandard deviation of point count,
* average point density per length,
* differential surface (sum and arithmetic mean),
* average differential surface per length and area,
* hausdorff distance.

## Data source
Vector layers are from [Geofabrik](http://download.geofabrik.de/).

## Exemplary results
Look for results in [/example](https://github.com/karomamczi/hexagonal_simplification/tree/master/example) directory.


## Requirements
ArcGIS for Desktop 10.x Advanced License

## Author
Karolina Mamczarz
