Author: Karolina Mamczarz

ArcToolbox: Hexagonal Simplification
Description: Implementation of line simplification methods based on hexagonal tessellation generated on bounding box and oriented rectangles. Implementation of automated tests on efficiency and quality of processed simplification.
	
Tool: Hexagonal Simplification
Description: Tool calculates simplified lines using hexagonal tessellation. You can choose,which method to use in order to orientate the tessellation in Cartesian coordinate system:
	•	bounding box (parallel to axes),
	•	minimum input data width rectangle or minimum area rectangle (oriented),
	•	furthest point of input data from line, which joins first and last point of the data(oriented).

Tool: Statistics
Description: Set of tools useful for testing the quality of simplification.Statistics cover:
	•	length difference between original and simplified line,
	•	percentage change in number of coordinates,
	•	percentage change instandard deviation of point count,
	•	average point density per length,
	•	differential surface (sum and arithmetic mean),
	•	average differential surface per length and area,
	•	hausdorff distance.
			
Look for results in /example directory.
			
Data source: vector layers http://download.geofabrik.de/

Requirements: ArcGIS for Desktop ^10.3 Advanced License 