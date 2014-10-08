geometry.py
1.0
======================

Overview
--------

Classes for building regular hexagons. Classes exist for defining points, hexagons, and hexagon grids (in this case "grid" implies only a collection of hexagons, not necessarily regular arrangement). Functions exist for manipulating geometries of all classes and for building various arrangements of hexagons within a HexGrid, including: tessalation, regular interval placement, random placement, and fit based on desired HexCell count in x and y dimensions.


New Features
------------

- All. Version 1.0


Deprecated Features
-------------------

- None. Version 1.0


Questions and feedback
----------------------

David Mangold, mangoldd
GIS Analyst
Clean Water Services
Phone: 210.287.8018
Email: mangoldd@cleanwaterservices.org


Required Software
-----------------

ArcGIS 10.0 or later (only required for writing out datasets)


Supported Systems
-----------------

Tested on Windows 7 Enterprise SP1


Example Usage
-------------

from geometry import *

in_fc = r'\\fileserv\USA\GIS\Code_Catalog\Python\PyPackages\hextools\test_input\test_in.shp'

out_ws = r'\\fileserv\USA\GIS\Code_Catalog\Python\PyPackages\hextools\test_output\test_out.gdb'

out_name = 'test_out'

refcell = HexCell(center=HexPoint(0,0), side=50.0, tag='TEST')

grid = HexGrid()

scur = arcpy.da.SearchCursor(in_fc, 'SHAPE@')

count = 0

for row in scur:

    print "processing feature: %s" % count

    geom = row[0]

    ext = geom.extent

    grid.set_extent(ext.XMin, ext.YMin, ext.XMax, ext.YMax)

    grid.tessalate(refcell)

    print "  %s cells" % len(grid.current_cells)

    print ""

    count += 1

del row, scur

grid.write_featclass(out_ws, out_name)