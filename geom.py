# created by: David Mangold, mangoldd
# last updated: 20140630
# purpose: tools for working with regular hexagons
# input parameters: NA, module
# dev notes:

import os
import math
from random import uniform
from numpy import *
import env
from copy import deepcopy
import matplotlib.pyplot as plt

# ----- CLASSES -----


class Result(object):
    """Generic result object definition."""

    def __init__(self):

        pass


class ResultFactory(object):
    """Generate Result instances."""

    def make_result(self):
        """Return result instance."""

        return Result()


class Extent(object, env.Spatial):
    """Extent (minimum bounding box) for geom objects."""

    def __init__(self, hex_geom_obj=None):
        """Initialize Extent instance..."""

        env.Spatial.__init__(self)

        self.parent_geom = hex_geom_obj
        self.hexpoint_factory = HexPointFactory()
        self.xmin = 0.0
        self.xmax = 0.0
        self.ymin = 0.0
        self.ymax = 0.0
        self.xrange = 0.0
        self.yrange = 0.0
        self.vertices = []
        self.__build_extent(self.parent_geom)

    def __build_extent(self, hex_geom_obj=None):
        """Build extent instance based on x and y limits."""

        if hex_geom_obj:

            if type(hex_geom_obj) is HexCell:

                x_vals = [v.x for v in hex_geom_obj.vertices]
                y_vals = [v.y for v in hex_geom_obj.vertices]

                xmin = min(x_vals)
                xmax = max(x_vals)
                ymin = min(y_vals)
                ymax = max(y_vals)

                self.set_limits(xmin, xmax, ymin, ymax)

            elif type(hex_geom_obj) is HexGroup:

                x_vals = []
                y_vals = []

                if len(hex_geom_obj.cells) == 0:

                    x_vals.append(0.0)
                    y_vals.append(0.0)

                for cell in hex_geom_obj.cells:

                    x_vals += [v.x for v in cell.vertices]
                    y_vals += [v.y for v in cell.vertices]

                xmin = min(x_vals)
                xmax = max(x_vals)
                ymin = min(y_vals)
                ymax = max(y_vals)

                self.set_limits(xmin, xmax, ymin, ymax)

            else:

                print("%s is not valid. Use a HexCell, or HexGroup instance instead." % hex_geom_obj)
                return

        else:

            self.set_limits()

        return self

    def __build_vertices(self, xmin, xmax, ymin, ymax):
        """Build vertices from extent limits."""

        vertices = [self.hexpoint_factory.make_hexpoint(x, y)
                    for x, y in [(xmin, ymin),
                                 (xmin, ymax),
                                 (xmax, ymax),
                                 (xmax, ymin)]]

        return vertices

    def update(self):
        """Update extent based on hexpy geom object."""

        self.__build_extent(self.parent_geom)

    def set_limits(self, xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0):
        """Set extent limits: xmin, ymin, xmax, ymax."""

        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.ymin = float(ymin)
        self.ymax = float(ymax)
        self.xrange = xmax - xmin
        self.yrange = ymax - ymin

        self.vertices = self.__build_vertices(xmin, xmax, ymin, ymax)

        if self.xrange < 0.0:

            raise Warning("xrange less than zero")

        if self.yrange < 0.0:

            raise Warning("yrange less than zero")

        return self.xmin, self.xmax, self.ymin, self.ymax

    def get_limits(self):
        """Return tuple of extent limits: xmin, ymin, xmax, ymax."""

        return self.xmin, self.xmax, self.ymin, self.ymax

    def copy(self):
        """Return a deep copy of the current instance."""

        return deepcopy(self)


class ExtentFactory(object):
    """Make extent instances."""

    def make_extent(self, hex_geom_obj=None):
        """Return extent for hex_geom_obj."""

        return Extent(hex_geom_obj)


class HexPoint(object, env.Spatial):
    """Represents point in cartesian two-space.

    Used as HexCell vertex and center, and as origin for transformations.
    """

    def __init__(self, x=None, y=None):
        """Initialize HexPoint instance with specified x, y coordinates."""

        env.Spatial.__init__(self)

        self.x = self.round_value(float(x))
        self.y = self.round_value(float(y))
        self.location = array([self.x, self.y])

    def shift(self, dx=0.0, dy=0.0):
        """Update HexPoint location from dx, dy offsets."""

        self.x += self.round_value(dx)
        self.y += self.round_value(dy)
        self.location = array([self.x, self.y])

        return self

    def place(self, x=None, y=None):
        """Update HexPoint location from x, y coordinates."""

        self.x = self.round_value(x)
        self.y = self.round_value(y)
        self.location = array([self.x, self.y])

        return self

    def distance(self, hexpoint=None):
        """Return cartesian distance to other HexPoint instance."""

        loc1, loc2 = self.location, hexpoint.location
        deltas = loc2 - loc1
        distance = (deltas**2).sum()**0.5

        return self.round_value(distance)

    def copy(self):
        """Return a deep copy of the current instance."""

        return deepcopy(self)

    def save(self, path=None, name=None):
        """Save plot of HexCell vertices connected by lines."""

        filepath = os.path.join(path, name)
        plt.savefig(filepath)
        plt.clf()


class HexPointFactory(object):
    """Make HexPoint instances."""

    def make_hexpoint(self, x=None, y=None):
        """Return HexPoint instance at x, y location."""

        return HexPoint(x, y)


class HexCell(object, env.Spatial):
    """Represents a single, regular, hexagonal cell."""

    def __init__(self, center=None, side=0.0, ident='', tag=''):
        """Initialize HexCell with specified center and side length.

        Center must be a HexPoint instance; ident and tag are optional.
        """

        env.Spatial.__init__(self)

        self.hexpoint_factory = HexPointFactory()
        self.extent_factory = ExtentFactory()
        self.ident = ident
        self.tag = tag        
        self.center = None        
        self.side = 0.0
        self.diagonal = 0.0
        self.width = 0.0
        self.perimeter = 0.0        
        self.area = 0.0        
        self.vertices = []
        self.extent = None
        self.__build_cell(center, side)

    def __build_cell(self, center=None, side=0.0):
        """Build HexCell instance based on center and side length."""

        try:

            side = self.round_value(float(side))
            
            if side == 0.0:

                print "Build cell failed. Side length must be greater than zero."

                return

            if side < 0.0:

                side = self.round_value(side * -1.0)

        except ValueError:

            print "Build cell failed. Side length must be a number."

            return

        try:

            if isinstance(center, HexPoint):

                self.center = center

            else:

                print "Build cell failed. Center point must be HexPoint instance."

                return

        except NameError:

            print "Build cell failed. Center point not defined."

            return
            
        self.side = side
        self.vertices = self.__build_vertices()
        self.diagonal = self.get_diagonal()
        self.width = self.get_width()
        self.perimeter = self.get_perimeter()
        self.area = self.get_area()
        self.extent = self.extent_factory.make_extent(self)

        return self

    def __build_vertices(self):
        """Build list of HexPoint instances representing HexCell vertices."""

        if self.side > 0.0:

            angles = [math.radians(deg) for deg in range(0, 360, 60)]
            vectors = [array([math.sin(a), math.cos(a)]) for a in angles] * array([self.side])

            try:

                cent_location = self.center.location

            except AttributeError:

                print "Failed to build vertices. Center location not defined."

                return

            coords = [cent_location + v for v in vectors]
            vertices = [self.hexpoint_factory.make_hexpoint(c[0], c[1]) for c in coords]

            return vertices

        else:

            print "Failed to build vertices. Side length must be greater than zero."

            return

    def update_properties(self):
        """Update calculated properties of HexCell instance from current vertices.

        Calculated properties include: side, diagonal, width, and area
        """

        points = self.get_points()

        if len(points) == 7:

            self.side = points[0].distance(points[1])
            self.diagonal = self.get_diagonal()
            self.width = self.get_width()
            self.perimeter = self.get_perimeter()
            self.area = self.get_area()
            self.extent.update()

            return vars(self)

        else:

            print "Failed to update properties. Not all points are defined."

            return

    def shift(self, dx=0.0, dy=0.0):
        """Update HexCell location by adding dx, dy offsets."""

        points = self.get_points()

        for pnt in points:

            pnt.shift(dx, dy)

        self.update_properties()

        return self

    def place(self, x=None, y=None):
        """Update HexCell location by setting x, y coordinates."""

        new_coords = array([x, y])
        old_coords = array(self.center.location)
        deltas = new_coords - old_coords
        dx, dy = deltas[0], deltas[1]
        points = self.get_points()

        for pnt in points:

            pnt.shift(dx, dy)

        self.update_properties()

        return self

    def rotate(self, angle=0.0, origin=None):
        """Rotate HexCell by the specified angle (radians) about the origin."""

        if not origin:

            origin = self.center

        if isinstance(origin, HexPoint):

            rotate_array = array([[math.cos(angle), -math.sin(angle)],
                                  [math.sin(angle), math.cos(angle)]])            
            origin_array = origin.location
            points = self.get_points()

            for pnt in points:

                pnt_array = pnt.location
                pnt_array = pnt_array - origin_array
                pnt_array = dot(pnt_array, rotate_array)
                pnt_array = pnt_array + origin_array
                pnt.place(pnt_array[0], pnt_array[1])

            self.update_properties()

            return self

        else:

            print "Failed to rotate. Origin must be HexPoint instance."

            return

    def scale(self, factor=1.0, origin=None):
        """Scale HexCell by the specified factor centered at origin."""

        if not origin:

            origin = self.center

        if isinstance(origin, HexPoint):

            scale_array = array([[factor, 0.0],
                                 [0.0, factor]])
            origin_array = origin.location
            points = self.get_points()

            for pnt in points:

                pnt_array = pnt.location
                pnt_array = pnt_array - origin_array
                pnt_array = dot(pnt_array, scale_array)
                pnt_array = pnt_array + origin_array
                pnt.place(pnt_array[0], pnt_array[1])

            self.update_properties()

            return self

        else:

            print "Failed to scale. Origin must be HexPoint instance."

            return

        pass

    def get_points(self):
        """Return list of all HexPoint instances included in HexCell.

        Returned points list includes all vertices and center point.
        """

        return [self.center] + self.vertices

    def get_side(self):
        """Calculate side length from current area.

        Side measured from any vertex to an adjacent vertex.
        """

        return self.round_value(self.vertices[0].distance(self.vertices[1]))

    def get_diagonal(self):
        """Calculate diagonal length from current side length.

        Diagonal measured from any vertex to the opposing vertex.
        """
        
        return self.round_value(2.0 * self.get_side())

    def get_width(self):
        """Calculate width from current side length.

        Width measured from any edge midpoint to the opposing edge midpoint.
        """
        
        return self.round_value(3.0**0.5 * self.get_side())

    def get_perimeter(self):
        """Calculate perimeter from current side length."""

        return self.round_value(6.0 * self.get_side())

    def get_area(self):
        """Calculate area from current side length."""

        return self.round_value(1.5 * 3.0**0.5 * self.get_side()**2.0)

    def copy(self):
        """Return a deep copy of the current instance."""

        return deepcopy(self)


class HexCellFactory(object):
    """Make HexCell instances."""

    def __init__(self):
        """Initialize the HexCellFactory with ident -1

        ident is incrementally assigned to each new HexCell instance,
        ensuring a unique id for each.
        """

        self.ident = '-1'

    def make_hexcell(self, center=None, side=0.0, tag=''):
        """Increment ident and Return HexCell instance."""

        self.ident = str(int(self.ident) + 1)
        ident = self.ident

        return HexCell(center, side, ident, tag)


class HexGroup(object, env.Spatial):
    """Collection of HexCell instances.

    Class includes fuctions for tessalation, random distribution, and
    placement at regular intervals.
    """

    def __init__(self):
        """Initialize HexGroup instance with default values."""

        env.Spatial.__init__(self)

        self.result_factory = ResultFactory()
        self.hexcell_factory = HexCellFactory()
        self.hexpoint_factory = HexPointFactory()
        self.extent_factory = ExtentFactory()
        self.cells = []
        self.extent = self.extent_factory.make_extent(self)

    def update_properties(self):
        """Update group properties.

        HexCell properties are automatically updated for all cell instance transformations.
        """

        self.extent.update()

        return vars(self)

    def add_cells(self, cells_list=[]):
        """Add one or more cells to the group."""

        valid = 1

        for cell in cells_list:

            if not isinstance(cell, HexCell):

                valid = 0

        if valid == 1:

            for cell in cells_list:

                self.cells.append(cell)

            self.update_properties()

        else:

            print("Input list contains invalid objects. All objects must be HexCell instances.")

    def remove_cells(self, cells_list=[]):
        """Remove one or more cells from the group."""

        valid = 1

        for cell in cells_list:

            if not isinstance(cell, HexCell):

                valid = 0

        if valid == 1:

            for cell in cells_list:

                self.cells.remove(cell)

            self.update_properties()

        else:

            print("Input list contains invalid objects. All objects must be HexCell instances.")

    def clear_cells(self):
        """Remove all cells from the group."""

        self.cells = []

    def rotate(self, angle=0.0, origin=None):
        """Rotate all cells in the grid by angle(radians) about origin."""

        if not origin:

            xmean = self.extent.xmin + self.extent.xrange / 2
            ymean = self.extent.ymin + self.extent.yrange / 2
            origin = HexPoint(xmean, ymean)

        for cell in self.cells:

            cell.rotate(angle, origin)

        self.update_properties()

        return self

    def shift(self, dx=0.0, dy=0.0):
        """Update HexCell locations by adding dx, dy offsets."""

        for cell in self.cells:

            points = cell.get_points()

            for pnt in points:

                pnt.shift(dx, dy)

        self.update_properties()

        return self

    def scale(self, factor=1.0, origin=None):
        """Scale all cells in the grid by factor centered at origin."""

        if not origin:

            xmean = self.extent.xmin + self.extent.xrange / 2
            ymean = self.extent.ymin + self.extent.yrange / 2
            origin = HexPoint(xmean, ymean)

        for cell in self.cells:

            cell.scale(factor, origin)

        self.update_properties()

        return self
    
    def tessalate(self, refcell=None):
        """Return list of HexCell instances tessalated across grid extent."""

        refcell = refcell.copy()

        result = self.result_factory.make_result()
        outcells = []
        tag = refcell.tag
        side = refcell.side
        width = refcell.width
        diagonal = refcell.diagonal
        count = 0
        y_offset = diagonal - side / 2
        x_offset = width
        
        for y in arange(self.extent.ymin, self.extent.ymax + diagonal / 2, y_offset):
            
            for x in arange(self.extent.xmin, self.extent.xmax + width / 2, x_offset):

                if count % 2 == 0:

                    center = self.hexpoint_factory.make_hexpoint(x, y)

                else:

                    center = HexPoint(x + width / 2, y)

                hexcell = self.hexcell_factory.make_hexcell(center, side, tag)
                outcells.append(hexcell)

            count += 1

        self.cells += outcells
        result.cells = outcells

        self.update_properties()

        return result
                                  
    def offsets(self, refcell=None, dx=0.0, dy=0.0):
        """Return list of HexCell instances placed at specified intervals across extent.

        Cell spacing is set by dx, dy offsets.
        """

        refcell = refcell.copy()

        result = self.result_factory.make_result()
        outcells = []
        tag = refcell.tag
        side = refcell.side

        for y in arange(self.extent.ymin, self.extent.ymax, dy):

            for x in arange(self.extent.xmin, self.extent.xmax, dx):

                center = HexPoint(x, y)
                hexcell = self.hexcell_factory.make_hexcell(center, side, tag)
                outcells.append(hexcell)

        self.cells += outcells
        result.cells = outcells

        self.update_properties()

        return result

    def random(self, refcell=None, count=1.0):
        """Return list of HexCell instances placed at random locations in HexGroup.

        Each cell position is the result of randomly and independently generated x, y coordinates.
        """

        refcell = refcell.copy()

        result = self.result_factory.make_result()
        outcells = []
        tag = refcell.tag
        side = refcell.side

        while count > 0:

            y = uniform(self.extent.ymin, self.extent.ymax)
            x = uniform(self.extent.xmin, self.extent.xmax)
            center = HexPoint(x, y)
            hexcell = self.hexcell_factory.make_hexcell(center, side, tag)
            outcells.append(hexcell)
            count -= 1

        self.cells += outcells
        result.cells = outcells

        self.update_properties()

        return result

    ## IN TESTING
    ## fit in yrange returns one too few rows of cells may be result of number precision
    ## consider increasing yrange by some fraction of y offset (diagonal - side / 2)
    def fit(self, refcell=None, axis=None, count=1.0):
        """Return list of HexCell instances uniformly scaled to fit inside HexGroup extent.

        Input axis is 'x' or 'y', count is the desired cell count along axis.
        """

        refcell = refcell.copy()

        result = self.result_factory.make_result()
        outcells = []
        count = float(count)
        fromside = refcell.side

        if axis.lower() == 'y':

            toside = self.extent.yrange / (1.5 * count)

        elif axis.lower() == 'x':

            toside = self.extent.xrange / (3.0**0.5 * count)

        else:

            print("Argument toside is required.")

            return

        factor = toside / fromside
        refcell.scale(factor)
        tag = refcell.tag
        side = refcell.side
        width = refcell.width
        diagonal = refcell.diagonal
        y_offset = diagonal - side / 2
        x_offset = width
        
        for y in arange(self.extent.ymin + side, self.extent.ymax - side, y_offset):
            
            for x in arange(self.extent.xmin + width / 2, self.extent.xmax - width / 2, x_offset):

                if count % 2 == 0:

                    center = self.hexpoint_factory.make_hexpoint(x, y)                    

                else:

                    center = HexPoint(x + width / 2, y)

                hexcell = self.hexcell_factory.make_hexcell(center, side, tag)
                outcells.append(hexcell)

            count += 1

        self.cells += outcells
        result.cells = outcells

        self.update_properties()

        return result

    def copy(self):
        """Return a deep copy of the current instance."""

        return deepcopy(self)

    def write(self, workspace=None, name=None, overwrite=False):
        """Write HexGroup data out to a file."""

        if workspace and name and overwrite:

            if os.path.splitext(workspace)[-1] in ('.gdb', '.sde', '.mdb'):

                self.__arcpy_write_features(workspace, name, overwrite)

            else:

                self.__osgeo_write_features(workspace, name, overwrite)

        else:

            print("Missing arguments: workspace: %s name: %s overwrite: %s" % (workspace, name, overwrite))

    def __arcpy_write_features(self, workspace=None, name=None, overwrite=False):
        """Create fgdb feature class containing one feature for each HexCell in HexGroup using arcpy library."""

        import arcpy

        if overwrite is True:

            arcpy.env.overwriteOutput = True

        if not arcpy.Exists(workspace):

            print("Workspace: %s" % workspace)

            answer = raw_input("Workspace does not exist. Would you like to create it (Y/N)? ")

            if answer.lower() in ('y', 'yes'):

                if os.path.splitext(workspace)[-1] != '.gdb':

                    out_folder_path, out_gdb_name = os.path.split(workspace)[0], os.path.split(workspace)[1]
                    arcpy.CreateFileGDB_management(out_folder_path, out_gdb_name, 'CURRENT')
                    print("Created workspace: %s" % workspace)

                elif os.path.splitext(workspace)[-1] != '.mdb':

                    out_folder_path, out_mdb_name = os.path.split(workspace)[0], os.path.split(workspace)[1]
                    arcpy.CreatePersonalGDB_management(out_folder_path, out_mdb_name, 'CURRENT')
                    print("Created workspace: %s" % workspace)

                elif os.path.splitext(workspace)[-1] == '':

                    os.mkdir(workspace)
                    print("Created workspace: %s" % workspace)

                else:

                    print("Unable to create workspace of type %s" % os.path.splitext(workspace)[-1])

                    return

            else:

                print "Write feature class operation aborted."

                return

        template = os.path.join(get_module_dir(), 'shp_template', 'hex_template.shp')
        spatial_reference = arcpy.Describe(template).spatialReference
        geometry_type = 'POLYGON'
        has_m = 'DISABLED'
        has_z = 'DISABLED'
        out_fc = arcpy.CreateFeatureclass_management(workspace, name, geometry_type,
                                                     template, has_m, has_z, spatial_reference)
        icur = arcpy.InsertCursor(out_fc)

        for cell in self.cells:

            pointarray = arcpy.Array([arcpy.Point(v.location[0], v.location[1]) for v in cell.vertices])
            geom = arcpy.Polygon(pointarray)
            row = icur.newRow()
            row.setValue('SHAPE', geom)
            row.setValue('HexId', cell.ident)
            row.setValue('Tag', cell.tag)
            icur.insertRow(row)

        if row:

            del row

        del icur

        arcpy.env.overwriteOutput = False

        return out_fc

    # DEPRECATED
    def write_shapefile(self, out_dir=None, out_name=None, overwrite=False):
        """Create shapefile containing one feature for each HexCell in HexGroup using arcpy library.

        out_dir must be a filesystem folder and out_name must include .shp extension.
        """

        import arcpy

        if overwrite is True:

            arcpy.env.overwriteOutput = True

        if not arcpy.Exists(out_dir):

            print "Path: %s" % out_dir
            answer = raw_input("Path does not exist. Would you like to create it (Y/N)? ")

            if answer.lower() in ('y', 'yes'):

                os.mkdir(out_dir)
                print "Created path: %s" % out_dir

            else:

                print "Write shapefile operation aborted."

                return

        template = os.path.join(get_module_dir(), 'shp_template', 'hex_template.shp')
        spatial_reference = arcpy.Describe(template).spatialReference
        geometry_type = 'POLYGON'
        has_m = 'DISABLED'
        has_z = 'DISABLED'
        shp = arcpy.CreateFeatureclass_management(out_dir, out_name, geometry_type,
                                                  template, has_m, has_z, spatial_reference)
        icur = arcpy.InsertCursor(shp)

        for cell in self.cells:

            pointarray = arcpy.Array([arcpy.Point(v.location[0], v.location[1]) for v in cell.vertices])
            geom = arcpy.Polygon(pointarray)
            row = icur.newRow()
            row.setValue('Shape', geom)
            row.setValue('HexId', cell.ident)
            row.setValue('Tag', cell.tag)
            icur.insertRow(row)

        if row:

            del row

        del icur

        arcpy.env.overwriteOutput = False

        return shp

    def __osgeo_write_features(self, workspace=None, name=None, overwrite=False):
        """Create shapefile containing one feature for each HexCell in HexGroup using osgeo library.

        out_dir must be a filesystem folder and out_name must include .shp extension.
        """

        from osgeo import ogr, osr

        if not os.path.isdir(workspace):

            print "Workspace: %s" % workspace
            answer = raw_input("Workspace does not exist. Would you like to create it (Y/N)? ")

            if answer.lower() in ('y', 'yes'):

                os.mkdir(workspace)
                print "Created workspace: %s" % workspace

            else:

                print "Write shapefile operation aborted."

                return

        driver = ogr.GetDriverByName('ESRI Shapefile')
        temp_path = os.path.join(get_module_dir(), 'shp_template', 'hex_template.shp')
        temp_data = driver.Open(temp_path)
        temp_layer = temp_data.GetLayer()        
        temp_sr = temp_layer.GetSpatialRef()
        temp_sr_wkt = temp_sr.ExportToWkt()
        temp_defn = temp_layer.GetLayerDefn()
        temp_field_count = temp_defn.GetFieldCount()
        spatial_reference = osr.SpatialReference()
        spatial_reference.ImportFromWkt(temp_sr_wkt)
        shp = os.path.join(workspace, name)

        if overwrite is True and os.path.exists(shp):

            driver.DeleteDataSource(shp)

        data = driver.CreateDataSource(shp)
        layer = data.CreateLayer('hexcells', spatial_reference, ogr.wkbPolygon)        

        for i in range(temp_field_count):
            
            temp_field = temp_defn.GetFieldDefn(i)
            field = ogr.FieldDefn(temp_field.GetName(), temp_field.GetType())            
            field.SetWidth(temp_field.GetWidth())            
            field.SetPrecision(temp_field.GetPrecision())
            layer.CreateField(field)

        layer_defn = layer.GetLayerDefn()
        this_fid = 0

        for cell in self.cells:

            ring = ogr.Geometry(ogr.wkbLinearRing)

            for vert in cell.vertices:

                ring.AddPoint(vert.location[0], vert.location[1])

            ring.AddPoint(cell.vertices[0].location[0], cell.vertices[0].location[1])
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            feature = ogr.Feature(layer_defn)
            feature.SetGeometry(poly)
            feature.SetFID(this_fid)
            feature.SetField('HexId', cell.ident)
            feature.SetField('Tag', cell.tag)
            this_fid += 1
            layer.CreateFeature(feature)

        data.Destroy()


def get_module_dir():

    return os.path.split(__file__)[0]
