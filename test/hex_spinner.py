from hexpy.geom import *
from hexpy.render import *

class Twister():
    
    def __init__(self, group, path):

        self.plt = Plotter()        
        self.count = 0
        self.group = group
        self.path = path

    def twist(self, angle):

        group1 = self.group
        group2 = group1.copy()
        group2.cells = []
        
        print("twisting...")

        self.plot_groups([group1, group2], self.path, "plot_%s" % self.count)
        print(self.count, len(group1.cells), len(group2.cells))
        self.count += 1

        while len(group1.cells) > 0:
            this_cell = group1.cells.pop()
            group2.cells.append(this_cell)
            group1.rotate(angle, this_cell.center)
            self.plot_groups([group1, group2], self.path, "plot_%s" % self.count)
            print(self.count, len(group1.cells), len(group2.cells))
            self.count += 1

        self.group.cells = group2.cells

    def plot_groups(self, groups, path, name):

        self.plt.fig.set_size_inches(6,6)
        self.plt.fig_sub.axis([-25.0, 75.0, -25.0, 75.0])

        for grp in groups:
            self.plt.plot(grp)
        self.plt.save(path, name)
        self.plt.clear()

## MY TESTS

path = r"\\fileserv\USA\GIS\Code_Catalog\Python\Packages\hexpy\test\pics5"
angle = math.pi/6

cell = HexCell(HexPoint(0,0), 1)
group = HexGroup()
group.extent.set_limits(0,10,0,10)
group.tessalate(cell)
twister = Twister(group, path)
twister.twist(angle)
twister.twist(-angle)
twister.twist(angle)
twister.twist(-angle)
