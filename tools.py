def distance(x1, y1, x2, y2):
    """Return cartesian distance between x, y coordiante pairs."""

    return ((x2 - x1)**2 + (y2 - y1)**2)**0.5


def side(area=1.0):
    """Calculate the side length for a regular hexagon from specified area.

    Side measued from any vertex to an adjacent vertex.
    """

    return ((2.0 * area) / (3.0 * 3.0**0.5))**0.5


def area(side=1.0):
    """Calculate the area of a regular hexagon from specified side length."""

    return 1.5 * 3.0**0.5 * side**2.0


def diagonal(side=1.0):
    """Calculate the diagonal length for a regular hexagon.

    Diagonal measured from any vertex to the opposing vertex.
    """
    
    return 2.0 * side


def width(side=1.0):
    """Calculate the width for a regular hexagon.

    Width measured from any edge midpoint to the opposing edge midpoint.
    """
    
    return 3.0**0.5 * side
