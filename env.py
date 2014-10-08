# created by: David Mangold, mangoldd
# last updated: 20140623
# purpose: hexpy environment settings
# input parameters: NA, module
# dev notes: set number precision


class Spatial():
    """Store spatial environment variables."""

    def __init__(self):
        self.precision = 10

    def round_value(self, number):
        """Round number to current precision."""

        return round(number, self.precision)