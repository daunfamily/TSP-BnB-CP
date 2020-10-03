from branchandbound import branchandbound
from cuttingplanes import cuttingplanes

tlim = 600
instance = "att48.tsp"

branchandbound(instance, tlim)
# cuttingplanes(instance, tlim)
