from pymoo.factory import get_problem, get_reference_directions
from pymoo.visualization.pcp import PCP

ref_dirs = get_reference_directions("das-dennis", 6, n_partitions=5) * [2, 4, 8, 16, 32, 64]
F = get_problem("dtlz1").pareto_front(ref_dirs)

PCP().add(F).show()
plot = PCP()
plot.set_axis_style(color="grey", alpha=0.5)
plot.add(F, color="grey", alpha=0.3)
plot.add(F[50], linewidth=5, color="red")
plot.add(F[75], linewidth=5, color="blue")
plot.show()

plot.reset()
plot.normalize_each_axis = False
plot.bounds = [[1, 1, 1, 2, 2, 5], [32, 32, 32, 32, 32, 32]]
plot.show()
