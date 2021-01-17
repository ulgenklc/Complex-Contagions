# Complex-Contagions

Complex Contagions is a python library where one can run a simple threshold model(Watts Threshold model) on noisy geometric networks. Then, one can run persistent homology on the distance matrix obtained from the first activation matrix, T(i-->j), with the `gudhi` python package and display the signal propagation patterns as a function of a topological invariant. [Documentation](https://complex-contagions.readthedocs.io/en/latest/index.html) is available!!

The complexity of the system may be increased by varying the step function
![alt text](https://github.com/ulgenklc/Complex-Contagions/blob/master/Excitations_python/figures/figg.png)

which simply corresponds to modifying the hill function and introducing stochasticity into the system. Also, by adding k-regular edge noise is another way of introducing stochasticty which is handled by the API (see the documentation) since the spatial distance of nodes are also reflected from the path distance between nodes in geometric networks.

Complex Contagions are deprecated by [Neuronal Cascades](https://github.com/ulgenklc/Neuronal_Cascades).


