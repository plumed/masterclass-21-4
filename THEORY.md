# PLUMED Masterclass 21.4: Metadynamics

## Background theory

In [metadynamics](https://www.pnas.org/doi/10.1073/pnas.202427399), an external history-dependent bias potential is constructed in the space of a few selected degrees of freedom $\vec{s}({q})$, generally called collective variables (CVs).  This potential is built as a sum of Gaussian kernels deposited along the trajectory in the CVs space:

$$
V(\vec{s},t) = \sum_{ k \tau < t} W(k \tau)
\exp\left(
-\sum_{i=1}^{d} \frac{(s_i-s_i({q}(k \tau)))^2}{2\sigma_i^2}
\right).
$$

where $\tau$ is the Gaussian deposition stride, $\sigma_i$ the width of the Gaussian for the $i$th CV, and $W(k \tau)$ is the height of the Gaussian. The effect of the metadynamics bias potential is to push the system away from local minima into visiting new regions of the phase space. Furthermore, in the long time limit, the bias potential converges to minus the free energy as a function of the CVs:

$$
V(\vec{s},t\rightarrow \infty) = -F(\vec{s}) + C.
$$

In standard metadynamics, Gaussian kernels of constant height are added for the entire course of a simulation. As a result, the system is eventually pushed to explore high free-energy regions and the estimate of the free energy calculated from the bias potential oscillates around the real value.  In [well-tempered metadynamics](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.020603), the height of the Gaussian is decreased with simulation time according to:

$$
 W (k \tau ) = W_0 \exp \left( -\frac{V(\vec{s}({q}(k \tau)),k \tau)}{k_B\Delta T} \right ),
$$

where $W_0$ is an initial Gaussian height, $\Delta T$ an input parameter with the dimension of a temperature, and $k_B$ is the Boltzmann constant.  With this rescaling of the Gaussian height, the bias potential smoothly converges in the long time limit, but it does not fully compensate the underlying free energy:

$$
V(\vec{s},t\rightarrow \infty) = -\frac{\Delta T}{T+\Delta T}F(\vec{s}) + C.
$$

where $T$ is the temperature of the system.  In the long time limit, the CVs thus sample an ensemble at a temperature $T+\Delta T$ which is higher than the system temperature $T$.  The parameter $\Delta T$ can be chosen to regulate the extent of free-energy exploration: $\Delta T = 0$ corresponds to standard MD, $\Delta T \rightarrow\infty$ to standard metadynamics. In well-tempered metadynamics literature and in PLUMED, you will often encounter the term "bias factor" which is the ratio between the temperature of the CVs ($ T+\Delta T$) and the system temperature ($T$):

$$
\gamma = \frac{T+\Delta T}{T}.
$$

The bias factor should thus be carefully chosen in order for the relevant free-energy barriers to be crossed efficiently in the time scale of the simulation.
