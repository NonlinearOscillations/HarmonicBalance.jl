
# Limit cycles {#limit_cycles_bg}

We explain how HarmonicBalance.jl uses a new technique to find limit cycles in systems of nonlinear ODEs. For a more in depth overview see Chapter 6 in [Jan Košata&#39;s PhD theses](https://www.doi.org/10.3929/ethz-b-000589190) or [del_Pino_2024](https://www.doi.org/10.1103/PhysRevResearch.6.03318). 

## Limit cycles from a Hopf bifurcation {#Limit-cycles-from-a-Hopf-bifurcation}

The end product of the [harmonic balance technique](/background/harmonic_balance#intro_hb) are what we call the harmonic equations, i.e., first-order ODEs for the harmonic variables $\mathbf{U}(T)$:

$$\frac{d \mathbf{U}(T)}{d T}=\overline{\mathbf{G}}(\mathbf{U})$$

These Odes have no explicit time-dependence - they are autonomous. We have mostly been searching for steady states, which likewise show no time dependence. However, time-dependent solutions to autonomous ODEs can also exist. One mechanism for their creation is a [Hopf bifurcation](https://en.wikipedia.org/wiki/Hopf_bifurcation) - a critical point where a stable solution transitions into an unstable one. For a stable solution, the associated eigenvalues $\lambda$ of the linearisation all satisfy $\operatorname{Re}(\lambda)<0$. When a Hopf bifurcation takes place, one complex-conjugate pair of eigenvalues crosses the real axis such that $\operatorname{Re}(\lambda)>0$. The state is then, strictly speaking, unstable. However, instead of evolving into another steady state, the system may assume a periodic orbit in phase space, giving a solution of the form

$$\mathbf{U}(T)=\mathbf{U}_0+\mathbf{U}_{\mathrm{lc}} \cos \left(\omega_{\mathrm{lc}} T+\phi\right)$$

which is an example of a limit cycle. We denote the originating steady state as Hopf-unstable.

We can continue to use harmonic balance as the solution still describes a harmonic response [Allwright (1977)](https://www.doi.org/10.1017/S0305004100054128). If we translate back to the the lab frame [variable $x(t)$], clearly, each frequency $\omega_j$ constituting our harmonic ansatz [$\mathbf{U}(T)$], we obtain frequencies $\omega_j$ as well as $\omega_j \pm \omega_{\text {lc }}$ in the lab frame. Furthermore, as multiple harmonics now co-exist in the system, frequency conversion may take place, spawning further pairs $\omega_j \pm k \omega_{\text {lc }}$ with integer $k$. Therefore, to construct a harmonic ansatz capturing limit cycles, we simply add an integer number $K$ of such pairs to our existing set of $M$ harmonics,

$$\left\{\omega_1, \ldots, \omega_M\right\} \rightarrow\left\{\omega_1, \omega_1 \pm \omega_{\mathrm{lc}}, \omega_1 \pm 2 \omega_{\mathrm{lc}}, \ldots, \omega_M \pm K \omega_{\mathrm{lc}}\right\}$$

## Ansatz

### Original ansatz {#Original-ansatz}

Having seen how limit cycles are formed, we now proceed to tackle a key problem: how to find their frequency $\omega_{\mathrm{lc}}$. We again demonstrate by considering a single variable $x(t)$. We may try the simplest ansatz for a system driven at frequency $\omega$,

$$x(t)=u_1(T) \cos (\omega t)+v_1(T) \sin (\omega t)$$

In this formulation, limit cycles may be obtained by solving the resulting harmonic equations with a Runge-Kutta type solver to obtain the time evolution of $u_1(T)$ and $v_1(T)$. See the [limit cycle tutorial](/tutorials/limit_cycles#limit_cycles) for an example.

### Extended ansatz {#Extended-ansatz}

Including newly-emergent pairs of harmonics is in principle straightforward. Suppose a limit cycle has formed in our system with a frequency $\omega_{\mathrm{lc}}$, prompting the ansatz

$$\begin{aligned}
& x(t)=u_1 \cos (\omega t)+v_1 \sin (\omega t) \\
& +\quad u_2 \cos \left[\left(\omega+\omega_{\mathrm{lc}}\right) t\right]+v_2 \sin \left[\left(\omega+\omega_{\mathrm{lc}}\right) t\right] \\
& \quad+u_3 \cos \left[\left(\omega-\omega_{\mathrm{lc}}\right) t\right]+v_3 \sin \left[\left(\omega-\omega_{\mathrm{lc}}\right) t\right]+\ldots
\end{aligned}$$

where each of the $\omega \pm k \omega_{\text {lc }}$ pairs contributes 4 harmonic variables. The limit cycle frequency $\omega_{\mathrm{lc}}$ is also a variable in this formulation, but does not contribute a harmonic equation, since $d \omega_{\mathrm{lc}} / d T=0$ by construction. We thus arrive at a total of $2+4 K$ harmonic equations in $2+4 K+1$ variables. To obtain steady states, we must thus solve an underdetermined system, which has an infinite number of solutions. Given that we expect the limit cycles to possess $U(1)$ gauge freedom, this is a sensible observation. We may still use iterative numerical procedures such as the Newton method to find solutions one by one, but homotopy continuation is not applicable. In this formulation, steady staes states are characterised by zero entries for $u_2, v_2, \ldots u_{2 K+1}, v_{2 K+1}$. The variable $\omega_{\text {lc }}$ is redundant and may take any value - the states therefore also appear infinitely degenerate, which, however, has no physical grounds. Oppositely, solutions may appear for which some of the limit cycle variables $u_2, v_2, \ldots u_{2 K+1}, v_{2 K+1}$ are nonzero, but $\omega_{\text {lc }}=0$. These violate our assumption of distinct harmonic variables corresponding to distinct frequencies and are therefore discarded.

### Gauge fixing {#gauge_fixing}

We now constrain the system to remove the $U(1)$ gauge freedom. This is best done by explicitly writing out the free phase. Recall that our solution must be symmetric under a time translation symmetry, that is, taking $t \rightarrow t+2 \pi / \omega$. Applying this $n$ times transforms $x(t)$ into

$$\begin{aligned}
x(t)=u_1 \cos (\omega t)+ & v_1 \sin (\omega t) \\
& \quad+u_2 \cos \left[\left(\omega+\omega_{\mathrm{lc}}\right) t+\phi\right]+v_2 \sin \left[\left(\omega+\omega_{\mathrm{lc}}\right) t+\phi\right] \\
& +u_3 \cos \left[\left(\omega-\omega_{\mathrm{lc}}\right) t-\phi\right]+v_3 \sin \left[\left(\omega-\omega_{\mathrm{lc}}\right) t-\phi\right]+\ldots
\end{aligned}$$

where we defined $\phi=2 \pi n \omega_{\text {lc }} / \omega$. Since $\phi$ is free, we can fix it to, for example,

$$\phi=-\arctan u_2 / v_2$$

which turns into

$$\begin{gathered}
x(t)=u_1 \cos (\omega t)+v_1 \sin (\omega t)+\left(v_2 \cos \phi-u_2 \sin \phi\right) \sin \left[\left(\omega+\omega_{\mathrm{lc}}\right) t\right] \\
+\left(u_3 \cos \phi-v_3 \sin \phi\right) \cos \left[\left(\omega-\omega_{\mathrm{lc}}\right) t\right]+\left(v_3 \cos \phi+u_3 \sin \phi\right)\left[\left(\omega-\omega_{\mathrm{lc}}\right) t\right]+\ldots
\end{gathered}$$

We see that fixing the free phase has effectively removed one of the variables, since $\cos \left[\left(\omega+\omega_{\text {lc }}\right) t\right]$ does not appear any more. Discarding $u_2$, we can therefore use $2+4 K$ variables as our harmonic ansatz, i.e.,

$$\mathbf{U}=\left(\begin{array}{c}
u_1 \\
v_1 \\
v_2 \\
\vdots \\
v_{2 K+1} \\
\omega_{\mathrm{lc}}
\end{array}\right)$$

to remove the infinite degeneracy. Note that $\phi$ is only defined modulo $\pi$, but its effect on the harmonic variables is not. Choosing $\phi=-\arctan u_2 / v_2+\pi$ would invert the signs of $v_2, u_3, v_3$. As a result, each solution is doubly degenerate. Combined with the sign ambiguity of $\omega_{\text {lc }}$, we conclude that under the new ansatz, a limit cycle solution appears as a fourfold-degenerate steady state.

The harmonic equations can now be solved using homotopy continuation to obtain all steady states. Compared to the single-harmonic ansatz however, we have significantly enlarged the polynomial system to be solved. As the number of solutions scales exponentially ([Bézout bound](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_theorem)), we expect vast numbers of solutions even for fairly small systems.
