# [Stability and linear response](@id linresp_background)

The core of the harmonic balance method is expressing the system's behaviour in terms of Fourier components or _harmonics_. For an $N$-coordinate system, we _choose_ a set of $M_i$ harmonics to describe each coordinate $x_i$ : 
```math
x_i(t) = \sum_{j=1}^{M_i} u_{i,j}  (T)  \cos(\omega_{i,j} t)+ v_{i,j} (T) \sin(\omega_{i,j} t) \:,
```
This means the system is now described using a discrete set of variables $u_{i,j}$ and $v_{i,j}$. Constructing a vector $\mathbf{u}(T) = (u_{1,1}(T), v_{1,1}(T), \ldots u_{N,M_N}(T), v_{N, M_N}(T))$, we may obtain the _harmonic equations_ (see [an example of this procedure](@ref Duffing_harmeq))
```math
\begin{equation} \label{eq:harmeq}
\frac{d\mathbf{u}(T)}{dT}  = \bar{\mathbf{F}} (\mathbf{u})
\end{equation}
```
where $\bar{\mathbf{F}}(\mathbf{u})$ is a nonlinear function. A steady state $\mathbf{u}_0$ is defined by $\bar{\mathbf{F}}(\mathbf{u}_0) = 0$.

### Stability

Let us assume that we found a steady state $\mathbf{u}_0$. When the system is in this state, it responds to small perturbations either by returning to $\mathbf{u}_0$ over some characteristic timescale (_stable state_) or by evolving away from $\mathbf{u}_0$ (_unstable state_). To analyze the stability of $\mathbf{u}_0$, we linearize Eq. \eqref{eq:harmeq} around $\mathbf{u}_0$ for a small perturbation $\delta \mathbf{u} = \mathbf{u} - \mathbf{u}_0$ to obtain
```math
\begin{equation} \label{eq:Jaceq}
\frac{d}{dT} \left[\delta \mathbf{u}(T)\right] =  J(\mathbf{u}_0) \delta \mathbf{u}(T) \,,
\end{equation}
```
where $J(\mathbf{u}_0)=\nabla_{\mathbf{u}}  \bar{\mathbf{F}}|_{\mathbf{u}=\mathbf{u}_0}$ is the _Jacobian matrix_ of the system evaluated at $\mathbf{u}=\mathbf{u}_0$.

Eq. \eqref{eq:Jaceq} is exactly solvable for $\delta \mathbf{u}(T)$ given an initial condition $\delta \mathbf{u}(T_0)$. The solution can be expanded in terms of the complex eigenvalues $\lambda_r$ and eigenvectors $\mathbf{v}_r$ of $J(\mathbf{u}_0)$, namely

```math
\begin{equation} \label{eq:fluct_evo}
    \delta \mathbf{u}(T) = \sum_{r}(\mathbf{v}_r\cdot \delta\mathbf{u}(T_0))\hspace{1mm}\mathbf{v}_r e^{\lambda_r T}.
\end{equation}
```

The dynamical behaviour near the steady states is thus governed by $e^{ \lambda_r T}$: if $\mathrm{Re}(\lambda_r)<0$ for all $\lambda_r$, the state $\mathbf{u}_0$ is stable. Conversely, if $\mathrm{Re}(\lambda_r)>0$ for at least one $\lambda_r$, the state is unstable - perturbations such as noise or a small applied drive will force the system away from $\mathbf{u}_0$. 


### Linear response (WIP)

The response of a stable steady state to an additional oscillatory force, caused by weak probes or noise, is often of interest. It can be calculated by solving for the perturbation $\delta \mathbf{u}$ in the presence of an additional drive term.

The linear response of the system in the state $\mathbf{u}_0$ is thus encoded in the complex eigenvalues and eigenvectors of $J(\mathbf{u}_0)$. 

[Check out this example](@ref linresp_ex) of the linear response module of HarmonicBalance.jl