
# Stability and linear response {#linresp_background}

The core of the harmonic balance method is expressing the system&#39;s behaviour in terms of Fourier components or _harmonics_. For an $N$-coordinate system, we _choose_ a set of $M_i$ harmonics to describe each coordinate $x_i$ : 

$$\begin{equation}
x_i(t) = \sum_{j=1}^{M_i} u_{i,j}  (T)  \cos(\omega_{i,j} t)+ v_{i,j} (T) \sin(\omega_{i,j} t) \:,
\end{equation}$$

This means the system is now described using a discrete set of variables $u_{i,j}$ and $v_{i,j}$. Constructing the vector 

$$\begin{equation}
\mathbf{u}(T) = (u_{1,1}(T), v_{1,1}(T), \ldots u_{N,M_N}(T), v_{N, M_N}(T))\,,
\end{equation}$$

we may obtain the _harmonic equations_ (see [an example of this procedure](/background/harmonic_balance#Duffing_harmeq))

$$\begin{equation}
\frac{d\mathbf{u}(T)}{dT}  = \bar{\mathbf{F}} (\mathbf{u})
\end{equation}$$

where $\bar{\mathbf{F}}(\mathbf{u})$ is a nonlinear function. A steady state $\mathbf{u}_0$ is defined by $\bar{\mathbf{F}}(\mathbf{u}_0) = 0$.

### Stability

Let us assume that we found a steady state $\mathbf{u}_0$. When the system is in this state, it responds to small perturbations either by returning to $\mathbf{u}_0$ over some characteristic timescale (_stable state_) or by evolving away from $\mathbf{u}_0$ (_unstable state_). To analyze the stability of $\mathbf{u}_0$, we linearize the equations of motion around $\mathbf{u}_0$ for a small perturbation $\delta \mathbf{u} = \mathbf{u} - \mathbf{u}_0$ to obtain

$$\begin{equation}
\frac{d}{dT} \left[\delta \mathbf{u}(T)\right] =  J(\mathbf{u}_0) \delta \mathbf{u}(T) \,,
\end{equation}$$

where $J(\mathbf{u}_0)=\nabla_{\mathbf{u}}  \bar{\mathbf{F}}|_{\mathbf{u}=\mathbf{u}_0}$ is the _Jacobian matrix_ of the system evaluated at $\mathbf{u}=\mathbf{u}_0$.

The linearised system is exactly solvable for $\delta \mathbf{u}(T)$ given an initial condition $\delta \mathbf{u}(T_0)$. The solution can be expanded in terms of the complex eigenvalues $\lambda_r$ and eigenvectors $\mathbf{v}_r$ of $J(\mathbf{u}_0)$, namely

$$\begin{equation} 
    \delta \mathbf{u}(T) = \sum_{r} c_r \hspace{1mm}\mathbf{v}_r e^{\lambda_r T}.
\end{equation}$$

The dynamical behaviour near the steady states is thus governed by $e^{ \lambda_r T}$: if $\mathrm{Re}(\lambda_r)<0$ for all $\lambda_r$, the state $\mathbf{u}_0$ is stable. Conversely, if $\mathrm{Re}(\lambda_r)>0$ for at least one $\lambda_r$, the state is unstable - perturbations such as noise or a small applied drive will force the system away from $\mathbf{u}_0$. 

### Linear response {#Linear-response}

The response of a stable steady state to an additional oscillatory force, caused by weak probes or noise, is often of interest. It can be calculated by solving for the perturbation $\delta \mathbf{u}(T)$ in the presence of an additional drive term.

$$\begin{equation} 
\frac{d}{dT} \left[\delta \mathbf{u}(T)\right] =  J(\mathbf{u}_0) \delta \mathbf{u}(T) + \boldsymbol{\xi} \,e^{i \Omega T}\,,
\end{equation}$$

Suppose we have found an eigenvector of $J(\mathbf{u}_0)$ such that $J(\mathbf{u}) \mathbf{v} = \lambda \mathbf{v}$. To solve the linearised equations of motion, we insert $\delta \mathbf{u}(T) = A(\Omega)\, \mathbf{v} e^{i \Omega T}$. Projecting each side onto $\mathbf{v}$ gives

$$A(\Omega) \left( i \Omega  - \lambda \right)  = \boldsymbol{\xi} \cdot \mathbf{v} \quad \implies \quad A(\Omega) = \frac{\boldsymbol{\xi} \cdot \mathbf{v} }{-\text{Re}[\lambda] + i \left( \Omega - \text{Im}[\lambda] \right)}$$

We see that each eigenvalue $\lambda$ results in a linear response that is a Lorentzian centered at $\Omega = \text{Im}[\lambda]$. Effectively, the linear response matches that of a harmonic oscillator with resonance frequency $\text{Im}[\lambda]$ and damping $\text{Re}[\lambda]$.

Knowing the response of the harmonic variables $\mathbf{u}(T)$, what is the corresponding behaviour of the &quot;natural&quot; variables $x_i(t)$? To find this out, we insert the perturbation back into the harmonic ansatz. Since we require real variables, let us use $\delta \mathbf{u}(T) = A(\Omega) \left( \mathbf{v} \, e^{i \Omega T} +   \mathbf{v}^* \, e^{-i \Omega T} \right)$. Plugging this into

$$\begin{equation}
\delta x_{i}(t) = \sum_{j=1}^{M_i} \delta{u}_{i,j}(t) \cos(\omega_{i,j} \,t) + \delta v_{i,j} (t) \sin(\omega_{i,j} \,t) 
\end{equation}$$

and multiplying out the sines and cosines gives

$$\begin{align}
\delta x_i(t) = \sum_{j=1}^{M_i} \bigg\{ \left( \text{Re}[\delta u_{i,j}] - \text{Im}[\delta v_{i,j}] \right) \,\cos[(\omega_{i,j} - \Omega) t]  \\
+ \left( \text{Im}[\delta u_{i,j}] + \text{Re}[\delta v_{i,j}] \right) \,\sin[(\omega_{i,j} - \Omega) t] \nonumber \\
+ \left( \text{Re}[\delta u_{i,j}] + \text{Im}[\delta v_{i,j}] \right) \,\cos[(\omega_{i,j} + \Omega) t] \nonumber \\
+ \left( -\text{Im}[\delta u_{i,j}] + \text{Re}[\delta v_{i,j}] \right) \,\sin[(\omega_{i,j} + \Omega) t] \bigg\} \nonumber
\end{align}$$

where $\delta u_{i,j}$ and $\delta v_{i,j}$ are the components of $\delta \mathbf{u}$ corresponding to the respective harmonics $\omega_{i,j}$. 

We see that a motion of the harmonic variables at frequency $\Omega$ appears as motion of $\delta x_i(t)$ at frequencies $\omega_{i,j}\pm \Omega$. 

To make sense of this, we normalize the vector $\delta \mathbf{u}$ and use normalised components $\delta \hat{u}_{i,j}$ and $\delta \hat{v}_{i,j}$. We also define the Lorentzian distribution

$$\begin{equation}
L(x)_{x_0, \gamma} = \frac{1}{\left( x - x_0 \right)^2 + \gamma^2 }
\end{equation}$$

We see that all components of $\delta x_i(t)$ are proportional to $L(\Omega)_{\text{Im}[\lambda], \text{Re}[\lambda]}$. The first and last two summands are Lorentzians centered at $\pm \Omega$ which oscillate at $\omega_{i,j} \pm \Omega$, respectively. From this, we can extract the linear response function in Fourier space, $\chi (\tilde{\omega})$

$$\begin{multline}
| \chi [\delta x _i](\tilde{\omega}) |^2 = \sum_{j=1}^{M_i} \bigg\{  \left[ \left( \text{Re}[\delta \hat{u}_{i,j}] - \text{Im}[\delta \hat{v}_{i,j}] \right)^2 + \left( \text{Im}[\delta \hat{u}_{i,j}] + \text{Re}[\delta \hat{v}_{i,j}] \right)^2 \right] L(\omega_{i,j} - \tilde{\omega})_{\text{Im}[\lambda], \text{Re}[\lambda]} \\
+ \left[ \left( \text{Re}[\delta \hat{u}_{i,j}] + \text{Im}[\delta \hat{v}_{i,j}] \right)^2 + \left(  \text{Re}[\delta \hat{v}_{i,j}] - \text{Im}[\delta \hat{u}_{i,j}] \right)^2 \right] L(\tilde{\omega} - \omega_{i,j})_{\text{Im}[\lambda], \text{Re}[\lambda]} \bigg\}
\end{multline}$$

Keeping in mind that $L(x)_{x_0, \gamma} = L(x + \Delta)_{x_0 + \Delta, \gamma}$ and the normalization $\delta \hat{u}_{i,j}^2 + \delta \hat{v}_{i,j}^2 = 1$, we can rewrite this as

$$\begin{equation}
|\chi [\delta x _i](\tilde{\omega})|^2 = \sum_{j=1}^{M_i} \left( 1 + \alpha_{i,j} \right) L(\tilde{\omega})_{\omega_{i,j} - \text{Im}[\lambda], \text{Re}[\lambda]}
+ \left( 1 - \alpha_{i,j} \right) L(\tilde{\omega})_{\omega_{i,j} + \text{Im}[\lambda], \text{Re}[\lambda]}
\end{equation}$$

where 

$$\alpha_{i,j} = 2\left( \text{Im}[\delta \hat{u}_{i,j}] \text{Re}[\delta \hat{v}_{i,j}] - \text{Re}[\delta \hat{u}_{i,j}] \text{Im}[\delta \hat{v}_{i,j}] \right)$$

The above solution applies to every eigenvalue $\lambda$ of the Jacobian. It is now clear that the linear response function $\chi [\delta x _i](\tilde{\omega})$ contains for each eigenvalue $\lambda_r$ and harmonic $\omega_{i,j}$ : 
- A Lorentzian centered at $\omega_{i,j}-\text{Im}[\lambda_r]$ with amplitude $1 + \alpha_{i,j}^{(r)}$ 
  
- A Lorentzian centered at $\omega_{i,j}+\text{Im}[\lambda_r]$ with amplitude $1 - \alpha_{i,j}^{(r)}$ 
  

_Sidenote:_ As $J$ a real matrix, there is an eigenvalue $\lambda_r^*$ for each $\lambda_r$. The maximum number of peaks in the linear response is thus equal to the dimensionality of $\mathbf{u}(T)$.

The linear response of the system in the state $\mathbf{u}_0$ is thus fully specified by the complex eigenvalues and eigenvectors of $J(\mathbf{u}_0)$. In HarmonicBalance.jl, the module [LinearResponse](/manual/linear_response#linresp_man) creates a set of plottable [`Lorentzian`](/manual/linear_response#HarmonicBalance.LinearResponse.Lorentzian-manual-linear_response) objects to represent this.

[Check out this example](/tutorials/linear_response#linresp_ex) of the linear response module of HarmonicBalance.jl
