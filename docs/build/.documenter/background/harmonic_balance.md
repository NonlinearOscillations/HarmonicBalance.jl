
# The method of harmonic balance {#intro_hb}

HarmonicBalance.jl focuses on systems with equations of motion that are subject to time-dependent terms of harmonic type.

## Harmonic generation in oscillating nonlinear systems {#prelude}

Let us take a general nonlinear system of $N$ second-order ODEs with real variables $x_i(t)$, $i = 1,2,\cdots,N$, and time $t$ as the independent variable,

$$\begin{equation}
  \ddot{\mathbf{x}}(t)+ \mathbf{F}(\mathbf{x}(t), t)=0\:.
\end{equation}$$

The vector $\mathbf{x}(t) = (x_1(t), ..., x_N(t))^{\text T}$ fully describes the state of the system.  In physics, the coordinates $\mathbf{x}(t)$ encode the amplitudes of oscillators, e.g., mechanical resonators, voltage oscillations in RLC circuits, oscillating electrical dipole moments, or standing modes of an optical cavity. 

As an example, let us first solve the harmonic oscillator in frequency space. The equation of motion is

$$\begin{equation}
\ddot{x}(t) + \gamma \dot{x}(t) + \omega_0^2 x(t) = F \cos(\omega_d t),
\end{equation}$$

where $\gamma$ is the damping coefficient and $\omega_0$ the natural frequency. Fourier-transforming both sides of this equation yields

$$\begin{equation}
(\omega_0^2 - \omega^2 + i \omega \gamma) \tilde{x}(\omega) = \frac{F}{2} \left[ \delta(\omega + \omega_d) + \delta(\omega - \omega_d) \right] \,.
\end{equation}$$

Evidently, $\tilde{x}(\omega)$ is only non-vanishing for $\omega = \pm \omega_d$. The system thus responds solely at the driving frequency, i.e., the time evolution can be captured by a single harmonic. This illustrates the general point that _linear systems are exactly solvable_ by transforming to Fourier space, where the equations are diagonal.

The situation becomes more complex if nonlinear terms are present, as these cause _frequency conversion_. Suppose we add a quadratic nonlinearity $\beta x^2(t)$ to the equations of motion; attempting to Fourier-transform, now, leads to

$$\begin{equation}
    \text{FT}[x^2](\omega) =  \int x^2(t) e^{-i\omega t} \: dt = \int_{-\infty}^{+\infty} \tilde{x}(\omega')\tilde{x}(\omega'') \delta(\omega''+\omega'-\omega) \: d\omega' \: d\omega'' \,,
\end{equation}$$

which couples all harmonics $\omega, \omega', \omega''$ such that $\omega + \omega' + \omega'' = 0$. To the lowest order, this means that the induced motion at the drive frequency generates a higher harmonic, $\omega_d \rightarrow 2\omega_d$. The frequency conversion couples the response at different frequencies and propagates through the spectrum, thus, _coupling an infinite number of harmonics_. Hence, the system is not easily solvable in Fourier space anymore!

## Harmonic ansatz &amp; harmonic equations {#Harmonic-ansatz-and-harmonic-equations}

Even though we need an infinity of Fourier components to describe our system exactly, some components are more important than others. The strategy of harmonic balance is to describe the motion of any variable $x_i(t)$ in a truncated Fourier space

$$\begin{equation}
x_i(t) = \sum_{j=1}^{M_i} u_{i,j}  (T)  \cos(\omega_{i,j} t)+ v_{i,j} (T) \sin(\omega_{i,j} t) \,.
\end{equation}$$

Within this space, the system is described by a finite-dimensional vector

$$\begin{equation}
\mathbf{u}(T) = (u_{1,1}(T), v_{1,1}(T), \ldots u_{N, M_N}(T), v_{N, M_N}(T))\,.
\end{equation}$$

Under the assumption that $\mathbf{u}(T)$ evolves at much slower timescales than the oscillatory harmonics, we may neglect all of its higher order time derivatives. Notice that once the ansatz is used in the equations of motion, all terms become oscillatory - each prefactor of $\cos(\omega_{i,j} t)$ and $\sin(\omega_{i,j} t)$ thus generates a separate equation. Collecting these, we obtain a set of coupled 1st order nonlinear ODEs,

$$\begin{equation}
\frac{d\mathbf{u}(T)}{dT}  = \bar{\mathbf{F}} (\mathbf{u})\,,
\end{equation}$$

which we call the _harmonic equations_. The main purpose of HarmonicBalance.jl is to obtain and solve these Harmonic equations. We are primarily interested in _steady states_ $\mathbf{u}_0$ defined by $\bar{\mathbf{F}}(\mathbf{u}_0) = 0$.

The process of obtaining the harmonic equations is best shown using the example below:

## Example: the Duffing oscillator {#Duffing_harmeq}

Here, we derive the harmonic equations for a single Duffing resonator, governed by the equation

$$\begin{equation}
    \ddot{x}(t) + \omega_0^2 x(t) + \alpha x^3(t) = F \cos(\omega_d t + \theta)\,.
\end{equation}$$

As explained [above](/background/harmonic_balance#prelude), for a periodic driving at frequency $\omega_d$ and a weak nonlinearity $\alpha$, we expect the response at frequency $\omega_d$ to dominate, followed by a response at $3\omega_d$ due to frequency conversion.

### Single-frequency ansatz {#Single-frequency-ansatz}

We first attempt to describe the steady states of the Duffing equations of motion using only a single harmonic, $\omega_d$, in the ansatz, i.e., our starting point is the following harmonic ansatz for $x$

$$\begin{equation}
	x(t) = u(T) \cos(\omega_d t) + v(T) \sin(\omega_d t)\:,
\end{equation}$$

with the harmonic variables $u$ and $v$. The _slow time_ $T$ is, for now, equivalent to $t$. Substituting this ansatz into the Duffing equation of motion results in

$$\begin{align} 
	\left[\ddot{u} + 2 \omega_d \dot{v} + u \left(\omega_0^2 - \omega_d^2 \right) +  \frac{3 \alpha \left(u^3 + uv^2\right)}{4} + F \cos{\theta}\right] &\cos(\omega_d t)& \\
	+ \left[\ddot{v} - 2 \omega_d \dot{u} + v \left(\omega_0^2 - \omega_d^2 \right)  +\frac{3 \alpha \left(v^3 + u^2 v\right)}{4} - F \sin{\theta}\right] &\sin(\omega_d t)& \nonumber \\
	+ \frac{\alpha \left(u^3 - 3 u v^2\right)}{4} \cos(3 \omega_d t) +  \frac{\alpha \left(3u^2 v - v^3\right)}{4} \sin(3 \omega_d t) &= 0. \nonumber
\end{align}$$

We see that the $x^3$ term has generated terms that oscillate at $3\omega_d$, describing the process of frequency upconversion. We now Fourier-transform both sides of the above equations with respect to $\omega_d$ to obtain the harmonic equations. This process is equivalent to extracting the respective coefficients of $\cos(\omega_d t)$ and $\sin(\omega_d t)$. Here, the distinction between $t$ and $T$ becomes important: since the evolution of $u(T)$ and $v(T)$ is assumed to be slow, they are treated as constant for the purpose of the Fourier transformation. Since we are interested in steady states, we drop the higher-order derivatives and rearrange the resulting equation to

$$\begin{equation}
	\frac{d}{dT} \begin{pmatrix} u \\ v  \end{pmatrix} = \frac{1}{8 \omega_d} \begin{pmatrix} 4 v \left(\omega_0^2-\omega_d^2 \right) + 3 \alpha \left(v^3 + u^2 v  \right) - 4 F \sin{\theta}  \\ 4 u \left(\omega_d^2-\omega_0^2 \right)  - 3 \alpha \left(u^3 + u v^2 \right) - 4 F \cos{\theta}  \end{pmatrix} \,.
\end{equation}$$

Steady states can now be found by setting the l.h.s. to zero, i.e., assuming $u(T)$ and $v(T)$ constant and neglecting any transient behaviour. This results in a set of 2 nonlinear polynomial equations of order 3, for which the maximum number of solutions set by [Bézout&#39;s theorem](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_theorem) is $3^2=9$. Depending on the parameters, the number of real solutions is known to be between 1 and 3 [[see proof](https://www.sciencedirect.com/science/article/pii/S0021782423001563?via%3Dihub)].

### Sidenote: perturbative approach {#Sidenote:-perturbative-approach}

The steady states describe a response that may be recast as $x_0(t) = X_0 \cos(\omega_d t + \phi)$, where $X_0=\sqrt{u^2+v^2}$ and $\phi=-\text{atan}(v/u)$. Frequency conversion from $\omega_d$ to $3 \omega_d$ can be found by setting $x(t) \equiv x_0(t) + \delta x(t)$ with $|\delta x(t)|\ll|x_0(t)|$ and expanding the equations of motion to first-order in $\delta x(t)$. The resulting equation

$$\begin{equation}
    \delta \ddot{x}(t) + \left[\omega_0^2 + \frac{3 \alpha X_0^2}{4} \right]\delta x(t) = - \frac{\alpha X_0^3}{4} \cos(3 \omega_d t + 3 \phi)\,,
\end{equation}$$

describes a simple harmonic oscillator, which is exactly soluble. Correspondingly, a response of $\delta x(t)$ at frequency $3 \omega_d$ is observed. Since this response is obtained &#39;on top of&#39; each steady state of the equations of motion, no previously-unknown solutions are generated in the process.

### Two-frequency ansatz {#Two-frequency-ansatz}

An approach in the spirit of harmonic balance is to use both harmonics $\omega_d$ and $3 \omega_d$ on the same footing, i.e., to insert the ansatz

$$\begin{equation}
	x(t) = u_1(T) \cos(\omega_d t) + v_1(T) \sin(\omega_d t) + u_2(T) \cos(3 \omega_d t) + v_2(T) \sin(3\omega_d t)\:,
\end{equation}$$

with $u_1, u_2, v_1, v_2$ being the harmonic variables. As before we substitute the ansatz into the equations of motion, drop second derivatives with respect to $T$ and Fourier-transform both sides. Now, the respective coefficients correspond to $\cos(\omega_d t)$, $\sin(\omega_d t)$, $\cos(3 \omega_d t)$ and $\sin(3\omega_d t)$. Rearranging, we obtain

$$	\begin{align}
	\begin{split}
	\frac{du_1}{dT} &=  \frac{1}{2\omega_d} \left[ \left({\omega_0}^{2} - \omega_d^2 \right) v_1 + \frac{3\alpha}{4} \left( v_1^3 + u_1^2 v_1 + u_1^2 v_2 - v_1^2 v_2 + 2 u_2^2 v_1 + 2 v_2^2 v_1 - 2 u_1 u_2 v_1\right)  + F \sin{\theta} \right],
	\\
	\frac{dv_1}{dT} &= \frac{1}{2\omega_d} \left[ \left({\omega_d}^{2} - \omega_0^2 \right) {u_1} - \frac{3 \alpha}{4} \left( u_1^3 + u_1^2 u_2 + v_1^2 u_1 - v_1^2 u_2+ 2 u_2^2 u_1 + 2 v_2^2 u_1  + 2 u_1 v_1 v_2\right) - F \cos{\theta} \right],
	\\
	\frac{d u_2}{dT} &= \frac{1}{6 \omega_d} \left[ \left(\omega_0^{2} - 9\omega_d^2 \right) {v_2} + \frac{\alpha}{4} \left( - v_1^3 + 3 v_2^3 + 3 u_1^2 v_1 + 6 u_1^2 v_2 + 3 u_2^2 v_2 + 6 v_1^2 v_2\right) \right],
	\\
	\frac{dv_2}{dT} &= \frac{1}{6 \omega_d} \left[ \left(9\omega_d^2 - \omega_0^2\right) {u_2} - \frac{\alpha}{4} \left( u_1^3 + 3 u_2^3 + 6 u_1^2 u_2 - 3 v_1^2 u_1 + 3 v_2^2 u_2 + 6 v_1^2 u_2\right) \right] \:.
	\end{split}
	\end{align}$$

In contrast to the single-frequency ansatz, we now have 4 equations of order 3, allowing up to $3^4=81$ solutions (the number of unique real ones is again generally far smaller). The larger number of solutions is explained by higher harmonics which cannot be captured perturbatively by the single-frequency ansatz. In particular, those where the $3 \omega_d$ component is significant. Such solutions appear, e.g., for $\omega_d \approx \omega_0 / 3$ where the generated $3 \omega_d$ harmonic is close to the natural resonant frequency. See the [examples](/tutorials/steady_states#Duffing) for numerical results.
