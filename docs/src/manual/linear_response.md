# [Linear response (WIP)](@id linresp_man)

This module currently has two goals. One is calculating the first-order Jacobian, used to obtain stability and approximate (but inexpensive) the linear response of steady states. The other is calculating the full response matrix as a function of frequency; this is more accurate but more expensive. 

The methodology used is explained in Jan's thesis.

## Stability

The Jacobian is used to evaluate stability of the solutions. It can be shown explicitly,

```@docs
HarmonicBalance.get_Jacobian
```

## Linear response

The response to white noise can be shown with `plot_linear_response`. Depending on the `order` argument, different methods are used. 

```@docs
HarmonicBalance.LinearResponse.plot_linear_response
```

### First order

The simplest way to extract the linear response of a steady state is to evaluate the Jacobian of the harmonic equations. Each of its eigenvalues $\lambda$ describes a Lorentzian peak in the response; $\text{Re}[\lambda]$ gives its center and $\text{Im}[\lambda]$ its width. Transforming the harmonic variables into the non-rotating frame (that is, inverting the harmonic ansatz) then gives the response as it would be observed in an experiment.

The advantage of this method is that for a given parameter set, only one matrix diagonalization is needed to fully describe the response spectrum. However, the method is inaccurate for response frequencies far from the frequencies used in the harmonic ansatz (it relies on the response oscillating slowly in the rotating frame). 

Behind the scenes, the spectra are stored using the dedicated structs `Lorentzian` and `JacobianSpectrum`.
```@docs
HarmonicBalance.LinearResponse.JacobianSpectrum
HarmonicBalance.LinearResponse.Lorentzian
```

### Higher orders

Setting `order > 1` increases the accuracy of the response spectra. However, unlike for the Jacobian, here we must perform a matrix inversion for each response frequency.  

```@docs
HarmonicBalance.LinearResponse.ResponseMatrix
HarmonicBalance.get_response
HarmonicBalance.LinearResponse.get_response_matrix
```
