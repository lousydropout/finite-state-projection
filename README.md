# Finite State Projection Method

## Summary

Many systems of chemical reactions, that are of interest to scientists and engineers, assumes the the number of molecules involved is large. This allows them to assume that said number of molecules can be represented as a continous variable, a real number, which then allows the system of chemical reactions to be represented as a system of ordinary differential equations (ODEs). This large number assumption, however, is not always valid. In these situations, the deterministic system of ODEs is replaced with a stochastic version.

There are two main approaches too solving the stochastic system of chemical reactions: a microscopic approach where you consider each possible series of events (out of potentially infinitely many possibilities) individually and a macroscopic approch where you track how the distribution of the various outcomes evolves over time. The Finite State Projection method is one such macroscopic approach, invented by Brian Munsky and Mustafa Khammash.

## Chemical Master Equation

Suppose that there are $M$ chemical reactions involving $N$ different chemical species $\mathbf{X} = (X_1, \dots, X_N)$. Let $a_\mu$ (known as the propensity function) be such that $a_\mu \Delta t$ represents the probability the $\mu^{th}$ chemical reaction occurs between time $t$ and $t + \Delta t$ for a given state $\mathbf{X}$ at time $t$. Further, let $\nu_\mu$ be stochiometric transition vector.

Then, without derivation, the probability that the system will be in state $\mathbf{X}$ at time $t + \Delta t$ will be given by
$$P(t + \Delta t, \mathbf{X}) = P(t, \mathbf{X}) \left(1 - \sum_{\mu=1}^M a_\mu(\mathbf{X}) \Delta t\right) + \sum_{mu=1}^M \left(P(t, \mathbf{X}) a_\mu (\mathbf{X} - \nu_\mu) \Delta t \right).$$

Since
$$
\dfrac{P(t + \Delta t, \mathbf{X}) - P(t, \mathbf{X})}{\Delta t} = - P(t, \mathbf{X}) \sum_{\mu=1}^M a_\mu(\mathbf{X}) + \sum_{mu=1}^M P(t, \mathbf{X}) a_\mu (\mathbf{X} - \nu_\mu),
$$
we let $\Delta t$ approach zero to arrive at the Chemical Master Equation:

$$
\dfrac{dP(t, \mathbf{X})}{d t} = - P(t, \mathbf{X}) \sum_{\mu=1}^M a_\mu(\mathbf{X}) + \sum_{mu=1}^M P(t, \mathbf{X}) a_\mu (\mathbf{X} - \nu_\mu).
$$

## The Finite State Projection method

The Finite State Projection Method is actually less a method and more a framework for solving the Chemical Master Equation directly. As it is, the Chemical Master Equation is likely easily solvable by many with a modicrum of knowledge regarding numerical differentiation methods. However, what's special about the Finite State Projection method is that it provides, thanks to the two theorems Munsky and Khammash derived in their paper but will not be shown here, a strict upper-bound for the error estimate (for any one-step method with the assumption that the initial condition for said state is perfectly known).

## The traditional approach and my approach

The traditional Finite State Projection approach to the Chemical Master Equation is to go directly from an initial time $t_i$ to a final time $t_f$ in one step, via an integral equation. The reason for this is simply because the Finite State Projection method assume that the initial state be perfectly known. If, instead, they proceeded with an intermediate timestep $t_{inter}$, then the error estimate provided by the Finite State Projection is no longer stricly applicable for the transition $t_{inter} \rightarrow t_f$.

Nevertheless, I have found this to be overly restrictive, if only because such an approach makes understanding how the probability distribution function evolves over time a much more cumbersome and expensive task. Also, I doubted the efficiency of the one-step approach, even if all we cared about is the state at the $t_f$, $P(t_f, \mathbf{X})$.

Instead, I apply a Runge-Kutta approach (just the adaptive first-order forward and backward Euler methods) and compare it to the traditional, integral approach. As a stand in for the exact solution, I use the fixed timestep Forward Euler method with incredibly small timesteps.

## Running the code

I have written the code twice, once in C++ and once in Ada 2012.

For the C++ programmers, the C++ codes lies in the `cpp/` folder with the driver code being `fsp.cpp`. The main section of the code is `line 311-370` and should be readily understood.

Further, a `makefile` has been provided and, so, compilation should be straight forward. Simply go inside the `cpp/` direction and type
```bash
make
```
and run the resulting binary.
