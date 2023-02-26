# Machine Component Design Optimization

Design projects for Boğaziçi University's **Machine Design II** (ME 424) class.

### Project 1: Gearing Design

Designing a planetary gear that will transmit the maximum amount of power with the given specifications and constraints. The design requires a bending safety factor at least 1.5 and a surface safety factor of at least 1.2 with 99 % reliability. Safe use of $1\times 10^{8}$ revolutions of the arm is required. 

<p align="center">
  <img src="https://github.com/edizferit/Machine_Component_Design_Optimization/blob/main/figures/p1.jpg?raw=true" width="40%">
</p>

### Project 2: Bolted Joint Design

Designing a bolted joint by choosing a suitable property class, diameter, and preload for the bolt such that the lowest of the safety factors for the possible failure modes is above 1.6 and the cost is minimum. Equilibrium and compatibility principles are used.

<p align="center">
  <img src="https://github.com/edizferit/Machine_Component_Design_Optimization/blob/main/figures/p2.jpg?raw=true" width="40%">
</p>

### Project 3: Spring Design

Designing a spring that has lowest cost and safe against fatigue and buckling failures with an infinite fatigue life.

<p align="center">
  <img src="https://github.com/edizferit/Machine_Component_Design_Optimization/blob/main/figures/p3.jpg?raw=true" width="50%">
</p>

**The governing equation:**

$$\frac{\partial^2 u(x,t)}{\partial x^2} = \frac{1}{\beta^2}\frac{\partial^2 u(x,t)}{\partial t^2}$$

**Axial force in the spring:**

$$F(x,t) = Lk\frac{\partial u(x,t)}{\partial x}$$

where $u$ is the displacement in the $x$ direction, $L$ is the length of the spring, $k$ is the spring rate, $m$ is the total mass of the spring, $\beta^2 = L^2k/m$



