# Julia Playground
random julia files, mainly chemical enginering

- ```water.jl```: Requires ```ForwardDiff```- a working but barebones iapws95(2018) implementation (water properties), with the auxiliary function ```pressure(rho,T)``` that computes the pressure, using the derivative of the helholtz function. all other properties can be defined in terms of the gradient and the hessian of the main function.
