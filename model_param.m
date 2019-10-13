classdef model_param < handle
    properties
       model_name
       rm1 % jam density of class 1
       rm2 % jam density of class 2
       vm1 % max speed of class 1
       vm2 % max speed of class 2
       vm 
       s1 
       s2
       alpha1
       alpha2
       lda
       c1
       c2
       tfinal % length of simulation
       len % length of roadway
       x 
       dx % discretization in space
       lambda
       dt % discretization in time, dt, dx and vm safisty CFL condition
       t
       M % total simulation timesteps
       N % total number of discretized cells
       d1l % upstream, class1
       d2l % upstream, class 2
       d1r % downstream, class 1
       d2r % downstream, class 2
       pce
    end
end
