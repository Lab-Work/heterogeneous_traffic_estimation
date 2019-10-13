classdef params_PF < handle

    properties
        Np % number of particles
        model_stdev % noise stdev of model
        meas_stdev % noise stdev of measurement
        init_stdev % noise stdev of initial condition
        bound_stdev % noise stdev on boundary condition
        meas_pt % sensor locations
    end
    methods
        function a = getMeasStdev(obj)
            a = obj.meas_stdev;
        end
    end
end

