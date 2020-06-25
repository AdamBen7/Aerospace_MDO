% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

classdef FlightCondition
    %FLIGHTCONDITIONS Summary of this class goes here
    %   All in SI Units
    
    properties
        a;          % speed of sound
        rho;        % air density
        velocity;
        AoA_deg;
        AoA_rad;
        Mach;
        dynamic_pressure;
        
    end
    
    methods
        function obj = FlightCondition(velocity, AoA_deg)
            %assumptions: based on SL condition. Can be a function of 
            %height or use a accepted higher fidelity model.
            obj.a = 340;
            obj.rho = 1.226;
            
            %inputs
            obj.velocity = velocity;
            obj.AoA_deg = AoA_deg; 
            
            %calculated
            obj.AoA_rad = AoA_deg*pi/180;
            obj.Mach = obj.velocity/obj.a;
            obj.dynamic_pressure = (1/2)*obj.rho*(obj.velocity^2);
        end
        
        function outputArg = method1(obj,inputArg)
            outputArg = obj.Property1 + inputArg;
        end
    end
end

