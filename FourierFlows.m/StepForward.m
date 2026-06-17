function prob = StepForward(prob, Nt)

arguments
    prob struct
    Nt   double {isinteger}
end

dt = prob.dt;
t = prob.t;

RHS = @(u, grid, t) prob.L .* u + prob.N(u, grid, t);

switch prob.stepper

    case 'Euler'

        for it = 1:Nt
        
            prob.u = prob.grid.filter .* (prob.u + dt * RHS(prob.u, prob.grid, t));
            t = t + dt;
        
        end

    case 'RK2'

        for it = 1:Nt

            K1 = RHS(prob.u, prob.grid, t);
            K2 = RHS(prob.u + dt * K1 / 2, prob.grid, t + dt / 2);
            prob.u = prob.grid.filter .* (prob.u + dt * K2);
            t = t + dt;

        end

    case 'RK4'

        for it = 1:Nt

            K1 = RHS(prob.u, prob.grid, t);
            K2 = RHS(prob.u + dt * K1 / 2, prob.grid, t + dt / 2);
            K3 = RHS(prob.u + dt * K2 / 2, prob.grid, t + dt / 2);
            K4 = RHS(prob.u + dt * K3, prob.grid, t + dt);
            prob.u = prob.grid.filter .* (prob.u + dt / 6 * (K1 + 2 * K2 + 2 * K3 + K4));
            t = t + dt;

        end

    otherwise

        error('Stepper:UnrecognisedValue', 'Stepper is not recognised.')

end

prob.iter = prob.iter + Nt;
prob.t = t;

end