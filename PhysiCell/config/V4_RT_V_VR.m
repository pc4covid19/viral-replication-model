function [err, timepoints, species_out, observables_out] = V4_RT_V_VR( timepoints, species_init, parameters, suppress_plot )
%V4_RT_V_VR Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'V4_RT_V_VR' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. V4_RT_V_VR returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = V4_RT_V_VR( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 12 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 6 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 60, 1e4, 200, 250, 1, 0 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 6  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 6].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 12  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 12].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,10,20+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'TIME_CONV', 'ACE2_0', 'r_rep_half', 'max_infected_half_max', 'apoptosis_hill', 'virus_fraction_released_at_death' };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-8,   ...
               'AbsTol',   0.0001,   ...
               'Stats',    'off',  ...
               'BDF',      'off',    ...
               'MaxOrder', 5   );


% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [~, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
    if(length(timepoints) ~= size(species_out,1))
        exception = MException('ODE15sError:MissingOutput','Not all timepoints output\n');
        throw(exception);
    end
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), 9 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'nV', 'R_EB', 'V', 'U', 'R', 'P', 'A', 'E', 'Cell' };

    % construct figure
    plot(timepoints,observables_out);
    title('V4_RT_V_VR observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end


%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%

% Define if function to allow nested if statements in user-defined functions
function [val] = if__fun (cond, valT, valF)
% IF__FUN Select between two possible return values depending on the boolean
% variable COND.
    if (cond)
        val = valT;
    else
        val = valF;
    end
end

% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,12);
    species_init(1) = 10;
    species_init(2) = params(2);
    species_init(3) = 1;
    species_init(4) = 1;
    species_init(5) = 0;
    species_init(6) = 0;
    species_init(7) = 0;
    species_init(8) = 0;
    species_init(9) = 0;
    species_init(10) = 0;
    species_init(11) = 0;
    species_init(12) = 0;

end


% user-defined functions
% function infected_apoptosis_rate
function [val] = infected_apoptosis_rate(expressions, observables)
    val = (expressions(15)*((observables(7)^expressions(17))/((observables(7)^expressions(17))+(expressions(16)^expressions(17)))));
end

% function rna_rep_rate
function [val] = rna_rep_rate(expressions, observables)
    val = (expressions(10)/(observables(5)+expressions(11)));
end

% function virus_export_rate
function [val] = virus_export_rate(expressions, observables)
    val = (observables(9)*expressions(14));
end

% function vout_day
function [val] = vout_day(expressions, observables)
    val = (((observables(9)*expressions(14))*observables(7))*24);
end

% function dead_cell_release
function [val] = dead_cell_release(expressions, observables)
    val = (expressions(18)*infected_apoptosis_rate(expressions,observables));
end




% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,18);
    expressions(1) = parameters(1);
    expressions(2) = (0.01*expressions(1));
    expressions(3) = (0.01*expressions(1));
    expressions(4) = (0.01*expressions(1));
    expressions(5) = (0.01*expressions(1));
    expressions(6) = parameters(2);
    expressions(7) = (0.01*expressions(1));
    expressions(8) = (0.01*expressions(1));
    expressions(9) = (0.01*expressions(1));
    expressions(10) = (3*expressions(1));
    expressions(11) = parameters(3);
    expressions(12) = (0.010*expressions(1));
    expressions(13) = (0.01*expressions(1));
    expressions(14) = (0.01*expressions(1));
    expressions(15) = (0.002*expressions(1));
    expressions(16) = parameters(4);
    expressions(17) = parameters(5);
    expressions(18) = parameters(6);
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,9);
    observables(1) = species(1);
    observables(2) = species(5);
    observables(3) = species(8);
    observables(4) = species(9);
    observables(5) = species(10);
    observables(6) = species(11);
    observables(7) = species(12);
    observables(8) = species(4);
    observables(9) = species(3);

end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,9);
    ratelaws(1) = (0.01*expressions(1))*species(2)*species(1);
    ratelaws(2) = infected_apoptosis_rate(expressions,observables)*species(3);
    ratelaws(3) = (0.01*expressions(1))*species(5);
    ratelaws(4) = (0.01*expressions(1))*species(6);
    ratelaws(5) = (0.01*expressions(1))*species(7);
    ratelaws(6) = (0.01*expressions(1))*species(8);
    ratelaws(7) = (0.01*expressions(1))*species(9);
    ratelaws(8) = (0.01*expressions(1))*species(10);
    ratelaws(9) = rna_rep_rate(expressions,observables)*species(10);
    ratelaws(10) = (0.010*expressions(1))*species(10);
    ratelaws(11) = (0.01*expressions(1))*species(11);
    ratelaws(12) = virus_export_rate(expressions,observables)*species(12);
    ratelaws(13) = dead_cell_release(expressions,observables)*species(12);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(12,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calc_ratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = -ratelaws(1);
    Dspecies(2) = -ratelaws(1) +ratelaws(5);
    Dspecies(3) = -ratelaws(2);
    Dspecies(4) = ratelaws(12) +ratelaws(13);
    Dspecies(5) = ratelaws(1) -ratelaws(3);
    Dspecies(6) = ratelaws(3) -ratelaws(4);
    Dspecies(7) = ratelaws(4) -ratelaws(5);
    Dspecies(8) = ratelaws(4) -ratelaws(6);
    Dspecies(9) = ratelaws(6) -ratelaws(7);
    Dspecies(10) = ratelaws(7) +ratelaws(9) -ratelaws(10);
    Dspecies(11) = ratelaws(8) -ratelaws(11);
    Dspecies(12) = ratelaws(11) -ratelaws(12) -ratelaws(13);

end


end
