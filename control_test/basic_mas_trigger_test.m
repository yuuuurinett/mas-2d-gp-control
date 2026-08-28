%% Basic first-order MAS: consensus dynamics and event-trigger diagnostics
% Model:
%   x_dot_i = u_i
%   u_i     = -Kappa_P * sum_j a_ij (x_i - x_j)
%
% The controller uses the current states. Events are recorded separately so
% that this script tests the trigger/termination logic without GP, DAC, AC,
% or time-trigger interactions.

clear; clc; close all;

%% Parameters
N = 6;
Kappa_P = 1.0;
dt = 1e-3;
T_end = 20;
t = 0:dt:T_end;

% Same undirected ring as the control task: 1--2--3--4--5--6--1
A = zeros(N);
for i = 1:N
    j = mod(i, N) + 1;
    A(i,j) = 1;
    A(j,i) = 1;
end
L = diag(sum(A,2)) - A;
edge_mask = triu(A > 0, 1);

x0 = [-2.0; 1.2; 3.0; -0.8; 2.1; -1.4];

% Event: agent i broadcasts when its change from the last broadcast reaches
% delta_trigger. A fixed deadband lets events stop near consensus.
delta_trigger = 2e-2;

% Terminate only when BOTH errors remain small and no new event occurs for
% settle_time seconds. The no-event requirement prevents the local change
% error from appearing converged merely because an event just reset it.
eps_change = delta_trigger;
eps_consensus = 1e-2;
settle_time = 0.25;
settle_steps = max(1, round(settle_time/dt));

%% Simulation storage
n_steps = numel(t);
x = zeros(N, n_steps);
x(:,1) = x0;
x_broadcast = x0;

change_error = zeros(1, n_steps);
consensus_error = zeros(1, n_steps);
trigger_log = false(N, n_steps);
trigger_log(:,1) = true;
trigger_times = cell(N,1);
for i = 1:N
    trigger_times{i} = 0;
end

change_error(1) = max(abs(x(:,1) - x_broadcast));
consensus_error(1) = max_edge_error(x(:,1), edge_mask);

settled_counter = 0;
stop_index = n_steps;
terminated = false;

%% Forward-Euler simulation
for k = 1:n_steps-1
    u = -Kappa_P * L * x(:,k);
    x(:,k+1) = x(:,k) + dt*u;

    % Local event condition: same-agent state change x_i(t)-x_i(t_k^i).
    local_change = abs(x(:,k+1) - x_broadcast);
    fired = local_change >= delta_trigger;
    if any(fired)
        x_broadcast(fired) = x(fired,k+1);
        trigger_log(fired,k+1) = true;
        fired_agents = find(fired);
        for q = 1:numel(fired_agents)
            i = fired_agents(q);
            trigger_times{i}(end+1) = t(k+1); %#ok<SAGROW>
        end
    end

    % Evaluate the two errors after processing events.
    change_error(k+1) = max(abs(x(:,k+1) - x_broadcast));
    consensus_error(k+1) = max_edge_error(x(:,k+1), edge_mask);

    both_small = change_error(k+1) <= eps_change && ...
                 consensus_error(k+1) <= eps_consensus;
    stable_without_event = both_small && ~any(fired);
    if stable_without_event
        settled_counter = settled_counter + 1;
    else
        settled_counter = 0;
    end

    if settled_counter >= settle_steps
        stop_index = k + 1;
        terminated = true;
        break;
    end
end

% Trim unused preallocated samples after termination.
t = t(1:stop_index);
x = x(:,1:stop_index);
change_error = change_error(1:stop_index);
consensus_error = consensus_error(1:stop_index);
trigger_log = trigger_log(:,1:stop_index);

%% Summary
trigger_counts = sum(trigger_log,2) - 1; % exclude the initial broadcast
fprintf('\nBasic MAS event-trigger test\n');
fprintf('Termination rule: change error <= %.3g AND consensus error <= %.3g\n', ...
    eps_change, eps_consensus);
fprintf('Terminated: %s, final time: %.3f s\n', string(terminated), t(end));
fprintf('Final change error: %.6g\n', change_error(end));
fprintf('Final consensus error: %.6g\n', consensus_error(end));
fprintf('Trigger counts (initial broadcast excluded):\n');
disp(trigger_counts.');

for i = 1:N
    intervals = diff(trigger_times{i});
    if isempty(intervals)
        fprintf('Agent %d: no post-initial triggers.\n', i);
    else
        fprintf(['Agent %d: first interval %.4f s, last interval %.4f s, ' ...
                 'mean interval %.4f s.\n'], ...
                 i, intervals(1), intervals(end), mean(intervals));
    end
end

%% Plots
figure('Color','w','Name','Basic MAS trigger test');

subplot(2,2,1);
plot(t, x.', 'LineWidth', 1.3);
grid on;
xlabel('Time (s)'); ylabel('x_i');
title('Agent states');
legend(compose('Agent %d',1:N), 'Location','best');

subplot(2,2,2);
semilogy(t, max(consensus_error,eps), 'LineWidth',1.4); hold on;
yline(eps_consensus,'--','Consensus threshold');
grid on;
xlabel('Time (s)'); ylabel('max_{(i,j) in E}|x_i-x_j|');
title('Current-state consensus error');

subplot(2,2,3);
semilogy(t, max(change_error,eps), 'LineWidth',1.4); hold on;
yline(eps_change,'--','Change threshold');
grid on;
xlabel('Time (s)'); ylabel('max_i|x_i-x_i(t_k^i)|');
title('Same-agent change error');

subplot(2,2,4); hold on;
for i = 1:N
    event_indices = find(trigger_log(i,:));
    plot(t(event_indices), i*ones(size(event_indices)), '|', ...
        'MarkerSize',8,'LineWidth',1.2);
end
grid on;
xlabel('Time (s)'); ylabel('Agent');
yticks(1:N); ylim([0.5, N+0.5]);
title('Trigger raster (including t=0)');

sgtitle(sprintf('First-order MAS, Kappa_P=%.2f, delta=%.3f', ...
    Kappa_P, delta_trigger));

%% Local helper
function value = max_edge_error(x_now, edge_mask)
    pairwise = abs(x_now - x_now.');
    edge_values = pairwise(edge_mask);
    if isempty(edge_values)
        value = 0;
    else
        value = max(edge_values);
    end
end
