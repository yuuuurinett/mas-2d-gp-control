# Control-task aggregation modes: logic to confirm with advisor

## Common online learning clock

- Physical integration step: `h = 0.01 s`.
- Each agent adds one local data pair and updates its local GP every
  `T_online = 0.1 s`.
- For PoE, at query/inducing point `z`, agent `i` contributes the information
  vector

  `P_i(z,k) = [mu_i(z,k)/sigma_i^2(z,k); 1/sigma_i^2(z,k)]`.

- The exact centralized information average used only as a diagnostic is

  `P_bar(z,k) = (1/N) sum_i P_i(z,k)`.

The PoE prediction is decoded from the aggregated numerator and denominator.

## IP-DAC (online dynamic-average tracker)

1. Evaluate every local GP on the same fixed set of `M` inducing points.
2. Keep the DAC internal state `Zeta_i` across all physical steps and all
   updates of `P_i`.
3. Form the local aggregation state `Xi_i = P_i - Zeta_i` and evolve `Zeta`
   with neighbor differences of the most recently broadcast `Xi_hat`.
4. At each communication opportunity, agent `i` checks the packaged trigger

   `||Xi_hat_i-Xi_i||_F^2/M >`
   `[sum_j a_ij ||Xi_hat_i-Xi_hat_j||_F^2/M + epsilon_i^2]/(4 d_i)`.

5. If triggered, agent `i` broadcasts its entire `p x M` state and this counts
   as one agent-level broadcast. There is no permanent terminal condition.

Purpose: track the changing `P_bar(k)` while reusing historical consensus
state and reducing communication.

## IP-AC-R (restarted static-average baseline)

At every `0.1 s` update of `P_i(k)`:

1. Restart `Xi_i^0(k) = P_i(k)` and `Xi_hat_i^0(k) = P_i(k)`.
2. Execute exactly `R` static consensus rounds

   `Xi^(l+1)(k) = Xi^l(k) - alpha L Xi_hat^l(k)`.

3. In each round, use the packaged agent-level trigger

   `||Xi_hat_i-Xi_i||_F^2/M > c_i ||sum_j a_ij(Xi_i-Xi_j)||_F^2/M`,

   where `c_i = sigma_i a(1-a N_i)/N_i`.
4. Use `Xi^R(k)` even if it has not reached an evaluation tolerance. The
   tolerance is diagnostic only and cannot terminate the fixed-round run.

`AC-10` corresponds directly to ten consensus opportunities in a `0.1 s`
online interval when the communication interval is `0.01 s`. `AC-20` is a
communication-budget sensitivity case; interpreting all 20 rounds as physical
communication inside `0.1 s` requires a `0.005 s` communication interval.

Purpose: compare repeatedly solving a static snapshot with continuously
tracking the moving average.

## TP-DAC

- Evaluate each local GP directly at each agent's current test/query state.
- For each query, treat the local information vectors as time-varying DAC
  inputs.
- Preserve the DAC state across physical steps; use neighbor communication.
- No inducing-point approximation is used.

Purpose: isolate DAC error from inducing-point projection error.

## TP-AC

- Evaluate local GPs at the current query state.
- Restart a static AC problem for that query snapshot.
- Use either a fixed communication budget (preferred for a fair DAC
  comparison) or a convergence tolerance (high-communication reference), but
  label the choice explicitly.

## CEN

- A central node has all agents' local GP predictions at the current query.
- It computes the selected PoE/BCM/RBCM aggregation exactly.
- There is no consensus approximation and no event-trigger communication.

Purpose: ideal centralized aggregation reference.

## NBR

- Agent `i` aggregates only its own local GP and the local GPs of its graph
  neighbors at its current query.
- There is no iterative global consensus and no event trigger in the current
  implementation.

Purpose: one-hop information-sharing baseline.

## Local

- Agent `i` uses only its own online local GP at its own current state.
- No model aggregation or inter-agent GP communication.

## Exact

- The controller is supplied with the true unknown dynamics used by the plant.
- No GP learning, approximation, consensus, or event-trigger communication.

Purpose: oracle control reference. It is not a GP aggregation method.

## Offline

- The local datasets are generated before the control run and remain fixed.
- No new online point is added during the run.
- If AC/DAC is applied to these fixed local models, label the consensus method
  separately; `offline` describes the dataset clock, not the aggregation rule.

## Metrics reported with AC-10, AC-20, and DAC

1. Consensus error to exact current information average:

   `E(k,l) = sqrt((1/(N M)) sum_i ||Xi_i^l(k)-P_bar(k)||_F^2)`.

2. Prediction error against true unknown dynamics.
3. Full-trajectory control tracking error.
4. Agent-level packaged broadcasts per agent and communication raster.

## Questions for advisor

1. Does `0.1 s / 0.01 s` mean exactly ten communication opportunities for
   each update of `P_i`?
2. Should AC be a restarted fixed-round static baseline, while DAC preserves
   its internal state?
3. Should AC-20 be interpreted as a `0.005 s` communication interval, or only
   as an algorithmic communication-budget sensitivity test?
4. Is one broadcast of the complete inducing-point matrix counted as one
   agent-level event?
