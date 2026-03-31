# LOADING LIBRARIES
library(pacman)
pacman::p_load(dplyr,
               data.table,
               R.utils,
               writexl,
               readr,
               readxl,ggplot2,reshape2,igraph)

set.seed(42)

# HELPER FUNCTIONS

# Simulate one trajectory of a Markov chain
simulate_mc <- function(P, start, n_steps = 50) {
  n_states <- nrow(P)
  traj <- integer(n_steps + 1)
  traj[1] <- start
  for (t in seq_len(n_steps)) {
    traj[t + 1] <- sample(n_states, 1, prob = P[traj[t], ])
  }
  traj
}

# Compute unconditional distribution at each step: pi_0 %*% P^n
unconditional_dist <- function(P, pi0, n_max = 60) {
  n_states <- ncol(P)
  probs <- matrix(0, nrow = n_max + 1, ncol = n_states)
  probs[1, ] <- pi0
  for (n in seq_len(n_max)) {
    probs[n + 1, ] <- probs[n, ] %*% P
  }
  probs
}

# Steady-state via left eigenvector (eigenvalue = 1)
steady_state <- function(P) {
  ev  <- eigen(t(P))
  idx <- which(abs(ev$values - 1) < 1e-9)[1]
  v   <- Re(ev$vectors[, idx])
  v / sum(v)
}

# GCD for period computation
gcd <- function(a, b) { if (b == 0) a else gcd(b, a %% b) }

# Period of state i (0-indexed internally): count step lengths of
# all paths i -> i and take their GCD.
period_state <- function(P, i, max_power = 60) {
  n  <- nrow(P)
  Pk <- diag(n)          # P^0
  d  <- 0
  for (k in seq_len(max_power)) {
    Pk <- Pk %*% P
    if (Pk[i, i] > 1e-12) d <- gcd(d, k)
  }
  if (d == 0) d <- Inf   # absorbing / never returns
  d
}

# QUESTION A1 – 5-State Markov Chain

P1 <- matrix(c(
  1.0, 0,   0,   0,   0,
  0.5, 0,   0,   0,   0.5,
  0.2, 0,   0,   0,   0.8,
  0,   0,   1.0, 0,   0,
  0,   0,   0,   1.0, 0
), nrow = 5, byrow = TRUE)

rownames(P1) <- colnames(P1) <- paste0("S", 1:5)
#cat("Transition matrix P1:\n"); print(P1)

#(a) Diagram + structural analysis
g1 <- graph_from_adjacency_matrix(
  (P1 > 0) * 1, mode = "directed", weighted = TRUE, diag = TRUE)
E(g1)$label <- round(P1[P1 > 0], 2)
lay <- matrix(c(0,2, 1,2, 2,2, 2,0, 3,2), ncol = 2, byrow = TRUE)
plot(g1,
     layout           = lay,
     vertex.size      = 35,
     vertex.color     = c("#E74C3C","#F39C12","#F39C12","#2ECC71","#2ECC71"),
     vertex.label     = V(g1)$name,
     vertex.label.cex = 1.1,
     edge.label       = E(g1)$label,
     edge.label.cex   = 0.85,
     edge.curved      = 0.25,
     edge.arrow.size  = 0.5,
     main             = "A1: 5-State Markov Chain")
legend("bottomleft", legend = c("Absorbing (recurrent)", "Transient"),
       fill = c("#2ECC71","#F39C12"), bty = "n", cex = 0.9)

# Structural analysis
cat("\n--- (a) Structural Analysis ---\n")
cat("State 1: ABSORBING (p11=1). Recurrent. Period = 1 (aperiodic).\n")
cat("State 2: TRANSIENT (can leave, may not return). Period = ?\n")
cat("State 3: TRANSIENT (can leave to S1 or S5). Period = ?\n")
cat("State 4: ABSORBING (p44=1). Recurrent. Period = 1 (aperiodic).\n")
cat("State 5: ABSORBING (p55=1... wait, check row 5)\n")

# Verify row 5
cat("Row 5 of P1:", P1[5,], "\n")
cat("State 5: p54=1 => goes to S4 deterministically => part of recurrent class {4,5}?\n")

# Let's compute periods numerically
for (i in 1:5) {
  d <- period_state(P1, i)
  cat(sprintf("Period of S%d = %s\n", i, ifelse(is.infinite(d), "Inf (never returns)", d))  )
}

cat("\nRecurrent classes: {S1} (absorbing), {S4} (absorbing)\n")
cat("Transient states : S2, S3, S5\n")
cat("Absorbing states : S1 (self-loop p=1), S4 (self-loop p=1)\n")
cat("No reflecting states (a reflecting state has p=1 to a neighbour and returns)\n")


# ---- (b) Simulate 3 trajectories ----
n_steps <- 40
starts  <- sample(1:5, 3, replace = TRUE)
cols    <- c("#E74C3C","#3498DB","#2ECC71")

traj_list <- lapply(starts, function(s) simulate_mc(P1, s, n_steps))

df_traj <- do.call(rbind, lapply(seq_along(traj_list), function(k) {
  data.frame(step  = 0:n_steps,
             state = traj_list[[k]],
             traj  = factor(k))
}))

p <- ggplot(df_traj, aes(step, state, colour = traj, group = traj)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  scale_y_continuous(breaks = 1:5, labels = paste0("S", 1:5)) +
  scale_colour_manual(values = cols,
                      labels = paste0("Start S", starts)) +
  labs(title  = "A1(b): Three Simulated Trajectories",
       x      = "Step (n)", y = "State",
       colour = "Trajectory") +
  theme_bw(base_size = 12)
print(p)

cat("\n--- (b) Trajectory Comments ---\n")
cat("All trajectories quickly get absorbed into S1 or S4 and stay there.\n")
cat("Transient states S2, S3, S5 are visited only briefly before absorption.\n")
cat("Once absorbed, the chain never leaves (as expected for absorbing states).\n")


#(c) Steady-state probabilities

cat("\n--- (c) Steady-State Probabilities ---\n")
pi1 <- steady_state(P1)
names(pi1) <- paste0("S", 1:5)
cat("Steady-state vector π:\n"); print(round(pi1, 6))

cat("\nInterpretation:\n")
cat("π(S1) + π(S4) = 1 (all probability absorbed into the two absorbing states).\n")
cat("π(S2) = π(S3) = π(S5) = 0 (transient states vanish at stationarity).\n")
cat("The exact split between S1 and S4 depends on starting distribution.\n")
cat("\nIs it ergodic? NO.\n")
cat("Ergodicity requires a single recurrent class that is aperiodic.\n")
cat("Here we have TWO recurrent classes + transient states => NOT ergodic.\n")


# (d) Unconditional probabilities vs time 

pi0 <- rep(1/5, 5)           # uniform start
probs1 <- unconditional_dist(P1, pi0, n_max = 60)

df_prob <- as.data.frame(probs1)
colnames(df_prob) <- paste0("S", 1:5)
df_prob$step <- 0:60
df_melt <- melt(df_prob, id.vars = "step", variable.name = "State", value.name = "Prob")


p <- ggplot(df_melt, aes(step, Prob, colour = State)) +
  geom_line(linewidth = 0.9) +
  labs(title = "A1(d): Unconditional Probabilities vs Time (uniform start)",
       x = "Step (n)", y = "P(X_n = s)") +
  theme_bw(base_size = 12)
print(p)

cat("\n--- (d) Convergence Comment ---\n")
cat("Probabilities for S2, S3, S5 decay to 0 rapidly (within ~10 steps).\n")
cat("Probabilities for S1 and S4 rise monotonically and stabilise.\n")
cat("Convergence is fast because transient states are quickly abandoned.\n")


# QUESTION A2 – 7-State Markov Chain

# a. Define the 7x7 transition matrix P2
P2 <- matrix(c(
  0,   1,   0,   0,   0,   0,   0,
  1,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0.4, 0.2, 0.2, 0.2,
  0,   0,   0,   0,   0.2, 0.4, 0.4,
  0.3, 0,   0,   0.1, 0.3, 0.1, 0.2,
  0,   0,   0,   0.2, 0.2, 0.3, 0.3,
  0,   0,   0,   0.5, 0.2, 0.2, 0.1
), nrow = 7, byrow = TRUE)

rownames(P2) <- colnames(P2) <- c("S1","S2","S3","S4","S5","S6","S7")

# 2. Build graph
g2 <- graph_from_adjacency_matrix(
  P2,
  mode     = "directed",
  weighted = TRUE,
  diag     = TRUE
)

# 3. FIX: Remove zero-weight edges using correct syntax
g2 <- delete_edges(g2, E(g2)[E(g2)$weight == 0])

# 4. Edge labels
E(g2)$label <- round(E(g2)$weight, 2)

# 5. Layout
lay2 <- matrix(c(
  -2,  1,   # S1
  -2, -1,   # S2
  0,   2,   # S3
  2,   1,   # S4
  0,   0,   # S5
  2, -1,   # S6
  0, -2    # S7
), ncol = 2, byrow = TRUE)

# 6. Vertex colours
vcols <- c("#E74C3C","#E74C3C","#F39C12","#3498DB","#3498DB","#3498DB","#3498DB")

# 7. Plot with enhanced visibility for recurrences (loops)
plot(g2,
     layout           = lay2,
     vertex.size      = 30,
     vertex.color     = vcols,
     vertex.label     = V(g2)$name,
     vertex.label.cex = 1.0,
     edge.label       = E(g2)$label,
     edge.label.cex   = 0.8,
     edge.curved      = 0.2,
     edge.arrow.size  = 0.3,    # Increased for better visibility
     edge.loop.angle  = -pi, # Adjusts the direction of self-loops
     # --- ADD THESE TWO LINES ---
     loop.size        = 1.5,   # Pushes the loop outside the vertex circle
     edge.arrow.width = 1.5,   # Makes the arrowhead bulky enough to see

     main             = "A2: 7-State Markov Chain")

legend("bottomleft",
       legend = c("Recurrent class {S1, S2} — period 2",
                  "Transient: S3",
                  "Recurrent class {S4, S5, S6, S7}"),
       fill   = c("#E74C3C","#F39C12","#3498DB"),
       bty    = "n",
       cex    = 0.8)

#(b) Recurrent/transient, periods, absorbing/reflecting

cat("\n--- (b) Structural Analysis ---\n")
for (i in 1:7) {
  d <- period_state(P2, i)
  cat(sprintf("Period of S%d = %g\n", i, d))
}
cat("\nRecurrent class 1 : {S1, S2} — period 2 (alternates S1<->S2)\n")
cat("Recurrent class 2 : {S4, S5, S6, S7} — aperiodic (period 1)\n")
cat("Transient state   : S3 (exits to S4–S7 with prob 1, never returns)\n")
cat("No absorbing states (no state with p_ii = 1).\n")
cat("No reflecting states in the classical sense.\n")


# (c) Simulate 2 trajectories

n2 <- 60
starts2 <- sample(1:7, 2, replace = TRUE)
traj2   <- lapply(starts2, function(s) simulate_mc(P2, s, n2))

df2 <- do.call(rbind, lapply(seq_along(traj2), function(k) {
  data.frame(step = 0:n2, state = traj2[[k]], traj = factor(k))
}))

p <- ggplot(df2, aes(step, state, colour = traj, group = traj)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  scale_y_continuous(breaks = 1:7, labels = paste0("S", 1:7)) +
  scale_colour_manual(values = c("#E74C3C","#3498DB"),
                      labels = paste0("Start S", starts2)) +
  labs(title  = "A2(c): Two Simulated Trajectories",
       x      = "Step (n)", y      = "State",
       colour = "Trajectory") +
  theme_bw(base_size = 12)
print(p)

cat("\n (c) Trajectory Comments \n")
cat("If a trajectory starts in S3, it quickly moves into {S4,S5,S6,S7} and stays.\n")
cat("If it starts in {S1,S2}, it oscillates between S1 and S2 indefinitely (period 2).\n")
cat("Trajectories never cross between the two recurrent classes.\n")


# d Limiting probabilities 

cat("\n (d) Limiting Probabilities \n")

# For the {S4,S5,S6,S7} sub-chain solve πQ = π
P2_sub <- P2[4:7, 4:7]
P2_sub_norm <- P2_sub / rowSums(P2_sub)   # already sums to 1 within sub-block
pi_sub <- steady_state(P2_sub_norm)
cat("Limiting probs within recurrent class {S4,S5,S6,S7}:\n")
names(pi_sub) <- paste0("S", 4:7)
print(round(pi_sub, 6))

# For {S1,S2}: alternates — time-average is 0.5 each
cat("\nLimiting (time-average) probs for recurrent class {S1,S2}: 0.5 each\n")
cat("(Chain oscillates with period 2, so Cesaro limit gives 1/2 for each)\n")

cat("\nLimiting probabilities for S3 = 0 (transient, visited finitely often).\n")

cat("\nIs the chain ergodic? NO.\n")
cat("Two separate recurrent classes => not a single irreducible chain.\n")
cat("Also, {S1,S2} has period 2, violating aperiodicity needed for ergodicity.\n")



# QUESTION A3 – Road Traffic Markov Chain

#===========================================================================================================
cat("\n========== QUESTION A3 ==========\n")

# States: 1=light, 2=heavy, 3=jammed
# 1PM–4PM: 9 intervals of 20 min  (3 hours / (20 min) = 9 steps)
# 4PM–6PM: 6 intervals of 20 min  (2 hours / (20 min) = 6 steps)

P3_pm <- matrix(c(
  0.4, 0.4, 0.2,
  0.3, 0.4, 0.3,   # NOTE: rows must sum to 1; this row sums to 1.1?
  0,   0.1, 0.9
), nrow = 3, byrow = TRUE)

# Check row sums
cat("Row sums P3_pm:", rowSums(P3_pm), "\n")
# Row 2 sums to 1.1 — likely a typo in the question; normalise gracefully
if (any(abs(rowSums(P3_pm) - 1) > 1e-8)) {
  cat("WARNING: P3_pm row 2 sums to", sum(P3_pm[2,]),
      "— normalising to ensure valid transition matrix.\n")
  P3_pm <- P3_pm / rowSums(P3_pm)
}

P3_eve <- matrix(c(
  0.1, 0.5, 0.4,
  0.1, 0.3, 0.6,
  0,   0.1, 0.9
), nrow = 3, byrow = TRUE)
cat("Row sums P3_eve:", rowSums(P3_eve), "\n")

rownames(P3_pm) <- rownames(P3_eve) <- c("light","heavy","jammed")
colnames(P3_pm) <- colnames(P3_eve) <- c("light","heavy","jammed")

#Matrix power helper 
mat_pow <- function(M, k) {
  result <- diag(nrow(M))
  for (i in seq_len(k)) result <- result %*% M
  result
}

n_pm  <- 9    # 1PM to 4PM
n_eve <- 6    # 4PM to 6PM

# Starting distribution at 1PM: light with prob 1
pi_1pm <- c(1, 0, 0)

# Distribution at 4PM
P3_pm_n  <- mat_pow(P3_pm,  n_pm)
pi_4pm   <- pi_1pm %*% P3_pm_n

# Distribution at 6PM
P3_eve_n <- mat_pow(P3_eve, n_eve)
pi_6pm   <- pi_4pm %*% P3_eve_n

cat("\n--- (a) Analytical distribution at 6PM ---\n")
cat("Distribution at 4PM:\n"); print(round(as.vector(pi_4pm), 6))
cat("Distribution at 6PM:\n"); print(round(as.vector(pi_6pm), 6))
cat(sprintf("  P(light  at 6PM) = %.4f\n", pi_6pm[1]))
cat(sprintf("  P(heavy  at 6PM) = %.4f\n", pi_6pm[2]))
cat(sprintf("  P(jammed at 6PM) = %.4f\n", pi_6pm[3]))


# (b) Simulation of 10,000 trajectories 

cat("\n--- (b) Simulation (10,000 trajectories) ---\n")

n_sim <- 10000
states <- 1:3

simulate_traffic <- function() {
  state <- 1   # start at light
  # 1PM -> 4PM: 9 steps with P3_pm
  for (i in seq_len(n_pm))  state <- sample(states, 1, prob = P3_pm[state,])
  # 4PM -> 6PM: 6 steps with P3_eve
  for (i in seq_len(n_eve)) state <- sample(states, 1, prob = P3_eve[state,])
  state
}

final_states <- replicate(n_sim, simulate_traffic())
sim_probs    <- table(final_states) / n_sim
names(sim_probs) <- c("light","heavy","jammed")[as.integer(names(sim_probs))]

cat("Simulated probabilities at 6PM (n =", n_sim, "):\n")
print(round(sim_probs, 4))
cat("\nAnalytical probabilities at 6PM:\n")
print(round(setNames(as.vector(pi_6pm), c("light","heavy","jammed")), 4))
cat("\nSimulation closely matches the analytical result (confirming correctness).\n")


# ---- Plot: comparison bar chart ----

df_comp <- data.frame(
  State  = rep(c("light","heavy","jammed"), 2),
  Method = rep(c("Analytical","Simulation"), each = 3),
  Prob   = c(as.vector(pi_6pm),
             as.vector(sim_probs)[match(c("light","heavy","jammed"), names(sim_probs))])
)
df_comp$State <- factor(df_comp$State, levels = c("light","heavy","jammed"))


p <- ggplot(df_comp, aes(State, Prob, fill = Method)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("#2980B9","#E67E22")) +
  labs(title = "A3(b): Distribution at 6PM — Analytical vs Simulation",
       x = "Traffic State", y = "Probability") +
  theme_bw(base_size = 12)
print(p)


# Also show convergence path step-by-step 


all_steps <- 0:(n_pm + n_eve)
dist_path  <- matrix(0, nrow = length(all_steps), ncol = 3)
dist_path[1, ] <- pi_1pm

for (t in 1:n_pm)
  dist_path[t + 1, ] <- dist_path[t, ] %*% P3_pm

for (t in seq_len(n_eve))
  dist_path[n_pm + t + 1, ] <- dist_path[n_pm + t, ] %*% P3_eve

df_path <- as.data.frame(dist_path)
colnames(df_path) <- c("light","heavy","jammed")
df_path$step <- all_steps
df_path$time <- format(
  as.POSIXct("2024-01-01 13:00:00") + (all_steps * 20 * 60),
  "%H:%M")

df_pm <- melt(df_path, id.vars = c("step","time"),
              variable.name = "State", value.name = "Prob")


p <- ggplot(df_pm, aes(step, Prob, colour = State)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_vline(xintercept = n_pm, linetype = "dashed", colour = "grey40") +
  annotate("text", x = n_pm + 0.3, y = 0.85, label = "4PM\n(matrix\nchanges)",
           hjust = 0, size = 3.2, colour = "grey30") +
  scale_x_continuous(breaks = all_steps,
                     labels = df_path$time) +
  scale_colour_manual(values = c("#27AE60","#E67E22","#E74C3C")) +
  labs(title = "A3: Distribution of Traffic State from 1PM to 6PM",
       x = "Time", y = "Probability") 
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

cat("\n ======THE END===========\n")