# Expert's evaluation
# H. Achicanoy
# CIAT, 2014

# =========================================================================== #
# PLS path modeling #
# =========================================================================== #

# Se asume que la opinión de un experto esta influyendo sobre la opinión de
# los otros (poco probable, ya que la opinión se considera como algo personal
# y por lo tanto independiente)

# --------------------------------------------------------------------------- #
# Avena
# --------------------------------------------------------------------------- #

EXP_1 = c(0,0)
EXP_2 = c(1,0)

exp_path = rbind(EXP_1,EXP_2)
innerplot(exp_path)

# blocks
exp_blocks = list(1:6, 7:12)

# modes
exp_mod = rep("A",2)

# scaling
exp_scale = list(c(rep("NUM",2),rep("ORD",4)), c(rep("NUM",2),rep("ORD",4)))

# modeling
exp_pls = plspm(x, exp_path, exp_blocks, modes=exp_mod)
plot(exp_pls)

plot(exp_pls, what="weights")
plot(exp_pls, what="loadings") # Correlación de la variable original con el índice

exp_pls$inner_model
exp_pls$inner_summary

exp_pls$unidim
exp_pls$scores

plot(exp_pls$scores)
abline(lm(exp_pls$scores[,2]~exp_pls$scores[,1]),col=2)
cor(exp_pls$scores)^2

rowMeans(exp_pls$scores)

mean(exp_pls$scores)

plot(x[,7],x[,1])

# --------------------------------------------------------------------------- #
# Cajanus
# --------------------------------------------------------------------------- #

EXP_1 = c(0,0,0)
EXP_2 = c(1,0,0)
EXP_3 = c(1,1,0)

exp_path = rbind(EXP_1,EXP_2,EXP_3)
innerplot(exp_path)

# blocks
exp_blocks = list(1:5, 6:11, 12:16)

# modes
exp_mod = rep("A",3)

# modeling
exp_pls = plspm(x, exp_path, exp_blocks, modes=exp_mod)
plot(exp_pls)

plot(exp_pls, what="weights")
plot(exp_pls, what="loadings") # Correlación de la variable original con el índice

# --------------------------------------------------------------------------- #
