Nb_pts = 30
range_pts = 1:Nb_pts
desired_clique_size = 10
range_desired_clique = sample(1:Nb_pts, desired_clique_size, replace=FALSE)

clique = expand.grid(range_desired_clique, range_desired_clique)

# Remove link to self : Point 1 --> Point 1
good_links = rownames(clique[which(!(clique[,1] == clique[, 2])),])
clique = clique[good_links,]

random_links = data.frame()

for (i in seq(1, (20*desired_clique_size), by=2))
{
  already_picked=list()
  rnd_id1 = sample(1:Nb_pts, 1)
  rnd_id2 = sample(1:Nb_pts, 1)
  while(is.element(rnd_id2, range_desired_clique) || rnd_id1==rnd_id2){
    rnd_id2 = sample(1:Nb_pts, 1)
  }
  already_picked= c(already_picked,rnd_id2)
  #print(already_picked)
  random_links[i, 1] = rnd_id1
  random_links[i, 2] = rnd_id2
  
  random_links[i + 1, 1] = rnd_id2
  random_links[i + 1, 2] = rnd_id1
}


colnames(random_links) = colnames(clique)

X = rbind(clique, random_links)

X = X[order(X["Var1"]), ]
X=unique(X)
rownames(X) = 1:length(X[, 1])

#head(X, 10)


df = data.frame(matrix(0, nrow=Nb_pts, ncol=Nb_pts))

for (i in 1:length(X[, 1]))
{
  line_axis = X[i, 1]
  column_axis  = X[i, 2]
  
  df[line_axis, column_axis] = 1
}

rownames(df) = colnames(df)

head(df, 10)


if (!isSymmetric(as.matrix(df)))
{
  print("There's a problem with the simulated link table : it's not symmetric, whereas it should be !")
} else
{
  print("Simulation is OK")
}


nmax = 5000



compute_error = function(G, links_table)
{
  H = 0
  for (i in 1:length(links_table[, 1]))
  {
    #print(links_table[i, 2])
    #print(G)
    if (G[links_table[i, 1]] == 1 && G[links_table[i, 2]] == 1)
    {
      H = H + 1
    }
  }
  return(H)
}  



transition1 = function(g, H, links_table, beta)
{
  # g est l'ensemble candidat clique maximum à K élements
  gnew = g
  clique_sz = sum(g)
  max_link = factorial(clique_sz)/factorial(clique_sz-2)
  # on tire une coordonnée au hasard dans g

  #clique_elements = g[find(g == 1)]
  clique_elements = which(g == 1)
  non_clique_elements = which(g == 0)
  k0 = sample(clique_elements, 1)
  k1 = sample(non_clique_elements, 1)
  gnew[k0] = 0
  gnew[k1] = 1
  
  # On calcul H (ou delta) : c'est la différence entre les anciennes et les nouvelles coordonnées
  delta = compute_error(gnew, links_table) - compute_error(g, links_table)
  alpha = max(min(exp(beta * (delta/max_link - 1)), 1), 0)
  
  u = runif(1)
  
  if (u < alpha)
  {
    # on accepte la transition
    result = c(gnew, compute_error(gnew, links_table))
  }
  else
  {
    # on rejette la transition
    result = c(g, H)
  }
  
  return(result)
}


transition2 = function(g, H, links_table, beta)
{
  gnew = g
  clique_sz = sum(g)
  max_link = factorial(clique_sz)/factorial(clique_sz-2)

  current_error = compute_error(g, links_table)
  
  clique_elements = which(g == 1)
  k0 = sample(clique_elements, 1)
  gnew[k0] = 0
  
  all_probas = rep(0, Nb_pts)
  for (i_j in 1:Nb_pts)
  {
    if (!(g[i_j] == 1))
    {
      gnew = g
      gnew[k0] = 0
      gnew[i_j] = 1
      delta_num = 0
      delta = compute_error(gnew, links_table) - current_error
      alpha = max(min(exp(beta * (delta/max_link - 1)), 1), 0)
      
      all_probas[i_j] = alpha     
    }
    else
    {
      all_probas[i_j] = 0     
    }
    
  }
  
  best_point = which.max(all_probas)
  
  gnew = g
  gnew[k0] = 0
  gnew[best_point] = 1
  
  Hnew = compute_error(gnew, links_table)
  
  return(c(gnew, Hnew))
}


compute_G = function(G, nmax_simul, num_method, links_table)
{
  transition_method = c(transition1, transition2)
  clique_sz = sum(G)
  max_link = factorial(clique_sz)/factorial(clique_sz-2)
  nb_iteration=0
  stopifnot(match(num_method, 1:4) > 0)
  
  H = rep(0, nmax_simul)
  H[1] = compute_error(G, links_table)
  
  beta_fixed = 1
  beta_0 = 0.01
  for (i in 1:nmax_simul)
  {
    if (is.element(num_method, c(3, 4)))
    {
      result = transition_method[[num_method - 2]](G, H[i], links_table, beta_0 * sqrt(i))
    }
    else
    {
      result = transition_method[[num_method]](G, H[i], links_table, beta_fixed)
    }
    
    #G = result[1:clique_sz]
    G = result[1:Nb_pts]
    H[i+1] = result[Nb_pts+1]
    nb_iteration=i
    
    if (H[i+1] == max_link)
    {
      H[i+2:length(H)] = H[i+1]
      
      break
    }
    # TEST

  }
  
  return(list("G"= G, "H"= H,"nb_iter"=nb_iteration))
}

old_warn_val = getOption("warn")
options(warn = -1)

start_clique_sz = 3

nb_max_clique1 = 0
nb_max_clique2 = 0
G_opt1_list = list()
H_opt1_list = list()
G_opt2_list = list()
H_opt2_list = list()
G_opt3_list = list()
H_opt3_list = list()
G_opt4_list = list()
H_opt4_list = list()
i_row = 1
clique1_is_def = FALSE
clique2_is_def = FALSE
clique3_is_def = FALSE
clique4_is_def = FALSE
nbr_iteration_list=list()
G0 = rep(0, Nb_pts)
#G1:start_clique_sz0[random_index] = 1
for (i_clique_sz in start_clique_sz:Nb_pts)
{
  # Méthode Transition 1
  # Initialissation aleatoire de G0
  #G0 = sample(1:Nb_pts,i_clique_sz,replace=FALSE)
  G0[1:i_clique_sz]=1
  
  #result = compute_G(G0, nmax, 1, X)
  result = compute_G(G0, nmax, 1, X)
  nb_iteration_transition1=result["nb_iter"]
  G_opt1_list[i_row] = result["G"]
  H_opt1_list[i_row] = result["H"]
  
  H_opt1 = unlist(H_opt1_list[i_row])
  clique_exists1 = H_opt1[nmax] == factorial(i_clique_sz)/factorial(i_clique_sz-2)
  
  plot(H_opt1,
       main="Transition 1",
       xlab="Nombre de Simulations",
       ylab="Erreur totale H")
  mtext(paste("Clique de taille ", i_clique_sz,
              " ==> ", "Clique détectée : ", clique_exists1))
  
  if (!clique_exists1 && !clique1_is_def)
  {
    nb_max_clique1 = i_clique_sz - 1
    clique1_is_def = TRUE
  }
  
  ## Méthode Transition 2
  ## on enlève les warning car on en déclenche certains
  old_warn_val = getOption("warn")
  options(warn = -1)
  
  #result = compute_G(G0, nmax, 2, X)
  result = compute_G(G0, nmax, 2, X)
  nb_iteration_transition2=result["nb_iter"]
  G_opt2_list[i_row] = result["G"]
  H_opt2_list[i_row] = result["H"]
  
  H_opt2 = unlist(H_opt2_list[i_row])
  clique_exists2 = H_opt2[nmax] == factorial(i_clique_sz)/factorial(i_clique_sz-2)
  
  plot(H_opt2,
       main="Transition 2",
       xlab="Nombre de Simulations",
       ylab="Erreur totale H")
  mtext(paste("Clique de taille ", i_clique_sz,
              " ==> ", "Clique détectée : ", clique_exists2))
  
  ## On remet les warnings
  options(warn = old_warn_val)
  
  if (!clique_exists2 && !clique2_is_def)
  {
    nb_max_clique2 = i_clique_sz - 1
    clique2_is_def = TRUE
  }
  
  # Méthode Transition 3
  #result = compute_G(G0, nmax, 3, X)
  
  result = compute_G(G0, nmax, 3, X)
  nb_iteration_transition3=result["nb_iter"]
  G_opt3_list[i_row] = result["G"]
  H_opt3_list[i_row] = result["H"]
  
  H_opt3 = unlist(H_opt3_list[i_row])
  clique_exists3 = H_opt3[nmax] == factorial(i_clique_sz)/factorial(i_clique_sz-2)
  
  plot(H_opt3,
       main="Transition 3",
       xlab="Nombre de Simulations",
       ylab="Erreur totale H")
  mtext(paste("Clique de taille ", i_clique_sz,
              " ==> ", "Clique détectée : ", clique_exists3))
  
  if (!clique_exists3 && !clique3_is_def)
  {
    nb_max_clique3 = i_clique_sz - 1
    clique3_is_def = TRUE
  }
  
  # Méthode Transition 4
  #result = compute_G(G0, nmax, 4, X)
  result = compute_G(G0, nmax, 4, X)
  nb_iteration_transition4=result["nb_iter"]
  G_opt4_list[i_row] = result["G"]
  H_opt4_list[i_row] = result["H"]
  
  H_opt4 = unlist(H_opt4_list[i_row])
  clique_exists4 = H_opt4[nmax] == factorial(i_clique_sz)/factorial(i_clique_sz-2)
  
  plot(H_opt4,
       main="Transition 4",
       xlab="Nombre de Simulations",
       ylab="Erreur totale H")
  mtext(paste("Clique de taille ", i_clique_sz,
              " ==> ", "Clique détectée : ", clique_exists4))
  
  if (!clique_exists4 && !clique4_is_def)
  {
    nb_max_clique4 = i_clique_sz - 1
    clique4_is_def = TRUE
  }
  
  nbr_iteration_list[[i_row]]=list(i_row+2,nb_iteration_transition2,nb_iteration_transition4)
  print(nb_iteration_transition2)
  print(nb_iteration_transition4)
  #print(nbr_iteration_list[i_row])
  if (clique1_is_def && clique2_is_def && clique3_is_def && clique4_is_def)
  {
    break
  }
  
  i_row = i_row + 1
}

install.packages("xtable")
library("xtable")
output=NULL
output = matrix(unlist(nbr_iteration_list), ncol = 5, byrow = TRUE)
colnames(output)=c(" ","Transition 1", "Transition 2","Transition 3","Transition 4")

options(warn = old_warn_val)