require(tidyverse)
require(tidygraph)
require(ggraph)

#Load summary csv of all response functions ############
RFxSum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

#Generate network elements based on tutorial https://www.jessesadler.com/post/network-analysis-with-r/

#Layer 1 type of experiment ############
  sources1 <- RFxSum %>% group_by(Data_type) %>% 
    summarize(Samp_size = n()) %>% 
    rename(label = Data_type)
  
  dest1 <- RFxSum %>% group_by(Experiment) %>% 
    summarize(Samp_size = n()) %>% 
    rename(label = Experiment)
  
  nodes1 <- bind_rows(sources1, dest1) %>% 
    rowid_to_column("id")
  
  edges1.1 <- RFxSum %>% group_by(Data_type, Experiment) %>% 
    summarise(weight = n()) %>% 
    ungroup()
  
  edges1.2 <- edges1.1 %>% 
    left_join(nodes1, by = c("Data_type" = "label")) %>% 
    rename(from = id)
  
  edges1 <- edges1.2 %>% 
    left_join(nodes1, by = c("Experiment" = "label")) %>% 
    rename(to = id) %>% 
    dplyr::select(from, to, weight)

#Layer 2 class of chemical  ##############
  sources2 <- nodes1 %>% filter(label %in% dest1$label)
  
  dest2 <- RFxSum %>% group_by(Class) %>% 
    summarize(Samp_size = n()) %>% 
    rename(label = Class) %>% 
    rowid_to_column("init_id") %>% 
    mutate(id = init_id + max(sources2$id)) %>% 
    dplyr::select(id, label, Samp_size)
  
  nodes2 <- bind_rows(sources2, dest2)
  
  edges2.1 <- RFxSum %>% group_by(Experiment, Class) %>% 
      summarise(weight = n()) %>% 
      ungroup()
    
  edges2.2 <- edges2.1 %>% 
      left_join(nodes2, by = c("Experiment" = "label")) %>% 
      rename(from = id)
    
  edges2 <- edges2.2 %>% 
      left_join(nodes2, by = c("Class" = "label")) %>% 
      rename(to = id) %>% 
      dplyr::select(from, to, weight)

#Layer 3 pathway affected  ##############
  sources3 <- nodes2 %>% filter(label %in% dest2$label)
  
  dest3 <- RFxSum %>% group_by(Pathway) %>% 
    summarise(Samp_size = n()) %>% 
    rename(label = Pathway) %>% 
    rowid_to_column("init_id") %>% 
    mutate(id = init_id + max(sources3$id)) %>% 
    dplyr::select(id, label, Samp_size)
  
  nodes3 <- bind_rows(sources3, dest3)
  
  edges3.1 <- RFxSum %>% group_by(Class, Pathway) %>% 
      summarise(weight = n()) %>% 
      ungroup()
    
  edges3.2 <- edges3.1 %>% 
      left_join(nodes3, by = c("Class" = "label")) %>% 
      rename(from = id)
    
  edges3 <- edges3.2 %>% 
      left_join(nodes3, by = c("Pathway" = "label")) %>% 
      rename(to = id) %>% 
      dplyr::select(from, to, weight)

#Layer 4 system affected  ##############
  sources4 <- nodes3 %>% filter(label %in% dest3$label)
  
  dest4 <- RFxSum %>% group_by(System) %>% 
    summarise(Samp_size = n()) %>% 
    rename(label = System) %>% 
    rowid_to_column("init_id") %>% 
    mutate(id = init_id + max(sources4$id)) %>% 
    dplyr::select(id, label, Samp_size)
  
  nodes4 <- bind_rows(sources4, dest4)
  
  edges4.1 <- RFxSum %>% group_by(Pathway, System) %>% 
      summarise(weight = n()) %>% 
      ungroup()
    
  edges4.2 <- edges4.1 %>% 
      left_join(nodes4, by = c("Pathway" = "label")) %>% 
      rename(from = id)
    
  edges4 <- edges4.2 %>% 
      left_join(nodes4, by = c("System" = "label")) %>% 
      rename(to = id) %>% 
      dplyr::select(from, to, weight)
    
#Plot network of first two layers 
  edges_l2 <- bind_rows(edges1, edges2)
  nodes_l2 <- bind_rows(nodes1, nodes2) %>% distinct()

  graph_l2 <- tbl_graph(nodes = nodes_l2, edges = edges_l2, directed = T)
 
  ggraph(graph_l2, layout = "graphopt") + 
    geom_node_point(aes(size = Samp_size)) +
    geom_edge_link(aes(width = weight, label = weight), alpha = 0.8) + 
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_text(aes(label = label), repel = TRUE) +
    #labs(edge_width = "Letters") +
    theme_graph()
  
  
#Plot network of first three layers 
  edges_l3 <- bind_rows(edges1, edges2, edges3)
  nodes_l3 <- bind_rows(nodes1, nodes2, nodes3) %>% distinct()

  graph_l3 <- tbl_graph(nodes = nodes_l3, edges = edges_l3, directed = T)
 
  ggraph(graph_l3, layout = "kk") + 
    geom_node_point(aes(size = Samp_size)) +
    geom_edge_link(aes(width = weight, label = weight), alpha = 0.8) + 
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_text(aes(label = label), repel = TRUE) +
    #labs(edge_width = "Letters") +
    theme_graph()
  
#Combine layers and create tidygraph object #######
  edges_fin <- bind_rows(edges1, edges2, edges3, edges4)
  nodes_fin <- bind_rows(nodes1, nodes2, nodes3, nodes4) %>% distinct()

  graph_fin <- tbl_graph(nodes = nodes_fin, edges = edges_fin, directed = T)
  
plot_fin <- ggraph(graph_fin, layout = "graphopt") + 
    geom_node_point(aes(size = Samp_size)) +
    geom_edge_link(aes(width = weight, label = weight), alpha = 0.8) + 
    scale_edge_width(range = c(0.5, 4)) +
    geom_node_text(aes(label = label), repel = TRUE) +
    #labs(edge_width = "Letters") +
    theme_graph()
  
plot_fin$data$x <- c(0, 
                     -15, 15,
                     45,0,-45,
                     60,30,-30,-60,
                     40,20,-20,-40)

plot_fin$data$y <- c(60, 
                     30, 30,
                     5,0,5,
                     -10,-15,-15,-10,
                     -40,-45,-45,-40)

plot_fin

#ggnet https://briatte.github.io/ggnet/#example-2-bipartite-network
network_fin <- network::network(edges_fin, vertex.attr = nodes_fin, 
                                matrix.type = "edgelist", ignore.eval = FALSE)

ggnet2(network_fin, size = "Samp_size", label = "vertex.names", edge.label = "weight")
