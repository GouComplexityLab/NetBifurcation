clear all
close all
clc

load('./brain_data/Aij_human_Streamlines_scale82.mat');
load('./brain_data/Dij_human_Streamlines_scale82.mat');

adjacent_matrix_brain = Aij./(Dij+eps);

adj_matrix = sparse(adjacent_matrix_brain);

save adj_matrix_brainnetwork adj_matrix

