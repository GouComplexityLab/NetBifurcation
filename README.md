This repository provides the Matlab codes that taste the bifurcation dynamics of the reaction-diffusion models defined on complex networks from our manuscript titled "Utilize bifurcation theory to scrutinize critical dynamics of the reaction-diffusion models defined on networks" by Wei Gou, Jian-Meng Cui, Yong-Li Song, Wei Lin* and Zhen Jin* submitted in 2025.

![](https://github.com/GouComplexityLab/NetBifurcation/blob/main/BifurcationDynamicsOnNetwork.png)
Figure: The Turing-Hopf bifurcation dynamics of the reaction-diffusion ratio-dependent predator-prey model defined on complex networks $G_{sw}(10,2,0.05)$.

# Guide for Code
These codes are partially developed based on the MatCont toolkit, which is a Matlab software project for the numerical continuation and bifurcation study of continuous and discrete parameterized dynamical systems and can be download from https://sourceforge.net/projects/matcont/.
To run our program, the user needs to place this folder in the folder where the pre downloaded MatCont is stored.
This folder includes several subfolders for different Figures in main text and Supplementary Information.

## network0 brain
This subfolder provides the adjacent matrix data of the weighted brain network in adj_matrix_brainnetwork.mat.mat and the nado names in brainregionnames.xlsx.

## network1 npbcws_n=10k=2p=0number=1
This subfolder provides the adjacent matrix data of the underlying network $G_{la}(10,1)$ in adj_matrix_nwnetwork_01.mat.

## network2 ws_n=10k=2p=0.05number=796
This subfolder provides the adjacent matrix data of the underlying network $G_{sw}(10,1,0.05)$ in adj_matrix_wsnetwork_796.mat.

## network3 ws_n=10k=4p=0.05number=69
This subfolder provides the adjacent matrix data of the underlying network $G_{sw}(10,2,0.05)$ in adj_matrix_wsnetwork_69.mat.


## gwPlotNetworks:
* step00_plot_network1_npbcwsn10k2p0number1.m Plot the underlying network $G_{la}(10,1)$ in Fig.1a.
* step00_plot_network2_wsn10k2p0dot05number796.m Plot the underlying network $G_{sw}(10,1,0.05)$ in Fig.1b.
* step00_plot_network3_wsn10k4p0dot05number69.m Plot the underlying network $G_{sw}(10,2,0.05)$ in Fig.1c.
* step00_plot_network0_brain_withoutnodename.m Plot the brain network without node names in the inset of Fig.5a.
* step00_plot_network0_brain_withnodename.m Plot the brain network with node names in the Fig.S11e.
* bubbleGraph_gou.m and bubbleGraph_BrainGou.m are modified based on the original code bubble digraph, which can be download from https://www.mathworks.com/matlabcentral/fileexchange/125140-bubble-digraph.
 
## gwRDPPmodelnetwork1
This subfolder provides a series of codes to simulation the Hopf, steady state and Turing-Hopf bifurcation dynamics of the reaction-diffusion ratio-dependent predator-prey (RDPP) model defined on networks $G_{la}(10,1)$, with step by step.
* RDPPmodelnetwork1_step01_BifurcationDiagramPlot.m Plot the bifurcation diagram of the RDPP model defined onnetworks $G_{la}(10,1)$ in $b − a$ plane in Fig.1d.
* RDPPmodelnetwork1_step02_calc_HBNF.m Claculate the Hopf bifurcation (HB) normal form of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step02_calc_SSBNF.m Claculate the steady state bifurcation (SSB) normal form of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step02_calc_SSB_criticalmode.m Plot the criticla Laplacian eigenvector (also called as critical mode) in the SSB of network $G_{la}(10,1)$.
* RDPPmodelnetwork1_step02_calc_THBNF.m Claculate the Turing-Hopf bifurcation (THB) normal form of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step02_THB_region.m Plot the bifurcation diagram in $µ_{1} − µ_{2}$ plane of the calculated THB normal form in Fig.3a, where the origin $(0,0)$ of $µ_{1} − µ_{2}$ plane correspons the THB point $(b_{THB},a_{THB})$ of $b − a$ plane, and plot choosed points $(µ_{1},µ_{2})$ in each divided region to simulate THB dynamics later.
* RDPPmodelnetwork1_step03_write_NetRDsystemtxt.m Write the the RDPP model defined onnetworks $G_{la}(10,1)$ in RDPPSystemOnNetwork1.txt.
* RDPPmodelnetwork1_step04_matcont_HB.m Using MatCont to obtain the data for Plotting bifurcation diagram around HB point of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step04_matcont_SSB.m Using MatCont to obtain the data for Plotting bifurcation diagram around SSB point of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step05_matcont_HB_plot3.m  Prepare (load data) to plot bifurcation diagram around HB point of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step05_matcont_HB_plot3_plot_by_nodes.m Plot bifurcation diagram, for each node, around HB point of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step05_matcont_HB_plot3_plot_all_nodes.m Plot bifurcation diagram, for all nodes, around HB point of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step05_matcont_SSB_plot3.m Prepare (load data) to plot bifurcation diagram around SSB point of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step05_matcont_SSB_plot3_plot_by_nodes.m Plot bifurcation diagram, for each node, around SSB point of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step05_matcont_SSB_plot3_plot_all_nodes.m Plot bifurcation diagram, for all nodes, around SSB point of the RDPP model defined onnetworks $G_{la}(10,1)$.
* RDPPmodelnetwork1_step06_sim_ode45_THB1.m Simulate the THB dynamics from four different initial conditions of the RDPP model defined onnetworks $G_{la}(10,1)$ with $(b,a)=(b_{THB}+µ_{1},a_{THB}+µ_{2})$ where $(µ_{1},µ_{2})$ lies in region one.
* RDPPmodelnetwork1_step06_sim_ode45_THB2.m Simulate the THB dynamics from four different initial conditions of the RDPP model defined onnetworks $G_{la}(10,1)$ with $(b,a)=(b_{THB}+µ_{1},a_{THB}+µ_{2})$ where $(µ_{1},µ_{2})$ lies in region two.
* RDPPmodelnetwork1_step06_sim_ode45_THB3.m Simulate the THB dynamics from four different initial conditions of the RDPP model defined onnetworks $G_{la}(10,1)$ with $(b,a)=(b_{THB}+µ_{1},a_{THB}+µ_{2})$ where $(µ_{1},µ_{2})$ lies in region three.
* RDPPmodelnetwork1_step06_sim_ode45_THB4.m Simulate the THB dynamics from four different initial conditions of the RDPP model defined onnetworks $G_{la}(10,1)$ with $(b,a)=(b_{THB}+µ_{1},a_{THB}+µ_{2})$ where $(µ_{1},µ_{2})$ lies in region four.
* RDPPmodelnetwork1_step06_sim_ode45_THB5.m Simulate the THB dynamics from four different initial conditions of the RDPP model defined onnetworks $G_{la}(10,1)$ with $(b,a)=(b_{THB}+µ_{1},a_{THB}+µ_{2})$ where $(µ_{1},µ_{2})$ lies in region five.
* RDPPmodelnetwork1_step06_sim_ode45_THB6.m Simulate the THB dynamics from four different initial conditions of the RDPP model defined onnetworks $G_{la}(10,1)$ with $(b,a)=(b_{THB}+µ_{1},a_{THB}+µ_{2})$ where $(µ_{1},µ_{2})$ lies in region six.
* RDPPmodelnetwork1_step07_quickplot_n6.m Load data and then plot the evolution of the THB dynamics, as well the coresponding very long state of node i (with different colors) in $u_i − v_i$ phase plane of the RDPP model defined on $G_{la}(10,1)$ with 
(µ1, µ2) locating in each region.

## gwRDPPmodelnetwork2
This subfolder provides a series of codes to simulation the Hopf, steady state and Turing-Hopf bifurcation dynamics of the reaction-diffusion RDPP model defined on networks $G_{sw}(10,1,0.05)$, with step by step. 
The function of each code is the same as the corresponding one in subfolder gwRDPPmodelnetwork1.

## gwRDPPmodelnetwork3
This subfolder provides a series of codes to simulation the Hopf, steady state and Turing-Hopf bifurcation dynamics of the reaction-diffusion RDPP model defined on networks $G_{sw}(10,2,0.05)$, with step by step. 
The function of each code is the same as the corresponding one in subfolder gwRDPPmodelnetwork1.


## gwRDPPmodelnetwork3
This subfolder provides a series of codes to simulation the Hopf, steady state and Turing-Hopf bifurcation dynamics of the reaction-diffusion FitzHugh-Nagumo (FHN) model defined on weighted brain networks, with step by step. 
The function of each code is the same as the corresponding one in subfolder gwRDPPmodelnetwork1.


## gwRDPPmodelnetwork3
This subfolder provides a series of codes to simulation the Hopf, steady state and Turing-Hopf bifurcation dynamics of the reaction-diffusion FHN model defined on networks $G_{la}(10,1)$, with step by step. 
The function of each code is the same as the corresponding one in subfolder gwRDPPmodelnetwork1.


