#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import numpy as np
#from scipy import sparse
#from math import sqrt
#import matplotlib.pyplot as plt
#import networkx as nx
#import random
#import gudhi


# In[ ]:


class geometric_network():
    """
    Geometric Network Object to run complex contagions on.
        
    Attributes
    ----------
    geometric_network.N: int
        size, number of nodes in the network.
    geometric_network.M: int
        total number of edges in the network.
    geometric_network.graph: a Networkx object
        Networkx graph corresponding to the Geometric object. Having this attribute, you can use all the networkx library.
    geometric_network.pos: dict
        A dictionary of nodes and their spatial location.
    geometric_network.A: A Scipy sparse matrix
        Adjacency matrix of the graph. Use .todense() or .toarray() to manipulate.
    geometric_network.text: str
        a simple description of the network.
    
    Parameters
    -----------
    network_type: str
        type of the network to be created. It can be '2D_lattice' or 'ring_lattice'
    **kwargs: 
        tiling: int
            should be provided if the network type is 2D_lattice.
            tiling of the 2d lattice. It can be 3,4,6 for now. this is the number of neighbors to be 
            connected.
    **kwargs:
        periodic: bool
            should be provided if the network type is 2D_lattice. if True, edges of the network 
            are going to be glued together.
    **kwargs:
        banded: bool
            should be provided if the network type is ring_lattice. If True, the closest band_length many 
            neigbors from right and left is going to be connected to every node, creating a banding.
    **kwargs:
        band_length: int
            should be provided if the network type is ring_lattice. geometric degre divided by 2.
            Note that geometric degree must be an even number.
    size: int
        size of the network to be initated. If 2D_lattice, there will be size**2 many total nodes.
    
    """
    def __init__(self, network_type, size, **kwargs):

        if network_type == '2D_lattice': 
            tiling = kwargs['tiling']
            per = kwargs['periodic']
            if tiling == 3: 
                self.graph = nx.triangular_lattice_graph(size, size, periodic = per, with_positions = True)
                self.pos = nx.get_node_attributes(self.graph,'pos')
                self.M = len(self.graph.edges())
                self.N = len(self.graph.nodes())
                

            elif tiling == 4: 
                self.graph = nx.grid_2d_graph(size, size, periodic = per)
                self.pos = dict( (n, n) for n in self.graph.nodes() )
                self.labels = dict( ((i, j), i * size + j) for i, j in self.graph.nodes() )
                self.M = len(self.graph.edges())
                self.N = len(self.graph.nodes())
                
            elif tiling == 6: 
                self.graph = nx.hexagonal_lattice_graph(size, size, periodic = per, with_positions = True)
                self.pos = nx.get_node_attributes(self.graph,'pos')
                self.M = len(self.graph.edges())
                self.N = len(self.graph.nodes())
            
            
        elif network_type == 'ring_lattice':# TODO: banding for every node
            self.graph = nx.cycle_graph(size)
            theta = (2*np.pi)/size
            self.pos = dict((i,(np.sin(theta*i),np.cos(theta*i))) for i in range(size))
            self.M = len(self.graph.edges())
            self.N = len(self.graph.nodes())
            self.text = 'Ring Lattice'
            if kwargs['banded']:
                if kwargs['band_length'] >= int(self.N/2)-1: 
                    raise ValueError('Band length cannot exceed the half of the size of the network')
                if kwargs['band_length'] <2: 
                    raise ValueError('Band length should be a positive integer greater 1 since the closest neighbors are already connected')
                for u in range(self.N):
                    for i in range(2,kwargs['band_length']+1):
                        # ranges from 2 to k+2 to avoid the closest node and start
                        ## banding from the second closest node
                        if u + i >= self.N: v = u + i - self.N
                        else: v = u + i
                        self.graph.add_edge(u, v)
                        if u - i < 0: v = self.N + u - i
                        else: v = u - i
                        self.graph.add_edge(u, v)
                self.text = self.text + ' w/ bandlength %d'%kwargs['band_length']
            else:self.text = self.text + ' w/ bandlength 0'
                        
        else: raise ValueError('network type can be a lattice or a ring')
            
        self.A = nx.adjacency_matrix(self.graph)
        ## use .toarray() or .todense() attributes to call the sparse matrix
              
    def display(self, n_size = 15, labels = True):
        """
        Method to pass parameters into nx.draw().
        
        Parameters
        -----------
        n_size: int
            node sizes.
        labels: bool
            node labels.
            
        Returns
        --------
        nx.draw(self.graph)  
            
        """
        
        try: return(nx.draw(self.graph, 
                            pos = self.pos, 
                            node_size = n_size, 
                            with_labels = labels, 
                            font_weight = 'bold'))
        
        except: return(nx.draw(self.graph, with_labels = labels, font_weight = 'bold'))
        
    def add_noise_to_geometric(self, noise_type, d2):
        """
        This method adds non-geometric edges to the network that are long range. Depending on the 'noise_type'
        the way we add these long range edges differ. If noise_type = ER_like, then there will be d2 many
        non geometric edges ON AVARAGE for every node. When the noise_type = k_regular, every node will 
        have exactly d2 many long range edges.
        
        Parameters
        ----------
        noise_type: str
            ER_like or k_regular
        d2: int
            degree to assign non-geometric edges to every node
            
        Returns
        -------
        None. Updates the geometric_network.A
        
        Raises
        ------
        ValueError
            if the geometric_network.N * d2 is an odd number.
            
        """
        
        if noise_type == 'ER_like':   
            
            M = int(self.N * d2)
        
            if M%2 == 1: raise ValueError('Try providing an even non-geometric degree')
                
            edges_built = 0
            
            while edges_built <= M/2:
                
                nodes_to_connect = [random.randint(0,self.N-1) for i in range(2)]
                
                if nodes_to_connect[0] == nodes_to_connect[1]: #no self-loops
                    edges_built = edges_built
                elif self.A.toarray()[nodes_to_connect[0]][nodes_to_connect[1]] == 1 or self.A.toarray()[nodes_to_connect[1]][nodes_to_connect[0]] == 1: # no existing loops
                    edges_built = edges_built
                else: 
                    self.graph.add_edge(nodes_to_connect[0],nodes_to_connect[1])
                    edges_built = edges_built + 1
            self.text = self.text + ' & ER-like noise w/ degree %d'%d2

                    
        elif noise_type == 'k_regular':
            
            self.text = self.text + ' & k-regular noise w/ degree %d'%d2
            
            M = int(self.N * d2)
        
            if M%2 == 1: raise ValueError('Try providing an even non-geometric degree')
            
            flag_2 = True
            
            while flag_2:
                flag_2 = False
                #build stubs
                stubs = np.zeros(M)
                for i in range(self.N):
                    index = (i*d2) + np.arange(d2)
                    stubs[index[0]:index[-1]+1] = (i) * np.ones(d2)
                    
                #build undirected link list
                link_list = np.zeros((int(M/2),2))
                for m in range(int(M/2)):
                    flag_1 = True # turn on flag to enter while loop
                    count = 0
                    while flag_1:
                        flag_1 = False #turn off flag to exit while loop
                        rand = [random.randint(0,len(stubs)-1) for i in range(2)]
                        
                        node_A = int(stubs[rand[0]])
                        node_B = int(stubs[rand[1]])
                                                
                        if node_A == node_B: flag_1 = True
                        
                        for n in range(m):
                            if link_list[n,0] == node_A and link_list[n,1] == node_B:
                                flag_1 = True
                            if link_list[n,0] == node_B and link_list[n,1] == node_A:
                                flag_1 = True
                            if self.A.toarray()[node_A][node_B] == 1 or self.A.toarray()[node_B][node_A] == 1:
                                flag_1 = True
                                
                        count = count +1
                        
                        if count > M: flag_2 = True ; break
                            
                    #make link       
                    link_list[m,0] = node_A
                    link_list[m,1] = node_B
                    
                    #remove stubs from list
                    stubs = np.delete(stubs,[rand[0],rand[1]])
            
            #build network
            for m in range(int(M/2)):
                self.graph.add_edge(link_list[m,0],link_list[m,1])
            
        else: raise InputError('Noise type can only be ER_like or k_regular')
        
        self.A = nx.adjacency_matrix(self.graph)
        
    def excitation(self, T, C, seed, threshold, refractory = False, ax = None, spy = False):
        """
        THE CORE FUNCTION OF THE NEURONAL CONTAGION MODEL. 
        In this model, a neuron fires if the ratio of it's excited neighbors to the total number of neighbors
        is greater than the threshold. Let's call the difference between this ratio and the threshold = F 
        so that if F is positive, neuron is going to fire and it doesn't fire when it's negative. We add some 
        stocasticity to the model by defining the sigmoid function so that the probability that the 
        neuron is going to fire is not a step function, but a sigmoid function.
        
        Parameters
        ----------
        T: int
            Number of time steps contagions is going to be iterated.
        C: int
            A positive constant for the sigmoid function, if C is too large(>100),
            jump from 0 to 1 is gonna be too quick i.e. model is going to be deterministic.
        seed: int 
            node id to start the contagion, in the first time step, we infect
            the neighbors of the seed with probablity 1 then enter the while loop below
        threshold: float
            threshold to compare for a neuron's neighbor input. threshold must be in (0,1).
        refractory: bool
            if TRUE, sets the refractory period of 1 time step i.e. neuron cannot fire for 1 time step 
            right after it fires. if FALSE, neuron stays active once its activated.
        ax, spy: matplotlib.axis, bool
            if spy is TRUE, there have to be an axis provided to plot the contagion spread.
            
        Returns
        -------
        activation_times: array
            An array of n x 1 keeping track of the first time step the corresponding
            node gets activated.
        size_of_contagion: array  
            An array of (T+1) x 1 keeping track of the number of active nodes at a
            given time(at t = 0, all neighbors of the seed is active)
        """
        
        n = self.N
        x_1 = sparse.lil_matrix(((n,T+1)), dtype = bool)
        activation_times = np.ones(n, dtype = int)*T
        
        excited_nodes = self.A[seed].nonzero()[1]
        
        deg = len(excited_nodes)
        
        size_of_contagion = np.zeros(T+1)
        
        vals = [True] * deg
        
        t = 0
        
        for i in range(deg): 
            x_1[excited_nodes[i],t] = vals[i]
            
        all_nodes = set([k for k in range(n)])
        size_of_contagion[t] = len(excited_nodes)
        
        while t<T and len(excited_nodes)>0 and np.any(activation_times==T):
            x_old = x_1.toarray()[:,t]
            excited_nodes = np.nonzero(x_old)[0]         
            activation_times[excited_nodes] = np.minimum(activation_times[excited_nodes], 
                                                         np.array([t]*len(excited_nodes)))
            activation_times = activation_times.astype(np.int32)
            t+=1
            for j in list(all_nodes-set(excited_nodes)): 
                ids = np.nonzero(self.A[j])[1]
                F = np.sum(x_old[ids])/len(ids) - threshold 
                Z = 1/(1+np.exp(-C*F)) 
                rand = random.uniform(0,1)
                x_1[j,t] = bool(rand<Z) 
                
            if not refractory: x_1[x_old,t] = True 
            
            size_of_contagion[t] = np.sum(x_1[:,t])
      
        if spy: ax.imshow(x_1.todense())
        return(activation_times, size_of_contagion)
    
    
    def run_excitation(self, Trials, T, C, seed, threshold, refractory, ax = None):
        """
        Helper function to run the excitation over several trials.
        
        Parameters
        ----------
        Trials: int
            Number of trials to run the contagion on the same network.
        T: int
            Number of time steps contagions is going to be iterated.
        C: int
            A positive constant for the sigmoid function, if C is too large(>100),
            jump from 0 to 1 is gonna be too quick i.e. model is going to be deterministic.
        seed: int 
            node id to start the contagion, in the first time step, we infect
            the neighbors of the seed with probablity 1 then enter the while loop below
        threshold: float
            threshold to compare for a neuron's neighbor input. threshold must be in (0,1).
        refractory: bool
            if TRUE, sets the refractory period of 1 time step i.e. neuron cannot fire for 1 time step 
            right after it fires. if FALSE, neuron stays active once its activated.
        ax: matplotlib.axis, optional
            if not provided, contagion will not be shown.
            
        Returns
        -------
        activation_times: array
            An array of n x Trials keeping track of the first time step the corresponding
            node gets activated for each trial.
        size_of_contagion: array  
            An array of (T+1) x Trials keeping track of the number of active nodes at a
            given time(at t = 0, all neighbors of the seed is active) for each trial.
        
        """
        
        n = self.N
        
        first_excitation_times = np.zeros((n,Trials))
        size_of_contagion = np.zeros((T+1,Trials))
        
        if ax is not None: 
            for i in range(Trials):
    
                first_exct, contagion_size = self.excitation(T = T, 
                                                             C = C, 
                                                             seed = seed, 
                                                             threshold = threshold, 
                                                             refractory = refractory,
                                                             ax = ax[i], 
                                                             spy = True)
                first_excitation_times[:,i] = first_exct
                size_of_contagion[:,i] = contagion_size
                ax[i].set_title('Trial %d for %s'%(i+1, self.text), fontsize = 25)
                ax[i].set_xlabel('Time', fontsize = 35)
                ax[i].set_ylabel('Node ID',fontsize = 35)
                ax[i].set_xticks([5*t for t in range(int(T/5)+1)])
                ax[i].set_yticks([10*t for t in range(int(n/10))])
                ax[i].tick_params('both', which = 'major', labelsize = 15)
            plt.tight_layout()
            
        else:
            for i in range(Trials):
    
                first_exct, contagion_size = self.excitation(T = T, 
                                                             C = C, 
                                                             seed = seed, 
                                                             threshold = threshold, 
                                                             refractory = refractory)
                first_excitation_times[:,i] = first_exct
                size_of_contagion[:,i] = contagion_size
            
        return(first_excitation_times, size_of_contagion)
    

    def spy_first_activation(self, first_activation_times):
        """
        Helper function to visualize the first activation times.
        
        Parameters
        ----------
        first_activation_times: array of size n x Trials
            First output of the run_excitation showing the first time step that the contagion 
            reaches to a given node
        
        """
        
        fig,ax = plt.subplots(1,1,figsize = (10,8))
        
        pos = ax.imshow(first_activation_times, 
                        interpolation = 'nearest', 
                        aspect = 'auto', 
                        cmap = 'viridis', 
                        extent = [0.5,Trials+0.5,0,n])
        ax.set_title('First Activation Times for %s'%self.text, fontsize = 20)
        ax.set_xlabel('Trials', fontsize = 15)
        ax.set_ylabel('Node ID', fontsize = 15)
        ax.set_xticks([i for i in range(1,Trials+1)])
        
        fig.colorbar(pos, ax = ax)
        
    def average_over_trials(self, matrix): 
        """
        Helper funtion to take the averages over trials of the given matrix.
        
        Parameters
        ----------
        matrix: array k x Trials
            Matrix to take the average over trials. Matrix have to be k x Trials. This can be 
            the size of the contagion or first activation times depending on what you need.
            
        Returns
        -------
        mean_matrix: array k x 1
            Mean matrix
        
        """
        return(np.mean(matrix, axis = 1))
    
    def make_distance_matrix(self, T, C, threshold, Trials, refractory, spy_distance = False):
        """
        A shortcut to run all of the above functions in one function. This creates a distance matrix by running the
        contagion on starting from every node and encoding the first activation times of each node. The (i,j) entry
        corresponds to the average time(over the trials) that a contagion reaches node j starting from node i.
        
        Parameters
        ----------
        T: int
            Number of time steps contagions is going to be iterated.
        C: int
            A positive constant for the sigmoid function, if C is too large(>100),
            jump from 0 to 1 is gonna be too quick i.e. model is going to be deterministic.
        threshold: float
            threshold to compare for a neuron's neighbor input. threshold must be in (0,1).
        Trials: int
            Number of trials to run the contagion on the same network.
        refractory: bool
            if TRUE, sets the refractory period of 1 time step i.e. neuron cannot fire for 1 time step 
            right after it fires. if FALSE, neuron stays active once its activated.
        spy_distance: bool
            Check True if you want to take a peak at the distance matrix.
            
        Returns
        -------
        D: n x n array
            Distance matrix
        Q: n x T+1 array
            Array carrying the size of the contagion at every time step.
        
        """
        n = self.N
        D = np.zeros((n,n))
        Q = np.zeros((n,T+1))
        
        for i in range(n):
            fat, qs = self.run_excitation(Trials, T, C, i, threshold, refractory)
            D[i], Q[i] = self.average_over_trials(fat), self.average_over_trials(qs)
        
        if spy_distance: 
            fig,ax = plt.subplots(1,1, figsize = (15*int(n/20),10*int(n/20)))
            pos = ax.imshow(D, 
                            origin = 'lower', 
                            interpolation = 'nearest', 
                            aspect = 'auto', 
                            cmap = 'viridis', 
                            extent = [-0.5,n-0.5,-0.5,n-0.5])
            ax.set_title('Distance Matrix for %s'%self.text, fontsize = 20)
            ax.set_ylabel('Excitation starting from node i', fontsize = 15)
            ax.set_xlabel('First time node j gets activated', fontsize = 15)
            ax.set_xticks([i*5 for i in range(int(n/5))])
            ax.set_yticks([i*5 for i in range(int(n/5))])

            fig.colorbar(pos, ax = ax)
        return(D,Q)
    
    def compute_persistence(self, distances, spy = False): 
        """
        Helper to compute persistent homology using the distance matrix by building a Rips filtration up to 
        dimension 2(topological features to be observed are going to be 1 dimensional at max).
        First normalizes the distances before the computation.
        
        Parameters
        ----------
        distances: n x n array
            distance matrix. First output of the make_distance_matrix.
        spy: bool, optional
            Take a peak at the persistence diagram
        Returns
        -------
        diag: list
            The birth and death times of the topological features in all dimensions.
        """
        rips_complex = gudhi.RipsComplex(distance_matrix = distances/np.max(distances), max_edge_length = 1)
        simplex_tree = rips_complex.create_simplex_tree(max_dimension = 2)
        diag = simplex_tree.persistence(min_persistence = 0.0)
        if spy: 
            gudhi.plot_persistence_diagram(diag)
        return(diag)
    
    def one_d_Delta(self, persistence):
        """
        Helper to compute the specific topological features.
        
        Parameters
        ----------
        persistences: list
             a list of birth and death times of the topological features or the output of the compute_persistence.
        Returns
        -------
        Delta_min: float
            The difference between the life times of the longest and the second longest 1-cycle.
        Delta_max: float
            The difference between the life times of the longest and the shorthes 1-cycle.
        Delta_avg: float
            The average lifetime of all the 1-cycles.
        """
        oned_holes = [(0,0)]
        for i in range(len(persistence)):
            if persistence[i][0] == 1:
                oned_holes.append(persistence[i][1])
        oned_holes = np.array(oned_holes)
        persistence_life_times = oned_holes[:,1]-oned_holes[:,0]
        Delta_min = np.sort(persistence_life_times)[-1]-np.sort(persistence_life_times)[-2]
        Delta_max = np.sort(persistence_life_times)[-1]-np.sort(persistence_life_times)[1]
        Delta_avg = np.mean(persistence_life_times[1:])
        return(Delta_min, Delta_max, Delta_avg)
    
    def display_comm_sizes(self, Q, labels):
        """
        Helper to visualize the size of the active nodes during the contagion. Shades are indicating the max 
        and min values of the spread starting from different nodes.
    
        Parameters
        ----------
        Q: list, [n x T+1 array]
            Output of the make_distance_matrix appended in a list
        labels: figure labels corresponding to every list element, threshold, network type, C etc...

        """
        fig,ax = plt.subplots(1,1, figsize =(20,10))
        for i in range(len(Q)):
            Q_mean = np.mean(Q[i], axis = 0)
            ax.plot(Q_mean[:int(np.max(np.nonzero(Q)[2])+1)], 
                    label = '%s'%labels[i], marker = 'v')
        
            X = np.linspace(0,int(np.max(np.nonzero(Q)[2])),int(np.max(np.nonzero(Q)[2])+1))
            ax.fill_between(X, 
                            np.max(Q[i], axis = 0)[:int(np.max(np.nonzero(Q)[2])+1)], 
                            np.min(Q[i], axis = 0)[:int(np.max(np.nonzero(Q)[2])+1)], 
                            alpha = 0.2)
        fig.suptitle('%s'%self.text)
        ax.set_title('Size of the Active Nodes', fontsize = 25)
        ax.set_xlabel('Time', fontsize = 20)
        ax.set_ylabel('Number of Nodes', fontsize = 20)
        ax.legend(fontsize = 'x-large')


# In[ ]:





# In[ ]:





# In[ ]:




