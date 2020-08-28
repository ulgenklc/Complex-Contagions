
%% Run Kanouchi-Kopelli dynamics

function excitationDYNAMICS()

clc;
clear;


   %% make noisy geometric ring network

   net.N=2000;
   D1=200; % geometic degree
   D2=1; % non-geometric degree
   net.A = spalloc(net.N,net.N,(D1+D2)*net.N);

   noise.type='k_regular'; % type of non-geometric links  
   net.type=['banded_ring_lattice']; 
 
   net = make_geometric_net(net.type,net.N,D1 );
   
   net.A_geo = net.A;%adjacency matrix for geometric network without noisy links
   net.A_nongeo = add_noise_to_geometric(net,noise.type,D2);%adjacency matrix for noisy links
   net.A = net.A_geo + net.A_nongeo;%noisy geometric network
   
   

   net.A = sparse(net.A);
   %spy(net.A)


   %% run dynamics

   p = 0.050;% excitation probability
   pstart = .5;
   seed = net.N/2;
   T = 50;      
   Trials = 1;
   
   
   Threshold = 1;
   first_activation_times1 = run_excitation_difference(net,seed,p,pstart,T,Threshold,Trials);
   Threshold = 2;
   first_activation_times2 = run_excitation_difference(net,seed,p,pstart,T,Threshold,Trials);
   Threshold = 3;
   first_activation_times3 = run_excitation_difference(net,seed,p,pstart,T,Threshold,Trials);
   Threshold = 4;
   first_activation_times4 = run_excitation_difference(net,seed,p,pstart,T,Threshold,Trials);

   
   
   figure;
   plot( (first_activation_times1')','r','displayname','threshold1');
   hold on;
   plot( (first_activation_times2')','b','displayname','threshold2');
   plot( (first_activation_times3')','k','displayname','threshold3');
   plot( (first_activation_times4')','y','displayname','threshold4');
   legend show
   
   save
   
end

function [first_activation_times] = run_excitation_difference(net,seed,p,pstart,T,Threshold,Trials)

   first_activation_times = zeros(net.N,Trials);

   for trial=1:Trials      
      x = spalloc(net.N,T+1,net.N);
      initialized_cluster = find(net.A(seed,:));
      %initial condition - stochastic initialization at a seed
      x(initialized_cluster,1) = rand(length(initialized_cluster),1) > (1-pstart);
      plot(x);



      for t=1:T
         x_old = x(:,t);      
         excited_nodes = find(x_old)';
         if sum(first_activation_times(:,trial)>0)<net.N
            for n = excited_nodes
               if first_activation_times(n,trial)==0
                  first_activation_times(n,trial)=t;
               end
            end
         end
         %x(t+1) = spalloc(net.N,1,1);
         for n=setdiff(1:net.N,excited_nodes)
            ids = find(net.A(n,:));
            x(n,t+1) = sum( (x_old(ids).*rand(length(ids),1)) > 1-p  ) >= Threshold;

         end

         %plot(x(:,t+1));


         

         %Sigma = x*x';
         %imagesc(cov(x'))
         imagesc(x)
         pause(.01)

      end
   end
end






