%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [A] = add_noise_to_geometric(met,D2)
%
% input: net - is a struct giving properties of the network
%           net.A - is the adjacency matrix
%           net.N - is the number of nodes
%           net.M is the number of links
%           net.geometry is the node locations on the manifold
%        type - how to generate the noisy links
%        param - parameters for the noisy links
%
% output: A - adjacency matrix of noisy links that doesn't overlap with the
%         given network with adjacency matrix net.A, i.e., A.*net.A = 0
% 
% DRT 1-7-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [A] = add_noise_to_geometric(net,type,param)

   switch type
      case 'ER_like'% distribution of D2*N/2 links so each node has expected D2 non-geometric links
       
         D2 = param(1);%non-geometric nodal degree
         
         net.N = net.N;
         M = net.N*D2; % 2 times number of links to distribute
         
         %build network
         A = spalloc(net.N,net.N,M);
         
         edges_build=0; % initialize
         
         while edges_build < M/2
           nodestoconnect=randi(net.N,2,1);
           if nodestoconnect(1)==nodestoconnect(2)
             % do not allow self loops
             edges_build=edges_build;
           elseif net.A(nodestoconnect(1),nodestoconnect(2))==1 | A(nodestoconnect(1),nodestoconnect(2))==1
             % do not create allready existing links
             edges_build=edges_build;
           else
             A(nodestoconnect(1),nodestoconnect(2))=1;
             A(nodestoconnect(2),nodestoconnect(1))=1;
             edges_build=edges_build+1;
         end
         
         %A = A + A';%make links bidirectional
      end
         
      case 'k_regular'%each node has exactly D2 non-geometric links

         D2 = param(1);%non-geometric nodal degree

         net.N = net.N;
         M = net.N*D2;

         flag2 = 1;

         while(flag2)
            flag2 = 0;
            %build stubs
            stubs = zeros(M,1);
            for n=1:net.N;
               stubs((n-1)*D2 + (1:D2)) = n*ones(D2,1);
            end

            %build undirected link list
            LinkList = zeros(M/2,2); 
            for m=1:(M/2)

               flag = 1;%set flag to enter while loop
               count =0;
               while(flag)
                  flag = 0;%turn off flag to leave loop
                  %choose 2 random stubs
                  RAND = randi(length(stubs),1,2);
                  node_A = stubs(RAND(1));
                  node_B = stubs(RAND(2));

                  %prevent self links
                  if node_A == node_B
                     flag = 1;%try again
                  end

                  %don't make a link that already exists
                  for n = 1:(m-1)
                     if LinkList(n,1) == node_A && LinkList(n,2) == node_B
                        flag = 1;%try again
                     end
                     if LinkList(n,1) == node_B && LinkList(n,2) == node_A
                        flag = 1;%try again
                     end
                     if net.A(node_A,node_B) || net.A(node_B,node_A)
                        flag = 1;%try again
                     end
                  end
                  count = count+1;
                  if count > (M)
                     flag2=1;
                     break;
                  end
               end

               %make link 
               LinkList(m,1) = node_A;
               LinkList(m,2) = node_B;

               %remove stubs from list
               stubs(RAND(:)) = [];

            end

         end

         %build network
         A = spalloc(net.N,net.N,M);
         for m=1:(M/2)
             A(LinkList(m,1),LinkList(m,2)) = 1;
         end
         A = A + A';%make links bidirectional   

   end
end


