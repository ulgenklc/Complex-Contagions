%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% net = make_geometric_net(type,SIZE,param)
%
% input: type - is the type of network to be made
%        SIZE - is the number of nodes
%        param - is a vector with the entries being used differently for each
%               type
%
% output: net - is a struct giving properties of the network
%           net.A - is the adjacency matrix
%           net.N - is the number of nodes
%           net.M is the number of links
%           net.geometry is the node locations on the manifold
% 
% DRT 1-7-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function net = make_geometric_net(type,SIZE,param)

net.N = SIZE;% network size
net.A = spalloc(net.N,net.N,10*net.N);
net.type = type;

 
switch net.type
      
   case '2D_lattice'
     
      %param(1) = number of neighbors: 4, 8, 12, 20
      %currently it is only coded up for 8
         
      switch param(1)
         case 8;      
            for n = 1:net.N               
               %handle four corners
               if n==1
                  ids = [n+1,n+sqrt(net.N),n+sqrt(net.N)+1];
               elseif n==sqrt(net.N)
                  ids = [n-1,n+sqrt(net.N)-1,n+sqrt(net.N)];
               elseif n==net.N-sqrt(net.N)+1
                  ids = [n-sqrt(net.N),n-sqrt(net.N)+1,n+1];
               elseif n==net.N 
                  ids = [n-sqrt(net.N)-1,n-sqrt(net.N),n-1];
                  
               %handle sides
               elseif n>1 && n<sqrt(net.N) %left
                  ids = [n-1,n+1,n+sqrt(net.N)-1,n+sqrt(net.N),n+sqrt(net.N)+1];
               elseif net.N>n && (net.N-n)<sqrt(net.N) %right
                  ids = [n-sqrt(net.N)-1,n-sqrt(net.N),n-sqrt(net.N)+1,n-1,n+1];
              
               elseif mod(n-1,sqrt(net.N))==0 %top
                  ids = [n-sqrt(net.N),n-sqrt(net.N)+1,n+1,n+sqrt(net.N),n+sqrt(net.N)+1];
               elseif mod(n,sqrt(net.N))==0%bottom
                  ids = [n-sqrt(net.N)-1,n-sqrt(net.N),n-1,n+sqrt(net.N)-1,n+sqrt(net.N)];               
               else
                  ids = n + [-sqrt(net.N)-1,-sqrt(net.N),-sqrt(net.N)+1,...
                           -1,1,sqrt(net.N)-1,sqrt(net.N),sqrt(net.N)+1];
               end
               net.A(n,ids) = ones(1,length(ids));               
            end
           % net.A=net.A+net.A';
      end               
      
      %node locations in space
      net.geometry = zeros(net.N,2);
      net.geometry(:,1) = repmat((1:sqrt(net.N))',sqrt(net.N),1);
      count = 0;
      for m=1:sqrt(net.N)
         net.geometry(count+(1:sqrt(net.N)),2) = ones(size(net.geometry(count+(1:sqrt(net.N)),2))) * m;
         count = count+sqrt(net.N);
      end
  
   
   case 'ring_lattice'
      
      for n=2:net.N
         net.A(n,n-1) = 1;
         net.A(n-1,n) = 1;
      end
      net.A(1,net.N) = 1;
      net.A(net.N,1) = 1;
      
      %node locations in space
      theta = (1:net.N)'/net.N*(2*pi);%angle      
      net.geometry(:,1) = [cos(theta),sin(theta)];
           
   
   case 'chain_lattice'
      for n=2:net.N
         net.A(n,n-1) = 1;
         net.A(n-1,n) = 1;
      end
      
      %node locations in space
      theta = (1:net.N)'/net.N*(2*pi);%angle      
      net.geometry(:,1) = [theta,zeros(size(theta))];

   case 'banded_chain_lattice'
      Neighbors = param(1);%must be a multiple of 2 
      for p = 1:Neighbors/2
         for n = (1+p):net.N
            net.A(n,n-p) = 1;
            net.A(n-p,n) = 1;
         end
      end
      
      %node locations in space
      theta = (1:net.N)'/net.N*(2*pi);%angle      
      net.geometry(:,1) = [theta,zeros(size(theta))];
      
   case 'banded_ring_lattice'
      Neighbors = param(1);%must be a multiple of 2 
      for p = 1:(Neighbors/2)
         for n = (1+p):net.N
            net.A(n,n-p) = 1;
            net.A(n-p,n) = 1;
         end
      end
      for p = 1:(Neighbors/2)
         for n = 1:((Neighbors/2)-p+1)
            net.A(p,net.N-n+1) = 1;
            net.A(net.N-n+1,p) = 1;
         end 
      end
      
      %node locations in space
      theta = (1:net.N)/net.N*(2*pi);%angle      
      net.geometry = [cos(theta'),sin(theta')];
      
end  

net.M = length(find(net.A));

%laplacian matrix
   %net.D = diag(sum(net.A,2));
   %net.L = net.D-net.A;
 
%order by degree   
   %d = sum(net.A);
   %[d,ids] = sort(d,'descend');
   %net.A = net.A(ids,ids);

end

