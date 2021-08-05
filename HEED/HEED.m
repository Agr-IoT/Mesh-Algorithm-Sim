
clear;

load('XR.mat'); % 1000 nodes set
load('YR.mat');

load('X.mat'); % 300 nodes set
load('Y.mat');

Extra = 0; %this is used when want to put extral nodes in a specific area
NodeNums = length(X)+Extra; %the num of node  
AreaR = 100*100 ;   %the area of sensing field 
Bx=50;  %The Postion of Base Station  
By=50;


MaxInteral = 200; % the simulate time (number of rounds)

for i=1:1:NodeNums-Extra
    Node.x(i)=X(i);
    Node.y(i)=Y(i);
end

% the init energy of all node 
InitEn=0.5;
% the transit Radius 
% transmission range  ------> determin how many neighbours a sensor can have 
Tr=25;  
a=1;
Elec=50 * 10^(-9); % used for energy consumption for transmission/receiving 
Efs=10*10^(-12);
Eamp=0.0013*10^(-12);
MaxEn=InitEn*(1+a); 
Efusion= 5 * 10^(-9);

%data packet size
DataKbit = 100 * 8;
BroadcastKbit = 25 * 8 *2;
DataHeader = 25 * 8;

%CH assembling rate
Gathingcoefficient=0.8;  

TDMA=5; 
Cprob=0.05; % the desired percentage of cluster heads
Pmin = 1*10.^(-4);

NON_CH					= 0;%			// non cluster head 
TENTATIVE_CH			= 1; %			// tentative cluster head				 
FINAL_CH				= 2;%				// final cluster head 

sym ClusterHeadNum ; 
ClusterHeadNum=linspace(0,0,MaxInteral); 

TOS_LOCAL_ADDRESS = -1;       % TOS_LOCAL_ADDRESS  must <=0 

 A=linspace(0,0,MaxInteral);  % Track network coverage
     
Node.IsClusterHeads=linspace(0,0,NodeNums); % NON_CH,TENTATIVE_CH,FINAL_CH 
Node.IsCovered=linspace(0,0,NodeNums);      % Have Been Covered by a cluster head 1:yes 0:No 
Node.IsCoveredByFinal= linspace(0,0,NodeNums);


Node.c=linspace(0,0,NodeNums);              % the Cluster head of node, node number
Node.chcost=linspace(0,0,NodeNums);         % the Cluster head of node 
Node.d=linspace(0,0,NodeNums);              % the distance between cluster head and node 
%Node.l=zeros(1,NodeNums)+Kbit;             % the length of node i transmit packet 
Node.EnNode=zeros(1,NodeNums)+InitEn;       % the init energy of all node 

Node_Energy = zeros(MaxInteral,NodeNums);

Node.StateNode=ones(1,NodeNums);            % the State of all node 1: alive 0:dead 
Node.Listothernode=zeros(NodeNums);         % if node is a cluster head,Listothernode save the id of node belong to this cluster        
Node.csize=linspace(0,0,NodeNums);          % cluser size ,each cluster node num 

Node.Nbr=zeros(NodeNums);                   % neighbor of node 
Node.NumNbr=linspace(0,0,NodeNums);         % the neighbor's num of node 
%Node.DistNbr=linspace(0,0,NodeNums);         % the neighbor's dist of node 

Node.CHprob=zeros(1,NodeNums)+Cprob;  
Node.InitCHprob=zeros(1,NodeNums); 

Node.tent_CH=zeros(1,NodeNums)+NON_CH;  %?
Node.tent_CH_Cost=ones(1,NodeNums)+9999; %?

Node.IsaddDummyRound=linspace(0,0,NodeNums); 
Node.n_finalCH=linspace(0,0,NodeNums);  
Node.ListfinalCH=zeros(NodeNums); 
Node.ListfinalCH_Cost=zeros(NodeNums)+9999; % the cost is set to the maximum
Node.n_tentCH=linspace(0,0,NodeNums);  
Node.ListtentCH=zeros(NodeNums); 
Node.ListtentCH_Cost=zeros(NodeNums)+9999; 
Node.my_finalCH=linspace(0,0,NodeNums); 
Node.my_tentCH=linspace(0,0,NodeNums); 
Node.my_final_CH_Cost=ones(1,NodeNums)+9999; 
Node.Isstop=ones(1,NodeNums);    % clustering is end ? 1:no,0:yes 



for i=1:(MaxInteral) 
    AliveNode(i)=0;%NodeNums; 
end

for i=1:NodeNums 
    count =0 ; 
    for j=1:NodeNums 
        
        if(j~=i)  
        dist = ((Node.x(i)-Node.x(j)).^2)+((Node.y(i)-Node.y(j)).^2);  % the distance.^2 
               if dist < (Tr.^2) % original dist < (Tr.^2)
                   count=count+1;  % count how many neighbours one node has
                   Node.Nbr(i,count)=j; 
                    
               end 
                   
         end 
         if j== NodeNums  
                Node.NumNbr(i) = count ; 
         end   
    end  
 end 
   
 syms filen strnumnode tpye strround; 
     strnumnode = int2str(NodeNums);   
 sym iteration; 

 
 for Rounds = 1:MaxInteral   
      
      % the Setup phase of cluster 
      % beginning of a round
      Node.CHprob=max(Cprob.*((Node.EnNode)./MaxEn), Pmin); 
      Node.InitCHprob=Node.CHprob; 
     
      Node.IsaddDummyRound=0.*(Node.IsaddDummyRound);  %clear 
      Node.n_finalCH=linspace(0,0,NodeNums); %Node.n_finalCH-Node.n_finalCH;  
      Node.IsCovered=linspace(0,0,NodeNums);   
     
      Node.ListfinalCH=Node.ListfinalCH-Node.ListfinalCH;
      Node.ListfinalCH_Cost=Node.ListfinalCH_Cost-Node.ListfinalCH_Cost+9999; 
      Node.my_finalCH=Node.my_finalCH-Node.my_finalCH; 
      Node.my_final_CH_Cost=(Node.my_final_CH_Cost-Node.my_final_CH_Cost)+9999; 
      Node.n_tentCH=Node.n_tentCH-Node.n_tentCH;  
      Node.ListtentCH=Node.ListtentCH-Node.ListtentCH; 
      Node.ListtentCH_Cost=Node.ListtentCH_Cost-Node.ListtentCH_Cost+9999; 
      Node.csize=Node.csize-Node.csize;
      Node.Isstop = Node.StateNode ;
      % Node.Isstop is different from Node.StateNode
      Node.tent_CH=Node.tent_CH-Node.tent_CH+NON_CH; 
      Node.tent_CH_Cost=Node.tent_CH_Cost-Node.tent_CH_Cost+9999; 
     
      Node.c=Node.c-Node.c; 
      Node.d=Node.d-Node.d; 
      %ClusterHeadNum=0; 
      Node.IsClusterHeads=Node.IsClusterHeads-Node.IsClusterHeads+NON_CH; 
      iteration(Rounds)=0; 
      
      %update alive nodes
      for i =1:NodeNums  
         if Node.StateNode(i)==1
             AliveNode(Rounds)=AliveNode(Rounds)+1;% calculate the number of alive nodes for this interation
             Node_Energy(Rounds,i) = Node.EnNode(i);
         end
          
      end
      
      %update neighbours
      Node.Nbr=zeros(NodeNums);                    
      Node.NumNbr=linspace(0,0,NodeNums);
      for i=1:NodeNums 
        count =0 ; 
        if(Node.StateNode(i)==1)
        for j=1:NodeNums 
            
            if(j~=i && Node.StateNode(j)==1)  
            dist = ((Node.x(i)-Node.x(j)).^2)+((Node.y(i)-Node.y(j)).^2);  % the distance.^2 
                   if dist < (Tr.^2) % original dist < (Tr.^2)
                       count=count+1; 
                       Node.Nbr(i,count)=j; 

                   end 

             end 
             if j== NodeNums  
                    Node.NumNbr(i) = count ; 
             end   
        end  
        end
      end 
     % while sum(Node.CHprob)~=NodeNums
      while sum(Node.Isstop)~=0
         iteration(Rounds)=iteration(Rounds)+1;    
         
      for i =1:NodeNums          
       if Node.Isstop(i)==1 
         if Node.IsCovered(i) == 1 
           for j=1:Node.NumNbr(i)    
                if Node.IsClusterHeads(Node.Nbr(i,j)) ~= NON_CH  
                  if Node.my_final_CH_Cost(i) > Node.NumNbr(Node.Nbr(i,j));
                     Node.my_finalCH(i)= Node.Nbr(i,j); 
                     Node.my_final_CH_Cost(i)=Node.NumNbr(Node.Nbr(i,j)); 
                  end
                end   
           end
           if Node.my_finalCH(i) == i
             
               if Node.CHprob(i)==1
                 Node.IsClusterHeads(i)= FINAL_CH;  
                 Node.my_final_CH_Cost(i)= Node.NumNbr(i);%computeDegree(i);   
                 ClusterHeadNum(Rounds)=ClusterHeadNum(Rounds)+1;
                 %broadcast 
                 dist =Tr.^2;
                 EntranPCH=EnTran(Elec,Eamp,BroadcastKbit,dist) ;
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                 if Node.EnNode(i) <= 0 
                    Node.StateNode(i)=0; 
                    Node.EnNode(i)=0; 
                 end 
                 for k=1:Node.NumNbr(i)
                    Node.IsCovered(Node.Nbr(i,k)) = 1;
                    Node.IsCoveredByFinal(Node.Nbr(i,k)) = 1;
                    EnRecP=EnRec(Elec,BroadcastKbit); 
                    Node.EnNode(Node.Nbr(i,k)) = Node.EnNode(Node.Nbr(i,k))-EnRecP;
                    if Node.EnNode(Node.Nbr(i,k)) <= 0 
                       Node.StateNode(Node.Nbr(i,k))=0; 
                       Node.EnNode(Node.Nbr(i,k))=0; 
                    end 
                 end 
   
              else
                 Node.IsClusterHeads(i)=TENTATIVE_CH;  
                 Node.my_finalCH(i) = i;
                 Node.my_final_CH_Cost(i) = Node.NumNbr(i);
                 Node.c(i) =TOS_LOCAL_ADDRESS;
                 Node.tent_CH(i)=TOS_LOCAL_ADDRESS; 
                 Node.tent_CH_Cost(i)=Node.NumNbr(i); 
                 
                 %broadcast 
                 dist =Tr.^2;
                 EntranPCH=EnTran(Elec,Eamp,BroadcastKbit,dist) ;
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                 if Node.EnNode(i) <= 0 
                       Node.StateNode(i)=0; 
                       Node.EnNode(i)=0; 
                 end 
                 for k=1:Node.NumNbr(i)
                    Node.IsCovered(Node.Nbr(i,k)) = 1;
                    EnRecP=EnRec(Elec,BroadcastKbit); 
                    Node.EnNode(Node.Nbr(i,k)) = Node.EnNode(Node.Nbr(i,k))-EnRecP;
                    if Node.EnNode(Node.Nbr(i,k)) <= 0 
                       Node.StateNode(Node.Nbr(i,k))=0;  
                       Node.EnNode(Node.Nbr(i,k))=0; 
                    end 
                 end
               end
           end
         elseif Node.CHprob(i)==1
                 Node.IsClusterHeads(i)= FINAL_CH; 
                 Node.my_finalCH(i)=i; 
                 Node.my_final_CH_Cost(i)= Node.NumNbr(i);%computeDegree(i);  
                 Node.d(i)=((Node.x(i)-Bx).^2)+((Node.y(i)-By).^2);  % the distance.^2 
                 ClusterHeadNum(Rounds)=ClusterHeadNum(Rounds)+1;
                 %broadcast 
                 dist =Tr.^2;
                 EntranPCH=EnTran(Elec,Eamp,BroadcastKbit,dist) ;
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                 if Node.EnNode(i) <= 0 
                   Node.StateNode(i)=0; 
                   Node.EnNode(i)=0; 
                 end 
                 for k=1:Node.NumNbr(i)
                    Node.IsCovered(Node.Nbr(i,k)) = 1;
                    Node.IsCoveredByFinal(Node.Nbr(i,k)) = 1;
                    EnRecP=EnRec(Elec,BroadcastKbit); 
                    Node.EnNode(Node.Nbr(i,k)) = Node.EnNode(Node.Nbr(i,k))-EnRecP;
                    if Node.EnNode(Node.Nbr(i,k)) <= 0 
                       Node.StateNode(Node.Nbr(i,k))=0;  
                       Node.EnNode(Node.Nbr(i,k))=0; 
                    end 
                 end
         elseif rand(1,1)<Node.CHprob(i)
                 Node.IsClusterHeads(i)=TENTATIVE_CH;  
                 Node.my_finalCH(i) = i;
                 Node.my_final_CH_Cost(i) = Node.NumNbr(i);
                 Node.c(i) =TOS_LOCAL_ADDRESS;
                 Node.tent_CH(i)=TOS_LOCAL_ADDRESS; 
                 Node.tent_CH_Cost(i)=Node.NumNbr(i);
                 
                 %broadcast 
                 dist =Tr.^2;
                 EntranPCH=EnTran(Elec,Eamp,BroadcastKbit,dist) ;
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                 if Node.EnNode(i) <= 0 
                   Node.StateNode(i)=0; 
                   Node.EnNode(i)=0; 
                 end 
                 for k=1:Node.NumNbr(i)
                    Node.IsCovered(Node.Nbr(i,k)) = 1;
                    EnRecP=EnRec(Elec,BroadcastKbit); 
                    Node.EnNode(Node.Nbr(i,k)) = Node.EnNode(Node.Nbr(i,k))-EnRecP;
                    if Node.EnNode(Node.Nbr(i,k)) <= 0 
                       Node.StateNode(Node.Nbr(i,k))=0; 
                       Node.EnNode(Node.Nbr(i,k))=0; 
                    end 
                 end 
         end
         CHprevious = Node.CHprob(i);
         Node.CHprob(i)=min(Node.CHprob(i).*2,1); 
         if CHprevious ==1
             Node.Isstop(i) =0;
         end
       end
      end
      
      
      
      end
      

     
      %join a cluster
      for i=1:NodeNums        
        if Node.StateNode(i)==1      
         if Node.IsClusterHeads(i)~= FINAL_CH
          if Node.IsCoveredByFinal(i) == 1
           for j=1:Node.NumNbr(i)    
                if Node.IsClusterHeads(Node.Nbr(i,j)) == FINAL_CH % covered by final CH
                  if Node.my_final_CH_Cost(i) > Node.my_final_CH_Cost(Node.Nbr(i,j));
                     Node.my_finalCH(i)= Node.Nbr(i,j); 
                     Node.my_final_CH_Cost(i)=Node.my_final_CH_Cost(Node.Nbr(i,j));
                  end 
                else
                end 
           end
           dist =Tr.^2;
           EntranPCH=EnTran(Elec,Eamp,BroadcastKbit,dist) ;
           Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
           if Node.EnNode(i) <= 0 
             Node.StateNode(i)=0; 
             Node.EnNode(i)=0; 
           end
           
           %corresponding CH cost
           Node.csize(Node.my_finalCH(i)) = Node.csize(Node.my_finalCH(i))+1;
           EnRecP=EnRec(Elec,BroadcastKbit);
           Node.EnNode(Node.my_finalCH(i)) = Node.EnNode(Node.my_finalCH(i))-EnRecP;
           Node.csize(Node.my_finalCH(i))= Node.csize(Node.my_finalCH(i)) +1;
           if Node.EnNode(Node.my_finalCH(i)) <= 0 
             Node.StateNode(Node.my_finalCH(i))=0; 
             Node.Isstop(Node.my_finalCH(i))=0; 
             Node.EnNode(Node.my_finalCH(i))=0; 
           end
        else 
                 Node.IsClusterHeads(i)= FINAL_CH; 
                 Node.my_finalCH(i)=i; 
                 Node.my_final_CH_Cost(i)= Node.NumNbr(i);%computeDegree(i); 
                 ClusterHeadNum(Rounds)=ClusterHeadNum(Rounds)+1;
                 for k=1:Node.NumNbr(i)
                    Node.IsCovered(Node.Nbr(i,k)) = 1; 
                    Node.IsCoveredByFinal(Node.Nbr(i,k)) = 1;
                 end
                 
                 dist =Tr.^2;
                 EntranPCH=EnTran(Elec,Eamp,BroadcastKbit,dist) ;
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                 if Node.EnNode(i) <= 0 
                   Node.StateNode(i)=0; 
                   Node.EnNode(i)=0; 
                 end
                 
                 for k=1:Node.NumNbr(i)
                 Node.IsCovered(Node.Nbr(i,k)) = 1;
                 Node.IsCoveredByFinal(Node.Nbr(i,k)) = 1;
                 EnRecP=EnRec(Elec,Kbit); 
                 Node.EnNode(Node.Nbr(i,k)) = Node.EnNode(Node.Nbr(i,k))-EnRecP;
                 if Node.EnNode(Node.Nbr(i,k)) <= 0 
                   Node.StateNode(Node.Nbr(i,k))=0; 
                   Node.EnNode(Node.Nbr(i,k))=0; 
                 end 
                 end
        end
       end
      end
      end 
      
    
     %% TDMA  and TPC 
     for i=1:NodeNums 
       if Node.StateNode(i) ~=0 
         if Node.IsClusterHeads(i) == NON_CH 
             dist = ((Node.x(i)-Node.x(Node.my_finalCH(i))).^2)+((Node.y(i)-Node.y(Node.my_finalCH(i))).^2);
             EntranPCH=EnTran(Elec,Eamp,DataKbit+DataHeader,dist) ; 
             Node.EnNode(i)=Node.EnNode(i)-(TDMA.*EntranPCH); 
             if Node.EnNode(i) <= 0 
                        Node.StateNode(i)=0; 
                        Node.EnNode(i)=0; 
                        
             end  
              EnRecP=EnRec(Elec,DataKbit+DataHeader); 
              Node.EnNode(Node.my_finalCH(i))=Node.EnNode(Node.my_finalCH(i))-TDMA.*EnRecP;  
                   if Node.EnNode(Node.my_finalCH(i)) <= 0 
                        Node.StateNode(Node.my_finalCH(i))=0; 
                        Node.EnNode(Node.my_finalCH(i))=0; 
                    end 
         else 
             dist = ((Node.x(i)-Bx).^2)+((Node.y(i)-By).^2);
             Ef = Efusion * DataKbit*Node.csize(i);
             EntranPCH=EnTran(Elec,Eamp,DataKbit+DataHeader,dist) ; %Node.csize(i)
             Node.EnNode(i)=Node.EnNode(i)-(TDMA.*EntranPCH); %problem the transmission power used from the CH to the BS
             %Node.EnNode(i)=Node.EnNode(i)-(TDMA.*EntranPCH);
             Node.EnNode(i)=Node.EnNode(i)-(TDMA.*Ef);
             % TDMA.*EntranPCH
             if(TDMA.*EntranPCH>1)
                 TDMA.*EntranPCH
             end
             if Node.EnNode(i) <= 0 
                        Node.StateNode(i)=0; 
                        Node.EnNode(i)=0; 
                        
             end  
         end 
      end  
     end
     
 %%%%%%%%%%%%%%%%EVALUATION METRIC: NETWORK COVERAGE%%%%%%%%%%%%%%%%%%%
 
 A(Rounds) = grid(Node);%calcuate the coverage of the network
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 end

Alive=0;
for i =1:NodeNums  
   if Node.StateNode(i)==1
      Alive=Alive+1;
   end
end
 



for i=1:NodeNums
    count =0 ; 
    if Node.StateNode(i) == 1
    for j=1:NodeNums 
        
         if(j~=i && Node.StateNode(j) == 1)  
         dist = ((Node.x(i)-Node.x(j)).^2)+((Node.y(i)-Node.y(j)).^2);  % the distance.^2 
               if dist < (Tr.^2) % original dist < (Tr.^2)
                   count=count+1; 
                   Node.Nbr(i,count)=j; 
                    
               end 
                   
         end 
         if j== NodeNums  
                Node.NumNbr(i) = count ; 
         end   
    end  
    end
end
 
  