function clusters=orthoclust_2species(A2,A3,orthologs23,kappa,method)
%Syntax: clusters=orthoclust_2species(A2,A3,orthologs23,kappa,method)
%the output contains how genes are mapped to different clusters
%there are two ways to optimize the cost function of OrthoClust in the
%configuration space of size (N2+N3)^q, where q is the number of
%spins/labels allowed. one is based on simulated annealing, which could be
%very slow. the other is based on the Louvian algorithm (Blondel et.al.
%Journal of Statistical Mechanics: Theory and Experiment 2008), which is
%essentially a heuristic using the steepest descent approach. set method=1
%for steepest descent, and method=2 for simulated annealing. 

addpath(genpath('./3rd_party_code'));%add path to 3rd party code
N2=length(A2);
N3=length(A3);

fprintf('setup orthology relationship\n');
[big_couple,couple23]=generate_coupling_matrix(orthologs23,N2,N3);

if nargin<5 
    method=2;
end

if method==1
    clusters=steepest_descent_Potts_net_2species(A2,A3,big_couple,kappa);
else
    clusters=anneal_Potts_net_2species(A2,A3,orthologs23,couple23,kappa);
end

end

function [big_couple,couple12]=generate_coupling_matrix(orthologs23,N2,N3)
%transform the orthology relationships in a matrix form
O12=sparse(orthologs23(:,1),orthologs23(:,2),ones(size(orthologs23,1),1),N2,N3);
d1=full(sum(O12,2));
d2=full(sum(O12))';
o1=d1(orthologs23(:,1));
o2=d2(orthologs23(:,2));
couple12=(1./o1+1./o2)/2;
couple=sparse(orthologs23(:,1),orthologs23(:,2),couple12,N2,N3);
big_couple=[sparse(N2,N2) couple
    couple' sparse(N3,N3)];
end

function clusters = steepest_descent_Potts_net_2species(A2,A3,big_couple,kappa)

fprintf('initialize the steepest descent\n');
if sum(sum(A2<0))==0 || sum(sum(A3<0))==0
    clusters=louvain_unsigned(A2,A3,big_couple,kappa);
else
    clusters=louvain_signed(A2,A3,big_couple,kappa);
end

end

function [c,a2c,b2c]=my_union(a,b)
%[c a2c b2c]=my_union(a,b)
%given two arrays a, b, find the combine list of unique elts, then return
%how a and b are mapped to the combine list...
[au,~,Ja]=unique(a);
[bu,~,Jb]=unique(b);
c=union(a,b);

[~,ya,za]=intersect(c,au);
[~,yb,zb]=intersect(c,bu);

a2c=ya(Ja);
b2c=yb(Jb);
end

function E=get_energy_signed(A2,A3,M2P,M2N,M3P,M3N,K3SP,K3SN,K2SN,K2SP,sigma2,sigma3,big_couple,kappa)
%for a single network,
%E=-sum_{i~=j}(A_{ij}-\gamma*p_{ij})\delta(\sigma_i,\sigma_j), where
%p_{ij}=k_i*k_k/2M and 2M=sum(sum(A))
%for gamma=1, the energy function can be re-written as -sum_{i,j\in edges}
%\delta(\sigma_i,\sigma_j)+\gamma/2M*\sum_{s=1}^q Ks.^2/2, the
%implementation here is based on this equation
%a signed network is essentially the superposition of 2 networks, one for positive edges and one for negative edges, 
%the energy is the weighted sum (lambda and 1-lambda) of the 2 terms.
%the energy of the system is equal to the energy term from the first
%network, plus the energy term from the second network, plus the energy
%term from the coupling, orthologs having the same spin have lower energy.
lambda=1/2;
N2=length(A2);
N3=length(A3);

E2P=sum(K2SP.^2)/(2*M2P);
E2N=sum(K2SN.^2)/(2*M2N);

[u,v]=find(A2>0);
E2P=E2P-sum(sigma2(u)==sigma2(v));
[u,v]=find(A2<0);
E2N=E2N-sum(sigma2(u)==sigma2(v));

E3P=sum(K3SP.^2)/(2*M3P);
E3N=sum(K3SN.^2)/(2*M3N);

[u,v]=find(A3>0);
E3P=E3P-sum(sigma3(u)==sigma3(v));
[u,v]=find(A3<0);
E3N=E3N-sum(sigma3(u)==sigma3(v));

[I, J, L]=find(big_couple(1:N2,N2+1:N3+N2));
O12=-sum(L(sigma2(I)==sigma3(J)));

E=lambda*E2P+lambda*E3P+O12*kappa-(1-lambda)*E2N-(1-lambda)*E3N;
end

function E=get_energy_unsigned(A1,A2,M1,M2,K1S,K2S,sigma1,sigma2,big_couple,kappa)
%for a single network,
%E=-sum_{i~=j}(A_{ij}-\gamma*p_{ij})\delta(\sigma_i,\sigma_j), where
%p_{ij}=k_i*k_k/2M and 2M=sum(sum(A))
%for gamma=1, the energy function can be re-written as -sum_{i,j\in edges}
%\delta(\sigma_i,\sigma_j)+\gamma/2M*\sum_{s=1}^q Ks.^2/2, the
%implementation here is based on this equation
%the energy of the system is equal to the energy term from the first
%network, plus the energy term from the second network, plus the energy
%term from the coupling, orthologs having the same spin have lower energy.
N1=length(A1);
%K1=full(sum(A1))';
E1=sum(K1S.^2)/(2*M1);
[u,v]=find(A1>0);
E1=E1-sum(sigma1(u)==sigma1(v));

N2=length(A2);
%K2=full(sum(A2))';
E2=sum(K2S.^2)/(2*M2);
[u,v]=find(A2>0);
E2=E2-sum(sigma2(u)==sigma2(v));

[I, J, L]=find(big_couple(1:N1,N1+1:N2+N1));
O12=-sum(L(sigma1(I)==sigma2(J)));

E=E1+E2+O12*kappa;
end

function [A2P_array,A2N_array,A3P_array,A3N_array]=setup(A2,A3)
A2P_array=cell(length(A2),1);
A2N_array=cell(length(A2),1);
A3P_array=cell(length(A3),1);
A3N_array=cell(length(A3),1);

for i=1:length(A2);
    A2P_array{i}=find(A2(i,:)>0);
    A2N_array{i}=find(A2(i,:)<0);
end
for i=1:length(A3);
    A3P_array{i}=find(A3(i,:)>0);
    A3N_array{i}=find(A3(i,:)<0);
end
end

function clusters = louvain_unsigned(A1,A2,big_couple,kappa)
N1=length(A1);
N2=length(A2);

N=N1+N2;
q=N;
sigma1=(1:N1)';
sigma2=N1+(1:N2)';
gamma=1;
Niter=0;

A1_array=cell(N1,1);
A2_array=cell(N2,1);
for i=1:N1;
    A1_array{i}=find(A1(:,i));
end
for i=1:N2;
    A2_array{i}=find(A2(:,i));
end

K1=full(sum(A1))';
M1=sum(K1)/2;%number of edges
Ks1=zeros(q,1);
K2=full(sum(A2))';
M2=sum(K2)/2;
Ks2=zeros(q,1);

for k=1:q
        Ks1(k)=sum(K1(sigma1==k));%%for each spin, find all nodes with the spin label, and sum their degrees
        Ks2(k)=sum(K2(sigma2==k));
end

E_init=get_energy_unsigned(A1,A2,M1,M2,Ks1,Ks2,sigma1,sigma2,big_couple,kappa);
gain=1;
E=E_init;

seq=randperm(N);
%seq=1:N;
while (gain==1)
    Cost = zeros(1,N);
    gain = 0;
    
    
    for j=1:N        
        x = seq(j);
        if x<=N1
            spin=sigma1(x);
            neighbors_spin=sigma1(A1_array{x});
            antiferro_delta=(Ks1-Ks1(spin))*K1(x)+K1(x).^2;%new-old;
            antiferro_delta(spin)=0;
            M=M1;            
        elseif x>N1          
            spin=sigma2(x-N1);
            neighbors_spin=sigma2(A2_array{x-N1});
            antiferro_delta=(Ks2-Ks2(spin))*K2(x-N1)+K2(x-N1).^2;%new-old;
            antiferro_delta(spin)=0;
            M=M2;
        end
        nq_hist=histc(neighbors_spin,1:q,1);%among the neighbors, how many of them have 1,2 ..q?
        neighbors_delta=(-nq_hist(spin)+nq_hist);%the change with respect to all possibilites...%new-old% to all spin value..indexed by spin value        
        E_delta=-neighbors_delta+gamma*antiferro_delta/(2*M);%at the beginning, flipping the spin to values at other network cannot change E
        
        i_ortho=find(big_couple(x,:));
        couple=full(big_couple(x,i_ortho));
        if x<=N1
            i_ortho=i_ortho-N1;
            sigma_tmp=sigma2;
        else
            sigma_tmp=sigma1;
        end
        ortho_delta=zeros(q,length(i_ortho));
        if ~isempty(i_ortho)
            for o=1:length(i_ortho)
                if spin~=sigma_tmp(i_ortho(o))
                    ortho_delta(sigma_tmp(i_ortho(o)),o)=couple(o);%deltaE chart is [0 0...1 0 0]', just means New-OLD, with OLD=0, and New=1 for a particular q
                else
                    ortho_delta(:,o)=-couple(o);
                    ortho_delta(sigma_tmp(i_ortho(o)),o)=0;%deltaE chart is [-1 -1...0 -1 -1]'; OLD=1, New=1 for a particular q, others=0
                end
            end
        end
        
        E_delta=E_delta-kappa*sum(ortho_delta,2);
                
        [Cost(j),id]=min(E_delta);%
        new_spin=id;%pick the most negative change
        
        if x<=N1
            sigma1(x)=new_spin;
            Ks1(spin)=Ks1(spin)-K1(x);
            Ks1(new_spin)=Ks1(new_spin)+K1(x);
        else
            sigma2(x-N1)=new_spin;
            Ks2(spin)=Ks2(spin)-K2(x-N1);
            Ks2(new_spin)=Ks2(new_spin)+K2(x-N1);            
        end
        
    end
    
    [Nco,sigma1,sigma2]=my_union(sigma1,sigma2);
    
    q=length(Nco);
    Ks1=Ks1(Nco);
    Ks2=Ks2(Nco);
    %sigma=[sigma1' sigma2']';
    %Sco=histc(sigma,1:length(Nco));
    %Nco2 = length(Sco(Sco>1));
    
    E=E+sum(Cost);
    if sum(Cost)<0
        gain=1;
    end
   
    fprintf('Iteration %d - Energy=%f %d clusters \n',Niter,E,length(Nco));
    Niter = Niter + 1;

end

Niter = Niter - 1;
clusters.species{1} = sigma1;
clusters.species{2}= sigma2;
clusters.energy = E;
clusters.Niter = Niter;

end

function clusters = louvain_signed(A2,A3,big_couple,kappa)
%initial setup
N2=length(A2);
N3=length(A3);
N=N2+N3;
q=N;
sigma2=(1:N2)';
sigma3=N2+(1:N3)';
gamma=1;
lambda=1/2;
Niter=0;
[A2P_array,A2N_array,A3P_array,A3N_array]=setup(A2,A3);

K2P=full(sum(A2>0))';
K2N=full(sum(A2<0));
M2P=sum(K2P)/2;
M2N=sum(K2N)/2;
K2SP=zeros(q,1);
for k=1:q
        K2SP(k)=sum(K2P(sigma2==k));
end
K2SN=zeros(q,1);
for k=1:q
        K2SN(k)=sum(K2N(sigma2==k));
end

K3P=full(sum(A3>0))';
K3N=full(sum(A3<0));
M3P=sum(K3P)/2;
M3N=sum(K3N)/2;
K3SP=zeros(q,1);
for k=1:q
        K3SP(k)=sum(K3P(sigma3==k));
end
K3SN=zeros(q,1);
for k=1:q
        K3SN(k)=sum(K3N(sigma3==k));
end

E_init=get_energy_signed(A2,A3,M2P,M2N,M3P,M3N,K3SP,K3SN,K2SN,K2SP,sigma2,sigma3,big_couple,kappa);
seq=randperm(N);
%seq=1:N;
gain=1;
E=E_init;

fprintf('start iteration\n');
end_aux=0;

while (gain==1)
    Cost = zeros(1,N);
    gain = 0;
    
    %go through all nodes once, for each node, flip its spin to lower
    %the overall energy
    for j=1:N
        x = seq(j);
        if x<=N2
            spin=sigma2(x);
            neighbors_spin_p=sigma2(A2P_array{x});
            antiferro_delta_p=(K2SP-K2SP(spin))*K2P(x)+K2P(x).^2;%new-old;
            antiferro_delta_p(spin)=0;
            neighbors_spin_n=sigma2(A2N_array{x}); %Creates list of spins of neighbor nodes
            antiferro_delta_n=(K2SN-K2SN(spin))*K2N(x)+K2N(x).^2;%new-old;
            antiferro_delta_n(spin)=0;
            Mp=M2P;Mn=M2N;
        elseif x>N2
            spin=sigma3(x-N2);
            neighbors_spin_p=sigma3(A3P_array{x-N2});
            antiferro_delta_p=(K3SP-K3SP(spin))*K3P(x-N2)+K3P(x-N2).^2;%new-old;
            antiferro_delta_p(spin)=0;
            neighbors_spin_n=sigma3(A3N_array{x-N2}); %Creates list of spins of neighbor nodes
            antiferro_delta_n=(K3SN-K3SN(spin))*K3N(x-N2)+K3N(x-N2).^2;%new-old;
            antiferro_delta_n(spin)=0;
            Mp=M3P;Mn=M3N;
        end
        
        nq_hist_p=histc(neighbors_spin_p,1:q,1);
        neighbors_delta_p=-nq_hist_p(spin)+nq_hist_p;%the change with respect to all possibilites...%new-old
        nq_hist_n=histc(neighbors_spin_n,1:q,1);
        neighbors_delta_n=-nq_hist_n(spin)+nq_hist_n;
        E_delta_n=neighbors_delta_n-gamma*antiferro_delta_n/(2*Mn);
        E_delta_p=-neighbors_delta_p+gamma*antiferro_delta_p/(2*Mp);
        
        i_ortho=find(big_couple(x,:));
        couple=full(big_couple(x,i_ortho));
        if x<=N2
            i_ortho=i_ortho-N2;
            sigma_tmp=sigma3;
        else
            sigma_tmp=sigma2;
        end
        ortho_delta=zeros(q,length(i_ortho));
        if ~isempty(i_ortho)
            for o=1:length(i_ortho)
                if spin~=sigma_tmp(i_ortho(o))
                    ortho_delta(sigma_tmp(i_ortho(o)),o)=couple(o);
                else
                    ortho_delta(:,o)=-couple(o);
                    ortho_delta(sigma_tmp(i_ortho(o)),o)=0;
                end
            end
        end
        
        E_delta=E_delta_n*(1-lambda)+lambda*E_delta_p-kappa*sum(ortho_delta,2);
        
        [Cost(j),id]=min(E_delta);%
        new_spin=id;%pick the most negative change
        
        if x<=N2
            sigma2(x)=new_spin;
            K2SP(spin)=K2SP(spin)-K2P(x);
            K2SP(new_spin)=K2SP(new_spin)+K2P(x);
            K2SN(spin)=K2SN(spin)-K2N(x);
            K2SN(new_spin)=K2SN(new_spin)+K2N(x);
        else
            sigma3(x-N2)=new_spin;
            K2SN(spin)=K2SN(spin)-K2N(x-N2);
            K2SN(new_spin)=K2SN(new_spin)+K2N(x-N2);
        end
        
    end
    
    [Nco,sigma2,sigma3]=my_union(sigma2,sigma3);
    q_old=q;
    q=length(Nco);
    if q==q_old;
        end_aux=end_aux+1;
    else
        end_aux=0;
    end
    K3SP=K3SP(Nco);
    K2SP=K2SP(Nco);
    K3SN=K3SN(Nco);
    K2SN=K2SN(Nco);
    %sigma=[sigma2' sigma3']';
    %Sco=histc(sigma,1:length(Nco));
    %Nco2 = length(Sco(Sco>1));
    E=E+sum(Cost);
    if sum(Cost)<0 && end_aux<4
        gain=1;%stop if energy cannot be further lowered
    end
    %fprintf('Iteration %d - Energy=%f %d clusters (%d non isolated)\n',Niter,E,length(Nco),Nco2);
    fprintf('Iteration %d - Energy=%f %d clusters \n',Niter,E,length(Nco));
    Niter = Niter + 1;
end

Niter = Niter - 1;
clusters.species{1} = sigma2;
clusters.species{2}= sigma3;
clusters.energy = E;
clusters.Niter = Niter;

end

function clusters=anneal_Potts_net_2species(A1,A2,orthologs12,couple12,kappa)

fprintf('initialize the annealing process\n');
N1=length(A1);
N2=length(A2);
if sum(sum(A1<0))==0 && sum(sum(A2<0))==0;
    fprintf('the simulated annealing function currently is for signed network as illustrated in the paper. considering using the steepest descent function for unsigned networks.\n');
    clusters=[];
    return;
end

[A1P_array,A1N_array,A2P_array,A2N_array]=setup(A1,A2);
gamma=1;
lambda=1/2;%relative contribution of positive and negative links
q=250;%max. number of modules allowed.
T_init=3;%initial temp, if it's chosen very high, it takes longer to cool down the system.
dd=0.9;%cooling rate
stps=250;%number of sweeping
stps_frac=0.8;%only this frac. of time series is used
equ_threshold=0.05;%a threshold to test whether the system is in equilibrium, at equ. E-mean(E) follows a Gaussian with mean zero, 
stop_threshold=0.01;%if average flipping rate is lower than this threshold, the process stops.
nT=100;%number of cooling steps
T=T_init*dd.^(0:nT);%cooling process

rng('shuffle')
sigma1=randi(q,N1,1);
sigma2=randi(q,N2,1);

mkdir('trajectory');
mkdir('ground_states');
next_temp=0;
DeltaE_sweep_cat=[];
flip_cat=[];
fprintf('temperature = %f\n',T(1));
while next_temp==0;
    beta=1./T(1);
    [sigma1,sigma2,DeltaE_sweep,flip] = Potts_net_2species(A1P_array,A1N_array,A2P_array,A2N_array,...
        sigma1,sigma2,orthologs12,couple12,q,beta,gamma,lambda,kappa,stps);
    flip_cat=[flip_cat flip'];
    DeltaE_sweep_cat=[DeltaE_sweep_cat DeltaE_sweep'];
    E_sweep_aux=fliplr(cumsum(DeltaE_sweep));%NB this is energy of the system up to a constant...
    E_sweep_aux=E_sweep_aux(1:stps*stps_frac);
    H=swtest(E_sweep_aux,equ_threshold);%%Null Hypothesis:X is normal with unspecified mean and variance.
    if H==0;%H=0 means Do not reject the null hypothesis
        next_temp=1;
    end
end
savefile=['./trajectory/config_' num2str(T(1)) '.mat'];
save(savefile,'DeltaE_sweep_cat','sigma1','sigma2','flip_cat');

for i=2:length(T)
    fprintf('temperature = %f\n',T(i));
    if mean(flip)>stop_threshold;
        beta=1./T(i);
        next_temp=0;
        DeltaE_sweep_cat=[];
        flip_cat=[];
        while next_temp==0;
		[sigma1,sigma2,DeltaE_sweep,flip] = Potts_net_2species(A1P_array,A1N_array,A2P_array,A2N_array,...
        		sigma1,sigma2,orthologs12,couple12,q,beta,gamma,lambda,kappa,stps);
            flip_cat=[flip_cat flip'];
            DeltaE_sweep_cat=[DeltaE_sweep_cat DeltaE_sweep'];
            E_sweep_aux=fliplr(cumsum(DeltaE_sweep));%NB this is energy of the system up to a constant...
            E_sweep_aux=E_sweep_aux(1:stps*stps_frac);
            H=swtest(E_sweep_aux,equ_threshold);
            if H==0;
                next_temp=1;
                i_ed=i;
            end
        end
        savefile=['./trajectory/config_' num2str(T(i)) '.mat'];
        save(savefile,'DeltaE_sweep_cat','sigma1','sigma2','flip_cat');
    end
end
savefile=['./ground_states/config_' num2str(T(i_ed)) '.mat'];
save(savefile,'DeltaE_sweep_cat','sigma1','sigma2','flip_cat');

clusters.species{1} = sigma1;
clusters.species{2}= sigma2;

end

function [sigma1,sigma2,DeltaE_sweep,flip] = Potts_net_2species(A1P_array,A1N_array,A2P_array,A2N_array,...
    sigma1,sigma2,orthologs12,couple12,q,beta,gamma,lambda,kappa,stps)

N1=length(A1P_array);
N2=length(A2P_array);
N=N1+N2;
flip=zeros(stps,1);
DeltaE_sweep=zeros(stps,1);

K1p=zeros(N1,1);K1n=zeros(N1,1);
for i=1:N1
        K1p(i)=length(A1P_array{i});
        K1n(i)=length(A1N_array{i});
end
M1p=sum(K1p)/2;
M1n=sum(K1n)/2;
K1sp=zeros(q,1);
for k=1:q
        K1sp(k)=sum(K1p(sigma1==k));
end
K1sn=zeros(q,1);
for k=1:q
        K1sn(k)=sum(K1n(sigma1==k));
end

K2p=zeros(N2,1);K2n=zeros(N2,1);
for i=1:N2
        K2p(i)=length(A2P_array{i});
        K2n(i)=length(A2N_array{i});
end
M2p=sum(K2p)/2;
M2n=sum(K2n)/2;
K2sp=zeros(q,1);
for k=1:q
        K2sp(k)=sum(K2p(sigma2==k));
end
K2sn=zeros(q,1);
for k=1:q
        K2sn(k)=sum(K2n(sigma2==k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ortholog_spinchart12 = zeros(length(orthologs12),4);
ortholog_spinchart12(:,1) = orthologs12(:,1); %First row is the labels for the nodes that have orthologs
ortholog_spinchart12(:,3) = orthologs12(:,2); %Third row is the labels for the orthologs in the first column
ortholog_spinchart12(:,2) = sigma1(orthologs12(:,1)); %Second column is the spins of the orthologs of the first network
ortholog_spinchart12(:,4) = sigma2(orthologs12(:,2)); %Fourth column is the spins of the orthologs of the second network
%E_ortho12=-kappa*sum((ortholog_spinchart12(:,2)==ortholog_spinchart12(:,4)).*couple12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DeltaE_series_tmp=zeros(N,1);

for i=1:stps
        seq=randperm(N);
        for j=1:N;
                x = seq(j);
                
                if x<=N1
                    N2adjust_a = 0;
                    N2adjust_b=0;
                    N2adjust_c=-1;
                    spin=sigma1(x); %Original spin
                    neighbors_spin_p=sigma1(A1P_array{x}); %Creates list of spins of neighbor nodes
                    antiferro_delta_p=(K1sp-K1sp(spin))*K1p(x)+K1p(x).^2;%new-old;
                    antiferro_delta_p(spin)=0;
                    neighbors_spin_n=sigma1(A1N_array{x}); %Creates list of spins of neighbor nodes
                    antiferro_delta_n=(K1sn-K1sn(spin))*K1n(x)+K1n(x).^2;%new-old;
                    antiferro_delta_n(spin)=0;
                    Mp=M1p;Mn=M1n;
                else
                    x=x-N1;
                    N2adjust_a =2;
                    N2adjust_b=-1;
                    N2adjust_c=0;
                    spin=sigma2(x);
                    neighbors_spin_p=sigma2(A2P_array{x});%%this line is the most time consuming...
                    antiferro_delta_p=(K2sp-K2sp(spin))*K2p(x)+K2p(x).^2;%new-old;
                    antiferro_delta_p(spin)=0;
                    neighbors_spin_n=sigma2(A2N_array{x});
                    antiferro_delta_n=(K2sn-K2sn(spin))*K2n(x)+K2n(x).^2;%new-old;
                    antiferro_delta_n(spin)=0;
                    Mp=M2p;Mn=M2n;
                end
                
                nq_hist_p=histc(neighbors_spin_p,1:q,1);
                neighbors_delta_p=-nq_hist_p(spin)+nq_hist_p;%the change with respect to all possibilites...%new-old
                nq_hist_n=histc(neighbors_spin_n,1:q,1);
                neighbors_delta_n=-nq_hist_n(spin)+nq_hist_n;
                E_delta_n=neighbors_delta_n-gamma*antiferro_delta_n/(2*Mn);
                E_delta_p=-neighbors_delta_p+gamma*antiferro_delta_p/(2*Mp);
                
                O12=[];
                if N2adjust_a>=0
                    O12 = find(ortholog_spinchart12(:,1+N2adjust_a)==x); %Check if ortholog exists
                    O12_tmp=O12(:,ones(q,1))';
                    couple12_tmp=couple12(O12_tmp);
                end
     
                if ~isempty(O12)
                        ortholog_delta_12 = zeros(q,length(O12));
                        for r=1:length(O12)%ortholog exists for node j
                                if spin~=ortholog_spinchart12(O12(r),4-N2adjust_a) %original spin is not equal to ortholog's spin
                                        ortholog_delta_12(ortholog_spinchart12(O12(r),4-N2adjust_a),r) = 1; %deltaE chart is [0 0...1 0 0]', just means New-OLD, with OLD=0, and New=1 for a particular q
                                elseif spin==ortholog_spinchart12(O12(r),4-N2adjust_a)
                                        ortholog_delta_12(:,r) =-1;
                                        ortholog_delta_12(ortholog_spinchart12(O12(r),4-N2adjust_a),r) = 0;  %deltaE chart is [-1 -1...0 -1 -1]'; OLD=1, New=1 for a particular q, others=0
                                end
                        end
                        ortholog_delta_12=ortholog_delta_12.*couple12_tmp;
                else
                        ortholog_delta_12= zeros(q,1);
                end
                
                E_delta=E_delta_n*(1-lambda)+lambda*E_delta_p-kappa*sum(ortholog_delta_12,2);

                %to avoid exp gives Inf, we do the following scaling
                E_delta_aux=E_delta;
                if sum(beta*E_delta_aux<-709)>0
                        D=-min(E_delta_aux(beta*E_delta_aux<-709))-709/beta;%709 is the largest number that exp makes sense...
                        E_delta_aux=E_delta_aux+D;
                end
                
                Z=sum(exp(-beta*E_delta_aux));%NB the pos
                P=exp(-beta*E_delta_aux)./Z;
                cP=cumsum(P);cP(end)=1;
                R=rand;
                id=find(R<=cP);
                new_spin=id(1);
                if new_spin~=spin;
                        flip(i)=flip(i)+1;
                end

                if N2adjust_c<0%x in net 1
                    sigma1(x)=new_spin;
                    K1sp(spin)=K1sp(spin)-K1p(x);
                    K1sp(new_spin)=K1sp(new_spin)+K1p(x);
                    K1sn(spin)=K1sn(spin)-K1n(x);
                    K1sn(new_spin)=K1sn(new_spin)+K1n(x);
                    
                elseif N2adjust_b<0%x in net2
                    sigma2(x)=new_spin;
                    K2sp(spin)=K2sp(spin)-K2p(x);
                    K2sp(new_spin)=K2sp(new_spin)+K2p(x);
                    K2sn(spin)=K2sn(spin)-K2n(x);
                    K2sn(new_spin)=K2sn(new_spin)+K2n(x);
                end
                
                if ~isempty(O12)
                    ortholog_spinchart12(O12,2+N2adjust_a) = new_spin; %Update ortholog_spinchart
                end
                DeltaE_series_tmp(j)=E_delta(new_spin);
        
        end
        
        DeltaE_sweep(i)=sum(DeltaE_series_tmp);
        DeltaE_series_tmp=zeros(N,1);
        
end

flip=flip/N;

end
