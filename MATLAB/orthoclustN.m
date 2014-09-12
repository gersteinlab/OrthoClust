%[assign,Q]=orthoclustN(A,O_pairs,kappa)
%A is a N-by-1 cell array storing all the networks, eg. A{1} is human network, A{2}
%is worm network, A{3} is fly network
%O_pairs is a N by N cell array storing all pairwise orthologous
%relationships...i.e. O_pairs(1,2) is a e-by-2 array matching the e
%orthologous pairs between species 1 and 2. only the upper N(N-1)/2
%elements are used.
%kappa is the coupling constant defined in our paper.
%the output assign is a N-by-1 cell array that maps genes in each species
%to module IDs.
function [assign,Q]=orthoclustN(A,O_pairs,kappa)
T=length(A);
gamma=1;
N=zeros(T,1);
for i=1:T
    N(i)=length(A{i});
end
Ntot=sum(N);
o1=[];o2=[];ow=[];
for s1=1:T-1
    for s2=s1+1:T
        w=get_weight(O_pairs{s1,s2},N(s1),N(s2));
        O_pairs_tmp=O_pairs{s1,s2};
        if s1>1
            o1_tmp=O_pairs_tmp(:,1)+sum(N(1:s1-1));
        else
            o1_tmp=O_pairs_tmp(:,1);
        end
        if s2>1
            o2_tmp=O_pairs_tmp(:,2)+sum(N(1:s2-1));
        else
            o2_tmp=O_pairs_tmp(:,2);
        end
        o1=[o1;o1_tmp];
        o2=[o2;o2_tmp];
        ow=[ow;w];
    end
end
coupling=sparse(o1,o2,ow,Ntot,Ntot);
coupling=coupling+coupling';

ii=[]; jj=[]; vv=[];st_tmp=0;
twomu=0;
node2species=zeros(Ntot,1);
for s=1:T
    indx=(1:N(s))'+st_tmp;
    [i,j,v]=find(A{s});
    ii=[ii;indx(i)]; jj=[jj;indx(j)]; vv=[vv;v];
    k=sum(A{s});
    kv=zeros(Ntot,1);
    twom=sum(k);
    twomu=twomu+twom;
    kv(indx)=k'/twom;
    kcell{s}=kv;
    st_tmp=indx(end);
    node2species(st_tmp+1:end)=node2species(st_tmp+1:end)+1;
end
node2species=node2species+1;

AA = sparse(ii,jj,vv,Ntot,Ntot);
clear ii jj vv st_tmp;
kvec = full(sum(AA));

AA = AA + kappa*coupling;
twomu=twomu+T*kappa*Ntot*(T-1);
B = @(i) AA(:,i) - gamma*kcell{node2species(i)}*kvec(i);

[S,Q] = genlouvain(B);%this code is written by Mucha and Porter that can be download via
%http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
assign=cell(T,1);
for s=1:T
    assign{s}=S(node2species==s);
end
Q=Q/twomu;
end


function Ow=get_weight(list_orthologs,N1,N2)
O12=sparse(list_orthologs(:,1),list_orthologs(:,2),ones(size(list_orthologs,1),1),N1,N2);
d1=full(sum(O12,2));
d2=full(sum(O12))';
o1=d1(list_orthologs(:,1));
o2=d2(list_orthologs(:,2));
Ow=(1./o1+1./o2)/2;
end