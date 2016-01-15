using MAT;
using Graphs;
using DataFrames;


function optimize_network_modularity_louvain(W,order)
#order=1, forward
#order=-1, backward
#order=0, random
	
	N=size(W,1);
    i_no_dark=find(sum(abs(W),1).>0);
    N_no_dark=length(i_no_dark);

    res=1;
    k=sum(W,2);sW=sum(W);
    E_W=k*k'/sW;

    B=full(W)-E_W*res;
    E_W=0;
    Bcompact=B[i_no_dark,i_no_dark];
    Wcompact=W[i_no_dark,i_no_dark];
    B=0;W=0;

    num_run=1;
    (assign,Q,Brenorm,Wrenorm)=iterate_network_modularity(Bcompact,Wcompact,order);
    transfer=sparse(collect(1:N_no_dark),assign,ones(size(assign)));

    keep_doing=1;

    while keep_doing==1
        (tmp_assign,tmp_Q,tmp_Brenorm,tmp_Wrenorm)=iterate_network_modularity(Brenorm,Wrenorm,order);
        tmp_transfer=sparse(collect(1:size(tmp_assign,1)),tmp_assign,ones(size(tmp_assign)));
        if isequal(tmp_transfer,speye(length(tmp_assign)))
        	keep_doing=0;
        	#Brenorm, Wrenorm are optimal.
        else
	        transfer=transfer*sparse(collect(1:size(tmp_assign,1)),tmp_assign,ones(size(tmp_assign)));
        	Brenorm=tmp_Brenorm+0;
        	Wrenorm=tmp_Wrenorm+0;
        	Q=tmp_Q+0;
        end
    end

    (u,v)=findn(transfer);
    iu=sortperm(u);
    tmp_assign=v[iu];
    final_assign=zeros(N);
    final_assign[i_no_dark]=tmp_assign;

    return final_assign, Q, Brenorm, Wrenorm;

end

#(sigma2, Q2, Brenorm, Wrenorm)=iternate_network_modularity(Bcompact,Wcompact,order);
#the method of re-normalizating the nodes makes sense only for conventional network modularity defined by w_i*w_j
function iterate_network_modularity(Bcompact,Wcompact,order);

    Nb=size(Bcompact,1);
    sigma=collect(1:Nb);
    sW=sum(Wcompact);

	if order==1
        u=collect(1:Nb);
    elseif order==-1
        u=flipdim(collect(1:Nb),1);
    elseif order==0
        u=collect(1:Nb);
        u=u[randperm(Nb)];
    end

    gain=1;Niter=1;

    while (gain==1)
        gain = 0;
        for j=1:Nb
            x=u[j];
            # spin=sigma[i];
            display(j);
            spin=sigma[x];
            c=Bcompact[x,:];
            c[x]=0;#this is important step to make sure the deltaQ is right 
            neighbors_spin=sigma;
            DeltaQ=-sum(c'.*(sigma.==spin))+full(sparse(neighbors_spin,[1 for dd=1:Nb],vec(c)));
            #the 2nd term sum over the components from each community in advance
            #1st term, the effect of getting rid of the original spin contribution..
            id=indmax(DeltaQ);#note the dim of DeltaQ is the number of communities
            new_spin=id;
            if (new_spin!=spin)&(DeltaQ[id].>0);
                gain=1;
                sigma[x]=new_spin;
            end
        end
        Q=compute_modularity(sigma,Bcompact,sW);
        @printf("iteration %d - sum = %f %d modules \n",Niter,Q,length(unique(sigma)))
        Niter = Niter + 1
    end

    Q=compute_modularity(sigma,Bcompact,sW);

    sigma2=relabel_communities(sigma);
    usigma=sort(unique(sigma));
    N_renorm=length(usigma);

    #this is the normalization required for this null model, to make sense the Q remains the same
    Wrenorm=zeros(N_renorm,N_renorm);

    for i=1:N_renorm
        for j=1:N_renorm
            Wrenorm[i,j]=sum(sum(Wcompact[sigma.==usigma[i],sigma.==usigma[j]]));
        end
    end
    w=sum(Wrenorm,2);

    E_Wrenorm=(w*w')/sW;

    Brenorm=Wrenorm-E_Wrenorm;

    Q2=compute_modularity(collect(1:N_renorm),Brenorm,sW);
    @printf("step - Q = %f %d modules \n",Q2,length(unique(sigma)))
    #println(Q2), verified to be the same as Q...
    return sigma2, Q2, Brenorm, Wrenorm;
end

#sigma_new=relabel_communities(sigma);
function relabel_communities(sigma)
    u=unique(sigma);
    u=sort(u);
    sigma_new=zeros(size(sigma));
    for i=1:length(u)
        iz=findin(sigma,u[i]);
        sigma_new[iz]=i;
    end
    sigma_new=round(Int64,sigma_new);
    return sigma_new;
end


#Q=compute_modularity(sigma,Brenorm,sW);
function compute_modularity(sigma,Brenorm,sW);
    #good for both corr, based or else
    COMu = unique(sigma);
    Q=0;
    for k=1:length(COMu)
        id = find(sigma.==COMu[k]);
        Q=Q+sum(Brenorm[id,id]);
    end
    Q=Q/sW;
    return Q;
end


function matrix_to_graph(Z);
    Z=Z-diagm(diag(Z));
    Z=triu(Z);
    nnodes=size(Z,1);
    nedges=sum(Z)/2
    g=Graphs.simple_graph(nnodes,is_directed=false);
    for i=1:nnodes;
        tg=find(Z[i,:]);
        for j=1:length(tg)
            Graphs.add_edge!(g,i,tg[j]);
        end
    end
    return g;
end

function generate_mapping_id(networks);

    Nw=size(networks,1);
    mapping_s2b=Array{Any}(Nw);
    N_networks=zeros(Int64,Nw);
    for i=1:Nw;
        N_networks[i]=size(networks[i],1);
    end
    mapping_b2s=zeros(Int64,(sum(N_networks),2));
    for i=1:Nw;
        if i.==1
            mapping_s2b[i]=collect(1:N_networks[i]);
            mapping_b2s[collect(1:N_networks[i]),1]=1;
            mapping_b2s[collect(1:N_networks[i]),2]=1:N_networks[i];
        else
            mapping_s2b[i]=collect(1:N_networks[i])+N_networks[i-1];
            mapping_b2s[collect(1:N_networks[i])+N_networks[i-1],1]=i;
            mapping_b2s[collect(1:N_networks[i])+N_networks[i-1],2]=1:N_networks[i];
        end
    end
    return N_networks,mapping_b2s, mapping_s2b;    

end

function expand_mat(networks,ortho)

    num_species=size(networks,1);
    aux=ones(Int,(num_species,num_species));
    aux=triu(aux-diagm(diag(aux)));
    (aux1,aux2,aux3)=findnz(aux);
    pairs=[aux1 aux2];

    N_networks,mapping_b2s, mapping_s2b=generate_mapping_id(networks);
    tot_nodes=sum(N_networks);
    for i=1:num_species
        net=networks[i];
        if net==net'
            net=triu(net);
        end
        (u,v,w)=findnz(net);
        id=mapping_s2b[i];
        u=id[u];v=id[v];        
        if i==1
            all_e1=u;
            all_e2=v;
        else
            all_e1=vcat(all_e1,u);
            all_e2=vcat(all_e2,v);
        end
    end

    e12_values=ones(Float64,size(all_e1));
    ortho_transform=zeros(Int,size(ortho));
    ortho_transform[:,1]=ortho[:,1];
    ortho_transform[:,2]=ortho[:,2];

    for j=1:num_species
        iz=find(ortho[:,1].==j);
        mapping=mapping_s2b[j];
        if ~isempty(iz)
            ortho_transform[iz,4]=mapping[ortho[iz,4]];
        end
        iz2=find(ortho[:,2].==j);
        if ~isempty(iz2)
            ortho_transform[iz2,5]=mapping[ortho[iz2,5]];
        end
    end

    couple_const=zeros(Int,size(pairs,1));
    ortho_couple=[];
    for i=1:size(pairs,1);
        i_pick=find((ortho[:,1].==pairs[1,1]).*ortho[:,2].==pairs[1,2]);
        couple_const[i]=unique(ortho[i_pick,3])[1];
        couple_mat=sparse(ortho[i_pick,4],ortho[i_pick,5],ones(size(i_pick)));
        d1=full(sum(couple_mat,2));
        d2=full(sum(couple_mat,1))';
        o1=d1[ortho[i_pick,4]];
        o2=d2[ortho[i_pick,5]];
        Ow=(1./o1+1./o2)/2;
        if i.==1;
            ortho_couple=Ow;
        else
            ortho_couple=vcat(ortho_couple,Ow*couple_const[i]);
        end
    end

    all_e1=vcat(all_e1,ortho_transform[:,4]);
    all_e2=vcat(all_e2,ortho_transform[:,5]);
    e12_values=vcat(e12_values,ortho_couple);

    big_mat=sparse(all_e1,all_e2,e12_values,tot_nodes,tot_nodes);
    tmp=spdiagm(diag(big_mat));
    big_mat_aux=big_mat-tmp;

    big_mat=big_mat_aux+big_mat_aux'+tmp;
        
    return big_mat;

end

function multiplex_louvain(networks,ortho);

    big_mat=expand_mat(networks,ortho);
    final_assign, Q, Brenorm, Wrenorm=optimize_network_modularity_louvain(big_mat,0);
    N_networks,mapping_b2s, mapping_s2b=generate_mapping_id(networks);
    output=[mapping_b2s final_assign];
    writedlm("./genes_to_clusters.out",output);

    return output;

end


#function get_coexp_net(C,d)
#C is a pairwise correlation matrix, d is a rank based cutoff such that 
#every gene is connected to the top d most correlated genes
#    N=length(C);
#    C=C-diagm(diag(C));
    #[~,I]=sort(abs(C));
    #[~,I]=sort(I);
    #A=sparse(I>N-d);
    #A=(A.*C);
    #A=(A+A')/2;
    #A=sign(A);

#return A;



