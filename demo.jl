include("./lib_orthoclust.jl");

vars = matread("./data/worm_fly_networks_info.mat");
Af=vars["Af"];
#Af=full(Af);
Aw=vars["Aw"];
#Aw=full(Aw);
ortho=vars["orthologs"];
fly_id=vars["fly_id"];
worm_id=vars["worm_id"];

networks=Array{Any}(2);
networks[1]=Aw;
networks[2]=Af;

output=multiplex_louvain(networks,ortho);
