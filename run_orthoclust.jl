include("./lib_orthoclust.jl");

files=[];
for x in ARGS;
	files=vcat(files,x)
end

networks=Array{Any}(length(files)-2);

for i=1:length(networks);
	X=readdlm(files[i],Int64);
	N=maximum(X);
	net=sparse(X[:,1],X[:,2],ones(size(X,1)),N,N);
	networks[i]=net;
end

ortho=readdlm(files[end-1],Int64);
gamma=float(files[end]);

multiplex_louvain(networks,ortho,gamma);