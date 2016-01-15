using HDF5; 
using JLD;
using MAT;
using Graphs;
using DataFrames;
#using Winston;
using CurveFit;
using Distributions;


function load_annotation(annot_file);
	vars=matread(annot_file)
	Gene=vars["Gene_name"];
	ENSG=vars["Gene_id"];
	strand=vars["strand"];
	gene_st=vars["st"];
	gene_ed=vars["ed"];
	chr=vars["chr"];
	TPM=vars["TPM"];
	TPM_pe=vars["TPM_pe"]
	celltype=vars["celltype"]
	return Gene, ENSG, strand, gene_st, gene_ed, chr, celltype, TPM, TPM_pe;
end

function define_chr_size();
	chr_length=zeros(25,1);
	chr_length[1]=249250621;
	chr_length[2]=243199373;
	chr_length[3]=198022430;
	chr_length[4]=191154276;
	chr_length[5]=180915260;
	chr_length[6]=171115067;
	chr_length[7]=159138663;
	chr_length[8]=146364022;
	chr_length[9]=141213431;
	chr_length[10]=135534747;
	chr_length[11]=135006516;
	chr_length[12]=133851895;
	chr_length[13]=115169878;
	chr_length[14]=107349540;
	chr_length[15]=102531392;
	chr_length[16]=90354753;
	chr_length[17]=81195210;
	chr_length[18]=78077248;
	chr_length[19]=59128983;
	chr_length[20]=63025520;
	chr_length[21]=48129895;
	chr_length[22]=51304566;
	chr_length[23]=155270560;#X
	chr_length[24]=59373566;#Y
	chr_length[25]=16571;
	return chr_length;
end

function generate_arbitrary_mapping_files(bin_size);
	chr2bins=zeros(2,25);
	chr_length=define_chr_size();
	chr_num_bins=int(floor(chr_length/bin_size))+1;
	chr2bins[2,:]=cumsum(chr_num_bins)'-1;
	chr2bins[1,1]=0;
	chr2bins[1,2:end]=chr2bins[2,1:end-1]+1;
	X=int(chr2bins+1);
	bin2loc=zeros(3,X[2,end]);
	for c=1:25
		bin2loc[1,X[1,c]:X[2,c]]=c-1;
		bin2loc[2,X[1,c]:X[2,c]]=int(collect(1:bin_size:chr_length[c]))';
		bin2loc[3,X[1,c]:X[2,c]]=[int(collect(bin_size:bin_size:chr_length[c]))' chr_length[c]];
	end
	return int(chr2bins),int(bin2loc);
	
end

function define_ENCODE3_HiC_files()
	all_HiC_files=Any["ENCODE3-A549C-HindIII-R1__hg19__hdf",
"ENCODE3-A549D-HindIII-R2__hg19__hdf",
"ENCODE3-CAK12B-HindIII-R2__hg19__hdf",
"ENCODE3-Caki2A-HindIII-R1__hg19__hdf",
"ENCODE3-G401A-HindIII-R1__hg19__hdf",
"ENCODE3-G401B-HindIII-R2__hg19__hdf",
"ENCODE3-LNCaPC-HindIII-R1__hg19__hdf",
"ENCODE3-LNCaP-HindII-R2__hg19__hdf",
"ENCODE3-NCIH460A-HindIII-R1__hg19__hdf",
"ENCODE3-NCIH460B-HindIII-R2__hg19__hdf",
"ENCODE3-PANC1B-HindIII-R1__hg19__hdf",
"ENCODE3-PANCIC-HindIII-R2__hg19__hdf",
"ENCODE3-RPMI7951C-HindIII-R1__hg19__hdf",
"ENCODE3-RPMI7951D-HindIII-R2__hg19__hdf",
"ENCODE3-SJCRH30A-HindIII-R1__hg19__hdf",
"ENCODE3-SJCRH30B-HindIII-R2__hg19__hdf",
"ENCODE3-SKMEL5A-HindIII-R1__hg19__hdf",
"ENCODE3-SKMEL5B-HindIII-R2__hg19__hdf",
"ENCODE3-SKNDZA-HindIII-R1__hg19__hdf",
"ENCODE3-SKNDZB-HindIII-R2__hg19__hdf",
"ENCODE3-SKNMCC-HindIII-R1__hg19__hdf",
"ENCODE3-SKNMCD-HindII-R2__hg19__hdf",
"ENCODE3-T470A-HindIII-R1__hg19__hdf",
"ENCODE3-T470B-HindIII-R2__hg19__hdf"];
	return all_HiC_files;
end


function read_ENCODE3_HiC_files(input_file_loc,output_file_loc,cell_type,bin_size,Gene,chr,gene_st,gene_ed);
	filename=input_file_loc*cell_type*"/C-"*bin_size*"/iced/"*chop(chop(chop(cell_type)))*"genome__C-"*bin_size*"-iced.hdf5";
	c=h5open(filename,"r") do file
		read(file);
	end
	A=c["interactions"];
	chr2bins=c["chr_bin_range"];
	bin2loc=c["bin_positions"];
	z=sum(isnan(A),1);
	L=size(A,1);
	dark_bins=find(z.==L);
	#A=0;

	G2B_dict=Dict{ASCIIString,Any}();
	for i=1:length(Gene)
		display(i)
		ch=chr[i];
		st=gene_st[i];
		ed=gene_ed[i];
		p1=find((bin2loc[1,:].==ch-1).*(bin2loc[2,:].<st))
		p2=find((bin2loc[1,:].==ch-1).*(bin2loc[3,:].>ed))
		iz=intersect(p1,p2);
		#iz=find((bin2loc[1,:].==ch-1).*(bin2loc[2,:].<st).*(bin2loc[3,:].>ed));
		if ~isempty(iz)
			G2B_dict[Gene[i]]=iz[1];
		elseif ~isempty(p1) && ~isempty(p2)
			G2B_dict[Gene[i]]=[p1[end],p2[1]];
		else
			G2B_dict[Gene[i]]=[];
		end
	end
	out_file=output_file_loc*cell_type*"_mapping_"*bin_size*".jld";
	save(out_file,"chr2bins",chr2bins,"bin2loc",bin2loc,"G2B_dict",G2B_dict,"dark_bins",dark_bins);	
	return chr2bins, bin2loc, G2B_dict,dark_bins;
end

#text files in form of 3 columns..
function read_HiCPro_outputs(input_file_loc,bin_size);
	file_loc=input_file_loc*"/"*bin_size;
	filename=file_loc*"/reads_"*bin_size*"_iced.matrix";
	table=readtable(filename,separator='\t',header=false);
	N=maximum(table[:x1]);
	A=sparse(table[:x1],table[:x2],table[:x3],N,N);
	tmp=A-spdiagm(diag(A));
	A=A+tmp';
	#A=full(A);
	out_name=file_loc*"/contact.jld";
	save(out_name,"interaction",A);

	return A;

end

#the input_file should be the raw_observed M(i,j), input_file_aux are the normalization vectors
#i and j are written as the left edge of the bin at a given resolution; for example, at 100 kb resolution, 
#the entry corresponding to the first row and tenth column of the matrix 
#would correspond to M_i,j, where i=0, j=900000).
#confirmed, in large bin matrices, we see i=0 at the beginning..
#for norm file, M_norm=M(i,j)/(n_vec(i)*n_vec(j));
#it seems that the norm files have always one line more than the size of W.
#I have checked the dark bins matches, it seems the last term of the file (0) is crap..
#for expected fiel, first line of the expected vector file is the expected number of contacts between two loci separated 
#by 0*RES base pairs, the second line of the file is the expected number of contacts between two loci
#separated by 1*RES base pairs, and so on. M(i,j) / the abs(i-j) term
#we look at the expected files stored under different chr. , they have the same wc -l. it seems that calculation
#are done in a genome-wide fashion. since the largest separation belongs to chromosome one, the number of lines should be the same as 
#the matrix dim of chromosome 1. BUT, the actual number of lines is 1 less than the dim of chr1. it doesn;t make sense because we count 
#from 0 to N-1 . may be we can add the last entry manually as the 2nd last
#cf the expected file in different chr., they are not identical, but differ by an overall constant close to 1...I don't understand why.
function read_Cell2014_matrices(input_file,chr_num,bin_size,input_file_aux=[]);
	
	chr_length=define_chr_size();
	X=readdlm(input_file,'\t');
	N=int(floor(chr_length[chr_num]/bin_size))+1;
	#Rao did consider the last bin, look at the raw files, the largest value (left bins)...
	M=sparse(int(X[:,1]/bin_size)+1,int(X[:,2]/bin_size)+1,X[:,3],N,N);
	tmp=sparse(diagm(diag(M)));
	M=M+M'-tmp;
	if isempty(input_file_aux)
		n_vec=[];
	else
		n_vec=readdlm(input_file_aux,'\t');
	end
	return M, n_vec;

end

#use chr_str, instead of integer..
function read_Cell2014_domains(rao_domains_file,chr_str,bin_size);
	df=readtable(rao_domains_file,separator='\t',header=true);
	iz=find((df[:chr1].==chr_str)&(df[:chr2].==chr_str));
	rao_domains_st=df[iz,:x1];
	rao_domains_ed=df[iz,:x2];

	id=sortperm(rao_domains_st);
	rao_domains_st=rao_domains_st[id];
	id=sortperm(rao_domains_ed);
	rao_domains_ed=rao_domains_ed[id];

	
	Rao_list=DataFrame(chr=ASCIIString[],domain_st=Int64[],domain_ed=Int64[],domain_st_bin=Int64[],domain_ed_bin=Int64[],idx=Int64[]);

	if chr_str=="X"
		chr_num=24;
	else
		chr_num=int(chr_str);
	end

	rao_domains_st_bin=int(floor(rao_domains_st/bin_size)+1)+chr2bins[:,chr_num][1];
	rao_domains_ed_bin=int(floor(rao_domains_ed/bin_size)+1)+chr2bins[:,chr_num][1];

	chr_string="chr"*chr_str;
    for i=1:size(rao_domains_st_bin,1)
        push!(Rao_list,[chr_string rao_domains_st[i] rao_domains_ed[i] rao_domains_st_bin[i] rao_domains_ed_bin[i] i]);
    end

    return Rao_list


end


function load_mapping_file(file);
	d=load(file);
	chr2bins=d["chr2bins"];
	bin2loc=d["bin2loc"];
	G2B_dict=d["G2B_dict"];
	dark_bins=d["dark_bins"];
	d=0;
	return chr2bins, bin2loc, G2B_dict, dark_bins;
end;

function get_G2B_info(Gene,G2B_dict,dark_bins)
	L=zeros(size(Gene));
	Ld=zeros(size(Gene));
	for i=1:length(Gene);
		m1=G2B_dict[Gene[i]];
		L[i]=m1[end]-m1[1]+1;
		x=[m1[1]:m1[end]];
		for k1=1:length(x)
			if x[k1] in dark_bins
				Ld[i]=Ld[i]+1;
			end
		end
	end
	return L, Ld;
end


function load_contact_matrix(file_loc,cell_type,bin_size);
	filename=file_loc*cell_type*"/C-"*bin_size*"/iced/"*chop(chop(chop(cell_type)))*"genome__C-"*bin_size*"-iced.hdf5";
	c=h5open(filename,"r") do file
		read(file);
	end
	A=c["interactions"]
	return A;
end

#NB. chr num starts with 1 to 25
function extract_chr(A,chr2bins,chr_num);
	st=1+chr2bins[1,chr_num];
	ed=1+chr2bins[2,chr_num];
	A_chr=A[st:ed,st:ed];
	return A_chr;
end

function generate_matrix_for_Plotter(contact_mat,A);
	f = open(contact_mat, "w");
	write(f, "This is the contact matrix.\n")
	writedlm(f,A);
	close(f);
end

#to display the matrix using HiCPlotter
#copy the file to /home/ky26/software/HiCPlotter/data/HiC/Human
#go to /home/ky26/software/HiCPlotter
#python2.7 HiCPlotter.py -f data/HiC/Human/my.matrix -n hES -chr chr10 -r 40000 -o myDixon -ptd 1 -pptd 1 -w 8 -tr 10 -pi 1
#python2.7 HiCPlotter.py -f data/HiC/Human/my.matrix -n hES -chr chr10 -r 40000 -o myDixon_chr10 -ptd 1 -pptd 1 -w 8 -tr 10

#-n, a list of labels
#-chr chr10, for the various bed files
#-r resolution in bp
#-o output file name
#-ptd 1 will display Rao's domains
#-pptd 1 will display Dixon's domains, both cell lines
#-fh, default is 1, the way we generated the files need no such option
#-s, -e zoom in to the bins in the matrix
#-w and -tr are used to calculate the insulation score, and domains by arrowhead..
#-pi, default is 0, for insulation score..
#-ptdb 1, domains in bar instead of triangle..

function generate_TADs_bed(TAD_list,out_file);
	X=[TAD_list[:chr] TAD_list[:domain_st] TAD_list[:domain_ed]];
	writedlm(out_file,X);
end


#we can run plotter at any working directory
#hic_file:
#if there are multiple files, kust seperate by a space, teh same for labels..
#pcustom, custom domain files
#muitple tracks separateed by ,
function generate_HiCPlotter_command(hic_files,out_file,labels,chr_string,bin_size,pst=0,ped=0,ptd=0,pdixon=0,pcustom="",hist_tracks="",hist_tracks_label="");
	
	cmd1="python /gpfs/scratch/fas/gerstein/ky26/HiCPlotter/HiCPlotter.py -f "*hic_file;
	cmd2=" -n "*labels*" -chr "*chr_string*" -r "*bin_size*" -o "*out_file;
	if pst.>0 && ped.>0
		cmd3=" -s "*string(pst)*" -e "*string(ped);
	else
		cmd3="";
	end
	if ptd.==1
		cmd4=" -ptd 1 -w 8 -tr 10 -pi 1";
	else
		cmd4="";
	end
	if pdixon.==1
		cmd5=" -pptd 1 ";
	else 
		cmd5="";
	end
	if pcustom=="";
		cmd6="";
	else
		cmd6=" -pcd 1 -pcdf "*pcustom;
	end
	if hist_tracks=="";
		cmd7="";
	else
		cmd7=" -hist "*hist_tracks*" -hl "*hist_tracks_label;
	end

	cmd=cmd1*cmd2*cmd3*cmd4*cmd5*cmd6*cmd7;

	return cmd;
end


function get_null_naive(W);
	W[isnan(W)]=0;
    p=sum(W,2)/sum(W);
    #p[i] is very much like, out of all the possible connections in the whole,
    #network, how many stubs are there from node i, essentially p_i/2M in
    #un-weighted network..
    #if position j has something to k[j] (sum(W,2)[j]) stubs to go out, the expectation is 
    #p_i*p_j*sum(W);
    #cf. this code with the next one, here. f_W=ones(size(W)) with setting dark bins to 0..

    expect=p*p'*sum(W);
    return expect;

end

function get_expect_vs_d(W);

	W[isnan(W)]=0;

	dark_bins=find(sum(W,1).==0);

	N=size(W,1);
	
	f_d=zeros(N);
	tt_d=zeros(N);

	#use the one with dark bins...
	#sum of all f_d are the reads..	
	#tt_d are the number of loci pairs separated by a distance d
	#f_d are the number of contacts between loci pairs separated by a distance d
	#f_d./tt_d is expectation..
	for d=0:N-1
		display(d);
		cd=diag(W,d);
		x=[1:N-d];
		y=[1+d:N];
		#d =1 means 1 vs 2, 2 vs 3.....
		is_okx=zeros(size(x));
		is_oky=zeros(size(y));
		for k=1:length(x)
			is_okx[k]=x[k] in dark_bins;
			is_oky[k]=y[k] in dark_bins;			
		end
		iz=find((1-is_okx).*(1-is_oky).>0);
		f_d[d+1]=sum(cd[iz]);
		tt_d[d+1]=length(iz);
	end

	expect_d=f_d./tt_d;
	expect_d[isnan(expect_d)]=0;
	#NB. we average out the coverage for all loci in this calculation..

	return expect_d;

end

#the fct is for probabilty fct only.. sum(y.*delta_x)=1; with constant delta x..
function log_bin(x,y,num_bins);
	x_tmp=logspace(log10(x[1]),log10(x[end]),2*num_bins+1);
	#x_intervals=logspace(log10(x[1]),log10(x[end]),num_bins+1);
	x_intervals=x_tmp[1:2:end];
	x_bin=x_tmp[2:2:end];
	aux=broadcast(-,x,x_intervals');
	aux=int(aux.<0);
	aux=diff(aux,2);#at this point, dim of aux is actually the number of bins
	aux[end,end]=1;
	x2bin=zeros(size(x));
	for i=1:length(x)
		x2bin[i]=find(aux[i,:].==1)[1];
	end
	y_bin=zeros(num_bins);
	for i=1:num_bins
		y_bin[i]=y[x2bin.==i][1];
	end
	x_width=diff(x_intervals);
	y_bin=y_bin./x_width;

	
	
	return x_bin, y_bin;
end

#a new fct for fitting expect_d. a general fct fitting everything may not work...
#y=Kx^-gamma;
function fit_expect_d(expect_d);

	x=[1:length(expect_d)];

	iz=find(expect_d.==0)[1];
	x1=x[1:iz-1];
	y1=expect_d[1:iz-1]

	#x1=x[expect_d.>0];
	#y1=expect_d[expect_d.>0];
	#y[y.==0]=eps();
	#fit=curve_fit(PowerFit,x,y);
	tmp=linear_fit(log10(x1),log10(y1));#this line is the same as powerfit	

	#we cannot do that..too large contribution from the leading points.
	#tmp=linear_fit(log10(x_bin),log10(y_bin));
	
	gamma=tmp[2];
	K=10^tmp[1];
	return gamma, K;

end

function get_f_d(W,expect_d)

	N=size(W,1);
	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);
	num_dark=length(dark_bins);
	N_eff=N-num_dark;
	f_W=zeros(size(W));

	for d=0:N-1
		f_W[1+d:N+1:end-d*N]=expect_d[d+1];
	end
	tmp=f_W-diagm(diag(f_W));
	f_W=f_W+tmp';
	#sum(f_W[1,:])=1 here..

	f_W[dark_bins,:]=0;
	f_W[:,dark_bins]=0;
	
	f_W=f_W/sum(f_W)*N_eff^2;
	return f_W;

end

function get_f_d_by_fitting(W,expect_d);

	N=size(W,1);
	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);
	num_dark=length(dark_bins);
	N_eff=N-num_dark;
	f_W=zeros(size(W));#what's f_W? it's a generation of ones(size(W));

	x=[1:length(expect_d)];
	gamma,K=fit_expect_d(expect_d);
	expect_d2=K*x.^gamma;
	
	#this step is added for Ren's data..
	#expect_d2[1]=expect_d[1];
	#miss this thought..

	for d=0:N-1
		f_W[1+d:N+1:end-d*N]=expect_d2[d+1];
	end
	tmp=f_W-diagm(diag(f_W));
	f_W=f_W+tmp';
	#sum(f_W[1,:])=1 here..

	f_W[dark_bins,:]=0;
	f_W[:,dark_bins]=0;
	f_W=f_W/sum(f_W)*N_eff.^2;

	return f_W;

end

function get_null_d_local(W);

	N=size(W,1);
	E_W_local=zeros(size(W));
	for i=1:round(Int,ceil(N/2));
		z=W[i,:];
		zf=flipdim(z,2);
		z1=[zf[1:N-2*i+1]' z];
		#z1=[zeros(1,N-2*int(i)+1) z];
		#z2=[fliplr(z) zeros(1,N-2*int(i)+1)];
		z2=[zf z[2*i:end]'];
		#z1[N-i+1]=z2[N-i+1];
		z_null=(z1+z2)/2;
		st=N-2*i+2;
		ed=2*N-2*i+1;
		z_null=z_null[st:ed];
		E_W_local[i,:]=z_null;
	end
	for i=round(Int,ceil(N/2))+1:N;
		z=W[i,:];
		zf=flipdim(z,2);
		#z1=[z zeros(1,2*int(i)-N-1)];
		z1=[z flipdim(z[1:2*i-N-1]',2)];
		#z2=[zeros(1,2*int(i)-N-1) fliplr(z)];
		z2=[z[1:2*i-N-1]' flipdim(z,2)];
		#z1[i]=z2[i];
		z_null=(z1+z2)/2;
		st=1;
		ed=N;
		z_null=z_null[st:ed];
		E_W_local[i,:]=z_null;
	end

	return E_W_local;
	#off set by trace(W)/2

end

#what should we used? f_W? or E_W?
#use correlation matrix?
function get_compartment_A_B(W,f_W);

	iz=find(sum(W,2).>0);
	izz=find(sum(W,2).==0);
	Wn=W[iz,iz]./f_W[iz,iz];
	C=cor(Wn);
	(U,V)=eigs(C);
	i_max=indmax(U);
	ev=V[:,i_max];
	ev_whole=zeros(size(W,1));
	ev_whole[iz]=ev;
	ev_whole[izz]=NaN;

	(loc,span)=get_chunks_v2(sign(ev_whole),1);#

	return ev_whole,loc,span;

end

function cor_spearman(x,y);
	r=cor(tiedrank(x), tiedrank(y));
	return r;
end


function bin_bigWig(bw_file,bin2loc,chr_num);

	ip=find(bin2loc[1,:].==chr_num-1);
	st=bin2loc[2,ip]';
	ed=bin2loc[3,ip]'+1;
	chr_string="chr"*string(chr_num);
	X=[repmat([chr_string],length(st),1) st ed];
	writedlm("aux.bed",X);
	run(`bwtool extract bed aux.bed $bw_file aux.out`);

	#Z=readtable(IOBuffer(readall(`cut -f 5 aux.out`)));
	display("finished bwtool");
	bin_val=zeros(size(st));
	
	f=open("aux.out");
	line=1;
	while !eof(f);
		x=readline(f);
		#println("$line $x")
		xx=split(x,"\t");
		z=xx[end];
		z=chomp(z);
		tmp=split(z,",");
		ic=find(tmp.=="");
		ic2=setdiff([1:length(tmp)],ic);
		tmp=tmp[ic2];
		ii=find(tmp.=="NA");
		ii2=setdiff([1:length(tmp)],ii);
		if ~isempty(ii2)
			bin_val[line]=mean(float(tmp[ii2]));
		else
			bin_val[line]=NaN;			
		end    	
    	line += 1
    	display(line);
 	end
	close(f);
	run(`rm aux.bed`);
	run(`rm aux.out`);

	#call the program bwtool to bin bigwig file
	#use chr2bins wtc to make a bed file..
	#bwtool agg 3:3 agg1.bed main.bigWig /dev/stdout
	#replace 3:3 by -u:u span the entrie bins..
	#https://github.com/CRG-Barcelona/bwtool/wiki/aggregate
	#bwtool agg -100:3 chr10_TADs.bed ENCFF000RSD.bigWig /de
	return bin_val;

end


function get_null_polymer(W,f_W,err_threshold);

	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);

	coverage_real=sum(W,2);
	invC=1./coverage_real;
	invC=invC[:,1];
	invC[isinf(invC)]=0;
	invC=diagm(invC);


	coverage_est=sum(W,2)/sqrt(sum(W));#sum(W) is 2N
	coverage_est[dark_bins]=0;

	#y=1./coverage_est;
	iy=coverage_est;
	
	tmp=f_W*iy;
	y_new=invC*tmp;
	coverage_est_new=1./y_new;
	coverage_est_new[isinf(coverage_est_new)]=0;
	#coverage_est_new=coverage_est_new/sum(coverage_est_new)*sqrt(sum(W)); this is not the right normalization..
	nm=sum((coverage_est_new*coverage_est_new').*f_W);
	coverage_est_new=coverage_est_new/sqrt(nm)*sqrt(sum(W));

	err=sum(abs(coverage_est_new-coverage_est));

	while err>err_threshold;
		display(err);
		coverage_est=coverage_est_new+0;
		iy=coverage_est;
		tmp=f_W*iy;
		y_new=invC*tmp;
		coverage_est_new=1./y_new;
		coverage_est_new[isinf(coverage_est_new)]=0;
		nm=sum((coverage_est_new*coverage_est_new').*f_W);
		coverage_est_new=coverage_est_new/sqrt(nm)*sqrt(sum(W));
		#coverage_est_new=coverage_est_new/sum(coverage_est_new)*sqrt(sum(W));
		err=sum(abs(coverage_est_new-coverage_est));
	end

	E_W=(coverage_est_new*coverage_est_new').*f_W;

	return coverage_est_new,E_W;
end

#it look at the enrichment of contacts over a polymer model.
#0 means no enrichment. value are -log10(P) based on a Poisson distribution.
function get_P_value_observed_vs_expected(W,E_polymer);
	all_P=ones(size(W));
	iz=find(sum(W,2).>0);
	for i=1:length(iz);
		display(i);
		for j=i:length(iz);
			if E_polymer[iz[i],iz[j]].>0 && W[iz[i],iz[j]].>E_polymer[iz[i],iz[j]]
				dd=Poisson(E_polymer[iz[i],iz[j]]);
				lw=cdf(dd,floor(W[iz[i],iz[j]]));
				uw=cdf(dd,ceil(W[iz[i],iz[j]]));
				delta=W[iz[i],iz[j]]-floor(W[iz[i],iz[j]]);
				cc=1-(lw+(uw-lw)*delta);
				all_P[iz[i],iz[j]]=cc;
				all_P[iz[j],iz[i]]=cc;
			end
		end
	end
	max_log_P=-floor(log10(minimum(all_P[all_P.>0])))+1;
	enrich=-log10(all_P);
	enrich[isinf(enrich)]=max_log_P;
	
	return enrich;
end

#Given A as the interaction matrix,
function construct_GPC(A,Gene,G2B_dict,dark_bins,file_loc,bin_size,cell_type);
	for i=1:length(Gene)
		display(i)
		m1=G2B_dict[Gene[i]];
		if length(m1).==1
			if m1 in dark_bins
				b1=[];
			else
				b1=[m1];
			end
		else
			b1_aux=Array(Bool,m1[2]-m1[1]+1,1);
			b1=[m1[1]:m1[2]];
			for k1=1:length(b1_aux)
				b1_aux[k1]=b1[k1] in dark_bins;
			end
			b1=b1[~squeeze(b1_aux',1)];
		end

		for j=i:length(Gene);
			m2=G2B_dict[Gene[j]];
			
			if length(m2).==1
				if m2 in dark_bins
					b2=[];
				else
					b2=[m2];
				end
			else
				b2_aux=Array(Bool,m2[2]-m2[1]+1,1);
				b2=[m2[1]:m2[2]];
				for k2=1:length(b2_aux)
					b2_aux[k2]=b2[k2] in dark_bins;
				end
				b2=b2[~squeeze(b2_aux',1)];
			end

			if isempty(b1) || isempty(b2)
				GPC[i,j]=NaN;
				GPC[j,i]=NaN;
			elseif length(b1)==1 && length(b2)==1
				e1=(A[b1,b2]+A[b2,b1])/2;
				e2=(A[b1,b2]+A[b2,b1])/2;
				GPC[i,j]=e1[1];
				GPC[j,i]=e2[1];
			else
				Atmp=(A[b1,b2]+A[b2,b1]')/2;
				GPC[i,j]=maximum(Atmp);
				GPC[j,i]=maximum(Atmp);
				#not sure..max is the closet distance, 
				#GPN[i,j]=mean(Atmp);
				#GPN[j,i]=mean(Atmp);
			end
		end
	end
	#the only case GPN is NaN is, either genes belong completely to the dark bins..
	out_file=file_loc*cell_type*"_GPC_"*bin_size*".jld";
	save(out_file,"GPC",GPC);
	return GPC;
end

#Genes in dark bins are gone..return the existing genes, and mm as the minimal counts
function get_GPN(GPC,Gene,file_loc,cell_type,bin_size);
	i_mg=find(sum(~isnan(GPC),1).>0);
	GPC=GPC[i_mg,i_mg];
	GPC=GPC-diagm(diag(GPC));
	mm=minimum(GPC[GPC.>0]);
	d=floor(log10(mm));
	GPN=log10(GPC+10^d);
	GPN[GPN.<=0]=0;
	out_file=file_loc*cell_type*"_GPN_"*bin_size*".jld";
	save(out_file,"i_mg",i_mg,"mm",mm,"GPN",GPN);
	return i_mg, mm, GPN;
end

function get_GPN2(GPC,Gene,file_loc,cell_type,bin_size);
	i_mg=find((sum(~isnan(GPC),1).>0)'.*(diag(GPC).>0));
	#filter genes with no count to itself
	GPC=GPC[i_mg,i_mg];
	D=diagm(1./sqrt(diag(GPC)));
	GPC=D*GPC;
	GPC=GPC*D;
	GPN=GPC-diagm(diag(GPC));
	return GPN, i_mg;
end

function proximty_expr_energy(x,CIN);
	x_tmp=CIN*x;
	E=sum(-x.*x_tmp)/2;
	return E;
end

function proximty_energy_simulation(x,CIN,iter);
	
	DeltaE_iter=zeros(iter);
	E_init=proximty_expr_energy(x,CIN);
	N=length(x);
	x_new=x+0;

	for k=2:iter;
		display(k);
		a=find(x.>0);
		a=a[randperm(length(a))];
		a=a[1];
		b=find(x.<0);
		b=b[randperm(length(b))];
		b=b[1];
		x_new[a]=-1*x[a];
		tmp=CIN[:,a];
		tmp[a]=0;
		DE1=-sum(tmp.*x)*(x_new[a]-x[a]);#this is based on the eneergy defined above. NB the factor 1/2 and -ve
		x_new[b]=-1*x[b];
		tmp=CIN[:,b];
		tmp[b]=0;
		DE2=-sum(tmp.*x)*(x_new[b]-x[b]);
		DeltaE_iter[k]=DE1+DE2;
		x=x_new+0;
	end
	E_iter=E_init+cumsum(DeltaE_iter);

	return E_iter,x_new;
end


function proximty_energy_biased_simulation(x,CIN,iter);
	
	DeltaE_iter=zeros(iter);
	E_init=proximty_expr_energy(x,CIN);
	N=length(x);
	x_aux=x+0;
	x_new=x+0;

	for k=2:iter;
		display(k);
		a=find(x.>0);
		a=a[randperm(length(a))];
		a=a[1];
		b=find(x.<0);
		b=b[randperm(length(b))];
		b=b[1];
		x_aux[a]=-1*x[a];
		tmp=CIN[:,a];
		tmp[a]=0;
		DE1=-sum(tmp.*x)*(x_aux[a]-x[a]);#this is based on the eneergy defined above. NB the factor 1/2 and -ve
		x_aux[b]=-1*x[b];
		tmp=CIN[:,b];
		tmp[b]=0;
		DE2=-sum(tmp.*x)*(x_aux[b]-x[b]);
		if DE1+DE2 <0
			x_new=x_aux+0;
			DeltaE_iter[k]=DE1+DE2;
			x=x_new+0;
		else
			x_aux=x+0;
		end
	end
	E_iter=E_init+cumsum(DeltaE_iter);
	return E_iter,x_new;
end

#function proximty_expr_energy_normalized(x,CIN);
#	Ksroot=diagm(vec(1./sqrt(sum(CIN,2)));
#	CIN_normalized=Ksroot*CIN;
#	CIN_normalized=CIN_normalized*Ksroot;
#	x_tmp=CIN*x;
#	E=sum(-x.*x_tmp)/2;
#	return E;
#end

function generate_net_file(CINr,cutoff,outfile)
	Net=CINr.*(CINr.>cutoff);
	Net=sparse(Net);
	a=find(triu(Net));
	(u,v)=ind2sub(size(Net),a);
	w=log10(full(Net[a]));
	df=DataFrame();
	df[:Source]=u;
	df[:Target]=v;
	w=squeeze(w,2);
	df[:Stength]=w;
	df[:Type]="undirected";
	writetable(outfile, df);
end

function get_diffusion_mode(CIN,n,i_mg,L);
	
	K=sum(CIN,1)';
	sqK=1./sqrt(K);
	sqK=squeeze(sqK,2);
	sqK=spdiagm(sqK);
	CIN=sparse(CIN);
	S=sqK*CIN*sqK;
	(ev,evec)=eigs(S,nev=n);
	for i=1:size(evec,2);
    	evec[:,i]=sqK*evec[:,i];
    	M=sum((evec[:,i].^2));
    	evec[:,i]=evec[:,i]./sqrt(M);
	end
	ev=real(ev);
	evec=real(evec);
	ipr=zeros(size(evec,2));
	evec_tmp=zeros(L,n);
	for i=1:size(evec,2);
		evec_tmp[:,i]=full(sparse(i_mg,int(ones(size(i_mg))),evec[:,i],L,1));
	end
	for i=1:length(ipr)
		ipr[i]=1/sum(evec_tmp[:,i].^4)
	end
	return ev, evec_tmp, ipr;

end

function change_chr(chr)
	chr2=cell(size(chr));
	for i=1:length(chr)
		if chr[i]==25;
			chr2[i]="chrM";
		elseif chr[i]==24;
			chr2[i]="chrY";
		elseif chr[i]==23;
			chr2[i]="chrX";
		elseif chr[i]==22;
			chr2[i]="chr22";
		elseif chr[i]==21;
			chr2[i]="chr21";
		elseif chr[i]==20;
			chr2[i]="chr20";
		elseif chr[i]==19;
			chr2[i]="chr19";
		elseif chr[i]==18;
			chr2[i]="chr18";
		elseif chr[i]==17;
			chr2[i]="chr17";
		elseif chr[i]==16;
			chr2[i]="chr16";
		elseif chr[i]==15;
			chr2[i]="chr15";
		elseif chr[i]==14;
			chr2[i]="chr14";
		elseif chr[i]==13;
			chr2[i]="chr13";
		elseif chr[i]==12;
			chr2[i]="chr12";
		elseif chr[i]==11;
			chr2[i]="chr11";
		elseif chr[i]==10;
			chr2[i]="chr10";
		elseif chr[i]==9;
			chr2[i]="chr9";
		elseif chr[i]==8;
			chr2[i]="chr8";
		elseif chr[i]==7;
			chr2[i]="chr7";
		elseif chr[i]==6;
			chr2[i]="chr6";
		elseif chr[i]==5;
			chr2[i]="chr5";
		elseif chr[i]==4;
			chr2[i]="chr4";
		elseif chr[i]==3;
			chr2[i]="chr3";
		elseif chr[i]==2;
			chr2[i]="chr2";
		elseif chr[i]==1;
			chr2[i]="chr1";
		end
	end
	return chr2;
end

function isorank_distance(A,B);
	R_AA,match_AA=isorank(A,A,0);
	R_BB,match_BB=isorank(B,B,0);
	R_AB,match_AB=isorank(A,B,0);
	D=abs(log10(diag(R_AB./R_AA)))+abs(log10(diag(R_AB./R_BB)))
	return D;
end


function isorank(A,B,alpha,maxiter=1000);
	
	VA=size(A,1);
	VB=size(B,2);
	N=VA*VB;
	A=sparse(A);
	B=sparse(B);
	TA=spdiagm(vec(1./sum(A,1)))*A;
	TB=spdiagm(vec(1./sum(B,1)))*B;
	# X=spzeros(N,N);
	# id=[1:N];#R matrix	
	# for k=1:length(id)
	# 	display(k)
	# 	(i,j)=ind2sub((VA,VB),k)
	# 	i_n=A[i,:]/sum(A[i,:]);
	# 	j_n=B[:,j]/sum(B[:,j]);
	# 	Z=j_n*i_n;
	# 	nz=find(j_n*i_n);#vertical, from 1,1, ones(1,N) is the left evec..with ev=1
	# 	for k2=1:length(nz);
	# 		(b,a)=ind2sub((VA,VB),nz[k2]);#a in A and b in B.. Note the order..
	# 		X[k,sub2ind((VA,VB),a,b)]=A[i,a]/sum(A[a,:])*B[j,b]/sum(B[b,:]);
	# 		#sub2ind((VA,VB),j,a);
	# 		#sub2ind((VA,VB),b,i);
	# 	end
	# end
	# the block above could be replaced by Kronecker product
	X=kron(TA,TB)'; #sum col of X is one...same as PNAS paper...

	R=ones(N,1);
	err=100;iter=0;
	while (err>1e-10)&(iter<maxiter)
		nextR=(1-alpha)*X*R+alpha;
		display(err)	
		err=sum((nextR-R).^2);
		R=nextR;
		iter=iter+1; 
	end
	#sum(R) is conserved...if no dead end in diffusion..

	# R=ones(N,1);R=R/norm(R);
	# err=100;
	# while err>1e-10
	# 	nextR=X*R/norm(X*R);
	# 	display(err)	
	#  	err=sum((nextR-R).^2);
	#  	R=nextR;
	# end

	R=reshape(R,(VA,VB));
	R_tmp=R+0;
	A_ind=[1:size(R,1)];
	B_ind=[1:size(R,1)];
	match=zeros(size(R,1),2);
	for i=1:size(R,1);
		(tmp1,tmp2)=ind2sub(size(R_tmp),find(R_tmp.==maximum(R_tmp)));
		tmp1=tmp1[1];
		tmp2=tmp2[1];
		match[i,1]=A_ind[tmp1];
		match[i,2]=B_ind[tmp2];
		R_tmp=R_tmp[[1:tmp1-1,tmp1+1:end],:]
		R_tmp=R_tmp[:,[1:tmp2-1,tmp2+1:end]];
		A_ind=deleteat!(A_ind,tmp1);
		B_ind=deleteat!(B_ind,tmp2);
	end

	return R,match;

end


function get_candidate_domains(W,E_W,res);

	i_dark=find(sum(W,2).==0);

	final_assign1, Q1, Brenorm1=optimize_TADs_modlouvain(W,E_W,res,1);
	final_assign2, Q2, Brenorm2=optimize_TADs_modlouvain(W,E_W,res,-1);
  
    Z=(broadcast(-,final_assign1,final_assign1').==0)+(broadcast(-,final_assign2,final_assign2').==0);
    Z[i_dark,i_dark]=0;
    Z=(Z.==2);
    g=matrix_to_graph(Z);
    all_modules_aux=Graphs.connected_components(g);
    bins2modules=zeros(size(B,1),1);
    ct=1;
    #this method is not perfect, if a few contin. seqments form  a big modules. onky the biggest
    #will be picked..on the other hand, some noise is gone..to improve, we should rewrite iterate_TADs_modlouvain
    for i=1:length(all_modules_aux)
    	tmp=all_modules_aux[i];
    	tmp=sort(tmp);
    	a=diff(tmp);
    	a = [NaN, a, NaN];
		b=[diff(a).!=0];
		c= diff(b);
		st=find(c.==-1);
		ed=find(c.==1)+1;
		if ~isempty(st) && ~isempty(ed);
			L=ed-st;
			mm=findmax(L);
			bins2modules[tmp[[st[mm[2]]:ed[mm[2]]]]]=ct;
			ct=ct+1;
        end
    end
    i_dark=find(sum(W,1).==0);
    bins2modules[i_dark]=0;

    return bins2modules;
end

#using the modified Louvain update, the TADs are easier to get
#final_assign=-1, will not be assigned to domain
#but this code further extract certain sig. compoentns...
#bins2modules=0 means not assigned..
function get_candidate_domains_v2(W,E_W,res,trial=1,accept=1);

	Z=zeros(size(W));
	for x=1:trial;
		final_assign_x, Q1, Brenorm1=optimize_TADs_modlouvain(W,E_W,res);
		i_undeter=find(final_assign_x.<0);
		Z_tmp=(broadcast(-,final_assign_x,final_assign_x').==0);
		Z_tmp[i_undeter,i_undeter]=0;
		Z=Z+Z_tmp;
	end

    Z=(Z.>=accept);
    Z=Z-diagm(diag(Z));
    g=matrix_to_graph(Z);
    all_modules_aux=Graphs.connected_components(g);
	#k=sum(Z,1);   
    bins2modules=zeros(size(W,1),1);
    ct=1;
    for i=1:length(all_modules_aux)
    	tmp=all_modules_aux[i];
    	if sum(Z[tmp,tmp])>0
    		bins2modules[tmp[1]:tmp[end]]=ct;ct=ct+1;
    	end
    end

    return bins2modules;
end


#W is the observed - 
#E_W is null NOT w_i*w_j/sum(W)
#but sum(W)=sum(E_W)
function optimize_TADs_modlouvain(W,E_W,res,order=0);
#B is the generalized modularity matrix, and sW is sum of all elts in W, ~ to 2M
#perform the step of louvain
#but the do not define renormalization modularity....
#order=1, forward
#order=-1, backward
#order=0, random
	
	N=size(W,1);
    i_no_dark=find(sum(abs(W),1).>0);
    N_no_dark=length(i_no_dark);

    B=W-E_W*res;
    Bcompact=B[i_no_dark,i_no_dark];
    Wcompact=W[i_no_dark,i_no_dark];
    sW=sum(W);

    #(assign,Q,Brenorm)=iterate_TADs_modlouvain(Bcompact,sW,order);
    (assign,Q,Brenorm)=iterate_TADs_modlouvain_v2(Bcompact,sW);
    transfer=sparse(collect(1:N_no_dark),assign,ones(size(assign)));

    keep_doing=1;

    while keep_doing==1
        #(tmp_assign,tmp_Q,tmp_Brenorm)=iterate_TADs_modlouvain(Brenorm,sW,order);
        (tmp_assign,tmp_Q,tmp_Brenorm)=iterate_TADs_modlouvain_v2(Brenorm,sW);
        tmp_transfer=sparse(collect(1:size(tmp_assign,1)),tmp_assign,ones(size(tmp_assign)));
        if isequal(tmp_transfer,speye(length(tmp_assign)))
        	keep_doing=0;
        	#Brenorm, Wrenorm are optimal.
        else
	        transfer=transfer*sparse(collect(1:size(tmp_assign,1)),tmp_assign,ones(size(tmp_assign)));
        	Brenorm=tmp_Brenorm+0;
        	Q=tmp_Q+0;
        end
    end

    (u,v)=findn(transfer);
    iu=sortperm(u);
    tmp_assign=v[iu];
    final_assign=zeros(N);
    final_assign[i_no_dark]=tmp_assign;
    #pick up the zeros....dark bins...where should they belong?
    #look at two sides, if both sides are the same, merge them, if not, keep dark...
    (loc,span)=get_chunks_v2(final_assign,1);#

    if final_assign[loc[1]].==0;
    	final_assign[loc[1]:loc[1]+span[1]-1]=-1;
    	loc=loc[2:end];
    	span=span[2:end];
   	end
   	if final_assign[loc[end]+span[end]-1].==0;
   		final_assign[loc[end]:loc[end]+span[end]-1]=-1;
   		loc=loc[1:end-1];
   		span=span[1:end-1];
   	end

    for i=1:length(loc)
    	if final_assign[loc[i]].==0
    		if final_assign[loc[i]-1]==final_assign[loc[i]+span[i]];
    			final_assign[loc[i]:loc[i]+span[i]-1]=final_assign[loc[i]-1];
    		else
    			final_assign[loc[i]:loc[i]+span[i]-1]=-1;
    		end
    	end
    end

    return final_assign, Q, Brenorm;

end

#do an improved version.
#update only the adjacent bins, use random order
#if adj. bin is dark, move further..because the values (lower for further away) 
#should take care of the problem..
#the same for renormalization...at the end, things should not be broken..
function iterate_TADs_modlouvain_v2(Bcompact,sW);

    Nb=size(Bcompact,1);
    sigma=collect(1:Nb);
	u=collect(1:Nb);
	u=u[randperm(Nb)];
   
    gain=1;Niter=1;

    while (gain==1)
        gain = 0;
        for j=1:Nb
            x=u[j];
            spin=sigma[x];
            if x==1
            	spin_f=sigma[x+1];
            	spin_b=sigma[x];
            elseif x==Nb;
            	spin_f=sigma[x];
            	spin_b=sigma[x-1];
            else 
            	spin_f=sigma[x+1];
            	spin_b=sigma[x-1];
            end
            c=Bcompact[x,:];
            c[x]=0;#this is important step to make sure the deltaQ is right 
            neighbors_spin=sigma;
            DeltaQ=-sum(c'.*(sigma.==spin))+full(sparse(neighbors_spin,[1 for dd=1:Nb],vec(c)));
            #the 2nd term sum over the components from each community in advance
            #1st term, the effect of getting rid of the original spin contribution..
            #note the dim of DeltaQ is the number of communities
            spin_choice=[spin_b spin spin_f];
            id=indmax(DeltaQ[spin_choice]);#choose 1 out of 3...
            id=spin_choice[id];
            new_spin=id;
            if (new_spin!=spin)&(DeltaQ[id].>0);
                gain=1;
                sigma[x]=new_spin;
            end
        end
        Q=compute_modularity(sigma,Bcompact,sW);
        @printf("iteration %d - sum = %f %d communities \n",Niter,Q,length(unique(sigma)))
        Niter = Niter + 1
    end

    Q=compute_modularity(sigma,Bcompact,sW);

    sigma2=relabel_communities(sigma);
    usigma=sort(unique(sigma));
    N_renorm=length(usigma);

    #we abandon the idea of renormalization..
  	Bcompact_renorm=zeros(N_renorm,N_renorm);

    for i=1:N_renorm
        for j=1:N_renorm
            Bcompact_renorm[i,j]=sum(sum(Bcompact[sigma.==usigma[i],sigma.==usigma[j]]));
        end
    end

    Q2=compute_modularity(collect(1:N_renorm),Bcompact_renorm,sW);
    @printf("step - Q = %f %d communities \n",Q2,length(unique(sigma)))
    #println(Q2), verified to be the same as Q...
    return sigma2, Q2, Bcompact_renorm;

end

function iterate_TADs_modlouvain(Bcompact,sW,order);

    Nb=size(Bcompact,1);
    sigma=collect(1:Nb);

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
        @printf("iteration %d - sum = %f %d communities \n",Niter,Q,length(unique(sigma)))
        Niter = Niter + 1
    end

    Q=compute_modularity(sigma,Bcompact,sW);

    sigma2=relabel_communities(sigma);
    usigma=sort(unique(sigma));
    N_renorm=length(usigma);

    #we abandon the idea of renormalization..
  	Bcompact_renorm=zeros(N_renorm,N_renorm);

    for i=1:N_renorm
        for j=1:N_renorm
            Bcompact_renorm[i,j]=sum(sum(Bcompact[sigma.==usigma[i],sigma.==usigma[j]]));
        end
    end

    Q2=compute_modularity(collect(1:N_renorm),Bcompact_renorm,sW);
    @printf("step - Q = %f %d communities \n",Q2,length(unique(sigma)))
    #println(Q2), verified to be the same as Q...
    return sigma2, Q2, Bcompact_renorm;

end



#W is the observed - 
#E_W is null in w_i*w_j/sum(W)
function optimize_network_modularity_louvain(W,E_W,res,order)
#B is the generalized modularity matrix, and sW is sum of all elts in W, ~ to 2M
#perform the step of louvain, returns renormalized networks..
#order=1, forward
#order=-1, backward
#order=0, random
	
	N=size(W,1);
    i_no_dark=find(sum(abs(W),1).>0);
    N_no_dark=length(i_no_dark);

    B=W-E_W*res;
    Bcompact=B[i_no_dark,i_no_dark];
    Wcompact=W[i_no_dark,i_no_dark];
    sW=sum(W);

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
        @printf("iteration %d - sum = %f %d communities \n",Niter,Q,length(unique(sigma)))
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
    @printf("step - Q = %f %d communities \n",Q2,length(unique(sigma)))
    #println(Q2), verified to be the same as Q...
    return sigma2, Q2, Brenorm, Wrenorm;
end

function get_C_corrected(W,E_polymer);

	N=size(W,1);
    i_no_dark=find(sum(abs(W),1).>0);
    N_no_dark=length(i_no_dark);
	Wcompact=W[i_no_dark,i_no_dark];
	E_polymer_compact=E_polymer[i_no_dark,i_no_dark];

    C=cor(Wcompact);
    C_null=cor(E_polymer_compact);
    Cnorm=sum(C);

    lam_null_max=eigmax(C_null);
    (U,V)=eig(C);
    i_pick=find(U.>lam_null_max);

    C_corrected=[];

    if ~isempty(i_pick)

	    C_corrected=zeros(size(C));
    	for k=1:length(i_pick);
    		i=i_pick[k];
    		C_corrected=U[i]*V[:,i]*V[:,i]'+C_corrected;
    	end
    
    end

    return C_corrected;

end


#W is the observed 
#E_polymer is null
#instead of looking at raw counts, we define interactions by correlation
function optimize_correlation_modularity_louvain(W,E_polymer,order)

	N=size(W,1);
    i_no_dark=find(sum(abs(W),1).>0);
    N_no_dark=length(i_no_dark);
	Wcompact=W[i_no_dark,i_no_dark];
	E_polymer_compact=E_polymer[i_no_dark,i_no_dark];

    C=cor(Wcompact);
    C_null=cor(E_polymer_compact);
    Cnorm=sum(C);

    lam_null_max=eigmax(C_null);
    (U,V)=eig(C);
    i_pick=find(U.>lam_null_max);

    if ~isempty(i_pick)

	    C_corrected=zeros(size(C));
    	for k=1:length(i_pick);
    		i=i_pick[k];
    		C_corrected=U[i]*V[:,i]*V[:,i]'+C_corrected;
    	end

	    num_run=1;
    	(assign,Q,C_corrected_renorm)=iterate_correlation_modularity(C_corrected,Cnorm,order);
    	transfer=sparse(collect(1:N_no_dark),assign,ones(size(assign)));

	    keep_doing=1;

	    while keep_doing==1
    	    (tmp_assign,tmp_Q,tmp_C_corrected_renorm)=iterate_correlation_modularity(C_corrected_renorm,Cnorm,order);
        	tmp_transfer=sparse(collect(1:size(tmp_assign,1)),tmp_assign,ones(size(tmp_assign)));
        	if isequal(tmp_transfer,speye(length(tmp_assign)))
        		keep_doing=0;
        		#Brenorm, Wrenorm are optimal.
        	else
	        	transfer=transfer*sparse(collect(1:size(tmp_assign,1)),tmp_assign,ones(size(tmp_assign)));
        		C_corrected_renorm=tmp_C_corrected_renorm+0;
        		Q=tmp_Q+0;
        	end
    	end

	    (u,v)=findn(transfer);
    	iu=sortperm(u);
    	tmp_assign=v[iu];
    	final_assign=zeros(N);
    	final_assign[i_no_dark]=tmp_assign;

    else
    	final_assign=[];
    	Q=[];
    	C_corrected_renorm=[];
    end

    return final_assign, Q, C_corrected_renorm;

end

function iterate_correlation_modularity(C_corrected,Cnorm,order)

	Nb=size(C_corrected,1);
	sigma=collect(1:Nb);

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
            spin=sigma[x];
            c=C_corrected[x,:];
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
        Q=compute_modularity(sigma,C_corrected,Cnorm);
        @printf("iteration %d - sum = %f %d communities \n",Niter,Q,length(unique(sigma)))
        Niter = Niter + 1
    end

	Q=compute_modularity(sigma,C_corrected,Cnorm);

    sigma2=relabel_communities(sigma);
    usigma=sort(unique(sigma));
    N_renorm=length(usigma);

    #this is the normalization required for this null model, to make sense the Q remains the same
    #Wrenorm=zeros(N_renorm,N_renorm);

    #for i=1:N_renorm
    #    for j=1:N_renorm
    #        Wrenorm[i,j]=sum(sum(Wcompact[sigma.==usigma[i],sigma.==usigma[j]]));
    #    end
    #end
    #w=sum(Wrenorm,2);

    #E_Wrenorm=(w*w')/sW;

    #Brenorm=Wrenorm-E_Wrenorm;

    #correlation based modularity use a different renormalization procedure..

    C_corrected_renorm=zeros(N_renorm,N_renorm);

    for i=1:N_renorm
        for j=1:N_renorm
            C_corrected_renorm[i,j]=sum(sum(C_corrected[sigma.==usigma[i],sigma.==usigma[j]]));
        end
    end
    #Cnorm is preserved

    Q2=compute_modularity(collect(1:N_renorm),C_corrected_renorm,Cnorm);
    #@printf("step - Q = %f %d communities \n",Q2,length(unique(sigma)))
    #println(Q2), verified to be the same as Q...
    return sigma2, Q2, C_corrected_renorm;
end

# We tried to use C_corrected to further deoompse, C_corrected[assign.==1,assign.==1] cannot be further decomposed though it looks many sub-division from image
function optimize_correlation_modularity_louvain_recursive(W,E_polymer,order);
  
    #N=size(W,1);
    i_no_dark=find(sum(abs(W),1).>0);
    N_no_dark=length(i_no_dark);
	Wcompact=W[i_no_dark,i_no_dark];
	E_polymer_compact=E_polymer[i_no_dark,i_no_dark];

    N=size(Wcompact,1);
    iz=collect(1:N);#iz are loci

    node_id=1;#root
    nodes_queue=[node_id]

    communities_tree=Dict{Int64,Int64}();
    communities=Dict{Int,Array}();
    communities[node_id]=[1:N];
      
    while ~isempty(nodes_queue)

        active_ind=nodes_queue[1];
        display(active_ind)
        iz=communities[active_ind];
        assign, Q, C_corrected_renorm=optimize_correlation_modularity_louvain(Wcompact[iz,iz],E_polymer_compact[iz,iz],1);
        u_s=unique(assign);
        u_s=u_s[u_s.>0];
    
        if length(u_s)>1
            for i=1:length(u_s)
                node_id=node_id+1;#node_id keeps iterating
                communities_tree[node_id]=active_ind;
                communities[node_id]=iz[assign.==i];
                push!(nodes_queue,node_id);
            end
        end

        nodes_queue=setdiff(nodes_queue,active_ind);
        display(nodes_queue);

    end

    return communities, communities_tree
        
end

function louvain_recursive_hierarchy(communities, communities_tree)

	N=length(communities);
    g = simple_graph(N);

    for key in keys(communities_tree);
    	add_edge!(g,key,communities_tree[key]);
    end

    g2 = simple_graph(N);

    for key in keys(communities_tree);
    	add_edge!(g2,communities_tree[key],key);
    end

    N=num_vertices(g);

    kin=[in_degree(x,g) for x in vertices(g)]

    leaves=find(kin.==0);

    comty=Dict{Int64,Any}();
    
    for j=1:length(leaves)
        comty[j]=communities[leaves[j]]
    end

    r = dijkstra_shortest_paths(g2, 1);#distance to all vertices from root
    level=r.dists;

    return comty,g,level;
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

#id is the starting loc of a chunk, and d is the length it spans..
function get_chunks_v2(a,singleton=0);
	# adopt from ammatlab code by Jiro Doke;
	 a                 = [NaN, a, NaN];
	 b                 = diff(a);
	 b1                = b;  # to be used in fullList (below)
	 ii                = trues(size(b));
	 ii[b.==0] = false;
	 b[ii]             = 1;
	 c                 = diff(b);
	 id                = find(c.==-1);

	 #Get single-element chunks also
	 if singleton.==1
	 	b1[id]          = 0;
	 	ii2             = find(b1[1:end-1]);
	 	d               = vcat(find(c.==1) - id + 1, ones(length(ii2)));
	 	id              = [id,ii2];
	 	
	 	v=sortperm(id);
	 	id=sort(id);
	 	#(id,tmp)        = sort(id);
	 	d               = d[v];
	 else
	 	d               = find(c.==1) - id + 1;
	 end

	 return id,d;
end


#function get_chunks(bins2modules)#this code doesmn't get chunks less than or equal to length 2...
#    a=diff(bins2modules);
#    a = [NaN, a, NaN];#
#	b=[diff(a).!=0];
#	c= diff(b);
#	st=find(c.==-1);
#	ed=find(c.==1)+1;
#	return st, ed
#end

#this code just report, no more filtering..should be identical to bins2modules..
#those 0 in bins2modules are not reported..
function report_domains(chr2bins,bin2loc,chr_num,bins2modules)

	u,v=get_chunks_v2(bins2modules,1);
	TAD_st=Int64[];
	TAD_ed=Int64[];

	for i=1:length(u);
		if bins2modules[u[i]]>0
			push!(TAD_st,u[i])
			push!(TAD_ed,u[i]+v[i]-1);
		end
	end
    
    TADs_list=DataFrame(chr=ASCIIString[],domain_st=Int64[],domain_ed=Int64[],domain_st_bin=Int64[],domain_ed_bin=Int64[],idx=Int64[]);

    st=chr2bins[:,chr_num][1]+1;#the bin count from zero in the stored file.
    ed=chr2bins[:,chr_num][2]+1;#we here shift it..
    chr2all_bin=[st:ed];
    #this will be the array mapped by elements in the chr of interest.

    if chr_num<=22
        chr_string=string("chr",string(chr_num));
    elseif chr_num==23
        chr_string=string("chr","X");
    elseif chr_num==24 
        chr_string=string("chr","Y");
    end;

    TAD_st_bin=chr2all_bin[TAD_st];
    TAD_ed_bin=chr2all_bin[TAD_ed];   

    for i=1:size(TAD_st_bin,1)
        loc1=bin2loc[2,TAD_st_bin[i]];
        loc2=bin2loc[3,TAD_ed_bin[i]];
        push!(TADs_list,[chr_string loc1 loc2 TAD_st_bin[i] TAD_ed_bin[i] i]);
    end

    return TADs_list
end

#from genome coordinate to bins...
function loc2bin(chr_num,loc,bin2loc)
    p1=find((bin2loc[1,:].==chr_num-1).*(bin2loc[2,:].<loc))
    p2=find((bin2loc[1,:].==chr_num-1).*(bin2loc[3,:].>=loc))
    bin_go=intersect(p1,p2);
    return bin_go;
end

function allTADs2Mat(TADs_list,chr_num,chr2bins,bin2loc);
    
    stc=chr2bins[:,chr_num][1]+1;#the bin count from zero in the stored file.
    edc=chr2bins[:,chr_num][2]+1;#we here shift it..
    s=zeros(edc-stc+1,1);
    size_chr=length(s);

    chr2all_bin=[stc:edc];
    #this will be the array mapped by elements in the chr of interest.

    for i=1:size(TADs_list,1);
        stm=find(chr2all_bin.==TADs_list[i,4])[1];
        edm=find(chr2all_bin.==TADs_list[i,5])[1];
        s[stm:edm]=TADs_list[i,6];
    end

    iz=find(s.>0);
    s2=s[iz];
   
    Z=(broadcast(-,s2,s2').==0)+0;
    
    allTADs=extend_mat(Z,iz,size_chr);
    #imagesc(Z)

    return allTADs;

end

function zoom_in(W,chr_num,loc_st,loc_ed,chr2bins,bin2loc);
    st=loc2bin(chr_num,loc_st,bin2loc)[1];
    ed=loc2bin(chr_num,loc_ed,bin2loc)[1];
    i_pick=[st:ed];
    if chr_num>1
        bd=chr2bins[2,chr_num-1];
    else
        bd=0;
    end
    i_pick=i_pick-bd;
    W_pick=W[i_pick,i_pick];
    return W_pick;
end

function extend_mat(Z,iz,L);
    (u,v)=ind2sub(size(Z),find(Z.!=0));
    w=Z[find(Z)];
    #w=nonzeros(Z);
    u=iz[u];
    v=iz[v];
    Z_extend=sparse(u,v,w,L,L);
    Z_extend=full(Z_extend);
    return Z_extend;
end

function display_matrix(W,chr,pic_st,pic_ed,fig_out);
    p=imagesc(W);
    setattr(p.x1,label=chr,draw_ticks=false,ticks=2,ticklabels=[string(pic_st),string(pic_ed)])
    setattr(p.y1,label="",draw_ticks=false,ticklabels=[]);
    #setattr(p.y2,label="chr22",draw_ticks=false,ticklabels=[string(ed),string(st)]);
    p;
    savefig(fig_out);
end

function output_TADs(file_name,TADs_list)
	writetable(file_name,TADs_list)
end

function report_boundaries(TADs_list)
    
    chr_string=TADs_list[:chr][1];
    TADs_boundaries=DataFrame(chr=ASCIIString[],Bst=Int64[],Bend=Int64[]);
    # forget the 1st ...
    for i=2:size(TADs_list,1)-1;
        if TADs_list[i+1,4]-TADs_list[i,5]<=1 #Dixon file organized as 0...
            push!(TADs_boundaries,[chr_string TADs_list[i+1,2] TADs_list[i+1,2]]);
        elseif TADs_list[i+1,4]-TADs_list[i,5]>1
            push!(TADs_boundaries,[chr_string TADs_list[i,3] TADs_list[i+1,2]]);
        end
    end
    return TADs_boundaries;

end

#bigWigToWig -chrom=chr22 -start=15620000 -end=16620000  

function get_peak_density_near_TADs_boundaries(TAD_list,TADs_boundaries,bed_file,boundary_length_thres,no_small=0);
    bin_size=40000;
    bb_lim=600000;
    df=readtable(bed_file,separator='\t',header=false);
    peak_loc=df[2]+round((df[3]-df[2])/2);
    peak_chr=df[1];
    peak_density=zeros(int(bb_lim/bin_size)*2,1);
    L_boundaries=TADs_boundaries[:Bend]-TADs_boundaries[:Bst];
    size_TAD=TAD_list[:domain_ed]-TAD_list[:domain_st];
    #all_b=zeros(size(TADs_boundaries,1),size(peak_density,1));
    size_TAD=size_TAD[2:end-1];
    if no_small.==1;
	    TADs_boundaries2=TADs_boundaries[size_TAD.>bb_lim,:];
    else
    	TADs_boundaries2=TADs_boundaries;
    end

    for bb=1:size(TADs_boundaries2,1)
    	if (L_boundaries[bb]>=boundary_length_thres);
	        bb_chr=TADs_boundaries2[bb,:][1][1];
    	    bb_loc1=TADs_boundaries2[bb,:][2][1];
        	bb_loc2=TADs_boundaries2[bb,:][3][1];
        	bb_loc=floor((bb_loc1+bb_loc2)/2);
        	iz=find(peak_chr.==bb_chr);
        	pl=peak_loc[iz];
        	pl=sort(pl);
        	bb_bins=[bb_loc-bb_lim:bin_size:bb_loc+bb_lim];
        	#how to put pl in bb_bins..
        	(a,b)=hist(pl,bb_bins);
        	peak_density=peak_density+b;
        	#all_b[bb,:]=b;
    	end
    end
    peak_density=squeeze(peak_density,2)/4;#as bin size=40kb, we want the unit as number per 10kb
    x=[-bb_lim+bin_size/2:bin_size:bb_lim-bin_size/2];
    return peak_density,x;
end


function form_big_mat(C_corrected_array,couple)

    T=size(C_corrected_array,3);
    N=size(C_corrected_array[:,:,1],1);
    tmp_big_mat=eye(N);
    for i=1:T-1
        tmp_big_mat=[tmp_big_mat eye(N)];
    end
    big_mat=tmp_big_mat;
    for i=1:T-1
        big_mat=vcat(big_mat,tmp_big_mat)
    end
    big_mat=big_mat*couple
    for t=1:T
        iz=[1:N]+(t-1)*N
        big_mat[iz,iz]=C_corrected_array[:,:,t]
    end

    return big_mat;
end


function multiplex_louvain(C_corrected_array,couple)

    big_mat=form_big_mat(C_corrected_array,couple)
    C_norm=sum(big_mat);
    (N1,N2,T)=size(C_corrected_array);
    (sigma, Q2, C_corrected_renorm) = iterate(big_mat,C_norm);
    louvain(C_corrected,C_norm)
    N_renorm=Array(Int,T,1)
    for j=1:T
        N_renorm[j]=length(unique(sigma[(1:N1)+(j-1)*N1]));
    end
    #Array(Float64,452,452,3);
    big_mat_renorm=zeros(sum(N_renorm),sum(N_renorm));
    usigma=Array(Int,T,1)

    for i=1:N_renorm
        for j=1:N_renorm
            C_corrected_renorm[i,j]=sum(sum(C_corrected[sigma.==usigma[i],sigma.==usigma[j]]));
        end
    end
end

function ice_corr(contact)

    id=find(sum(contact,1).>0);
    contact2=contact[id,id];

    T=contact2;
    B_init=sum(T,2);
    B_init=B_init/mean(B_init);
    B=B_init;
    err=100;
    while err>1e-5        
        B_old=B;
        tmp=broadcast(/,T,B)
        #tmp=bsxfun(@rdivide,T,B);
        T=broadcast(/,tmp',B)
        #T=bsxfun(@rdivide,tmp',B)';
        B=sum(T,2);B=B/mean(B);
        err=maximum(abs(B-B_old)./B)
        println(err)
    end

    T=(T+T')/2;
    s=full(sum(T,2));
    T=broadcast(/,T,s);
    T=(T+T')/2;
    Bias=diag(contact2)./diag(T);
    Bias=sqrt(Bias);

    return T, Bias
end

#this code is not necessary, we can use merely bins2moudules.
#a1, 0 means no TADs..
function TADs2track(TADs_list,chr_num,chr2bins);
	z=chr2bins[:,chr_num];
	N=z[2]-z[1]+1;
	a1=zeros(N);
	for i=1:size(TADs_list,1);
		st=TADs_list[:domain_st_bin][i];
		ed=TADs_list[:domain_ed_bin][i];
		st=st-z[1];
		ed=ed-z[1];
        a1[st:ed]=TADs_list[:idx][i];
    end
    return a1;
end

function blocks2track(st,ed,chr_num,chr2bins,bin_size);
	z=chr2bins[:,chr_num];
	N=z[2]-z[1]+1;
	a1=zeros(N);
	st=ceil(st/bin_size);
	ed=floor(ed/bin_size);
    for i=1:size(st,1);
        a1[st[i]:ed[i]]=i;
    end
    return a1;
end

#input: the track, wig file in an array
function bin_wig_track(track,chr2bins,bin2loc,chr_num)
	iz=find(bin2loc[1,:].==chr_num-1);
	track_bin=zeros(length(iz));
	for i=1:length(iz)
		st=bin2loc[2,iz[i]];
		ed=bin2loc[3,iz[i]];
		track_bin[i]=mean(track[st:ed]);
	end
	return track_bin;
end


function mix_track(a1);
	a1_perm=zeros(size(a1));
	e=[1:maximum(a1)];
	e=e-.5;
	e=vcat(e,e[end]+1);
	(u,v)=hist(a1,e);
	u=int(u-.5);
	u=u[2:end];
	u_perm=u[randperm(length(u))];
	st=1;
	for i=1:length(u)
		a1_perm[st:st+v[u_perm[i]]-1]=u_perm[i];
		st=st+v[u_perm[i]];
	end
	return a1_perm;
end

function MI_two_partitions(a1,a2);

    m1=maximum(a1);
    m2=maximum(a2);
    (u1,v1)=hist(a1,0:m1)
    (u2,v2)=hist(a2,0:m2)
    (u12,v12,c12)=hist2d([a1 a2],0:m1,0:m2)
    P1=v1/sum(v1)
    P2=v2/sum(v2)
    P12=c12/sum(c12)
    
    H1=P1.*log2(P1);
    H1[isnan(H1)]=0;
    H1=-sum(H1)

    H2=P2.*log2(P2);
    H2[isnan(H2)]=0;
    H2=-sum(H2)

    Z=P12./(P1*P2');
    S=P12.*log2(Z);
    S[isnan(S)]=0;

    MI=sum(S);
    MI_norm=2*MI/(H1+H2);

    return MI, MI_norm

end

#have dataframe all_TADs_list
function map_Genes_to_TADs(G2B_dict,all_TADs_list);
	#it's ideal if st and ed of a gene fit inside a TAD
	#st belongs to one TAD and ed belongs to another
	#st belongs to one TAD, and ed belongs to no TAD, may be a bit off... 
	#there are cases the bins belong to no, set as none.
	G2TAD_dict=Dict{ASCIIString,ASCIIString}();
	for key in keys(G2B_dict);
		bin2=G2B_dict[key];
		i_match=find((all_TADs_list[:domain_st_bin].<=bin2[1])&(all_TADs_list[:domain_ed_bin].>=bin2[end]));
		if length(i_match)==1
			TAD=all_TADs_list[i_match,:chr]*"_"*string(all_TADs_list[i_match,:idx][1]);
			G2TAD_dict[key]=TAD[1];
		else
			if length(bin2)>1
				i_match1=find((all_TADs_list[:domain_st_bin].<=bin2[1])&(all_TADs_list[:domain_ed_bin].>=bin2[1]));
				i_match2=find((all_TADs_list[:domain_st_bin].<=bin2[end])&(all_TADs_list[:domain_ed_bin].>=bin2[end]));
				if isempty(i_match1)&isempty(i_match2);
					G2TAD_dict[key]="none";
				elseif ~isempty(i_match1) & isempty(i_match2);
					TAD1=all_TADs_list[i_match1,:chr]*"_"*string(all_TADs_list[i_match1,:idx][1]);
					G2TAD_dict[key]=TAD[1];
				elseif isempty(i_match1) & ~isempty(i_match2);
					TAD2=all_TADs_list[i_match2,:chr]*"_"*string(all_TADs_list[i_match2,:idx][1]);
					G2TAD_dict[key]=TAD[1];
				elseif ~isempty(i_match1) & ~isempty(i_match2);
					TAD1=all_TADs_list[i_match1,:chr]*"_"*string(all_TADs_list[i_match1,:idx][1]);
					TAD2=all_TADs_list[i_match2,:chr]*"_"*string(all_TADs_list[i_match2,:idx][1]);
					G2TAD_dict[key]=TAD1[1]*","*TAD2[1]; 
				end
			else 
				G2TAD_dict[key]="none";
			end
		end
	end
	return G2TAD_dict;
end

function map_B2G_dict(G2B_dict)
	B2G_dict=Dict{Int,ASCIIString}();
	for key in keys(G2B_dict);
		bin2=G2B_dict[key];
		if length(bin2)==1;
			B2G_dict[bin2]=key;
		else
			for s=bin2[1]:bin2[end];
				B2G_dict[s]=key;
			end
		end
	end
	return B2G_dict;
end



#read a wiggle file, for a histone mark, a particular chr..
function read_wig(wig_file,chr2bins,chr_num);
	N=chr2bins[2,chr_num]-chr2bins[1,chr_num]+1;
	track=zeros(N*40000);
	f=open(wig_file);
	lines=readlines(f);
	for i=1:length(lines)
		if ~contains(lines[i],"variable")
			z=split(lines[i],r"\t|\n");
			track[int(z[1]):int(z[1])+25-1]=float(z[2]);
		end
	end
	#wiggle=readtable(wig_file,separator='\t');
	return track;
end

#given a TAD, st and end bins, return all genes within..
function TAD_to_Genes(B2G_dict,st_bin,ed_bin);
	Genes_list=ASCIIString[];
	for j=st_bin:ed_bin;
		if j in keys(B2G_dict)
			g=B2G_dict[j];
			push!(Genes_list,g);
		end
	end
	Genes_list=unique(Genes_list)
	return Genes_list;
end

function evec_distance(x,y);
	d1=sum((x-y).^2);
	d2=sum((x+y).^2);
	if d1<d2
		d=d1;
	else 
		d=d2;
	end
	return sqrt(d);
end

#in 20277 dim
function randomize_exp_within_TADs(expr_bin,all_TADs_list,B2G_dict,Gene);
	#something is wrong
	#check the effect of a gene to multiple TADs..
	#x=zeros(size(all_TADs_list,1));
	#y=zeros(size(all_TADs_list,1));
	#X=zeros(size(all_TADs_list,1));
	#expr_rand2=expr_bin+0;
	expr_rand1=zeros(size(expr_bin));
	for i=1:size(all_TADs_list,1);
		#display(i)
		st=all_TADs_list[i,4];
		ed=all_TADs_list[i,5];
		gg=TAD_to_Genes(B2G_dict,st,ed);
		if length(gg)>=1
			gg_id=zeros(size(gg));
			for k=1:length(gg);
				gg_id[k]=find(Gene.==gg[k])[1];
			end
			gg_id=int(gg_id);
			z=expr_bin[gg_id];
			z_rand=z[randperm(length(z))];				
			expr_rand1[gg_id]=z_rand+0;
	#		expr_rand2[gg_id]=z_rand;
			#y[i]=sum(expr_rand.>0)
		#X[i]=isequal(sum(expr_rand[gg_id].>0),sum(expr_bin[gg_id].>0));
	
		#if sum(expr_rand.>0)!=sum(expr_bin.>0)
		#	break;
		#end
		#x[i]=sum(expr_rand[gg_id]);
		#y[i]=sum(expr_bin[gg_id]);
		end
	end

	iz=find(expr_rand1.==0);
	expr_rand1[iz]=expr_bin[iz];

	return expr_rand1;

end

function read_Dixon_TADs(input_file,bin2loc,chr2bins,chr_num)
	all_Ren_TADs=readtable(input_file,header=false,separator='\t');
	#info=matread("/gpfs/scratch/fas/gerstein/ky26/Hi-C_processed/hg18_bin_info.mat");
	#chr2bins=int(info["chr2bins"]);
	#bin2loc=info["bin2loc"];
	chr="chr"*string(chr_num);
	i_ren=find(all_Ren_TADs[:1].==chr);
	Ren_TADs_chr=all_Ren_TADs[i_ren,:];
	tmp=chr2bins[2,chr_num-1];
	Ren_TADs_chr[:x4]=tmp+floor(Ren_TADs_chr[:x2]/40000);
	Ren_TADs_chr[:x5]=tmp+floor(Ren_TADs_chr[:x3]/40000);
	Ren_TADs_chr[:x6]=[1:size(Ren_TADs_chr,1)];
	rename!(Ren_TADs_chr,:x1,:chr)
	rename!(Ren_TADs_chr,:x4,:domain_st_bin)
	rename!(Ren_TADs_chr,:x5,:domain_ed_bin)
	rename!(Ren_TADs_chr,:x6,:idx);
	rename!(Ren_TADs_chr,:x2,:domain_st)
	rename!(Ren_TADs_chr,:x3,:domain_ed)
	return Ren_TADs_chr;
end





