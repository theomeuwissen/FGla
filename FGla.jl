# using FGla algorithm:
# find first segregation indicators S
# next calculate Gla
#                                           
# usage: julia --threads 20 FGla.jl pedfile=\"lamip.ped\" plinkfilstem=\"plink\" lstfile=\"byr2022.A22.srt\" lstFfile=\"byr2022.A22.srt\"
# note:                                     ^^here \" means a literal " sign, which denotes the beginning and end of the filename 
#
# Theo Meuwissen April 2025
              
using Mmap,DelimitedFiles,LinearAlgebra,Statistics, Base.Threads
pedfile="FGla.ped"          #defaul pedigree file
plinkfilstem="plink"        #default  stem of plink files .bed, .bim, .fam 
lstfile="" #list of animals for which Gla relationships are needed (DEFAULT: no relationships are needed: lstfile="")
lstFfile=""               #list of animals for which inbreeding F is needed (DEFAULT: no F coefficients  needed: lstFfile="")
if(length(ARGS)>0)
  fargs=open("FGla.args","w")
  map(x->println(fargs,x),ARGS)
  close(fargs)
  include("FGla.args")  #read command line arguments
end




#code: reading data
#reading plink files
fbed=open(plinkfilstem*".bed")
plinkfam=readdlm(plinkfilstem*".fam",Int64)
plinkbim=readdlm(plinkfilstem*".bim")
Ngen=size(plinkfam,1)
Nsnp=size(plinkbim,1)
recomb=(plinkbim[end,4]-plinkbim[1,4])/Nsnp/1e8
gen_err=0.05
genotIDS=plinkfam[:,2]
          
nbytes=ceil(Int,Ngen/4)
a=zeros(UInt8,nbytes)
a[1:3]=read(fbed,3)
bit_mask=UInt8.([3; 12; 48; 192])
gen01=ones(Int32,Nsnp,2Ngen);  #genotypes coded "1 1", "1 2", "2 2"
for i=1:Nsnp
    global a=read(fbed,nbytes)
    ianim=0
    for j=1:nbytes
        for k=1:4
            ianim+=1
            if(ianim<=Ngen)
                igenot=count_ones(a[j]&bit_mask[k])
                if(igenot==1)
                    gen01[i,ianim*2]=2
                elseif (igenot==2)
                    gen01[i,ianim*2-1:ianim*2].=2
                end
                (ianim==1 && i<20) ? println(igenot) : nothing
             end
         end
     end
end

# read pedfile; open output file
ped=readdlm(pedfile,Int)
Nanim=size(ped,1)
iotxt="S.Int32$(Nsnp)x$(2Ngen).mmap"
println("S coefficients (1/2 & 0 for unknown) being written to ",iotxt)          
io=open(iotxt,"w+")  #mmap output
S=mmap(io,Array{Int32,2},(Nsnp,2Nanim))	


#code: processing data
	volgnr2genid=zeros(Int,Nanim);
	volgnr2genid[genotIDS]=collect(1:size(genotIDS,1))

              
# sort pat/mat inheritance of heterozygous genotypes
# code 1,2 and -1,-2 means known paternal/maternal inheritance
	for i=1:Ngen
	    (is,id)=ped[genotIDS[i],2:3];
            iss=0; idd=0;
	    (is>0) ? iss=volgnr2genid[is] : iss=0
	    (id>0) ? idd=volgnr2genid[id] : idd=0
            if(iss>0)    
	       for j=1:Nsnp
	           if (gen01[j,2i-1]!=gen01[j,2i]) && (gen01[j,2iss-1]==gen01[j,2iss])
                       if(idd>0)&&(gen01[j,2idd-1]==gen01[j,2idd]) && (gen01[j,2iss]==gen01[j,2idd])
                       else                           
	                   abs(gen01[j,2iss])!=abs(gen01[j,2i-1]) ? gen01[j,2i-1:2i]=-abs.(gen01[j,2i:-1:2i-1]) : gen01[j,2i-1:2i]=-abs.(gen01[j,2i-1:2i])
                       end
	          end
	        end
	    end
	    if(idd>0)
	      for j=1:Nsnp
	        if (gen01[j,2i-1]!=gen01[j,2i]) && (gen01[j,2idd-1]==gen01[j,2idd])
                       if(iss>0)&&(gen01[j,2iss-1]==gen01[j,2iss]) && (gen01[j,2iss]==gen01[j,2idd])
                       else                           
	                   abs(gen01[j,2idd])!=abs(gen01[j,2i]) ? gen01[j,2i-1:2i]=-abs.(gen01[j,2i:-1:2i-1]) : gen01[j,2i-1:2i]=-abs.(gen01[j,2i-1:2i])
                       end    
	        end
	      end
	    end
	end

# setup S
	       for i=1:size(genotIDS,1)
	           (is,id)=ped[genotIDS[i],2:3]; ii=genotIDS[i]
                   iss=0; idd=0
	           (is>0) ? iss=volgnr2genid[is] : iss=0
	           (id>0) ? idd=volgnr2genid[id] : idd=0
	           if(iss>0)
	             for j=1:Nsnp
	                if(gen01[j,2iss-1]!=gen01[j,2iss])&&(gen01[j,2iss]<0)
	                   if(gen01[j,2i-1]==gen01[j,2i])||(gen01[j,2i-1]<0)
	                     abs(gen01[j,2i-1])==abs(gen01[j,2iss-1]) ? S[j,2ii-1]=1 : S[j,2ii-1]=2
	                   end
	                end
	             end
	            end
	            if(idd>0)
	              for j=1:Nsnp
	                if(gen01[j,2idd-1]!=gen01[j,2idd])&&(gen01[j,2idd]<0)
	                  if(gen01[j,2i-1]==gen01[j,2i])||(gen01[j,2i]<0)
	                     abs(gen01[j,2i])==abs(gen01[j,2idd-1]) ? S[j,2ii]=1 : S[j,2ii]=2
	                  end
	                end
	              end
	            end
	        end

# interpolate S
include("viterbi.jl")
	for ii=1:Nanim*2
	    if(any(S[:,ii].>0))
#                println(ii," #Pos ",sum(S[:,ii].>0))
                S[:,ii]=viterbi(S[:,ii],recomb,gen_err)
                println(ii," #Pos#swaps ",sum(S[2:end,ii].!=S[1:end-1,ii]))
                sum(S[2:end,ii].!=S[1:end-1,ii])>10 ? S[:,ii].=0 : nothing
	   end
	end
println(" Calculation of S: done")


          
# calculate Gla          
if(length(lstfile)+length(lstFfile)>0)
if(length(lstfile)>0)          
  lst=readdlm(lstfile,Int) #list of animals for which relationships are needed
else
  lst=Int.([])
end
nlst=size(lst,1)
if(length(lstFfile)>0)          
  Flst=readdlm(lstFfile,Int) #list of animals for which relationships are needed
else
  Flst=Int.([])
end
nFlst=size(Flst,1)          
println("  No animals in pedigree           =  ",Nanim)
println("  Gla needed for number of animals =  ",nlst)
println("  F needed for number of animals   =  ",nFlst)

if(nlst>0)          
   iostxt="gla.Float32.$(nlst).mmap"
   println(" Gla matrix written to = ", iostxt)          
   ios=open(iostxt,"w+") #memorymapped matrix of relationship coefficients
   A22 = mmap(ios, Array{Float32,2}, (nlst,nlst))
end          
A22_= zeros(Float32,nlst,nlst)
F22_=zeros(Float32,nFlst)
          
# setup gametic_pedigree and list
n=Nanim*2 #no. of gametes
PED=zeros(Int32,n,2)
for i=1:Nanim
    if(ped[i,2]>0)
        PED[2i-1,1:2]=[2*ped[i,2]-1 2*ped[i,2]]
    end
    if(ped[i,3]>0)
        PED[2i,1:2]=[2*ped[i,3]-1 2*ped[i,3]]
    end
end
lstall=sort(unique(vcat(lst,Flst)))
lst22=zeros(Int,size(lstall,1)*2)
for i=1:size(lstall,1)
    lst22[2i-1]=lstall[i]*2-1
    lst22[2i]=lstall[i]*2
end

# set up founder-alleles
iox=open("foundall.Int32.mmap","w+")
foundall=mmap(iox, Array{Int32,2}, (n,Nsnp))
foundall.=0
nfoundall=0
for i=1:n
    if(PED[i,1]==0)
        global nfoundall+=1
        foundall[i,:].=nfoundall
    end
end

# sample founder alleles for gametes in lists and all positions
for j=1:Nsnp
    for i=nlst*2:-1:1
        i_=lst22[i]
        if(foundall[i_,j]==0)
            ancest_=[i_]  #ancestor list
            ii_=deepcopy(i_)
            ifoundall_=0
            while (ifoundall_==0)
              if(S[j,ii_]==0)
                Si_=rand(1:2) #sample at random ancestor
              else
                  Si_=S[j,ii_]
              end
              if(PED[ii_,Si_]>0)  
                ii_=PED[ii_,Si_]  #goto parent
                if(foundall[ii_,j]>0) #founder allele found
                    ifoundall_=foundall[ii_,j]
                    foundall[ancest_,j].=ifoundall_
                else
                  push!(ancest_,ii_)  #add parent to ancestor list
                end
              else
                  println("ERROR: founder-allele not found ",i_)
              end
            end  #while
        end #if
    end #i
end #j
println(" finished sampling of founder alleles")

# set up F and Gla relationships given founder alleles

if(nFlst>0)          
@threads for i=1:nFlst
    i_=[2*lst[i]-1, 2*lst[i]]
    F22_[i]=sum(foundall[i_[1],1:Nsnp].==foundall[i_[2],1:Nsnp])/Nsnp
end
writedlm("FGla.F",F22_)
println("  F coefficients written to FGla.F")                    
end #if


          
if(nlst>0)          
@threads for i=1:nlst
  for j=i:nlst
    i_=[2*lst[i]-1, 2*lst[i]]
    j_=[2*lst[j]-1, 2*lst[j]]
    if(i==j) #diagonal
        A22_[i,i]=1+sum(foundall[i_[1],1:Nsnp].==foundall[i_[2],1:Nsnp])/Nsnp
    else
        A22_[i,j]=sum(foundall[i_[1],1:Nsnp].==foundall[j_[1],1:Nsnp])+sum(foundall[i_[1],1:Nsnp].==foundall[j_[2],1:Nsnp])
        A22_[i,j]+=sum(foundall[i_[2],1:Nsnp].==foundall[j_[1],1:Nsnp])+sum(foundall[i_[2],1:Nsnp].==foundall[j_[2],1:Nsnp])
        A22_[i,j]/=(Nsnp*2)
        A22_[j,i]=A22_[i,j]
    end
  end
end
A22.=A22_
          close(ios)
println(" calculation of Gla: done")          
end #if

          
end #if              
exit()
