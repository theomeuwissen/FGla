# FGla
## Calculates IBD based relationship matrix Gla using linkage analysis

The pedigree relationship matrix A is based on loops in the pedigree and calculates the chance of receiving the same allele twice through the pedigree loops. For this it
assumes 50/50 probabilities for paternal/maternal inheritances of alleles. Gla does the same but uses marker data to determine paternal/maternal inheritances using linkage analysis.
For this FGla needs at least 3 generations of genotyped and pedigree data: 2 generation to determine which allele is paternal/maternal and an extra generation to see whether 
the maternal or paternal allele is transmitted. Since these inheritances differ per SNP position, the Gla matrix is averaged over all SNP loci provided. SNP loci may be provided for a region, or entire chromosome but not for several chromosomes (FGla will be confused by the large number of recombinations). The Gla matrices of several chromosome may be 
calculated in parallel, and averaged (weighted by the chromosome size or number of SNPs).            

For a detailed description and reference for the FGla method see Meuwissen, Yu and Berg (Gen. Sel. Evol. 2025)        



## Usage:
julia --threads 20 FGla.jl pedfile=\"example.ped\" plinkfilstem=\"example\" lstfile=\"example.lst22\" Flstfile=\"example.lst22\"     

Here the program uses 20 threads.     

pedfile is file with pedigree data. The pedigreed animals are numbered from 1,2,3,..,Nped, where animals are numbered from old till young. The pedigree file has 3 columns: 

'ID sireID damID, where for unknown sires/dams '0' is used as their ID number.    

plinkstemfile denotes the stem of the plink .bed, .bim and .fam files. Not all pedigreed animals need genotypes, but if an animal is in included in the .fam file, all
genotypes are assumed known (for now).    

lstfile (optional) contains a list of animal IDs for which the Gla relationship matrix is needed.    

Flstfile (optional) contains a list of animal IDs for which the inbreeding coefficients are needed.    

Note: the input file names are between " signs which need to be 'protected' by a \\, resulting in the \\" signs. Dont include additional spaces. 
If an additional space is included this will be considered the next commandline argument, i.e. filenames should not contain spaces.    

## Output files
FGla.F: contains the inbreeding coefficients in the order provided by the Flstfile.

gla.Float32.*.mmap (e.g. gla.Float32.10311.mmap): a .mmap file containing the Gla of the 10,311 IDs included in lstfile (in this order) that is to be read by julia:    

`    
using Mmap;    
io=open("gla.Float32.10311.mmap","r+");    
Gla=mmap(io,Array{Float32,2},(10311,10311));  #note 10311 = number of IDs in lstfile;   
`    
After these commands the 10311x10311 Gla matrix is available for use in Julia. Writing out the matrix in .txt form could be by:    
`    
using DelimitedFiles;   
writedlm("Gla.txt",Gla);  #assuming its smallish;    
#or;    
iout=open("Gla.txt","w"); #assuming its large;    
for i=1:10311, j=1:i;    
   println(iout,Gla[i,j]);    
end;   
`










