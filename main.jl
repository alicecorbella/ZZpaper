# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/"

# Loading the packages
include(string(CODEwd, "env.jl"))


Dim = 2

μ = zeros(Dim)
Σ = Diagonal(ones(Dim))+zeros( Dim,Dim)

# minus log density
function U(x::Vector, m=μ, s=Σ)
    -(-2/2*log(2π)-0.5*log(det(s))-
    0.5*(transpose(x-m)*inv(s)*(x-m))
     )
end

gradUx(x; Ux=U) = ForwardDiff.gradient(Ux, x)
