# ZZpaper
Suppement code to Corbella, Spencer and Roberts, 2022, Automatic Zig-Zag sampling in practice 
This repository contains the code used to obtain results on the performance of the Automatic Zig-Zag sampler.

Most of the code can be run in Julia. 
Specifically, the file’s content is outlined below:

<div>
  <ul>
    <li>Functions.jl – all the functions used to perform the sampling given a potential U(x). The tasks performed comprise, among others: gradient evaluation; rate computation; rate optimization; sample drawings via Automatic Zig-Zag; super-efficient sample drawing via sub-sampling</li>
    <li>Plotting.jl – plotting routines for, among others: the tuning of t_max, the skeleton points, the comparison of efficiency and the comparison of robustness</li>
    <li>Main.jl – main file for the performance analysis</li>
    <li>ExampleDugongs.jl – main file for the analysis of section 5.1
    <li>BDA_preprocessing.jl – pre-processing of the Simulacrum dataset for the analyses in sections 5.2 and 6.4  </li>
    <li>BDA4par.jl – main file for the analysis of section 5.2 </li>
    <li>BDA4ALLDATA.jl – main file for the analysis of section 6.4  </li>
    <li>StanDugongs.jl, RstanDugongs.R, Dugongs.stan – file to perform the analysis of section 5.1 via stan.  </li>
  </ul>
</div>
