#import "../lib.typ": *

Here is one band unfolding method to transfer from the deeph-format Hamitonian data to a unfolded band structure. 

Here are two methods:

- Method1 is the general method, which is siutable for almost all cases. 

- Method2 supports for two-dimensional Moire honeycomb lattice now, which is my original purpose. Actually, for other cases, only minor modifications are needed, you can try it by yourself.

= Band unfolding method
Hamilton basis: 
  $ angle.l bold(r)|alpha, bold(k) angle.r = 1/sqrt(N) sum_bold(R) e^(i bold(k) dot bold(R)) ϕ(bold(r)-bold(R)-bold(tau)_alpha) $
  Bloch wfc: 
  $ 
  |n,bold(k) angle.r &= sum_α phi_(n alpha)(bold(k))|alpha,bold(k) angle.r \ 
  angle.l bold(r)|n,bold(k) angle.r &= 1/sqrt(N) sum_(alpha,bold(R)) e^(i bold(k) dot bold(R)) ϕ(bold(r)-bold(R)-bold(tau)_alpha) phi_(n alpha)(bold(k))
  $
  Define 
  $ ϕ(bold(p)) = integral d^3 bold(r) e^(-i bold(p) dot (bold(r)-bold(R)-bold(tau)_alpha))ϕ(bold(r)-bold(R)-bold(tau)_alpha) $ as the Fourier transform of $ϕ(bold(r))$, then
  $
  angle.l bold(p)|n,bold(k) angle.r &= 1/N sum_(alpha,R) integral d^3 bold(r) e^(-i bold(p) dot bold(r)) e^(i bold(k) dot bold(R)) ϕ(bold(r)-bold(R)-bold(tau)_alpha) phi_(n alpha)(bold(k)) \ 
  &= 1/N sum_(alpha,bold(R)) e^(i(bold(k)-bold(p))dot bold(R))e^(-i bold(p)dot bold(tau)_alpha) ϕ(bold(p)) phi_(n alpha)(bold(k)) \ 
  &= sum_(alpha,bold(G)) delta_(bold(p)_(| |),bold(k+G)) e^(-i bold(p)dot bold(tau)_alpha) ϕ(bold(p)) phi_(n alpha)(bold(k))
  $
  For ARPES, the intensity is given by

  $ I(bold(p),E) &∝ sum_(n,bold(k)) |angle.l bold(p)|bold(A)dot hat(bold(v))|n,bold(k) angle.r|^2 delta(E-epsilon_(n bold(k))) \ 
    &= sum_(n,bold(k))|bold(A)dot sum_(alpha,beta,bold(k'))angle.l bold(p)|alpha,bold(k')angle.r angle.l alpha,bold(k')|hat(bold(v))|beta,bold(k')angle.r angle.l beta,bold(k')|n,bold(k)angle.r|^2 delta(E-epsilon_(n bold(k)))\ 
    &= |ϕ(bold(p))|^2sum_(n,bold(k))|bold(A)dot sum_(alpha,beta,bold(G))  delta_(bold(p)_(| |),bold(k+G)) e^(-i bold(p)dot bold(tau)_alpha) angle.l alpha,bold(k)|partial_(bold(k))hat(H)|beta,bold(k)angle.r  phi_(n beta)(bold(k))|^2 delta(E-epsilon_(n bold(k)))\
    o r &= |ϕ(bold(p))|^2 sum_(n,bold(k),alpha ) |(bold(A)dot bold(p))^2 delta_(bold(p)_(| |),bold(k+G)) e^(-i bold(p)dot bold(tau)_alpha) phi_(n beta)(bold(k))|^2 delta(E-epsilon_(n bold(k)))
  $
  To simplify the problem, we regard the real spherical harmonic basis as the delta function.

  For unfolding, similarly, we have  

    $
    I(bold(p),E) &∝ sum_(g) sum_(N,bold(K)) |angle.l bold(p+g)|N,bold(K) angle.r|^2 delta(E-epsilon_(N bold(K))) \ 
    &= |ϕ(bold(p))|^2 sum_(g,N,bold(K),alpha,G) |delta_(bold(p+g),bold(K+G)) e^(-i bold((p+g))dot bold(tau)_alpha) phi_(N alpha)(bold(K))|^2 delta(E-epsilon_(N bold(K)))
  $

  where capital letters for moire cell and small letters for primitive cell.

= Unfolding Method1

Using eq.6, for a certain $p$ in space, we can find the corresponding $K$ and $G$, ignore the sumover for $G$, then we can calculate the intensity $I(p)$, with three loops of sumover: 1. $sum_(N K)$; 2. $sum_g$; 3. $sum_alpha$.


= Unfolding Method2

Back to eq.6, in fact, we project the Bloch wave function using a set of plane waves. From the perspective of band folding, the bands at det(M) G-points in a pBZ are all folded into an mBZ. Conversely, we only need to distinguish which G-point each band comes from, that is, to look at the weight with the largest projection component.

#figure(image("fig/3.png",width:50%))
*Step1*  For each $p$ in the KPATH, find $p = K + G_(0)$, where K in the (0,0)'s mBZ. 

*Step2*  For the $K$ given in the Step1, find all $k = K + G_(m),$ where $m=1,2,...,det(M)$. 

*Step3* For each $k$ given in the Step2, calculate weight $w_(m)$ and save 
maximum $max(w_m)$.

*Step4* Compare $max(w_m)$ with $w_(0)$ in Step1, if $w_(0)$ is maximum, then accept this unfolding action.




// #definition(
//   "Stokes' theorem",
//   footer: "Information extracted from a well-known public encyclopedia"
// )[
//   Let $Sigma$ be a smooth oriented surface in $RR^3$ with boundary $diff Sigma
//   equiv Gamma$. If a vector field $iboxed(bold(F)(x,y,z))=(F_x (x,y,z), F_y (x,y,z),
//   F_z (x,y,z))$ is defined and has continuous first order partial derivatives
//   in a region containing $Sigma$, then

//   $ integral.double_Sigma (bold(nabla) times bold(F)) dot bold(Sigma) =
//   dboxed(integral.cont_(diff Sigma) bold(F) dot dif bold(Gamma)) $ 
// ] <stokes>

