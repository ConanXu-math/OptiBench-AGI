# A relative-error inexact ADMM splitting algorithm for convex optimization with inertial effects

**arXiv ID:** 2409.10311v1

**Authors:** M. Marques Alves, M. Geremia

**Abstract:** We propose a new relative-error inexact version of the alternating direction method of multipliers (ADMM) for convex optimization. We prove the asymptotic convergence of our main algorithm as well as pointwise and ergodic iteration-complexities for residuals. We also justify the effectiveness of the proposed algorithm through some preliminary numerical experiments on regression problems.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:2409.10311v1  [math.OC]  16 Sep 2024
A relative-error inexact ADMM splitting algorithm for conv ex
optimization with inertial eﬀects
M. Marques Alves ∗ Marina Geremia †
Dedicated to the memory of Professor Hedy Attouch
Abstract
We propose a new relative-error inexact version of the alter nating direction method of mul-
tipliers (ADMM) for convex optimization. We prove the asymp totic convergence of our main
algorithm as well as pointwise and ergodic iteration-compl exities for residuals. We also justify
the eﬀectiveness of the proposed algorithm through some pre liminary numerical experiments on
regression problems.
2000 Mathematics Subject Classiﬁcation: 90C25, 90C06, 49J 52.
Key words: Convex optimization, ADMM, inexact, relative-e rror, inertial algorithms.
1 Introduction
Consider the convex problem
minimize
x∈H
f (x) + g(Lx), (1)
where f : H → (−∞ , ∞ ] and g : G → (−∞ , ∞ ] are lower semicontinuous proper convex functions and
L : H → G is a linear operator ( H and G denote ﬁnite-dimensional inner product spaces). Problem
(1) appears in diﬀerent contexts in applied mathematics, in cluding optimization, inverse problems,
machine learning, among others.
One of the most popular numerical algorithms for solving (1) is the alternating direction method
of multipliers (ADMM) [26, 28, 29], which has now attracted a lot of attention from the numerical
optimization community (see, e.g., [9, 11, 12, 13, 14, 16, 22 , 23, 30, 32, 33, 35, 37, 42, 46]).
In this paper we propose and study a new inexact version of the ADMM allowing relative-error
criteria for the solution of the second subproblem (which wi ll appear in the formulation of the proposed
algorithm) and promoting inertial eﬀects on iterations.
Organization of the paper. The material is organized as follows. In Section 2 we motivat e
the deﬁnition of our main algorithm, review some related wor ks and discuss the main contributions
∗Departamento de Matemática, Universidade Federal de Santa Catarina, Florianópolis, Brazil, 88040-900
(maicon.alves@ufsc.br). The work of this author was partially supported by CNPq gra nt 308036/2021-2.
†Departamento de Matemática, Universidade Federal de Santa Catarina, Florianópolis, Brazil, 88040-900. Depar-
tamento de Ensino, Pesquisa e Extensão, Instituto Federal d e Santa Catarina (IFSC) ( marina.geremia@ifsc.edu.br).
1

<!-- page 2 -->
of this paper. In Section 3, we present our main algorithm (Al gorithm 1) and some preliminary
results that will be needed to study its convergence and iter ation-complexity. In Section 4, we study
the asymptotic behavior of Algorithm 1 under diﬀerent assum ptions under the inertial parameters
involved in the formulation of the method. The main results a re Theorems 4.2 and 4.4. In Section 5
we study the iteration-complexities of Algorithm 1. The mai n results in this section are Theorems 5.2
and 5.4. Numerical experiments will be presented in Section 6. Appendix A contains some auxiliary
results.
2 Motivation, related works and contributions
Motivation. We ﬁrst note that (1) is clearly equivalent to the separable p roblem
minimize f (x) + g(y),
subject to Lx − y = 0. (2)
An iteration of the standard ADMM [22] for solving (2) can be d escribed as follows: given a starting
point (y0, z 0) ∈ G 2 and a regularization parameter γ > 0, iterate for k ≥ 0:
xk+1 ∈ argmin
x∈H
{
f (x) + ⟨zk |Lx − yk⟩ + γ
2 ∥Lx − yk∥2
}
, (3)
yk+1 ∈ argmin
y∈G
{
g(y) + ⟨zk |Lxk+1 − y⟩ + γ
2 ∥Lxk+1 − y∥2
}
, (4)
zk+1 = zk + γ (Lxk+1 − yk+1) . (5)
We consider here the case in which (3) can be solved exactly an d, on the other hand, (4) is supposed
to be solved only approximately by some other (inner) algori thm, like, for instance, CG or BFGS,
depending on the particular structure of the function g(·) in (2).
With this in mind, we will introduce a notion of relative-err or approximate solution for (4) (more
details will be given on Section 3). To this end, ﬁrst note tha t (4) is an instance of the general family
of minimization problems
minimize
y∈G
{
g(y) + ⟨z |Lx − y⟩ + γ
2 ∥Lx − y∥2
}
, (6)
where x ∈ H , z ∈ G and γ > 0 are given (in the case of (4), we have (x, z ) = ( xk+1, z k)). Moreover,
since the function g(·) is convex, we have that (6) is also equivalent to the inclusio n/equation system
for the pair (y, v ):



v ∈ ∂g (y),
v − z + γ(y − Lx) = 0 .
(7)
A formal deﬁnition of approximate (inexact) solution of (7) (or, equivalently, (6)) will be given
in Deﬁnition 3.1 in Section 3; such a notion of approximate so lution will allow for errors in both the
inclusion and the equation in (7).
2

<!-- page 3 -->
The extended-solution set. The Fenchel dual of (1) is
maximize
z∈G
− f ∗(− L∗z) − g∗(z), (8)
where f ∗ : H → (−∞ , ∞ ] and g∗ : G → (−∞ , ∞ ] denote the Fenchel conjugates of f and g, respec-
tively, and L∗ : G → H denotes the adjoint operator of L. Under standard regularity conditions [10] on
f, g and L it is well-known that (1) and (8) are, respectively, equival ent to the (monotone) inclusions
0 ∈ ∂f (x) + L∗∂g (Lx), (9)
and
0 ∈ − L∂f ∗(− L∗z) + ∂g ∗(z). (10)
We make the blanket assumption:
Assumption 2.1. For the function f and the operator L as in (1), the following holds:
∂(f ∗ ◦ − L∗) = − L ◦ ∂f ∗ ◦ − L∗.
Several suﬃcient conditions for Assumption 2.1 to hold true can be found, e.g., in [10]. We will
also consider an extended-solution set S, attached to the pair of inclusions (9)–(10), deﬁned as
S = {(z, w ) ∈ G 2 | − w ∈ ∂(f ∗ ◦ − L∗)(z) and w ∈ ∂g ∗(z)}. (11)
Under Assumption 2.1, it is easy to check that if (z, w ) ∈ S , then it follows that there exists
x ∈ H such that x ∈ ∂f ∗(− L∗z), w = Lx and x and z are solutions of (9) and (10), respectively.
Throughout this work we will assume the following.
Assumption 2.2. We assume the extended solution set S as in (11) is nonempty.
Inertial algorithms. Iterative algorithms with inertial eﬀects for monotone inc lusions (and related
topics in optimization, saddle-point, equilibrium proble ms, etc) were ﬁrst proposed in the seminal
paper [2] and subsequently developed in various directions of research by diﬀerent authors and research
groups (see, e.g., [3, 5, 6, 7, 8, 12, 17] and references there in). Basically, the main idea consists in at
a current iterate, say pk, produce an “inertial eﬀect” by a simple extrapolation:
ˆpk = pk + α k(pk − pk− 1),
where α k ≥ 0, and then generate the next iterate pk+1 from ˆpk instead of pk (see (18)–(19) below).
Our main algorithm, namely Algorithm 1, will beneﬁt from ine rtial eﬀects on the iteration; see the
comments and remarks following Algorithm 1 for more discuss ions regarding the eﬀects of inertia.
Main contributions. We present a theoretical (asymptotic and iteration-comple xity analysis) and
computational study of a partially inexact ADMM splitting a lgorithm for solving (1). Our main
algorithm, namely Algorithm 1 below, beneﬁts from the addit ion of inertial eﬀects; see (18) and
(19). The convergence analysis is presented in Theorem 4.2, to which the proof incorporates some
elements of [3] and [46]. We also obtained iteration-comple xities for the proposed algorithm by
showing pointwise O(1/
√
k) and ergodic O(1/k ) global convergence rates for residuals; see Theorems
3

<!-- page 4 -->
5.2 and 5.4 below. We justify the eﬀectiveness of our main alg orithm through the realization of
numerical experiments on the LASSO problem (see Section 6).
Related works. A partially inexact ADMM splitting algorithm was recently p roposed and studied
in [46]. Paper [1] proposes a partially inexact ADMM for whic h the ﬁrst subproblem is supposed
to be solved inexactly. The analysis of the main algorithm in [1] is performed by viewing it as a
special instance of a non-Euclidean version of the hybrid pr oximal extragradient method [36]. In
contrast to this, analogously to [46], our main algorithm (A lgorithm 1 below) assumes the second
subproblem is solved inexactly. Moreover, since (22)–(24) below also allows for errors in ∂g , the
error criterion we propose here is potentially more ﬂexible than the corresponding one in [1]. Other
relative-error inexact versions of ADMM were also previous ly studied in [49, 50], but we notice that
the convergence results were restricted to the analysis of t he dual sequences. We also mention that
the relative-error inexact variants of the ADMM from [3, 23] only apply to (1) in the particular case
of L = I, and, additionally, these variants assume the ﬁrst subprob lem to be solved inexactly with the
error condition veriﬁed only a-posteriori, that is, only af ter the computation of second subproblem’s
solution.
General notation. We denote by H and G ﬁnite-dimensional inner product spaces with inner
product and induced norm denoted, respectively, by ⟨ | ⟩and ∥·∥=
√
⟨· | ·⟩. For any set X we denote
by X n the n-product X × · · · × X . In G2, we will consider the inner product and induced norm
deﬁned, respectively, by
⟨p |p′⟩γ := 1
γ ⟨z |z′⟩ + γ⟨w |w′⟩ and ∥p∥2
γ := ⟨p |p⟩γ , (12)
where p = (z, w ), p ′ = (z′, w ′) ∈ G 2 and γ > 0. More precisely, for p = (z, w ) ∈ G 2, the norm of p is
∥p∥2
γ = 1
γ ∥z∥2 + γ∥w∥2. (13)
An extended-real valued function f : H → (−∞ , ∞ ] is said to be convex whenever f (λx + (1 −
λ)y) ≤ λf (x) + (1 − λ)f (y) for all x, y ∈ H and λ ∈ (0, 1), and f is proper if its eﬀective domain ,
denoted by dom f , is nonempty. The Fenchel conjugate of a proper function f : H → (−∞ , ∞ ] is
f ∗ : H → (−∞ , ∞ ], deﬁned at any u ∈ H by f (u) = sup x∈H {⟨x |u⟩ − f (x)}. The ε-subdiﬀerential
and the subdiﬀerential of a convex function g : H → (−∞ , ∞ ] at x ∈ H are deﬁned as ∂εg(x) := {u ∈
H | g(y) ≥ g(x) + ⟨u |y − x⟩ − ε ∀y ∈ H} and ∂g (x) := ∂0g(x), respectively. For additional details
on standard notations and deﬁnitions of convex analysis we r efer the reader to the references [10, 41].
3 The main algorithm and some preliminary results
Consider the minimization problem (1), i.e.,
minimize
x∈H
f (x) + g(Lx), (14)
where f : H → (−∞ , ∞ ] and g : G → (−∞ , ∞ ] are lower semicontinuous proper convex functions
and L : H → G is a linear operator between ﬁnite-dimensional inner produ ct spaces H and G.
4

<!-- page 5 -->
In this section we present our main algorithm, namely Algori thm 1 below. This is a partially
inexact (the second block is allowed to be solved inexactly) ADMM with relative-error criterion for
the second subproblem. Recall the extended solution set S as in (11) and Assumptions 2.1 and 2.2.
The three technical lemmas 3.2, 3.3, 3.4 and 3.5 will be used i n the subsequent section.
Before presenting our main algorithm, as we discussed in the I ntroduction, we have to formalize
the notion of inexact solution that will be used to compute ap proximate solution for the second sub-
problem. Recall that the second subproblem of the standard A DMM (see (4)) belongs to the general
family of minimization problems (6), which is, in particula r, equivalent to the inclusion/equation
system (7) for the pair (y, v ), i.e.,



v ∈ ∂g (y),
v − z + γ(y − Lx) = 0 .
(15)
Deﬁnition 3.1 (σ -approximate solution of (6)) . For x ∈ H , (ˆz, ˆy) ∈ G 2 and γ > 0, a triple (˜y, v, ε ) ∈
G×G× R+ is said to be a σ -approximate solution of (6) (or, equivalently, of (7)) at (x, ˆz, ˆy) if σ ∈ [0, 1)
and









v ∈ ∂εg(˜y),
v − ˆz + γ(˜y − Lx) =: e
∥e∥2 + 2γε ≤ σ 2 min
{
γ2∥Lx − ˆy∥2, ∥v − ˆz∥2}
.
(16)
We will also write
˜y
σ
≈ argmin
y∈G
{
g(y) + ⟨ˆz |Lx − y⟩ + γ
2 ∥Lx − y∥2
}
meaning that there exists (v, ε ) such that (˜y, v, ε ) satisﬁes (16).
We now make some remarks regarding Deﬁnition 3.1:
(i) Note that if σ = 0 in (16), then it follows that e = 0 and ε = 0 , which is to say that the pair
(˜y, v ) satisﬁes the inclusion/equation system (7) (recall that ∂0g = ∂g ) and, in particular, ˜y is
an exact solution of (6).
(ii) The error criterion for (7) as in (16) belongs to the clas s of relative-error criteria for proximal-
type algorithms. Diﬀerent variants of such error condition s have been employed for computing
approximate solution for (sub) problems for a wide range of a lgorithms in monotone inclusions,
convex optimization, saddle-point problems, etc (see, e.g ., [3, 23, 36, 43, 44, 45]).
(iii) The error criterion (16) will be used to compute approx imate solutions in step 3 of our main
algorithm, namely Algorithm 1 below (see (22)–(24)).
(iv) As an illustrative example, consider the special case o f the LASSO problem [47]
min
x∈ Rd
{ 1
2 ∥Ax − b∥2 + ν∥x∥1
}
, (17)
5

<!-- page 6 -->
where A ∈ Rn× d, b ∈ Rn and ν > 0. Problem (17) is clearly a special instance of (1) in which
L := I, f (x) := ν∥x∥1 and g(x) := (1 / 2)∥Ax − b∥2 (see also Section 6 below). In this case, our
inclusion/equation system (7) clearly reduces to
v = A∗(Ay − b), v − z + γ(y − x) = 0 ,
or, in other words, in this special case, (7) is equivalent to the linear system (operator equation)
(
A∗A + γI
)
y = A∗b + z + γx.
The latter linear system can be solved by the CG algorithm [39 ], where e as in (16) will simply
denote the residual of the system and the inequality in (16) c an be used as a stopping criterion
for CG.
Next is our main algorithm.
Algorithm 1. An inexact inertial ADMM algorithm for solving (1)
(0) Let (z0, y 0) = ( z− 1, y − 1) ∈ G 2 and let α, σ ∈ [0, 1), τ ∈ (0, 1) and γ > 0. Set k = 0.
(1) Choose α k ∈ [0, α ] and let
ˆzk = zk + α k(zk − zk− 1), (18)
ˆyk = yk + α k(yk − yk− 1). (19)
(2) Compute
xk ∈ argmin
x∈H
{
f (x) + ⟨ˆzk |Lx − ˆyk⟩ + γ
2 ∥Lx − ˆyk∥2
}
. (20)
(3) Compute
˜yk
σ
≈ argmin
y∈G
{
g(y) + ⟨ˆzk |Lxk − y⟩ + γ
2 ∥Lxk − y∥2
}
(21)
at (xk, ˆzk, ˆyk) in the sense of Deﬁnition 3.1, i.e., compute (˜yk, v k, ε k) ∈ G × G × R+ such that
vk ∈ ∂εk g(˜yk), (22)
vk − ˆzk + γ(˜yk − Lxk) =: ek, (23)
∥ek∥2 + 2γεk ≤ σ 2 min
{
γ2∥Lxk − ˆyk∥2, ∥vk − ˆzk∥2}
. (24)
(4) Set
zk+1 = ˆzk + τ γ(Lxk − ˜yk), (25)
yk+1 = (1 − τ)ˆyk + τ
γ (ˆzk + γLx k − vk), (26)
k = k + 1 and go to step 1.
6

<!-- page 7 -->
We now make some remarks concerning Algorithm 1:
(i) Algorithm 1 is specially designed for instances of (1) in which (20) has a closed-form solution,
i.e., for problems in which (20) is easy to solve. In this rega rd, one example of interest is
when f (·) = ∥·∥1 and L = I, in which case (20) has a unique solution given explicitly by
xk = proxγ −1∥·∥1(ˆyk − γ− 1ˆzk). On the other hand, we assume that the computation of ˜yk as in
(21) demands the use of an (inner) algorithm, which the choic e of depends on the particular
structure of the function g, and, in this case, one can use (22)–(24) as a stopping criter ion for
the inner algorithm of choice.
(ii) Recall that we discussed in the Introduction (see “Rela ted works”) other ADMM-type algorithms
related to Algorithm 1.
(iii) The main results on the convergence and iteration-com plexity of Algorithm 1 are Theorems 4.2,
4.4, 5.2 and 5.4 below. Numerical experiments will be presen ted and discussed in Section 6.
(iv) The role of the parameter 0 < τ < 1 is to introduce (under) relaxation in the iterative process ;
see (25) and (26).
(v) We will also need the sequences (˘zk) and (˘yk), where, for all k ≥ 0,
˘zk := ˆzk + γ(Lxk − ˜yk), ˘yk := 1
γ (ˆzk + γLx k − vk). (27)
Note that ˘zk = zk+1 and ˘yk = yk+1 if we set τ = 1 in (25) and (26), respectively.
Next we present four technical lemmas – Lemmas 3.2, 3.3, 3.4 a nd 3.5 –, which will be useful in the
subsequent sections.
Lemma 3.2. Consider the sequences evolved by Algorithm 1 , let S be as in (11) and let (˘zk) and
(˘yk) be as in (27). Deﬁne
pk = (zk, y k), ˆpk = ( ˆzk, ˆyk), ˘pk = (˘zk, ˘yk) and ˜pk = (vk, Lx k) ∀k ≥ 0. (28)
(a) For all k ≥ 0,
pk+1 = (1 − τ)ˆpk + τ ˘pk.
(b) For all k ≥ 0,
− Lxk ∈ ∂ (f ∗ ◦ − L∗) (z′
k),
where
z′
k := ˆzk + γ (Lxk − ˆyk) . (29)
(c) For all k ≥ 0,
˘pk − ˆpk =
(
γ(Lxk − ˜yk), 1
γ
(
z′
k − vk
) )
.
7

<!-- page 8 -->
(d) For all k ≥ 0 and p = (z, w ) ∈ S ,
⟨˘pk − ˆpk |p − ˜pk⟩γ ≥ − εk.
(e) For all k ≥ 0,
ˆpk = pk + α k(pk − pk− 1).
Proof. (a) This result is a direct consequence of the deﬁnitions of pk, ˆpk and ˘pk as in (28) – see also
(27) – combined with (25) and (26).
(b) First note that from (20) and (29), we obtain 0 ∈ ∂f (xk) + L∗z′
k, or, equivalently, − L∗z′
k ∈
∂f (xk). As (∂f )− 1 = ∂f ∗, the latter inclusion is also equivalent to xk ∈ ∂f ∗(− L∗z′
k), which in turn
yields − Lxk ∈ − L∂f ∗(− L∗z′
k), which by Assumption 2.1 gives item (b).
(c) This follows easily from (27) – (29) and some simple algeb raic manipulations.
(d) As (∂εk g)− 1 = ∂εk g∗ – see, e.g., [51, p. 85, Theorem 2.4.4(iv)] –, we have that (22 ) is equivalent
to the inclusion
˜yk ∈ ∂εk g∗(vk). (30)
As p = ( z, w ) ∈ S , according to the deﬁnition of S in (11), we have − w ∈ ∂ (f ∗ ◦ − L∗) (z)
and w ∈ ∂g ∗(z). The latter inclusions combined with item (b) and (30) and th e monotonicity of
∂(f ∗ ◦ − L∗) and ∂g ∗ yield
⟨z′
k − z |w − Lxk⟩ ≥ 0 and ⟨z − vk |w − ˜yk⟩ ≥ − εk. (31)
Now using (28), item (c), (31) and the deﬁnition of ⟨· | ·⟩γ as in (12) we ﬁnd
⟨˘pk − ˆpk |p − ˜pk⟩γ = 1
γ ⟨γ(Lxk − ˜yk) |z − vk⟩ + γ
⟨ 1
γ
(
z′
k − vk
)
|w − Lxk
⟩
= ⟨Lxk − ˜yk |z − vk⟩ + ⟨z′
k − z |w − Lxk⟩ + ⟨z − vk |w − Lxk⟩
= ⟨z − vk |w − ˜yk⟩ + ⟨z′
k − z |w − Lxk⟩
≥ − εk.
(e) This follows directly from (18), (19) and the deﬁnition o f ˆpk as in (28).
Lemma 3.3. Consider the sequences evolved by Algorithm 1 and let (˘pk), (ˆpk) and (˜pk) be as in
(28). For all k ≥ 0,
(a) ∥˜pk − ˘pk∥2
γ = 1
γ
(
∥ek∥2 + ∥vk − ˆzk∥2
)
.
(b) ∥˘pk − ˆpk∥γ ≤ 2∥˜pk − ˆpk∥γ .
8

<!-- page 9 -->
Proof. (a) Direct use of (13), (23), (27) and (28) gives
∥˜pk − ˘pk∥2
γ = 1
γ ∥vk − ˘zk∥2 + γ∥Lxk − ˘yk∥2
= 1
γ ∥vk − ˆzk + γ(˜yk − Lxk)  
ek
∥2 + γ∥Lxk −
[
γ− 1(ˆzk + γLx k − vk)
]
∥2
= 1
γ
(
∥ek∥2 + ∥vk − ˆzk∥2
)
.
(b) In view of (24), (28) and item (a),
∥˜pk − ˘pk∥2
γ = 1
γ
(
∥ek∥2 + ∥vk − ˆzk∥2)
≤ 1
γ
(
σ 2∥γ(Lxk − ˆyk)∥2 + ∥vk − ˆzk∥2)
≤ 1
γ ∥vk − ˆzk∥2 + γ∥Lxk − ˆyk∥2
= ∥˜pk − ˆpk∥2
γ .
Hence, using the triangle inequality,
∥˘pk − ˆpk∥γ ≤ ∥ ˜pk − ˘pk∥γ + ∥˜pk − ˆpk∥γ ≤ 2∥˜pk − ˆpk∥γ .
Lemma 3.4. Consider the sequences evolved by Algorithm 1 and let (ˆpk), (˘pk) and (˜pk) be as in
(28).
(a) For all k ≥ 0 and p ∈ G 2,
∥p − ˆpk∥2
γ − ∥ p − ˘pk∥2
γ = ∥˜pk − ˆpk∥2
γ − ∥ ˜pk − ˘pk∥2
γ + 2⟨˘pk − ˆpk |p − ˜pk⟩γ .
(b) For all k ≥ 0 and p ∈ G 2,
∥p − ˆpk∥2
γ − ∥ p − ˘pk∥2
γ ≥ γ(1 − σ 2)∥Lxk − ˆyk∥2 + 2 [εk + ⟨˘pk − ˆpk |p − ˜pk⟩γ ] .
(c) For all k ≥ 0 and p = (z, w ) ∈ S ,
∥p − ˆpk∥2
γ − ∥ p − ˘pk∥2
γ ≥ γ(1 − σ 2)∥Lxk − ˆyk∥2.
(d) For all k ≥ 0 and p = (z, w ) ∈ S ,
∥p − ˆpk∥2
γ − ∥ p − pk+1∥2
γ ≥ τ γ(1 − σ 2)∥Lxk − ˆyk∥2 + (1 − τ)τ∥˘pk − ˆpk∥2
γ .
9

<!-- page 10 -->
(e) For all k ≥ 0 and p = (z, w ) ∈ S ,
∥p − ˆpk∥2
γ − ∥ p − pk+1∥2
γ ≥ τ(1 − τ)(1 − σ )2∥˜pk − ˆpk∥2
γ
≥ (1 − τ)(1 − σ )2
4τ ∥pk+1 − ˆpk∥2
γ .
Proof. (a) The desired result follows directly from the well-known identity ∥a − b∥2
γ − ∥ a − c∥2
γ =
∥d − b∥2
γ − ∥ d − c∥2
γ + 2⟨c − b |a − d⟩γ with a = p, b = ˆpk, c = ˘pk and d = ˜pk.
(b) Using (13) and the deﬁnitions of (˜pk) and (ˆpk) as in (28) we get
∥˜pk − ˆpk∥2
γ = 1
γ ∥vk − ˆzk∥2 + γ∥Lxk − ˆyk∥2,
which in turn combined with Lemma 3.3(a) yields
∥˜pk − ˆpk∥2
γ − ∥ ˜pk − ˘pk∥2
γ = γ∥Lxk − ˆyk∥2 − 1
γ ∥ek∥2.
From (24), item (a), the latter identity and some algebraic m anipulations,
∥p − ˆpk∥2
γ − ∥ p − ˘pk∥2
γ = γ∥Lxk − ˆyk∥2 − 1
γ ∥ek∥2 + 2⟨˘pk − ˆpk |p − ˜pk⟩γ
= γ∥Lxk − ˆyk∥2 − 1
γ
(
∥ek∥2 + 2γεk
)
+ 2 [εk + ⟨˘pk − ˆpk |p − ˜pk⟩γ ]
≥ γ(1 − σ 2)∥Lxk − ˆyk∥2 + 2 [εk + ⟨˘pk − ˆpk |p − ˜pk⟩γ ] ,
which ﬁnishes the proof of (b).
(c) This is a direct consequence of Lemma 3.2(d) and item (b) a bove.
(d) Using Lemma 3.2(a) and the identity ∥(1 − τ)a + τ b∥2
γ = (1 − τ)∥a∥2
γ +τ∥b∥2
γ − (1− τ)τ∥a − b∥2
γ
with a = p − ˆpk and b = p − ˘pk, we obtain
∥p − pk+1∥2
γ = ∥(1 − τ)(p − ˆpk) + τ(p − ˘pk)∥2
γ
= (1 − τ)∥p − ˆpk∥2
γ + τ∥p − ˘pk∥2
γ − (1 − τ)τ∥˘pk − ˆpk∥2
γ .
Now by multiplying the inequality in item (c) by τ > 0, using the latter identity and some simple
algebraic manipulations, we ﬁnd the desired result.
(e) Note ﬁrst that using (24) and the triangle inequality, we ﬁnd
∥vk − ˆzk∥ ≤ ∥ vk − ˆzk + γ(˜yk − Lxk)
 
ek
∥ + ∥γ(˜yk − Lxk)∥
≤ σ ∥vk − ˆzk∥ + ∥γ(˜yk − Lxk)∥,
so that
∥γ(˜yk − Lxk)∥ ≥ (1 − σ )∥vk − ˆzk∥,
10

<!-- page 11 -->
which in turn combined with (13), (27) and (28) yields
∥˘pk − ˆpk∥2
γ = 1
γ ∥˘zk − ˆzk∥2 + γ∥˘yk − ˆyk∥2
≥ 1
γ ∥˘zk − ˆzk∥2
= 1
γ ∥γ(˜yk − Lxk)∥2
≥ 1
γ (1 − σ )2∥vk − ˆzk∥2. (32)
Now using (12), (28), item (d) above and (32),
∥p − ˆpk∥2
γ − ∥ p − pk+1∥2
γ ≥ τ γ(1 − σ 2)∥Lxk − ˆyk∥2 + (1 − τ)τ∥˘pk − ˆpk∥2
γ
≥ τ γ(1 − σ 2)∥Lxk − ˆyk∥2 + (1 − τ)τ 1
γ (1 − σ )2∥vk − ˆzk∥2
≥ τ(1 − τ)(1 − σ )2
[ 1
γ ∥vk − ˆzk∥2 + γ∥Lxk − ˆyk∥2
]
= τ(1 − τ)(1 − σ )2∥˜pk − ˆpk∥2
γ . (33)
To prove the second inequality, one can use Lemmas 3.2(a) and 3.3(b) to conclude that ∥pk+1 − ˆpk∥γ ≤
2τ∥˜pk − ˆpk∥γ and then apply it in (33).
Lemma 3.5. Consider the sequences evolved by Algorithm 1 and let (pk) and (ˆpk) be as in (28).
Then, for all k ≥ 0,
∥ˆpk − p∥2
γ = (1 + α k)∥pk − p∥2
γ − α k∥pk− 1 − p∥2
γ + α k(1 + α k)∥pk − pk− 1∥2
γ ∀p ∈ G 2.
Proof. Recall that from Lemma 3.2(e) we have
ˆpk = pk + α k(pk − pk− 1), (34)
which is clearly equivalent to
pk − p = 1
1 + α k
(ˆpk − p) + α k
1 + α k
(pk− 1 − p) .
Now using the well-known identity ∥tx + (1 − t)y∥2
γ = t∥x∥2
γ + (1 − t)∥y∥2
γ − t(1 − t)∥x − y∥2
γ with
t = 1/ (1 + α k), x = ˆpk − p and y = pk− 1 − p, we ﬁnd
∥pk − p∥2
γ = 1
1 + α k
∥ˆpk − p∥2
γ + α k
1 + α k
∥pk− 1 − p∥2
γ − α k
(1 + α k)2 ∥ˆpk − pk− 1∥2
γ ,
which, in turn, when combined with the fact that ˆpk − pk− 1 = (1 + α k)(pk − pk− 1) – see (34) – and
after some simple algebraic manipulations it yields
∥ˆpk − p∥2
γ = (1 + α k)∥pk − p∥2
γ − α k∥pk− 1 − p∥2
γ + α k(1 + α k)∥pk − pk− 1∥2
γ .
11

<!-- page 12 -->
4 Asymptotic convergence of Algorithm 1
In this section, we study the asymptotic convergence of Algo rithm 1. The main results are Theorems
4.2 and 4.4.
Lemma 4.1. Consider the sequences evolved by Algorithm 1 and, for an arbitrary p = ( z, w ) ∈ S ,
deﬁne
hk = ∥pk − p∥2
γ ∀k ≥ − 1. (35)
Then h0 = h− 1 and, for all k ≥ 0,
hk+1 − hk − α k(hk − hk− 1) + τ(1 − τ)(1 − σ )2∥˜pk − ˆpk∥2
γ ≤ α k(1 + α k)∥pk − pk− 1∥2
γ ,
i.e., (hk) satisﬁes the assumptions of Lemma A.1 below, where, for all k ≥ 0,
sk+1 := τ(1 − τ)(1 − σ )2∥˜pk − ˆpk∥2
γ , (36)
δk := α k(1 + α k)∥pk − pk− 1∥2
γ . (37)
Proof. The fact that h0 = h− 1 follows directly from the fact that p0 = p− 1 (see step 0 in Algorithm
1 and the deﬁnition of pk as in (28)). On the other hand, from Lemma 3.5 and the deﬁnitio n of hk
as in (35),
∥ˆpk − p∥2
γ = (1 + α k) ∥pk − p∥2
γ
 
hk
− α k ∥pk− 1 − p∥2
γ
 
hk−1
+α k(1 + α k)∥pk − pk− 1∥2
γ .
The desired result now follows from the above displayed equa tion, Lemma 3.4(e) and the deﬁnition
of hk as in (35).
Next we present our ﬁrst result on the convergence of Algorit hm 1 when k → +∞ .
Theorem 4.2 (First result on the asymptotic convergence of Algorithm 1) . Consider the sequences
evolved by Algorithm 1 and let ∅ ̸= S be as in (11). Assume that
∞∑
k=0
α k∥pk − pk− 1∥2
γ < ∞ , (38)
where (pk) is as in (28). Then there exists (z∞ , w ∞ ) ∈ S such that
zk → z∞ and yk → w∞ . (39)
Additionally, we also have
vk → z∞ , Lx k → w∞ and ˜yk → w∞ . (40)
12

<!-- page 13 -->
Proof. We start by making a few remarks. First, from (38) and the fact that α k(1 + α k) ≤ 2α k
(because 0 ≤ α k < 1), we conclude that ∑∞
k=0 δk < ∞ , where δk is as in (37), which, in turn,
combined with Lemmas 4.1 and A.1 (below) gives
lim
k→∞
hk exists and
∞∑
k=1
sk < ∞ , (41)
where hk and sk+1 are as in (35) and (36), respectively. Using (36) and the seco nd statement in (41)
we also obtain ∥˜pk − ˆpk∥2
γ → 0, which in turn when combined with the deﬁnitions of ˜pk and ˆpk – as
in (28) –, (13), (23) and (24) yields
vk − ˆzk → 0, Lx k − ˆyk → 0, ˜yk − Lxk → 0 and εk → 0. (42)
Second, from (38) and the fact that α 2
k ≤ α k, we obtain
lim
k→∞
α k∥zk − zk− 1∥ = lim
k→∞
α k∥yk − yk− 1∥ = 0,
which, in turn, when combined with the deﬁnitions of ˆzk and ˆyk as in (18) and (19) yields
ˆzk − zk → 0 and ˆyk − yk → 0. (43)
Now, let pk = (zk, y k) be as in (28). Note that using the ﬁrst statement in (41), the d eﬁnition of
hk as in (35) and Lemma A.2 below, it follows that to prove the con vergence of (pk) to some element
in S – and hence the statement in (39) – it suﬃces to show that every cluster point of (pk) belongs
to S. To this end, let p∞ = ( z∞ , y ∞ ) ∈ G 2 be a cluster point of (pk) (we know from (41) and (35)
that (pk) is bounded), i.e., let z∞ and y∞ be cluster points of (zk) and (yk), respectively. Then let
also (kj ) be an increasing sequence of indexes such that
zkj → z∞ and ykj → y∞ . (44)
In view of (43) and (44), we have
ˆzkj → z∞ and ˆykj → y∞ , (45)
which, in particular, when combined with (42) gives
vkj → z∞ , Lx kj → y∞ , ˜ykj → y∞ and εkj → 0. (46)
From (29), the second statement in (42) (with k = kj) and the ﬁrst statement in (45) we also obtain
z′
kj → z∞ , which combined with Lemma 3.2(b) (with k = kj), the fact that the graph of ∂(f ∗ ◦− L∗) is
closed and the second statement in (46) yields − y∞ ∈ ∂(f ∗ ◦ − L∗)(z∞ ). As a consequence, according
to the deﬁnition of S as in (11), to prove that (z∞ , y ∞ ) ∈ S , it remains to verify that y∞ ∈ ∂g ∗(z∞ ).
To this end, recall ﬁrst that from (22) (with k = kj) we know that vkj ∈ ∂εkj
g(˜ykj ), which is equivalent
to ˜ykj ∈ ∂εkj
g∗(vkj ). Combining the latter inclusion with the ﬁrst, third and fou rth statements in
(46) as well as with the closedness of the graph of the ε-subdiﬀerential of g∗, we obtain the desired
result, namely y∞ ∈ ∂g ∗(z∞ ).
Altogether, we have proved that every cluster point of (pk) belongs to S and so, as we explained
above, it guarantees that (pk) converges to some element in S, i.e., here we ﬁnish the proof of (39).
Finally, the proof of (40) follows trivially from (39), (42) and (43).
13

<!-- page 14 -->
Remark: As we discussed in the Introduction (following Assumption 2 .1), under standard regularity
conditions on (1), the result on (z∞ , w ∞ ) as in Theorem 4.2, gives that there exists x∞ ∈ H such
that x∞ ∈ ∂f ∗(− L∗z∞ ), w∞ = Lx∞ and x∞ and z∞ are solutions of (9) and (10), respectively.
Moreover, the second statements in (39) and (40) give, in par ticular, that Lxk − yk → 0.
We will consider the following two suﬃcient conditions on th e sequences (α k) and/or (pk) to
ensure (38) holds – see (12) and the deﬁnition of pk as in (28) – :
Assumption A: for some 0 < θ < 1 and k0 ≥ 1,
α k ≤ min
{
α, θk
γ− 1∥zk − zk− 1∥2 + γ∥yk − yk− 1∥2
}
, ∀k ≥ k0; (47)
here we adopt the convention 1/ 0 = ∞ .
Assumption B: (α, σ, τ ) ∈ [0, 1) × [0, 1) × (0, 1) and the sequence (α k) satisfy
0 ≤ α k ≤ α k+1 ≤ α < β < 1 ∀k ≥ 0, (48)
where
β := 2η
1 + 2η + √ 1 + 8η (49)
and
η := (1 − τ)(1 − σ )2
4τ . (50)
Lemma 4.3. Under the Assumption B on Algorithm 1 , deﬁne the quadratic real function
q(t) := ( η − 1)t2 − (1 + 2η)t + η ∀t ∈ R. (51)
Then, q(α ) > 0 and, for every p = (z, w ) ∈ S ,
k∑
j=0
∥pj − pj− 1∥2
γ ≤ 2 ∥p0 − p∥2
γ
(1 − α )q(α ) ∀k ≥ 1. (52)
Proof. Note ﬁrst that combining Lemmas 4.1 and 3.4(e) (second inequ ality) and (50) we obtain, for
all k ≥ 0,
hk+1 − hk − α k(hk − hk− 1) + η∥pk+1 − ˆpk∥2
γ ≤ α k(1 + α k)∥pk − pk− 1∥2
γ , (53)
where hk is as in (35). On the other hand, using Lemma 3.2(e), the Cauch y-Schwarz inequality and
the Young inequality 2ab ≤ a2 + b2 with a := ∥pk+1 − pk∥γ and b := ∥pk − pk− 1∥γ we ﬁnd
∥pk+1 − ˆpk∥2
γ = ∥pk+1 − pk∥2
γ + α 2
k∥pk − pk− 1∥2
γ − 2α k⟨pk+1 − pk |pk − pk− 1⟩γ
≥ ∥ pk+1 − pk∥2
γ + α 2
k∥pk − pk− 1∥2
γ − α k (2∥pk+1 − pk∥γ ∥pk − pk− 1∥γ )
≥ (1 − α k)∥pk+1 − pk∥2
γ − α k(1 − α k)∥pk − pk− 1∥2
γ . (54)
14

<!-- page 15 -->
Using (53), (54) and some simple algebraic manipulations we ﬁnd
hk+1 − hk − α k(hk − hk− 1) + η(1 − α k)∥pk+1 − pk∥2
γ ≤ γk∥pk − pk− 1∥2
γ , (55)
where, for all k ≥ 0,
γk := (1 − η)α 2
k + (1 + η)α k. (56)
Deﬁne
µ 0 := (1 − α 0)h0 ≥ 0 and µ k := hk − α k− 1hk− 1 + γk∥pk − pk− 1∥2
γ ∀k ≥ 1, (57)
where hk is as in (35). Using (51), the assumption that (α k) is nondecreasing – see (48) – and
(55)–(57) we obtain, for all k ≥ 1,
µ k − µ k− 1 ≤
[
hk − hk− 1 − α k− 1(hk− 1 − hk− 2) − γk− 1∥pk− 1 − pk− 2∥2
γ
]
+ γk∥pk − pk− 1∥2
γ
≤ [γk − η(1 − α k)] ∥pk − pk− 1∥2
γ
= −
[
(η − 1)α 2
k − (1 + 2η)α k + η
]
∥pk − pk− 1∥2
γ
= − q(α k)∥pk − pk− 1∥2
γ . (58)
Note now that 0 < β < 1 as in (49) is either the smallest or the largest root of the qua dratic function
q(·). Hence, from (48), for all k ≥ 0,
q(α k) ≥ q(α ) > q (β ) = 0 .
The above inequalities combined with (58) yield
∥pk − pk− 1∥2
γ ≤ 1
q(α ) (µ k− 1 − µ k), ∀k ≥ 1, (59)
which combined with (48) and the deﬁnition of µ k as in (57) gives
k∑
j=0
∥pj − pj− 1∥2
γ ≤ 1
q(α ) (µ 0 − µ k),
≤ 1
q(α ) (µ 0 + αh k− 1) ∀k ≥ 1. (60)
Note now that using (48), (57) and (59) we also ﬁnd
µ 0 ≥ . . . ≥ µ k =hk − α k− 1hk− 1 + γk∥pk − pk− 1∥2
γ
≥ hk − αh k− 1, ∀k ≥ 1,
and so
hk ≤ α kh0 + µ 0
1 − α ≤ h0 + µ 0
1 − α ∀k ≥ 0. (61)
Hence, (52) follows directly from (60), (61), the deﬁnition of µ 0 as in (57) and the deﬁnition of h0 as
in (35).
Theorem 4.4 (Second result on the asymptotic convergence of Algorithm 1 ). Under the assumptions
A or B on the sequence (α k), all the conclusions of Theorem 4.2 hold true.
Proof. The proof follows form Theorem 4.2 (see (38)), Assumptions A and B above and Lemma
4.3.
15

<!-- page 16 -->
5 Global convergence rates of Algorithm 1
In this section, we study global convergence rates for Algor ithm 1. We obtain (global) pointwise
O(1/
√
k) and ergodic O(1/k ) rates for residuals; see Theorems 5.2 and 5.4 below.
Lemma 5.1. Consider the sequences evolved by Algorithm 1 and assume that
Assumption B holds.
Let (pk), (˜pk) and (ˆpk) be as in (28) and let also q(·) be as in (51). Then, for all p = (z, w ) ∈ S ,
∥pk − p∥2
γ + τ(1 − τ)(1 − σ )2
k∑
j=0
∥˜pj − ˆpj∥2
γ ≤
(
1 + 2α (1 + α )
(1 − α )2q(α )
)
∥p0 − p∥2
γ .
Proof. From Lemmas 4.1 and A.1(a) (below),
hk +
k∑
j=1
sj ≤ h0 + 1
1 − α
k− 1∑
j=0
δj,
where hk, sk and δk are as in (35), (36) and (37), respectively. Then, in view of ( 52),
∥pk − p∥2
γ + τ(1 − τ)(1 − σ )2
k∑
j=0
∥˜pj − ˆpj∥2
γ ≤ ∥ p0 − p∥2
γ + 1
1 − α
k− 1∑
j=0
α j(1 + α j)∥pj − pj− 1∥2
γ
≤
(
1 + 2α (1 + α )
(1 − α )2q(α )
)
∥p0 − p∥2
γ .
Theorem 5.2 (Pointwise global convergence rates of Algorithm 1) . Consider the sequences evolved
by Algorithm 1 and assume that
Assumption B holds.
Let (z′
k) be as in (29) and let d0 denote the distance of p0 = (z0, y 0) to the solution set S as in (11).
Then, for every k ≥ 0, there exists 0 ≤ i ≤ k such that

















− Lxi ∈ ∂ (f ∗ ◦ − L∗) (z′
i), ˜yi ∈ ∂εig∗(vi),
γ∥Lxi − ˜yi∥2 + 1
γ ∥z′
i − vi∥2 ≤ 2Cd 2
0
k ,
εi ≤ σ 2Cd 2
0
2k ,
(62)
where
C := 1
τ(1 − τ)(1 − σ )2
(
1 + 2α (1 + α )
(1 − α )2q(α )
)
. (63)
16

<!-- page 17 -->
Proof. Let p∗ = (z∗, w ∗) ∈ S be such that d0 = ∥p0 − p∗∥γ . From Lemma 5.1 (with p = p∗) and the
deﬁnition of C > 0 as in (63),
k∑
j=0
∥˜pj − ˆpj∥2
γ ≤ Cd 2
0. (64)
From (13) and Lemmas 3.2(c) and 3.3(b),
γ
2 ∥Lxk − ˜yk∥2 + 1
2γ ∥z′
k − vk∥2 = 1
2 ∥˘pk − ˆpk∥2
γ ≤ ∥ ˜pk − ˆpk∥2
γ . (65)
Due to (13), (24) and the deﬁnitions of ˜pk and ˆpk as in (28) we also ﬁnd
2εk
σ 2 ≤ γ∥Lxk − ˆyk∥2 + 1
γ ∥vk − ˆzk∥2 = ∥˜pk − ˆpk∥2
γ . (66)
Hence, from (64) – (66),
k∑
j=0
∆j ≤ Cd 2
0, (67)
where
∆j := max
{ γ
2 ∥Lxj − ˜yj∥2 + 1
2γ ∥z′
j − vj∥2, 2εj
σ 2
}
, j = 0, . . . , k. (68)
The two inequalities in (62) follow by choosing i ∈ { 0, . . . , k } such that ∆i ≤ ∆j for all j = 0, . . . , k
and using (67) and the deﬁnition of ∆i as in (68). To ﬁnish the proof of the theorem, note that
the inclusions in (62) follow directly from (22) (combined w ith the fact that (∂εk g)− 1 = ∂εk g∗) and
Lemma 3.2(b).
For the sequences generated by Algorithm 1 and (z′
k) as in (29), deﬁne the ergodic means
xa
k := 1
k + 1
k∑
j=0
xj, ˜ya
k := 1
k + 1
k∑
j=0
˜yj,
z′a
k := 1
k + 1
k∑
j=0
z′
j, v a
k := 1
k + 1
k∑
j=0
vj.
(69)
Deﬁne also, for all k ≥ 0,
δa
k := 1
k + 1
k∑
j=0
⟨z′
j |L(xa
k − xj)⟩,
εa
k := 1
k + 1
k∑
j=0
[εj + ⟨˜yj |vj − va
k⟩] .
(70)
17

<!-- page 18 -->
Lemma 5.3. Consider the sequences evolved by Algorithm 1 , let (xa
k), (˜ya
k), (z′a
k ) and (va
k) be as in
(69) and let (δa
k ) and (εa
k) be as in (70). Let also (˘pk), (ˆpk) and (˜pk) be as in (28). For all k ≥ 0,
(a) δa
k , ε a
k ≥ 0 and − Lxa
k ∈ ∂δa
k (f ∗ ◦ − L∗) (z′a
k ), ˜ya
k ∈ ∂εa
k g∗(va
k ).
(b) δa
k + εa
k = 1
k + 1
∑k
j=0 [εj + ⟨˘pj − ˆpj |˜pa
k − ˜pj⟩γ ], where
˜pa
k := 1
k + 1
k∑
j=0
˜pj = (va
k, Lx a
k). (71)
Proof. (a) The desired result follows from [36, Theorem 2.3] and the inclusions in (22) and in Lemma
3.2(b).
(b) In view of Lemma 3.2(c), for j = 0, . . . , k , we have ˘pj − ˆpj =
(
γ(Lxj − ˜yj), 1
γ
(
z′
j − vj
))
and
so by using the deﬁnition of (˜pj) and (71) we get
k∑
j=0
⟨˘pj − ˆpj |˜pa
k − ˜pj⟩γ =
k∑
j=0
[ 1
γ ⟨γ(Lxj − ˜yj) |va
k − vj⟩ + γ⟨ 1
γ (z′
j − vj) |Lxa
k − Lxj⟩
]
=
k∑
j=0
[
⟨Lxj − ˜yj |va
k − vj⟩ + ⟨z′
j − vj |L(xa
k − xj)⟩
]
=
k∑
j=0
[
⟨˜yj |vj − va
k⟩ + ⟨Lxj |va
k⟩ + ⟨z′
j |L(xa
k − xj)⟩ − ⟨ vj |Lxa
k⟩
]
=
k∑
j=0
[
⟨˜yj |vj − va
k⟩ + ⟨z′
j |L(xa
k − xj)⟩
]
.
The desired result now follows by adding the two equations in (70) and using the latter identity.
Theorem 5.4 (Ergodic global convergence rates of Algorithm 1) . Consider the sequences evolved by
Algorithm 1 , let (xa
k), (˜ya
k), (z′a
k ) and (va
k) be as in (69) and let (δa
k ) and (εa
k) be as in (70). Let also
d0 denote the distance of p0 = (z0, y 0) to the solution set S as in (11). Assume that α k ≡ α and that
Assumption B holds.
Then, for all k ≥ 0, δa
k, ε a
k ≥ 0 and

















− Lxa
k ∈ ∂δa
k (f ∗ ◦ − L∗) (z′a
k ), ˜ya
k ∈ ∂εa
k g∗(va
k),
γ∥Lxa
k − ˜ya
k∥2
γ + 1
γ ∥z′a
k − va
k∥2
γ ≤ D2d2
0
k2 ,
δa
k + εa
k ≤ 1
k
( α (1 + α )
τ(1 − α )q(α ) + D(1 + 2
√
3)
√
C
)
d2
0,
(72)
18

<!-- page 19 -->
where C > 0 is as in (63) and
D := 1 + α
τ
(
1 +
√
1 + 2α (1 + α )
(1 − α )2q(α )
)
. (73)
Proof. Note ﬁrst that the inclusions in (72) follow from Lemma 5.3(a ). Now let p∗ = ( z∗, w ∗) ∈ S
be such that d0 = ∥p0 − p∗∥γ . Using Lemma 3.2[(a) and (e)] and the assumption α k ≡ α we ﬁnd
τ(ˆpk − ˘pk) = pk − pk+1 + α (pk − pk− 1), for all k ≥ 0, and so (recall that p0 = p− 1)
τ






k∑
j=0
(ˆpj − ˘pj)






γ
= ∥p0 − pk+1 + α (pk − p0)∥γ ≤ ∥ p0 − pk+1∥γ + α ∥pk − p0∥γ . (74)
In view of Lemma 5.1 (with p = p∗), the deﬁnitions of d0 and D > 0, and the triangle inequality,
∥pk − p0∥γ ≤ ∥ pk − p∗∥γ + ∥p∗ − p0∥γ
≤
(
1 +
√
1 + 2α (1 + α )
(1 − α )2q(α )
)
d0
= τ Dd0
1 + α ∀k ≥ 0, (75)
which combined with (73) and (74) yields






k∑
j=0
(ˆpj − ˘pj)






γ
≤ Dd0. (76)
Recall that from Lemma 3.2(c) we have
˘pk − ˆpk =
(
γ(Lxk − ˜yk), 1
γ
(
z′
k − vk
) )
,
and so from the deﬁnitions of ergodic means as in (69), we ﬁnd
1
k + 1
k∑
j=0
(˘pj − ˆpj) =
(
γ(Lxa
k − ˜ya
k), 1
γ (z′a
k − va
k)
)
.
Hence, in view of (13) and (76),
γ∥Lxa
k − ˜ya
k∥2
γ + 1
γ ∥z′a
k − va
k∥2
γ = 1
(k + 1)2






k∑
j=0
(˘pj − ˆpj)






2
γ
≤ D2d2
0
k2 ,
which gives the ﬁrst inequality in (72).
Now let’s prove the second inequality in (72). To this end, le t p = (z, w ) ∈ G 2 and ﬁrst note that
from Lemma 3.4(b),
∥p − ˆpk∥2
γ − ∥ p − ˘pk∥2
γ ≥ 2 [εk + ⟨˘pk − ˆpk |p − ˜pk⟩γ ] . (77)
19

<!-- page 20 -->
By Lemma 3.2(a) and the identity ∥(1 − τ)a + τ b∥2
γ = (1 − τ)∥a∥2
γ + τ∥b∥2
γ − (1 − τ)τ∥a − b∥2
γ with
a = p − ˆpk and b = p − ˘pk,
∥p − pk+1∥2
γ = (1 − τ)∥p − ˆpk∥2
γ + τ∥p − ˘pk∥2
γ − (1 − τ)τ∥˘pk − ˆpk∥2
γ .
Now multiplying (77) by τ > 0 and using the latter identity we get
∥p − ˆpk∥2
γ − ∥ p − pk+1∥2
γ ≥ 2τ [εk + ⟨˘pk − ˆpk |p − ˜pk⟩γ ] + (1 − τ)τ∥˘pk − ˆpk∥2
γ
≥ 2τ [εk + ⟨˘pk − ˆpk |p − ˜pk⟩γ ] . (78)
Note that from Lemma 3.5 the assumption α k ≡ α we obtain
∥p − ˆpk∥2
γ = (1 + α )∥p − pk∥2
γ − α ∥p − pk− 1∥2
γ + α (1 + α )∥pk − pk− 1∥2
γ . (79)
Making the substitution of (79) into (78) and after some simp le algebra, we ﬁnd (now replacing the
index k ≥ 0 by j ≥ 0),
∥p − pj∥2
γ − ∥ p − pj+1∥2
γ + α (1 + α )∥pj − pj− 1∥2
γ ≥ 2τ [εj + ⟨˘pj − ˆpj |p − ˜pj⟩γ ]
− α
[
∥p − pj∥2
γ − ∥ p − pj− 1∥2
γ
]
.
Summing up the latter inequality from j = 0, . . . , k and with p = ˜pa
k — see (71) –, and using Lemma
5.3(b),
∥˜pa
k − p0∥2
γ − ∥ ˜pa
k − pk+1∥2
γ + α (1 + α )
k∑
j=0
∥pj − pj− 1∥2
γ
≥ 2τ
k∑
j=0
[εj + ⟨˘pj − ˆpj |˜pa
k − ˜pj⟩γ ] − α
[
∥˜pa
k − pk∥2
γ − ∥ ˜pa
k − p0∥2
γ
]
= 2τ(k + 1)(δa
k + εa
k) − α
[
∥˜pa
k − pk∥2
γ − ∥ ˜pa
k − p0∥2
γ
]
,
which combined with (52) (with p = p∗) and the deﬁnition of d0 yields
2τ(k + 1)(δa
k + εa
k) − 2α (1 + α )d2
0
(1 − α )q(α ) ≤
[
∥˜pa
k − p0∥2
γ − ∥ ˜pa
k − pk+1∥2
γ
]
+ α
[
∥˜pa
k − pk∥2
γ − ∥ ˜pa
k − p0∥2
γ
]
.
Now using the inequality ∥a∥2
γ − ∥ b∥2
γ ≤ 2∥a∥γ ∥a − b∥γ (in both terms in the right-hand side of the
latter inequality) and (75) we ﬁnd
2τ(k + 1)(δa
k + εa
k) − 2α (1 + α ) d2
0
(1 − α )q(α ) ≤ 2∥˜pa
k − p0∥γ ∥pk+1 − p0∥γ + 2α ∥˜pa
k − pk∥γ ∥pk − p0∥γ
≤ 2τ Dd0
1 + α
(
∥˜pa
k − p0∥γ + α ∥˜pa
k − pk∥γ
)
≤ 2τ Dd0 max{∥˜pa
k − p0∥γ , ∥˜pa
k − pk∥γ }. (80)
20

<!-- page 21 -->
Now deﬁne, for all k ≥ 0,
˘pa
k := 1
k + 1
k∑
j=0
˘pj. (81)
Using (24), (28), Lemma 3.3(a), (71), (81) as well as the conv exity of ∥·∥2
γ we ﬁnd
∥˜pa
k − ˘pa
k∥2
γ ≤ 1
k + 1
k∑
j=0
∥˜pj − ˘pj∥2
γ
= 1
k + 1
k∑
j=0
[ 1
γ ∥ej ∥2 + 1
γ ∥vj − ˆzj∥2
]
≤ 1
k + 1
k∑
j=0
[
γ∥Lxj − ˆyj∥2 + 1
γ ∥vj − ˆzj∥2
]
= 1
k + 1
k∑
j=0
∥˜pj − ˆpj∥2
γ
≤ Cd 2
0
k + 1
≤ Cd 2
0, (82)
where we used Lemma 5.1 (with p = p∗) and the deﬁnition of C > 0 as in (63). Similarly, using
(81), the triangle inequality and the well-known inequalit y ∥a + b∥2
γ ≤ 2
(
∥a∥2
γ + ∥b∥2
γ
)
, we get, for
all ℓ ≥ 1,
∥pℓ − ˘pa
k∥2
γ =






1
k + 1
k∑
j=0
(pℓ − ˘pj)






2
γ
≤ 1
k + 1
k∑
j=0
∥pℓ − ˘pj∥2
γ
≤ 2
k + 1
k∑
j=0
(
∥pℓ − pj+1∥2
γ + ∥pj+1 − ˘pj∥2
γ
)
= 2
k + 1
k∑
j=0
∥pℓ − pj+1∥2
γ + 2
k + 1
k∑
j=0
∥pj+1 − ˘pj∥2
γ . (83)
Using (again) Lemma 5.1 (with p = p∗) as well as the deﬁnition of C > 0 as in (63), we ﬁnd, for all
ℓ, j ≥ 0,
∥pℓ − pj+1∥2
γ ≤ 2
(
∥pℓ − p∗∥2
γ + ∥pj+1 − p∗∥2
γ
)
≤ 2τ(1 − τ)(1 − σ )2Cd 2
0. (84)
21

<!-- page 22 -->
Recall also that from Lemma 3.2(b) we have pj+1 − ˘pj = (1 − τ)(ˆpj − ˘pj) (0 ≤ j ≤ k) and so from
Lemma 3.3(a), for all j ≥ 0,
∥pj+1 − ˘pj∥γ = (1 − τ)∥ˆpj − ˘pj∥γ ≤ 2(1 − τ)∥˜pj − ˆpj∥γ ,
which then yields (by Lemma 5.1 (with p = p∗) and the deﬁnition of C > 0)
k∑
j=0
∥pj+1 − ˘pj∥2
γ = 4(1 − τ)2
k∑
j=0
∥˜pj − ˆpj∥2
γ ≤ 4(1 − τ)2Cd 2
0. (85)
Putting it all together, from (83) – (85) we obtain, for all ℓ ≥ 0,
∥pℓ − ˘pa
k∥2
γ ≤ 4τ(1 − τ)(1 − σ )2Cd 2
0 + 8
k + 1 (1 − τ)2Cd 2
0
= 4(1 − τ)
(
τ(1 − σ 2) + 2
k + 1 (1 − τ)
)
Cd 2
0
≤ 12Cd 2
0. (86)
Using now the triangle inequality, (82) and (86) we ﬁnd, for a ll ℓ ≥ 0,
∥pℓ − ˜pa
k∥γ ≤ ∥ pℓ − ˘pa
k∥γ + ∥˘pa
k − ˜pa
k∥γ ≤ (1 + 2
√
3)
√
Cd 0. (87)
Finally, using (80) and (87) (with ℓ = 0 and ℓ = k),
2τ(k + 1)(δa
k + εa
k) − 2α (1 + α ) d2
0
(1 − α )q(α ) ≤ 2τ D(1 + 2
√
3)
√
Cd2
0
so that
δa
k + εa
k ≤ 1
k
( α (1 + α )
τ(1 − α )q(α ) + D(1 + 2
√
3)
√
C
)
d2
0.
6 Numerical Experiments
This section presents some numerical experiments on the LAS SO problem, which is an instance of
the minimization problem (1). We compared Algorithm 1 from t his paper with and without inertial
eﬀects; they are called Inexact ADMM and Inexact inertial ADMM , respectively. We implemented
both algorithms in Matlab R2021a and, for both algorithms an d all problem classes, used the same
stopping criterion, namely
dist∞
(
0, ∂f (xk) + ∂g (xk)
)
≤ ε, (88)
where dist ∞ (0, S ) := inf {∥s∥∞ |s ∈ S} and ε is a tolerance parameter set to 10− 6.
The inertial parameter α k (as in step 1 of Algorithm 1) is updated according to the rule ( 47) with
θ = 0. 99 and k0 = 1. More precisely, we choose α k as
α k = min
{
α, θk
γ− 1∥zk − zk− 1∥2 + γ∥yk − yk− 1∥2
}
, ∀k ≥ 1,
22

<!-- page 23 -->
where 0 ≤ α < 1. The source codes are available under request (marina.gere mia@ifsc.edu.br).
The LASSO problem. We perform numerical experiments on the LASSO problem (as al ready
discussed in (17)), namely
min
x∈ Rd
{ 1
2 ∥Ax − b∥2 + ν∥x∥1
}
, (89)
where A ∈ Rn× d, b ∈ Rn and ν > 0. For the data matrix A and the vector b, we used ﬁve categories of
non-artiﬁcial datasets (available at the UCI Machine Learn ing Repository, https://archive.ics.uci.edu):
BlogFeedback: This category consists of one standard microarray dataset s that contain features
extracted from a blog post from [15]. This problem is called blogFeedback (with n = 60021 and
d = 280).
Breast Cancer Wisconsin (Prognostic) : This category consists of one prognostic Wisconsin
breast cancer database from [48], which has a dense matrix A. Each row represents follow-up
data for one breast cancer case. This problem is called Wisconsin (with n = 198 and d = 33).
DrivFace: This category comprises a single standard microarray data set containing image se-
quences of individuals driving in real-world scenarios fro m [34]. This problem is called DrivFace
(with n = 606 and d = 6400) and has a dense matrix A.
Gene expression : This category consists of six standard cancer DNA microarr ay datasets
from [19], which have dense and wide matrices A, with the number of rows n ∈ [42, 102]
and the number of columns d ∈ [2000, 6033]. These problems are called brain (with n = 42 and
d = 5597), colon (with n = 62 and d = 2000), leukemia (with n = 72 and d = 3571), lymphoma
(with n = 62 and d = 4026), prostate (with n = 102 and d = 6033) and srbct (with n = 63 and
d = 2308).
Single-Pixel camera: This category consists of four compressed image sensing da tasets from [21],
which have dense and wide matrices A, with n ∈ { 410, 1638} and d ∈ { 1024, 4096}]. These
problems are called Ball64_singlepixcam (with n = 1638 and d = 4096), Logo64_singlepixcam
(with n = 1638 and d = 4096 ), Mug32_singlepixcam (with n = 410 and d = 1024 ) and
Mug128_singlepixcam (with n = 410 and d = 1024).
We implemented both algorithms Inexact ADMM and Inexact inertial ADMM in Matlab R2021a,
combined with a CG procedure to approximately solve the subp roblems (21); see also the fourth
remark following Deﬁnition 3.1. As usual (see, e.g., [3]), w e solved the (easy) subproblem (20) by
using the standard-soft thresholding operator. We also set (α, σ, τ, γ ) = (0 . 33, 0. 99, 0. 999, 1) and
(σ, τ, γ ) = (0 . 99, 0. 999, 1) for Inexact inertial ADMM and Inexact ADMM , respectively. Moreover, as
in [13], we set the regularization parameter ν as 0. 1∥AT b∥∞ , and scaled the vector b and the columns
of matrix A to have ℓ2 unit norm.
Table 1 shows the number of outer iterations required by each algorithm on each problem instance,
the cumulative total number of inner iterations required by the CG algorithm for solving (21) and
runtimes in seconds demanded by each algorithm to achieve th e prescribed tolerance as in (88). From
Table 1, we see that Inexact inertial ADMM outperforms Inexact ADMM on average by about 30%,
25% and 25% on “Outer iterations”, “Total inner iterations” and “Runti mes”, respectively.
23

<!-- page 24 -->
Table 1: Comparison of performance in the LASSO problem
Inexact ADMM Inexact inertial ADMM
Outer
iterations
(outer1)
Total inner
iterations
(inner1)
Time
(time1)
Outer
iterations
(outer2)
Total inner
iterations
(inner2)
Time
(time2)
outer2
outer1
inner2
inner1
time2
time1
Brain 2923 14007 3.97 2120 11112 3.11 0.7253 0.7933 0.7834
Colon 505 2818 0.39 347 1866 0.27 0.6871 0.6622 0.6923
Leukemia 764 3695 0.83 544 2861 0.62 0.7121 0.7743 0.7469
Lymphoma 1101 6091 1.45 862 5119 1.13 0.7829 0.8404 0.7793
Prostate 2006 8328 4.13 1414 6561 3.46 0.7049 0.7878 0.8378
Srbct 511 3554 0.55 346 2293 0.38 0.6771 0.6452 0.6909
Ball64 313 490 8.77 216 325 5.93 0.6901 0.6633 0.6762
Logo64 316 495 8.47 221 359 5.98 0.6994 0.7253 0.7061
Mug32 134 282 0.09 88 197 0.06 0.6567 0.6986 0.6667
Mug128 955 1163 248.04 822 986 212.21 0.8607 0.8478 0.8555
DrivFace 2682 18727 165.47 1803 13284 115.36 0.6723 0.7094 0.6972
Wisconsin 285 605 0.07 182 456 0.05 0.6386 0.7537 0.7143
blogFeedback 386 1108 46.79 317 938 38.84 0.8212 0.8466 0.8301
Geometric mean 662.68 2187.47 3.21 473.79 1633.23 2.38 0.7149 0.7466 0.7414
A Auxiliary results
The following lemma was essentially proved by Alvarez and At touch in [2, Theorem 2.1] (see also [4,
Lemma A.4]).
Lemma A.1. Let the sequences (hk), (sk), (α k) and (δk) in [0, ∞ ) and α ∈ R be such that h0 = h− 1,
0 ≤ α k ≤ α < 1 and
hk+1 − hk + sk+1 ≤ α k(hk − hk− 1) + δk ∀k ≥ 0. (90)
The following hold:
(a) For all k ≥ 1,
hk +
k∑
j=1
sj ≤ h0 + 1
1 − α
k− 1∑
j=0
δj. (91)
(b) If ∑∞
k=0 δk < ∞ , then limk→∞ hk exists, i.e., the sequence (hk) converges to some element in
[0, ∞ ).
Lemma A.2 (Opial [40]) . Let H be a ﬁnite dimensional inner product space, let ∅ ̸ = S ⊂ H and
let {pk} be a sequence in H such that every cluster point of {pk} belongs to S and limk→∞ ∥pk − p∥
exists for every p ∈ S . Then {pk} converges to a point in S.
References
[1] V. A. Adona, M. L. N. Gonçalves, and J. G. Melo. A partially inexact proximal alternating
direction method of multipliers and its iteration-complex ity analysis. J. Optim. Theory Appl. ,
182(2):640–666, 2019.
24

<!-- page 25 -->
[2] F. Alvarez and H. Attouch. An inertial proximal method fo r maximal monotone operators via
discretization of a nonlinear oscillator with damping. Set-Valued Anal., 9(1-2):3–11, 2001.
[3] M. M. Alves, J. Eckstein, M. Geremia, and J.G. Melo. Relat ive-error inertial-relaxed inexact ver-
sions of Douglas-Rachford and ADMM splitting algorithms. Comput. Optim. Appl. , 75(2):389–
422, 2020.
[4] M. M. Alves and R. T. Marcavillaca. On inexact relative-e rror hybrid proximal extragradient,
forward-backward and Tseng’s modiﬁed forward-backward me thods with inertial eﬀects. Set-
Valued Var. Anal., 28(2):301–325, 2020.
[5] H. Attouch. Fast inertial proximal ADMM algorithms for c onvex structured optimization with
linear constraint. Minimax Theory Appl. , 6(1):1–24, 2021.
[6] H. Attouch and A. Cabot. Convergence of a relaxed inertia l proximal algorithm for maximally
monotone operators. Math. Program., 184(1-2, Ser. A):243–287, 2020.
[7] H. Attouch, A. Cabot, Z. Chbani, and H. Riahi. Inertial fo rward-backward algorithms with
perturbations: application to Tikhonov regularization. J. Optim. Theory Appl. , 179(1):1–36,
2018.
[8] H. Attouch and J. Peypouquet. Convergence of inertial dy namics and proximal algorithms
governed by maximally monotone operators. Math. Program., 174(1-2, Ser. B):391–432, 2019.
[9] H. Attouch and M. Soueycatt. Augmented Lagrangian and pr oximal alternating direction meth-
ods of multipliers in Hilbert spaces. Applications to games , PDE’s and control. Pac. J. Optim. ,
5(1):17–37, 2009.
[10] H. H. Bauschke and P. L. Combettes. Convex analysis and monotone operator theory in Hilbert
spaces. CMS Books in Mathematics/Ouvrages de Mathématiques de la SM C. Springer, New
York, 2011.
[11] M. Benning, F. Knoll, C.-B. Schönlieb, and T. Valkonen. Pr econditioned ADMM with nonlin-
ear operator constraint. In Lorena Bociu, Jean-Antoine Dési déri, and Abderrahmane Habbal,
editors, System Modeling and Optimization , pages 117–126, Cham, 2016. Springer International
Publishing.
[12] R. I. Boţ and E. R. Csetnek. ADMM for monotone operators: c onvergence analysis and rates.
Adv. Comput. Math. , 45(1):327–359, 2019.
[13] S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein. Dis tributed optimization and statistical
learning via the alternating direction method of multiplie rs. Found. Trends Mach. Learn., 3(1):1–
122, 2011.
[14] K. Bredies and H. Sun. A proximal point analysis of the pre conditioned alternating direction
method of multipliers. J. Optim. Theory Appl. , 173(3):878–907, 2017.
[15] K. Buza. BlogFeedback. UCI Machine Learning Repository , 2014.
[16] L. Chen, X. Li, D. Sun, and K.-C. Toh. On the equivalence o f inexact proximal ALM and ADMM
for a class of convex composite programming. Math. Program., 185(1-2, Ser. A):111–161, 2021.
25

<!-- page 26 -->
[17] P. L. Combettes and L. E. Glaudin. Quasi-nonexpansive i terations on the aﬃne hull of orbits:
from Mann’s mean value algorithm to inertial methods. SIAM J. Optim. , 27(4):2356–2380, 2017.
[18] W. Deng and W. Yin. On the global and linear convergence o f the generalized alternating
direction method of multipliers. J. Sci. Comput. , 66(3):889–916, 2016.
[19] M. Dettling and P. Bühlmann. Finding predictive gene gro ups from microarray data. J. Multi-
variate Anal. , 90(1):106–131, 2004.
[20] D. Dua and C. Graﬀ. UCI machine learning repository, 201 7.
[21] M. F. Duarte, M. A. Davenport, Dharmpal T., J. N. Laska, T ing S., K. F. Kelly, and R. G.
Baraniuk. Single-pixel imaging via compressive sampling: Bu ilding simpler, smaller, and less-
expensive digital cameras. IEEE Signal Processing Magazine , 25(2):83–91, 2008.
[22] J. Eckstein and D. P. Bertsekas. On the Douglas-Rachford splitting method and the proximal
point algorithm for maximal monotone operators. Math. Program., 55(3):293–318, 1992.
[23] J. Eckstein and W. Yao. Relative-error approximate ver sions of Douglas-Rachford splitting and
special cases of the ADMM. Math. Program., 170(2, Ser. A):417–444, 2018.
[24] J. Eckstein and W. Yao. Relative-error approximate ver sions of Douglas-Rachford splitting and
special cases of the ADMM. Math. Program., 170(2):417–444, 2018.
[25] R.-E. Fan, P.-H. Chen, and C.-J. Lin. Working set select ion using second order information for
training support vector machines. Journal of Machine Learning Research , 6:1889–1918, 12 2005.
[26] M. Fortin and R. Glowinski. On decomposition-coordina tion methods using an augmented La-
grangian. In M. Fortin and R. Glowinski, editors, Augmented Lagrangian methods: Applications
to the numerical solution of boundary-value problems , volume 15 of Studies in Mathematics and
its Applications , pages 97–146. North-Holland, Amsterdam, 1983.
[27] J. Franklin. The elements of statistical learning: dat a mining, inference and prediction. Math.
Intelligencer, 27(2):83–85, 2005.
[28] D. Gabay. Applications of the method of multipliers to v ariational inequalities. In M. Fortin and
R. Glowinski, editors, Augmented Lagrangian methods: Applications to the numeric al solution
of boundary-value problems , volume 15 of Studies in Mathematics and its Applications , pages
299–331. North-Holland, Amsterdam, 1983.
[29] R. Glowinski and A. Marroco. Sur l’approximation, par é léments ﬁnis d’ordre 1, et la résolution,
par pénalisation-dualité, d’une classe de problèmes de Dir ichlet non linéaires. C. R. Acad. Sci.
Paris Sér. A , 278:1649–1652, 1974.
[30] R. Glowinski, S. J. Osher, and W. Yin, editors. Splitting methods in communication, imaging,
science, and engineering . Scientiﬁc Computation. Springer, Cham, 2016.
[31] I. Guyon, S. Gunn, A. Ben-Hur, and G. Dror. Result analysi s of the nips 2003 feature selection
challenge. volume 17, 01 2004.
26

<!-- page 27 -->
[32] W. W. Hager and H. Zhang. Convergence rates for an inexac t ADMM applied to separable
convex optimization. Comput. Optim. Appl. , 77(3):729–754, 2020.
[33] B. He and X. Yuan. On the O(1/n ) convergence rate of the Douglas-Rachford alternating
direction method. SIAM J. Numer. Anal. , 50(2):700–709, 2012.
[34] A. Hernndez-Sabat, A. Lpez, and K. Diaz-Chito. DrivFac e. UCI Machine Learning Repository ,
2016.
[35] D. A. Lorenz and Q. Tran-Dinh. Non-stationary Douglas- Rachford and alternating direction
method of multipliers: adaptive step-sizes and convergenc e. Comput. Optim. Appl. , 74(1):67–92,
2019.
[36] R. D. C. Monteiro and B. F. Svaiter. On the complexity of th e hybrid proximal extragradient
method for the iterates and the ergodic mean. SIAM J. Optim. , 20(6):2755–2787, 2010.
[37] R. D. C. Monteiro and B. F. Svaiter. Iteration-complexit y of block-decomposition algorithms
and the alternating direction method of multipliers. SIAM J. Optim. , 23(1):475–507, 2013.
[38] A. Y. Ng. Feature selection, L1 vs. L2 regularization, and rotational invariance. In Proceedings,
Twenty-First International Conference on Machine Learnin g, ICML , pages 615–622, 2004.
[39] J. Nocedal and S. J. Wright. Numerical optimization . Springer Series in Operations Research
and Financial Engineering. Springer, New York, second edit ion, 2006.
[40] Z. Opial. Weak convergence of the sequence of successiv e approximations for nonexpansive
mappings. Bull. Amer. Math. Soc. , 73:591–597, 1967.
[41] R. T. Rockafellar. Convex Analysis . Princeton University Press, Princeton, NJ, 1970.
[42] R. Sheﬁ and M. Teboulle. Rate of convergence analysis of decomposition methods based on the
proximal method of multipliers for convex minimization. SIAM J. Optim. , 24(1):269–297, 2014.
[43] M. V. Solodov and B. F. Svaiter. A hybrid approximate extr agradient-proximal point algorithm
using the enlargement of a maximal monotone operator. Set-Valued Anal., 7(4):323–345, 1999.
[44] M. V. Solodov and B. F. Svaiter. A hybrid projection-prox imal point algorithm. J. Convex
Anal., 6(1):59–70, 1999.
[45] M. V. Solodov and B. F. Svaiter. A uniﬁed framework for som e inexact proximal point algorithms.
Numer. Funct. Anal. Optim. , 22(7-8):1013–1035, 2001.
[46] B. F. Svaiter. A partially inexact ADMM with o(1/n ) asymptotic convergence rate, O(1/n )
complexity, and immediate relative error tolerance. Optimization, 70(10):2061–2080, 2021.
[47] R. Tibshirani. Regression shrinkage and selection via the lasso. J. Roy. Statist. Soc. Ser. B ,
58(1):267–288, 1996.
[48] W. Wolberg, W. Street, and O. Mangasarian. Breast Cancer Wisconsin (Prognostic). UCI
Machine Learning Repository , 1995.
27

<!-- page 28 -->
[49] J. Xie. On inexact ADMMs with relative error criteria. Comput. Optim. Appl. , 71(3):743–765,
2018.
[50] J. Xie, A. Liao, and X. Yang. An inexact alternating dire ction method of multipliers with relative
error criteria. Optim. Lett., 11(3):583–596, 2017.
[51] C. Zălinescu. Convex analysis in general vector spaces . World Scientiﬁc Publishing Co., Inc.,
River Edge, NJ, 2002.
28