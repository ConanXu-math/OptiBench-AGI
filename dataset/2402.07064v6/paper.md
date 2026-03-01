# Piecewise SOS-Convex Moment Optimization and Applications via Exact Semi-Definite Programs

**arXiv ID:** 2402.07064v6

**Authors:** Queenie Yingkun Huang, Vaithilingam Jeyakumar, Guoyin Li

**Abstract:** This paper presents exact Semi-Definite Program (SDP) reformulations for infinite-dimensional moment optimization problems involving a new class of piecewise Sum-of-Squares (SOS)-convex functions and projected spectrahedral support sets. These reformulations show that solving a single SDP finds the optimal value and an optimal probability measure of the original moment problem. This is done by establishing an SOS representation for the non-negativity of a piecewise SOS-convex function over a projected spectrahedron. Finally, as an application and a proof-of-concept illustration, the paper presents numerical results for the Newsvendor and revenue maximization problems with higher-order moments by solving their equivalent SDP reformulations. These reformulations promise a flexible and efficient approach to solving these models. The main novelty of the present work in relation to the recent research lies in finding the solution to moment problems, for the first time, with piecewise SOS-convex functions from their numerically tractable exact SDP reformulations.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
Piecewise SOS-Convex Moment Optimization and
Applications via Exact Semi-Definite Programs
Q.Y. Huang∗, V. Jeyakumar † and G. Li ‡
July 3, 2024
Abstract
This paper presents exact Semi-Definite Program (SDP) reformulations for infinite-dimensional
moment optimization problems involving a new class of piecewise Sum-of-Squares (SOS)-
convex functions and projected spectrahedral support sets. These reformulations show that
solving a single SDP finds the optimal value and an optimal probability measure of the
original moment problem. This is done by establishing an SOS representation for the non-
negativity of a piecewise SOS-convex function over a projected spectrahedron. Finally, as an
application and a proof-of-concept illustration, the paper also presents numerical results for
the Newsvendor and revenue maximization problems with higher-order moments by solving
their equivalent SDP reformulations. These reformulations promise a flexible and efficient
approach to solving these models. The main novelty of the present work in relation to the
recent research lies in finding the solution to moment problems, for the first time, with
piecewise SOS-convex functions from their numerically tractable exact SDP reformulations.
Keywords. Moment Optimization; Sum-of-Squares Convex Polynomials; Piecewise Functions;
Generalized moment problems; Semi-Definite Programming
1 Introduction
Consider the generalized moment problem
min
µ∈PΩ
n
EΩ
µ

min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (ω)

: EΩ
µ[hj(ω)] ≤ γj, j = 1, . . . , J
o
(P)
where gk
ℓ , ℓ = 1 , . . . , L, k = 1 , . . . , r, and hj, j = 1 , . . . , J, are Sum-of-Squares (SOS)-convex
polynomials. The SOS-convexity for polynomials is a numerically tractable relaxation of convexity
∗Corresponding author. Department of Applied Mathematics, University of New South Wales, Sydney 2052,
Australia. Email: yingkun.huang@unsw.edu.au. The research of Ms Queenie Yingkun Huang was partially
supported by a grant from the Australian Research Council.
†Department of Applied Mathematics, University of New South Wales, Sydney 2052, Australia. Email: v.
jeyakumar@unsw.edu.au. The research of Prof Vaithilingam Jeyakumar was supported by a grant from the
Australian Research Council.
‡Department of Applied Mathematics, University of New South Wales, Sydney 2052, Australia. Email: g.li@
unsw.edu.au. Research of Prof Guoyin Li was supported by a grant from the Australian Research Council
1
arXiv:2402.07064v6  [math.OC]  2 Jul 2024

<!-- page 2 -->
since checking whether a polynomial is SOS-convex or not can be achieved by solving an SDP
[1, 15]. Recent applications of SOS-convexity in optimization may be found in [25, 18]. The set
Ω ⊂ Rm is a convex compact projected spectrahedron and PΩ is the set of probability measures
supported on Ω. The expectation of a random variable ω with respect to the probability measure
µ ∈ PΩ is denoted by EΩ
µ[ω].
Our model (P) is closely related to generalized moment problems [24] where gk
ℓ , ℓ = 1 , . . . , L,
k = 1, . . . , r, and hj, j = 1, . . . , J, are arbitrary real-valued continuous functions and the support
set Ω is an arbitrary compact subset of Rm. The dominant technique to reformulate a general-
ized moment problem as a semi-infinite optimization problem is via the duality theory [32, 36].
However, in general, a semi-infinite optimization problem can be computationally intractable.
The classical moment problem can be recovered from (P) if L = 1, r = 1, the random variable is
univariate, and the constraints are the first J moments [24]. The multivariate classical moment
problem extends to countably many moments, but numerically tractable forms are limited. When
L = 1, r = 1, g1
1 and hj, j = 1, . . . , J, are linear functions and the support Ω is a spectrahedron,
(P) admits an exact SDP reformulation [17]. When L = 1, gk
1, k = 1, . . . , r, are linear functions
with mean and variance constraints, the support Ω is an ellipsoid, it has been shown that (P)
admits an exact SDP reformulation [8].
When L = 1, r = 1, g1
1 and hj, j = 1 , . . . , J, are any polynomials and the support Ω is a
compact semi-algebraic set, the optimal solution for (P) can be found by solving a convergent
hierarchy of SDPs [22]. Further, polynomial optimization techniques were employed in [6] to
obtain exact SOS reformulations when P contains probability measures whose density functions
are SOS polynomials. A review of generalized moment problems and their applications is available
in [7].
The model problem (P) is of great interest in distributionally robust optimization [36], where
(P) appears as an uncertainty quantification problem to mitigate risks and uncertainties from the
input data [36, 13]. For example, in the Newsvendor model of [33], the cost of ordering represented
as a piecewise linear loss function is minimized to meet the uncertain demand. In the portfolio
management problem of [8], the profit of the portfolio represented as a piecewise linear utility
function is maximized amid random investment returns.
Motivated by the importance of such piecewise functions in optimization and the usefulness of
the minimax functions of the form mink=1,...,r maxℓ=1,...,L gk
ℓ in applications across diverse domains,
we introduce a new class of functions, called the piecewise SOS-convex functions. To the best of
our knowledge, this notion of piecewise SOS-convex function has so far not been studied in the
literature.
A function f on Rn is piecewise SOS-convex if there exist SOS-convex polynomialsgk
ℓ , ℓ = 1, . . . , L,
k = 1, . . . , r, on Rn such that f(v) = min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v) for all v ∈ Rn.
Note that a piecewise SOS-convex function is not necessarily a convex or a differentiable function
and it covers a broad range of functions that appear in applications across several domains, as
outlined below. Details of their piecewise representations are provided in Appendix A.
• The truncated ℓ1-norm, f(v) = min{1, ε|v|}, ε > 0, which is non-convex and non-smooth, is
used extensively in machine-learning regularization [26].
• The piecewise linear function f(v) = min
k=1,...,r
p⊤
k v + qk, pk ∈ Rn, qk ∈ R, k = 1, . . . , r, appears
in Newsvendor-type [3, 10], conditional value-at-risk [30], and portfolio selection [8] models,
2

<!-- page 3 -->
and it is a useful approximation for various important utility functions [8].
• The class of max-SOS-convex functions f(v) = max
ℓ=1,...,L
gℓ(v), where gℓ’s are SOS-convex
polynomials, is studied in [21]. Further, convex quadratics and convex separable polynomials
are SOS-convex [1, 19].
• The class of difference-of-convex functions of the form f(v) = max
ℓ=1,...,L
fℓ(v)− max
k=1,...,r
(p⊤
k v+qk),
where fℓ’s are SOS-convex polynomials, can be expressed as f(v) = min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v) for
gk
ℓ (v) = fℓ(v) − (p⊤
k v + qk), ℓ = 1 , . . . , L, k = 1 , . . . , r. These functions have attracted
applications in feature selection models of machine-learning [18].
• The piecewise convex quadratic function
f(v) =
(
a(v − b)2 + c, if 0 ≤ v ≤ b,
c, if v > b.
with a, b ≥ 0, c ∈ R, is used to model offer prices from customers [13].
• The Huber loss function with parameter ε > 0,
f(v) =
(
1
2 v2, if |v| ≤ ε,
ε|v| − 1
2 ε2, otherwise,
is studied in robust statistics [33].
The present study of moment optimization is motivated by two aspects. Firstly, moment prob-
lems are numerically challenging due to the presence of infinite-dimensional distributions and
multi-dimensional integrals. In response, we examine equivalent Semi-Definite Program (SDP) re-
formulations of these problems which can be efficiently solved numerically by commonly available
software via interior point methods.
Secondly, the study of piecewise SOS-convexity in (P) was stimulated by the computationally
attractive features of SOS-convex polynomial optimization and its applications in many areas. It
is known that SOS-convex polynomial optimization problems admit numerically favourable SDP
reformulations [19, 20, 21, 25].
The main goal of this paper is to provide numerically tractable reformulations for the moment
problem involving piecewise SOS-convex functions and show how the solution can be recovered
from its SDP reformulation.
Our contributions. Our main contributions are as follows.
(i) Firstly, we introduce the notion of piecewise SOS-convex functions and establish a new SOS
representation of the non-negativity of a piecewise SOS-convex function over a compact pro-
jected spectrahedron. This result not only generalizes the known representation result of
a non-negative SOS-convex polynomial over a spectrahedron [19] but also provides numeri-
cally tractable representations of non-negativity for a broad class of functions that are not
necessarily convex.
(ii) Secondly, exploiting the representation result, we derive an equivalent numerically tractable
SOS optimization reformulation for the generalized moment problem (P). Moreover, under
3

<!-- page 4 -->
suitable conditions, we show how to recover an optimal probability measure for (P) by
solving a single SDP. This extends the corresponding result of Lasserre [23] on SOS-convex
polynomial optimization.
(iii) Finally, as an application and a proof-of-concept illustration, we present numerical results
for two important practical models, i.e., the Newsvendor problem [12] with higher-order
moments and the revenue maximization problem [13] with a higher-order utility function,
by solving their equivalent SDP reformulations.
The novelty of the present work in relation to recent research in SOS-convexity and moment
optimization is the derivation of numerically tractable exact SDP reformulations for moment
problems involving, for the first time, piecewise SOS-convex functions . The main innovation is
the combined utilization of powerful tools from real algebraic geometry, convex analysis, and SOS
polynomials to produce the exact reformulations. We employ the key computationally favourable
features of SOS-convexity within an infinite-dimensional setting and exploit the geometry of the
projected spectrahedra to facilitate the SDP reformulations.
Our approach makes use of the piecewise structure of the functions for transforming infinite-
dimensional problems into finite-dimensional numerically tractable problems. Moreover, the pro-
jected spectrahedron covers a broad class of semi-algebraic sets [25], such as the spectrahedra,
ellipsoids, and boxes, used in robust optimization [2].
The paper is organized as follows. Section 2 provides SOS representation results for piecewise
SOS-convex functions. Section 3 presents the main theorems for the SOS reformulation of (P)
and the optimal solution recovery of (P). Section 4 describes applications to the Newsvendor
and revenue maximization models. Section 5 concludes with discussions on potential future work.
The appendices provide technical details related to explicit piecewise representations for some
piecewise SOS-convex functions (Appendix A), infinite-dimensional conic duality (Appendix B),
and SDP representations of SOS optimization problems (Appendix C).
2 Non-Negativity of Piecewise SOS-Convex Functions
In this section, we outline SOS representations for a class of piecewise functions involving SOS-
convex polynomials. This will play a vital role in reformulating the generalized moment problem
(P) as a numerically tractable SOS optimization problem.
We begin with some preliminaries on polynomials. Denote Rm the Euclidean space of dimension
m and Rm
+ the non-negative orthant of Rm. The standard inner product is a⊤b for a, b ∈ Rm. Let
R[v] be the space of polynomials with real coefficients over v ∈ Rm. A polynomial f ∈ R[v] is
called a Sum-of-Squares (SOS) polynomial if there exist polynomials fj ∈ R[v], j = 1, . . . , J, such
that f =PJ
j=1 f 2
j . For a polynomial f ∈ R[v], we use deg f to denote its degree. We also use
Σ2
d(v) to denote the set of all SOS polynomials f of degree at most d with respect to the variable
v ∈ Rm. Next, we recall the definition of SOS-convex polynomials.
Definition 2.1 (SOS-convex polynomial [1, 16]). A polynomial f ∈ R[v] is SOS-convex if its
Hessian H(v) is an SOS matrix polynomial, that is, if there exists an (s × m) polynomial matrix
P (v) for some s ∈ N such that H(v) = P (v)⊤P (v).
Several equivalent conditions for SOS-convexity are available in [1]. For instance, f ∈ R[v] is
4

<!-- page 5 -->
SOS-convex whenever the polynomial g(v1, v2) = f(v1) − f(v2) − ∇f(v2)⊤(v1 − v2) on Rm × Rm
is an SOS polynomial with respect to the variable ( v1, v2).
Any SOS-convex polynomials are convex polynomials, but the converse is not true [1]. In other
words, the class of SOS-convex polynomials is a proper subclass of convex polynomials. More-
over, the class of SOS-convex polynomials covers affine functions, convex quadratics, and convex
separable polynomials, while SOS-convex polynomials can be non-quadratic or non-separable in
general [19].
Now, we formally define the notion of a piecewise SOS-convex function.
Definition 2.2 (Piecewise SOS-convex function). We call a function f on Rm piecewise
SOS-convex if there exist SOS-convex polynomials gk
ℓ , ℓ = 1, . . . , L, k = 1, . . . , r, on Rm such that
f(v) = min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v) for all v ∈ Rm.
The following proposition links non-negative SOS-convex polynomials to SOS polynomials which
will be deployed in the tractable representation of piecewise SOS-convex functions (Theorem 2.4).
Proposition 2.3. (see [16] and [21, Corollary 2.1]). Let f ∈ R[v] be a non-negative SOS-convex
polynomial. Then, f is an SOS polynomial.
Denote ei ∈ Rm the i-th standard basis vector of Rm, i = 1, . . . , m. The simplex in RL is defined as
∆ = {δ ∈ RL
+ :PL
ℓ=1 δℓ = 1}. Let Sν be the set of ( ν × ν) symmetric matrices. The trace product
on Sν is tr(M N) for any M, N ∈ Sν. A matrix M ∈ Sν is positive semi-definite, denoted as M ⪰ 0
(resp. positive definite, denoted as M ≻ 0), if x⊤M x ≥ 0 for all x ∈ Rν (resp. x⊤M x > 0 for all
x ∈ Rν, x ̸= 0). Let Sν
+ be the cone of symmetric ( ν × ν) positive semi-definite matrices.
A projected spectrahedron [15], or a spectrahedral shadow [29], takes the form
Ω =
(
v ∈ Rm : ∃ξ ∈ RN , F0 +
mX
i=1
viFi +
NX
t=1
ξtMt ⪰ 0
)
, (LMI)
for some Fi ∈ Sν, i = 0 , 1, . . . , m, and Mt ∈ Sν, t = 1 , . . . , N. This set will later be used as a
support set for probability measures. This set is versatile as any semi-algebraic set {v ∈ Rm :
fi(v) ≤ 0, i = 1, . . . , n} described by SOS-convex polynomials fi, i = 1, . . . , n, can be represented
as a projected spectrahedron [15]. Moreover, spectrahedra, ellipsoids, and boxes are special forms
of projected spectrahedra.
The following representation of a piecewise SOS-convex function extends the result in [19] for an
SOS-convex polynomial over a spectrahedron and in [21] for a max-SOS-convex polynomial over
a subclass of projected spectrahedron.
Theorem 2.4 (SOS representation for non-negative piecewise SOS-convex functions ).
Let Ω ⊂ Rm be given as in (LMI), gk
ℓ , ℓ = 1, . . . , L, k= 1, . . . , r, be SOS-convex polynomials on Rm.
Assume that Ω is compact and there exist ¯v ∈ Rm and ¯ξ ∈ RN such that F0+
mX
i=1
¯viFi+
NX
t=1
¯ξtMt ≻ 0.
Then, the following statements are equivalent.
(i) v ∈ Ω = ⇒ min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v) ≥ 0.
5

<!-- page 6 -->
(ii) For each k = 1, . . . , r, there exist Zk ∈ Sν
+ and δk = (δk
1 , . . . , δk
L) ∈ ∆ such that



LX
ℓ=1
δk
ℓ gk
ℓ (v) − tr(ZkF0) −
mX
i=1
vitr(ZkFi) ∈ Σ2
d(v), k = 1, . . . , r,
tr(ZkMt) = 0, t = 1, . . . , N, k = 1, . . . , r.
where d is the smallest even integer larger than max
k=1,...,r
max
ℓ=1,...,L
deg gk
ℓ .
Proof. [(ii) =⇒ (i)]. Assume that there exist Zk ∈ Sν
+ and δk ∈ ∆, k = 1, . . . , r, such that



LX
ℓ=1
δk
ℓ gk
ℓ (v) − tr(ZkF0) −
mX
i=1
vitr(ZkFi) ∈ Σ2
d(v), k = 1, . . . , r,
tr(ZkMt) = 0, t = 1, . . . , N, k = 1, . . . , r.
Fix k ∈ {1, . . . , r}. Since an SOS polynomial always takes non-negative values, it follows that
LX
ℓ=1
δk
ℓ gk
ℓ (v) − tr(ZkF0) −
mX
i=1
vitr(ZkFi) ≥ 0, for all v ∈ Rm.
Further, tr(ZkF0) +Pm
i=1 vitr(ZkFi) = tr(Zk(F0 +Pm
i=1 viFi +PN
t=1 ξtMt)) for any ξ ∈ RN. This
implies that, for all ( v, ξ) ∈ Rm × RN,
LX
ℓ=1
δk
ℓ gk
ℓ (v) ≥ tr

Zk

F0 +
mX
i=1
viFi +
NX
t=1
ξtMt

. (1)
Since ( δk
1 , . . . , δk
L) ∈ ∆, we have max
ℓ=1,...,L
gk
ℓ (v) ≥
LX
ℓ=1
δk
ℓ gk
ℓ (v). Further, for an arbitrary v ∈ Ω,
there exists ξ ∈ RN such that F0 +
mX
i=1
viFi +
NX
t=1
ξtMt ∈ Sν
+. Also, since Zk ∈ Sν
+, we have
tr

Zk

F0 +
mX
i=1
viFi +
NX
t=1
ξtMt

≥ 0. It follows from Eqn. (1) that max
ℓ=1,...,L
gk
ℓ (v) ≥ 0 for any
v ∈ Ω. This holds for any k ∈ {1, . . . , r}, so min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v) ≥ 0 for any v ∈ Ω.
[(i) =⇒ (ii)]. Assume Statement (i) is true. Fix k ∈ {1, . . . , r}. We have
min
v∈Ω
max
ℓ=1,...,L
gk
ℓ (v) = min
v∈Ω
max
δ∈∆
LX
ℓ=1
δℓgk
ℓ (v) ≥ 0. (2)
Define a bifunction Pk : Ω × ∆ → R by Pk(v, δ) = PL
ℓ=1 δℓgk
ℓ (v). Direct verification shows that
Pk(v, ·) is linear for each fixed v ∈ Ω and Pk(·, δ) is convex continuous for each fixed δ ∈ ∆. Note
that both Ω ⊆ Rm and ∆ ⊆ RL are convex compact sets. The Convex-Concave Minimax Theorem
6

<!-- page 7 -->
[cf. 35, Theorem 2.10.2] gives us that
min
v∈Ω
max
δ∈∆
Pk(v, δ) = max
δ∈∆
min
v∈Ω
Pk(v, δ).
Thus, there exists δk ∈ ∆ such that max
δ∈∆
min
v∈Ω
Pk(v, δ) = min
v∈Ω
Pk(v, δk) = min
v∈Ω
LX
ℓ=1
δk
ℓ gk
ℓ (v).
Observe that
min
v∈Ω
LX
ℓ=1
δk
ℓ gk
ℓ (v) = min
v∈Rm
( LX
ℓ=1
δk
ℓ gk
ℓ (v) : F0 +
mX
i=1
viFi +
NX
t=1
ξtMt ⪰ 0 for some ξ ∈ RN
)
= min
(v,ξ)∈Rm×RN
( LX
ℓ=1
δk
ℓ gk
ℓ (v) : F0 +
mX
i=1
viFi +
NX
t=1
ξtMt ⪰ 0
)
. (3)
By the Lagrangian Duality under the strict feasibility assumption on Ω [35, Theorem 2.9.2], there
exists Zk ∈ Sν
+ such that
min (3) = min
(v,ξ)∈Rm×RN
( LX
ℓ=1
δk
ℓ gk
ℓ (v) − tr(ZkF0) −
mX
i=1
vitr(ZkFi) −
NX
t=1
ξttr(ZkMt)
)
.
Hence, Eqn. (2) implies that there exist Zk ∈ Sν
+ and δk ∈ ∆ such that
LX
ℓ=1
δk
ℓ gk
ℓ (v) − tr(ZkF0) −
mX
i=1
vitr(ZkFi) −
NX
t=1
ξttr(ZkMt) ≥ 0, (4)
for all ( v, ξ) ∈ Rm × RN. Letting ξ = 0 ∈ RN in Eqn. (4) gives
σk(v) :=
LX
ℓ=1
δk
ℓ gk
ℓ (v) − tr(ZkF0) −
mX
i=1
vitr(ZkFi) ≥ 0, for all v ∈ Rm.
Since σk is a sum of finitely many SOS-convex polynomials, σk is a non-negative SOS-convex
polynomial. By Proposition 2.3, σk is thus an SOS polynomial. Moreover, by fixing v = bv for
somebv ∈ Rm, Eqn. (4) further implies that σk(bv) −PN
t=1 ξttr(ZkMt) ≥ 0 for any ξ ∈ RN where
σk(bv) ≥ 0 is fixed. This forces tr( ZkMt) = 0 for all t = 1 , . . . , N. The conclusion now follows
because k is chosen arbitrarily from {1, . . . , r}.
3 Piecewise SOS-Convex Moment Problems
In this section, we focus on reformulating a generalized moment optimization problem with piece-
wise SOS-convex functions as an SDP.
We begin by fixing some terminologies. Let Y and Y ′ be real topological vector spaces. They are
paired if a bilinear form ⟨·, ·⟩ : Y ′ × Y → R is defined. The cone generated by any set S ⊆ Y and
the interior of S are denoted as cone( S) and int( S) respectively. The (positive) polar of a cone
7

<!-- page 8 -->
S ⊆ Y is defined as S+ = {y′ ∈ Y ′ : ⟨y′, y⟩ ≥ 0, ∀y ∈ S}. See Appendix B for details on strong
conic duality.
For a non-empty compact set Ω ⊂ Rm, CΩ denotes the vector space of all continuous real-valued
functions on Ω equipped with the supremum norm. Denote B the Borel σ-algebra of Ω. From this
section onwards, we let X be the vector space of (finite signed regular) Borel measures on (Ω, B).
It is known that X can be identified as the topological dual space of CΩ, i.e., X = C∗
Ω. We equip
X with the weak ∗ topology. For any point v ∈ Ω, the Dirac measure 1 v which takes a mass of 1
at the point v ∈ Ω and 0 otherwise belongs to X [32].
Consider the following generalized moment problem with piecewise SOS-convex functions:
min
µ∈X
EΩ
µ
h
min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (ω)
i
(P)
s.t. EΩ
µ [hj(ω)] ≤ γj, j = 1, . . . , J,
EΩ
µ[1] = 1, µ ⪰B 0,
where γj ∈ R, j = 1 , . . . , J, and hj, j = 1 , . . . , J, and gk
ℓ , ℓ = 1 , . . . , L, k = 1 , . . . , r, are
SOS-convex polynomials on Rm. Define g : Rm → R, g(v) = min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v) the piecewise
SOS-convex function, and note that g ∈ C Ω and hj ∈ C Ω, j = 1, . . . , J. We let X ′ be the vector
space generated by g, hj, j = 1, . . . , J, and the constant function of 1. Note that X ′ ⊆ C Ω, and
X ′ is equipped with the supremum norm. Following [32], X and X ′ are paired spaces, and the
continuous bilinear form between them is given by
⟨f, µ⟩ = EΩ
µ[f(ω)] :=
Z
Ω
f(ω) dµ(ω), for all f ∈ X ′ and µ ∈ X,
where EΩ
µ[f(ω)] refers to the expectation of the random variable f(ω) with respect to the measure
µ supported on Ω. Denote the set of probability measures by PΩ = {µ ∈ X : EΩ
µ[1] = 1, µ ⪰B 0},
where µ ⪰B 0 means that µ(A) = µ(ω ∈ A) ≥ 0 for all B-measurable sets A.
The moment problem (P) is known as an uncertainty quantification problem in the distributionally
robust optimization literature [14], and the feasible set {µ ∈ PΩ : EΩ
µ[hj(ω)] ≤ γj, j = 1, . . . , J} is
called a moment ambiguity set.
We assume that the support set Ω is a (convex) projected spectrahedron:
Ω =
n
v ∈ Rm : ∃ξ ∈ RN , F 0 +
mX
i=1
viFi +
NX
t=1
ξtMt ⪰ 0
o
, (LMI)
for some Fi ∈ Sν, i = 0, 1, . . . , m, and Mt ∈ Sν, t = 1, . . . , N. We associate (P) with the following
8

<!-- page 9 -->
SOS optimization problem
max
λ∈RJ
+×R
Zk∈Sν
+,δk∈∆
k=1,...,r
−
JX
j=1
λjγj − λJ+1 (D)
s.t.
LX
ℓ=1
δk
ℓ gk
ℓ (v) +
JX
j=1
λjhj(v) + λJ+1 − tr(ZkF0) −
mX
i=1
vitr(ZkFi) ∈ Σ2
d(v), k = 1, . . . , r,
tr(ZkMt) = 0, t = 1, . . . , N, k = 1, . . . , r,
where δk = (δk
1 , . . . , δk
L) ∈ ∆ = {δ ∈ RL
+ :PL
ℓ=1 δℓ = 1}, Σ 2
d(v) is the set of SOS polynomials of
degree at most d with respect to the variable v ∈ Rm, and d is the smallest even integer with
d ≥ max

max
k=1,...,r
max
ℓ=1,...,L
deg gk
ℓ , max
j=1,...,J
deg hj

.
The SOS program (D) can be equivalently rewritten as a Semi-Definite Program (SDP), and thus
can be efficiently solved via existing SDP software. We refer the readers to Appendix C for the
procedure.
The following useful characterization states that the polar cone of the cone of non-negative mea-
sures is the convex cone of functions in X ′ that are non-negative on Ω. This characterization is
known (e.g., [4, Examples 2.37, 2.122], [17, Lemma 3.1]), but we provide the proof here for the
sake of self-containment.
Note that, for the convex cone of non-negative measures C = {µ ∈ X : µ ⪰B 0}, we denote
C+ = {f ∈ X ′ : ⟨f, µ⟩ ≥ 0, ∀µ ∈ C}. Here, C+ is the intersection of the topological polar cones
of C with the subspace X ′ of the vector space of continuous real-valued functions CΩ.
Lemma 3.1 (Characterization of cone of non-negative measures). Let Ω ̸= ∅ be any
convex compact subset of Rm, C = {µ ∈ X : µ ⪰B 0}, and f ∈ X ′. Then, f ∈ C+ if and only if
minv∈Ω f(v) ≥ 0.
Proof. Suppose that f ∈ C+. Then, ⟨f, µ⟩ ≥ 0 for all µ ∈ C. Note that the Dirac measure 1 v at
any point v ∈ Ω belongs to C. This gives f(v) = ⟨f, 1 v⟩ ≥ 0 for any v ∈ Ω. Hence, min
v∈Ω
f(v) ≥ 0.
Conversely, suppose that min
v∈Ω
f(v) ≥ 0. Now, for any µ ∈ PΩ = {µ ∈ X :
R
Ω 1dµ(ω) = 1, µ ⪰B 0},
we have
Z
Ω
f(ω) dµ(ω) ≥
Z
Ω

min
v∈Ω
f(v)

dµ(ω) =

min
v∈Ω
f(v)
Z
Ω
1dµ(ω) = min
v∈Ω
f(v).
This implies ⟨f, µ⟩ ≥ 0 for any µ ∈ PΩ. Since C = cone(PΩ), we see f ∈ C+.
Now, we show that (D) is an exact SOS reformulation for (P) in the sense that min (P) = max (D).
Assumption 3.2 (Interior point condition for moment problem).
(γ1, . . . , γJ , 1) ∈ int

{(⟨h1, µ⟩, . . . ,⟨hJ , µ⟩, ⟨1, µ⟩) : µ ∈ X, µ ⪰B 0} + RJ
+ × {0}

.
9

<!-- page 10 -->
This interior point condition is commonly used in the moment problem literature and can be
satisfied in many situations. See [32, 34] for discussions.
Assumption 3.3. The projected spectrahedronΩ as in (LMI) is compact and there exist ¯v ∈ Rm
and ¯ξ ∈ RN such that F0 +Pm
i=1 ¯viFi +PN
t=1
¯ξtMt ≻ 0. For the problem (P), the polynomials gk
ℓ ,
ℓ = 1, . . . , L, k = 1, . . . , r, and hj, j = 1, . . . , J, are SOS-convex polynomials on Rm.
Theorem 3.4 (Exact SOS program for piecewise SOS-convex moment optimization ).
Assume that (P) admits a minimizer and Assumptions 3.2-3.3 are satisfied. Then, min (P) =
max (D).
Proof. [Conic duality] . Let g(v) = min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v), b = ( γ1, . . . , γJ , 1), K = RJ
+ × {0},
and C = {µ ∈ X : µ ⪰B 0}. Note that g ∈ X ′, and C and K are convex cones that are
closed in the respective topologies. Define a continuous linear map A : X → RJ+1 by Aµ =
(⟨h1, µ⟩, . . . ,⟨hJ , µ⟩, ⟨1, µ⟩). Then,
min (P) = min
µ∈C
{⟨g, µ⟩ : −Aµ + b ∈ K}.
The adjoint mapping A∗ : RJ+1 7→ X ′ is given by A∗λ =PJ
j=1 λjhj + λJ+1 · 1, and the (positive)
polar cone of K is K+ = RJ
+ × R. Recall that C+ = {f ∈ X ′ : ⟨f, µ⟩ ≥ 0, ∀µ ∈ C}. Then,
Assumption 3.2 can be equivalently rewritten as b ∈ int(A(C) + K). It follows from strong conic
linear duality (Corollary B.2) that
min (P) = max
(
−
JX
j=1
λjγj − λJ+1 : g +
JX
j=1
λjhj + λJ+1 · 1 ∈ C+, λ ∈ K+
)
, (5)
and the maximum of the problem in Eqn. (5) is attained.
[Equivalent SOS representation of conic constraint]. This means, there exists λ ∈ RJ
+ × R
with g +
JX
j=1
λjhj + λJ+1 · 1 ∈ C+ such that min (P) = −
JX
j=1
λjγj − λJ+1. By Lemma 3.1,
g +PJ
j=1 λjhj + λJ+1 · 1 ∈ C+ if and only if
min
v∈Ω
(
min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v) +
JX
j=1
λjhj(v) + λJ+1
)
≥ 0. (6)
The polynomial gk
ℓ +PJ
j=1 λjhj + λJ+1 is SOS-convex for each ℓ = 1 , . . . , L, k = 1 , . . . , r. By
Theorem 2.4, Eqn. (6) is equivalent to, for each k = 1, . . . , r,



LX
ℓ=1
δk
ℓ gk
ℓ (v) +
JX
j=1
λjhj(v) + λJ+1 − tr(ZkF0) −
mX
i=1
vitr(ZkFi) ∈ Σ2
d(v),
tr(ZkMt) = 0, t = 1, . . . , N,
for some Zk ∈ Sν
+ and δk ∈ ∆. Therefore, min (P) ≤ max (D).
[Weak duality]. Take any feasible point λ ∈ RJ
+ × R, Zk ∈ Sν
+, δk ∈ ∆, k = 1, . . . , r, for (D).
10

<!-- page 11 -->
Theorem 2.4 gives min
k=1,...,r
max
ℓ=1,...,L
n
gk
ℓ (v) +
JX
j=1
λjhj(v) + λJ+1
o
≥ 0 for all v ∈ Ω, and thus,
min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v) ≥ −
JX
j=1
λjhj(v) − λJ+1, for all v ∈ Ω.
Take any feasible point µ for (P), then, EΩ
µ
h
min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (ω)
i
≥ EΩ
µ
h
−
JX
j=1
λjhj(ω) − λJ+1
i
.
Further, −
JX
j=1
λjγj − λJ+1 ≤ EΩ
µ
h
−
JX
j=1
λjhj(ω) − λJ+1
i
. Thus, EΩ
µ
h
min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (ω)
i
≥
−
JX
j=1
λjγj − λJ+1, implying min (P) ≥ max (D). Together, min (P) = max (D).
The SOS optimization problem (D) can be expressed equivalently as the following SDP (Proposi-
tion C.1),
max
λ∈RJ
+×R,δk∈∆
Zk∈Sν
+,Qk∈Ss(m,d/2)
+
k=1,...,r
−
JX
j=1
λjγj − λJ+1 (R)
s.t.
LX
ℓ=1
δk
ℓ (gk
ℓ )α +
JX
j=1
λj(hj)α + λJ+1(1)α − tr(ZkFα) = tr(QkBα),
for all α ∈ N , k = 1, . . . , r,
tr(ZkMt) = 0, t = 1, . . . , N, k = 1, . . . , r,
where s(m, d) =
 m+d
d

, N = {α = (α1, . . . , αm) ∈ ({0} ∪ N)m :Pm
i=1 αi ≤ d}, and Bα ∈ Ss(m,d/2),
α ∈ N , are matrices given in Appendix C. The notation (f)α of any polynomial f ∈ R[v] of degree
at most d refers to the α-th coefficient of f. In particular, (1) α refers to the α-th coefficient of
the constant function 1, that is, (1) α = 1 when α = 0, and (1) α = 0 otherwise. The matrices are
F0 := F0 ∈ Sν, Fei := Fi ∈ Sν, i = 1, . . . , m, and Fα, α ∈ N \ { 0, e1, . . . , em}, is the zero matrix.
11

<!-- page 12 -->
The dual problem of (R), which is also an SDP, is given by
min
yk∈Rs(m,d),ξk∈RN
zk∈R,k=1,...,r
rX
k=1
zk (S)
s.t.
rX
k=1
X
α∈N
yk
α(hj)α ≤ γj, j = 1, . . . , J, (7)
X
α∈N
yk
α(gk
ℓ )α ≤ zk, ℓ = 1, . . . , L, k = 1, . . . , r, (8)
yk
0 F0 +
mX
i=1
yk
eiFi +
NX
t=1
ξk
t Mt ⪰ 0, k = 1, . . . , r, (9)
rX
k=1
yk
0 = 1,
X
α∈N
yk
αBα ⪰ 0, k = 1, . . . , r, (10)
Denote the minimum value of the problem in (S) as min (S). Sufficient conditions for the strong
duality between the conic programs (R) and (S) are standard [cf. 4, 35].
We will make use of (S) to recover an optimal probability measure for (P) under suitable sign
assumptions. The relationship between these problems is depicted in Figure 1.
Moment problem (P)
Exact SOS reformulation (D) Equivalent SDP representation (R)
Conic program (S)
Theorem 3.4
Proposition C.1
Conic duality
Theorem 3.6
Figure 1: Recovering an optimal probability measure for (P) from (S).
We make use of the following generalized Jensen’s inequality [23, Theorem 2.6].
Proposition 3.5. Let f ∈ R[v] be an SOS-convex polynomial of an even degree d, y = (yα)α∈N
satisfy y0 = 1 and
X
α∈N
yαBα ⪰ 0, and Ly : R[v] → R be the linear functional Ly(f) =
X
α∈N
(f)αyα.
Then, Ly(f) ≥ f(Ly(v)), where Ly(v) = (Ly(v1), . . . , Ly(vm)), and vi denotes the polynomial that
maps a vector v ∈ Rm to its i-th coordinate vi, i = 1, . . . , m.
Theorem 3.6 (Recovering an optimal solution ). Suppose that (P) admits a minimizer, As-
sumptions 3.2-3.3 are satisfied, and strong duality between (R) and (S) holds, i.e., max (R) =
min (S). Let (¯yk, ¯ξk, ¯zk) ∈ Rs(m,d) × RN × R, ¯yk = (¯yk
α)α∈N , be a minimizer for (S) with ¯yk
0 ≥ 0
for all k = 1, . . . , r. Denote K := {k ∈ {1, . . . , r} : yk
0 > 0} and buk := 1
¯yk
0
(¯yk
e1, . . . ,¯yk
em) ∈ Rm for
k ∈ K. Suppose that, for all k /∈ K,
zk ≥ 0 and
X
α∈N
yk
α(hj)α ≥ 0 for all j = 1, . . . , J.
12

<!-- page 13 -->
Then, bµ :=
X
k∈K
¯yk
0 1 buk, the linear combination of Dirac measures at points buk ∈ Rm, k ∈ K, is a
minimizer for (P).
Proof. The equalities min (P) = max (D), max (D) = max (R), and max (R) = min (S) follow
from Theorem 3.4, Proposition C.1, and the strong duality assumption, respectively. Hence,
min (P) = min (S). Let (¯yk, ¯ξk, ¯zk) ∈ Rs(m,d) × RN × R, k = 1, . . . , r, be a solution for problem (S)
with ¯yk
0 ≥ 0 for all k = 1, . . . , r, so it satisfies Eqns. (7)-(10).
For k ∈ K, definebyk := 1
¯yk
0
¯yk ∈ Rs(m,d), which satisfies byk
0 = 1 and by Eqn. (10),
X
α∈N
byk
αBα = 1
¯yk
0
X
α∈N
¯yk
αBα ⪰ 0.
Letbuk = (buk
1, . . . ,buk
m) = (byk
e1, . . . ,byk
em), k ∈ K. Clearly, buk ∈ Ω for all k ∈ K because we can take
bξk := 1
¯yk
0
¯ξk ∈ RN such that, by Eqn. (9),
F0 +
mX
i=1
buk
i Fi +
NX
t=1
bξk
t Mt = 1
¯yk
0

¯yk
0 F0 +
mX
i=1
¯yk
eiFi +
NX
t=1
¯ξk
t Mt

⪰ 0, k ∈ K.
Definebµ :=P
k∈K ¯yk
0 1 buk. We note that for any f ∈ X ′,
EΩ
bµ[f(ω)] = ⟨f,bµ⟩ = ⟨f,
X
k∈K
¯yk
0 1 buk⟩ =
X
k∈K
¯yk
0 ⟨f, 1 buk⟩ =
X
k∈K
¯yk
0 f(buk).
By Eqn. (10), bµ is a probability measure supported on Ω since EΩ
bµ[1] =P
k∈K ¯yk
0 =Pr
k=1 ¯yk
0 = 1.
Now,buk = Lbyk(v), giving hj(buk) = hj(Lbyk(v)) for j = 1, . . . , J, k ∈ K. Using Proposition 3.5 and
the fact that hj’s are SOS-convex polynomials, we have
hj(buk) ≤ Lbyk(hj) =
X
α∈N
byk
α(hj)α = 1
¯yk
0
X
α∈N
¯yk
α(hj)α, j = 1, . . . , J, k ∈ K.
By the assumption that
X
k /∈K
X
α∈N
¯yk
α(hj)α ≥ 0 and Eqn. (7), it follows that
EΩ
bµ[hj(ω)] =
X
k∈K
¯yk
0 hj(buk) ≤
X
k∈K
X
α∈N
¯yk
α(hj)α ≤
rX
k=1
X
α∈N
¯yk
α(hj)α ≤ γj, j = 1, . . . , J,
sobµ is feasible for (P).
Similarly, one can derive that ¯ yk
0 gk
ℓ (buk) ≤
X
α∈N
¯yk
α(gk
ℓ )α, ℓ = 1 , . . . , L, k ∈ K. From Eqn. (8),
for each k ∈ K,
X
α∈N
¯yk
α(gk
ℓ )α ≤ ¯zk. This shows that ¯ yk
0 gk
ℓ (buk) ≤ ¯zk for all ℓ = 1 , . . . , L, and so,
max
ℓ=1,...,L
¯yk
0 gk
ℓ (buk) ≤ ¯zk. Let g(v) = min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (v). Recall that g ∈ X ′ and, for each k = 1, . . . , r,
13

<!-- page 14 -->
g(v) ≤ max
ℓ=1,...,L
gk
ℓ (v) for any v ∈ Ω. Then, since ¯yk
0 > 0 for all k ∈ K, one has
EΩ
bµ[g(ω)] =
X
k∈K
¯yk
0 g(buk) ≤
X
k∈K
¯yk
0

max
ℓ=1,...,L
gk
ℓ (buk)

=
X
k∈K

max
ℓ=1,...,L
¯yk
0 gk
ℓ (buk)

≤
X
k∈K
¯zk.
Using the assumption that
X
k /∈K
¯zk ≥ 0, one has EΩ
bµ[g(ω)] ≤
X
k∈K
¯zk ≤
rX
k=1
¯zk = min (S), giving
min (P) ≤ min (S).
Conversely, min (P) = min (S) ensures that min (P) = EΩ
bµ
h
min
k=1,...,r
max
ℓ=1,...,L
gk
ℓ (ω)
i
. Therefore, bµ is
an optimal solution for the original moment problem (P).
The optimal probability measure bµ :=P
k∈K ¯yk
0 1 buk refers to the discrete distribution with prob-
abilities P(ω = buk) = ¯yk
0, k ∈ K, i.e., the distribution with masses of ¯ yk
0 at points buk ∈ Ω,
k ∈ K.
We present a recovery result under a stronger condition that is easy to verify.
Corollary 3.7 (Recovery under a simple stronger condition ). Suppose that (P) admits a
minimizer, Assumptions 3.2-3.3 are satisfied, and max (R) = min (S). Let (¯yk, ¯ξk, ¯zk) ∈ Rs(m,d) ×
RN × R, ¯yk = (¯yk
α)α∈N , be a minimizer for (S) with ¯yk
0 > 0 for all k = 1 , . . . , r. Denote buk :=
1
¯yk
0
(¯yk
e1, . . . ,¯yk
em) ∈ Rm, k = 1, . . . , r. Then, bµ :=Pr
k=1 ¯yk
0 1 buk is a minimizer for (P).
Proof. This follows from Theorem 3.6 by setting K = {1, . . . , r}.
As an application to polynomial optimization, our recovery result reduces to the known Lasserre’s
result in [23]. More explicitly, consider
min
v∈Ω
g(v), (11)
where Ω ⊂ Rm is given as in (LMI) and g is an SOS-convex polynomial. Notice that, in [23,
Theorem 3.3], the feasible set Ω is a compact semi-algebraic set formed by SOS-convex polynomial
inequalities, and thus can be represented as a projected spectrahedron [15].
We associate (11) with the following conic programs
max
λ∈R,Z∈Sν
+
n
− λ : g(v) + λ − tr(ZF0) −
mX
i=1
vitr(ZF i) ∈ Σ2
d(v), tr(ZM t) = 0, t = 1, . . . , N
o
, (12)
and its dual problem
min
y∈Rs(m,d),ξ∈RN
nX
α∈N
yα(g)α : F0 +
mX
i=1
yeiFi +
NX
t=1
ξtMt ⪰ 0, y 0 = 1,
X
α∈N
yαBα ⪰ 0
o
. (13)
The following Corollary shows how an optimal solution of (11) can be obtained by solving the
SDP program (13).
14

<!-- page 15 -->
Corollary 3.8. For the problem (11), suppose that Ω is compact, there exist ¯v ∈ Rm and ¯ξ ∈ RN
such that F0 +Pm
i=1 ¯viFi +PN
t=1
¯ξtMt ≻ 0, and max (12) = min (13). Let (¯y, ¯ξ) ∈ Rs(m,d) × RN be
a minimizer for (13). Then, (¯ye1, . . . ,¯yem) ∈ Rm is a minimizer for (11).
Proof. Consider the moment problems
min
µ∈PΩ
EΩ
µ[g(ω)] (14)
and
min
µ∈D
EΩ
µ[g(ω)], (15)
where D = {1 v ∈ X : v ∈ Ω} ⊂ P Ω is the set of all Dirac measures supported on Ω, giving
min (14) ≤ min (15). Note that min (15) = min (11) because EΩ
1 v[g(ω)] = g(v) for all v ∈ Ω.
The interior point condition (Assumption 3.2) for the moment problem (14) requires that 1 ∈
int{EΩ
µ[1] : µ ∈ C}, where C = {µ ∈ X : µ ⪰B 0}. Since C = cone(PΩ) and EΩ
P [1] = 1 for any
P ∈ P Ω, a direct verification shows that {EΩ
µ[1] : µ ∈ C} = R+. So, Assumption 3.2 is satisfied
automatically.
By Corollary 3.7 with ¯ y0 = 1, we have that bµ := 1 (¯ye1 ,...,¯yem) is a minimizer for (14) and
(¯ye1, . . . ,¯yem) ∈ Ω. As EΩ
bµ[g(ω)] = min (14) ≤ min (15) and bµ ∈ D, so bµ is also a minimizer
for (15). Noting that EΩ
bµ[g(ω)] = g(¯ye1, . . . ,¯yem), one conclude that (¯ye1, . . . ,¯yem) is a minimizer
for the polynomial optimization problem (11).
4 Applications to Newsvendor and Revenue Maximization
We devote this section to presenting SDP reformulations and numerical results for the Newsvendor
and revenue maximization problems.
Generalized Newsvendor Problems . A company orders n goods, where n ≥ 1. For each of
the i-th goods, the unit cost of upfront ordering is $ ci with 0 < c i < 1. For a fixed order quantity
0 ≤ xi ≤ Ri, where Ri refers to the capacity, the company wishes to estimate a worst-case (upper)
bound for the cost of ordering to meet the random demands. If the demand ωi exceeds the
upfront order xi, the company incurs a unit back-ordering cost of $1, or a total back-ordering cost
of $( ωi − xi), otherwise, the cost of back-ordering is $0. This results in a total ordering cost of
cixi + max{ωi − xi, 0}.
Assume that the demand for each goods is independent and the support is Qn
i=1 Ωi, where Ω i =
[ωi, ωi] with 0 ≤ ωi < ωi, i = 1 , . . . , n. Then, the multi-product Newsvendor model reduces
to single-product Newsvendor models, and the worst-case ordering cost can be found by solving
max
µ∈Pi
{cixi + EΩi
µ [max{ωi − xi, 0}]} for all i = 1 , . . . , n, where Pi consists of various moment
constraints.
We focus on the single-product Newsvendor model and suppress the index i ∈ { 1, . . . , n}. For
a fixed order quantity 0 ≤ x ≤ R, the worst-case ordering cost can be found from the following
15

<!-- page 16 -->
moment problem,
max
µ∈P
cx + EΩ
µ[max{ω − x, 0}], (16)
where Ω = [ω, ω] = {v ∈ R : F0 + vF1 ⪰ 0}, ω < ω, with F0 =

−ω 0
0 ω

and F1 =

1 0
0 −1

.
Another interpretation of (16) at c = 0 comes from pricing a European call option [3, 22]. The
random variable ω represents the price of the underlying stock, and the fixed value x is regarded
as the strike price. This gives the expected price of EΩ
µ[max{ω − x, 0}], and thus a sharp upper
bound max
µ∈P
EΩ
µ[max{ω − x, 0}] of the price.
Numerically tractable reformulations for (16) are available for a limited choice ofP. For example,
[27] provides a closed-form formula when P specifies known values for the mean and variance;
recently, [12] presents a closed-form formula when P specifies known values for mean plus one
other moment of any order; [3] reformulates (16) as an SDP when P contains the first J mo-
ments; and [33] reformulates (16) as a conic program when the random demand variable follows a
distribution defined by robust statistics. For other related models with tractable reformulations,
see [8].
We study the problem (16) when P1 = {µ ∈ PΩ : EΩ
µ[ω] ≤ γ1, EΩ
µ[ω2] ≤ γ2} and P2 = P1 ∩ {µ ∈
PΩ : EΩ
µ[ω4] ≤ γ3}, where PΩ is the set of probability measures on Ω. The set P1 contains upper
bounds for the mean and the second-order moment which is related to the variance, while P2
contains an additional bound for the fourth-order moment which is related to the kurtosis of a
distribution. The mean, variance, and kurtosis of the random demand distribution are usually
accessible in real-world applications.
The standard minimization form for the Newsvendor problem is
min (17) := min
µ∈P
EΩ
µ
h
min
k=1,2
gk
1(ω)
i
, g 1
1(v) = x − v, g 2
1(v) = 0, (17)
and the worst-case expected ordering cost is equal to [ cx − min (17)].
The associated SOS optimization problem for (17) when P = P1 is
max (18) := max
λ∈R2
+×R
Z,Z ′∈S2
+
− λ1γ1 − λ2γ2 − λ3 (18)
s.t. [x + λ3 − tr(ZF0)] + [−1 + λ1 − tr(ZF1)]v + λ2v2 ∈ Σ2
2(v),
[λ3 − tr(Z ′F0)] + [λ1 − tr(Z ′F1)]v + λ2v2 ∈ Σ2
2(v).
Its equivalent SDP representation is
max (19) := max
λ∈R2
+×R
Z,Z ′∈S2
+
− λ1γ1 − λ2γ2 − λ3 (19)
s.t.

x + λ3 − tr(ZF0) 1
2[−1 + λ1 − tr(ZF1)]
1
2[−1 + λ1 − tr(ZF1)] λ2

⪰ 0,

λ3 − tr(Z ′F0) 1
2[λ1 − tr(Z ′F1)]
1
2[λ1 − tr(Z ′F1)] λ2

⪰ 0.
16

<!-- page 17 -->
Proposition 4.1. Assume that (17) admits a minimizer for P = P1 and (γ1, γ2, 1) ∈
int{{(EΩ
µ[ω], EΩ
µ[ω2], EΩ
µ[1]) : µ ∈ X, µ ⪰B 0} + R2
+ × {0}}. Then, min (17) = max (19).
Proof. By Theorem 3.4, min (17) = max (18). Next, following the approach in Proposition C.1,
the constraint [x+λ3 −tr(ZF0)]+ v[−1+ λ1 −tr(ZF1)]+ λ2v2 ∈ Σ2
2(v) is equivalent to the existence
of Q ∈ S2
+ such that tr(QB0) = x+ λ3 −tr(ZF0), tr(QB1) = −1+ λ1 −tr(ZF1), and tr(QB2) = λ2,
where B0 =

1 0
0 0

, B1 =

0 1
1 0

, and B2 =

0 0
0 1

. This gives
Q =

x + λ3 − tr(ZF0) 1
2[−1 + λ1 − tr(ZF1)]
1
2[−1 + λ1 − tr(ZF1)] λ2

∈ S2
+.
Similar arguments apply for the other constraint, and hence max (18) = max (19).
Now, we illustrate how one can recover an optimal measure for (17) via the following SDP
min (20) := min
y1,y2∈R3
z1,z2∈R
z1 + z2 (20)
s.t. y 1
1 + y2
1 ≤ γ1, y 1
2 + y2
2 ≤ γ2, y 1
0 + y2
0 = 1,
xy1
0 − y1
1 ≤ z1, 0 ≤ z2,
yk
0 F0 + yk
1 F1 ⪰ 0,

yk
0 yk
1
yk
1 yk
2

⪰ 0, k = 1, 2.
Proposition 4.2. Assume the same conditions as in Proposition 4.1 and max (19) = min (20).
Let (¯yk, ¯zk) ∈ R3 × R with ¯yk = (¯yk
0 , ¯yk
1 , ¯yk
2) be a minimizer for (20) with ¯yk
0 > 0, k = 1, 2. Denote
buk := 1
¯yk
0
¯yk
1 ∈ R, k = 1, 2. Then, ¯y1
01 bu1 + ¯y2
01 bu2 is a minimizer for (17) for P = P1.
Proof. This follows from Corollary 3.7.
Note that the probability measure optimal for the standard minimization form (17) will also be
optimal for the original Newsvendor maximization problem (16).
In a similar vein, we examine the Newsvendor model with P = P2. Its associated SOS optimiza-
tion problem is given by
max (21) := max
λ∈R3
+×R
Z,Z ′∈S2
+
− λ1γ1 − λ2γ2 − λ3γ3 − λ4 (21)
s.t. [x + λ4 − tr(ZF0)] + [−1 + λ1 − tr(ZF1)]v + λ2v2 + λ3v4 ∈ Σ2
4(v),
[λ4 − tr(Z ′F0)] + [λ1 − tr(Z ′F1)]v + λ2v2 + λ3v4 ∈ Σ2
4(v).
17

<!-- page 18 -->
Its equivalent SDP representation is
max (22) := max
λ∈R3
+×R
Z,Z ′∈S2
+,Q,Q′∈S3
+
− λ1γ1 − λ2γ2 − λ3γ3 − λ4 (22)
s.t. tr(QB0) = x + λ4 − tr(ZF0), tr(Q′B0) = λ4 − tr(Z ′F0),
tr(QB1) = −1 + λ1 − tr(ZF1), tr(Q′B1) = λ1 − tr(Z ′F1),
tr(QB2) = tr(Q′B2) = λ2,
tr(QB3) = tr(Q′B3) = 0,
tr(QB4) = tr(Q′B4) = λ3,
where
B0 =


1 0 0
0 0 0
0 0 0

 , B1 =


0 1 0
1 0 0
0 0 0

 , B2 =


0 0 1
0 1 0
1 0 0

 , B3 =


0 0 0
0 0 1
0 1 0

 , B4 =


0 0 0
0 0 0
0 0 1

 .
(23)
Proposition 4.3. Assume that (17) admits a minimizer for P = P2 and (γ1, γ2, γ3, 1) ∈
int{{(EΩ
µ[ω], EΩ
µ[ω2], EΩ
µ[ω4], EΩ
µ[1]) : µ ∈ X, µ ⪰B 0} + R3
+ × {0}}. Then, min (17) = max (22).
Proof. The proof is similar to Proposition 4.1.
We recover an optimal solution for (17) using the following SDP,
min (24) := min
y1,y2∈R5
z1,z2∈R
z1 + z2 (24)
s.t. y 1
1 + y2
1 ≤ γ1, y 1
2 + y2
2 ≤ γ2, y 1
4 + y2
4 ≤ γ3, y 1
0 + y2
0 = 1,
xy1
0 − y1
1 ≤ z1, 0 ≤ z2,
yk
0 F0 + yk
1 F1 ⪰ 0,


yk
0 yk
1 yk
2
yk
1 yk
2 yk
3
yk
2 yk
3 yk
4

 ⪰ 0, k = 1, 2.
Proposition 4.4. Assume the same conditions as in Proposition 4.3 and max (22) = min (24).
Let (¯yk, ¯zk) ∈ R5 × R with ¯yk = (¯yk
0 , ¯yk
1 , ¯yk
2 , ¯yk
3 , ¯yk
4) be a minimizer for (24) with ¯yk
0 > 0, k = 1, 2.
Denotebuk := 1
¯yk
0
¯yk
1 ∈ R, k = 1, 2. Then, ¯y1
01 bu1 + ¯y2
01 bu2 is a minimizer for (17) for P = P2.
Proof. This follows from Corollary 3.7.
We illustrate how the worst-case costs of ordering [ cx − max (19)] and [cx − max (22)] change by
varying the order quantity x ∈ [0, R] = [0, 10] for P1 and P2. The SDPs are modelled by CVX
MATLAB [11] and solved by Mosek [28]. Set Ω = [ ω, ω] = [0, 100] for the support, c = 0.1 for the
unit cost of upfront ordering, and γ1 = γ2 = γ3 = 1 for the moment bounds. Figure 2 illustrates
the cost of ordering ( $) for P = P1 (blue) and P = P2 (orange).
18

<!-- page 19 -->
Figure 2: Costs of ordering ( $) for P = P1 (blue) and P = P2 (orange).
The curves in Figure 2 exhibit “U” shapes. Before the turning points ( x ∈ [0, 1.5811] blue,
x ∈ [0, 1.3337] orange), the more upfront orders at cost c = $0.1 the company makes, the more it
saves from back-ordering at cost $1. For larger upfront ordering quantities ( x ∈ [1.5811, 10] blue,
x ∈ [1.3337, 10] orange) that can largely cover the demand, back-ordering is not needed. The
slope of the curves in these regions is approximately c = 0.1.
Notice that P2 is a subset of P1, Figure 2 confirms that the costs of ordering for P1 (blue) are
at least the costs for P2 (orange), and this leads to a higher minimum cost for P1. This could be
justified that, by knowing more information about the demand, i.e., P = P2, the future demand
can be better predicted, and thus the ordering strategy could be better determined.
The minimum costs of ordering are achieved at the turning points of the curves. Specifically, the
minimum costs are $0.3162 and $0.1778 if the company orders x = 1.5811 and x = 1.3337 units
of goods under the ambiguity constraints P = P1 (blue) and P = P2 (orange), respectively.
For P1 with order quantity x = 1.5811, an optimal probability measure is 0 .11 3.1623 + 0.91 0.0003.
It refers to the discrete distribution with P(ω = 3.1623) = 0.1 and P(ω = 0.0003) = 0.9. For P2
with x = 1.3337, an optimal probability measure is 0 .11 1.7782 + 0.91 0.0200.
Revenue Maximization. A merchant supplies a random quantity ω ∈ [0, R] ⊆ R of goods and
sells them to n customers. Each customer offers a different price based on the supply quantity.
The goods can be sold to one customer exclusively at a time, and the merchant wishes to maximize
revenue by selling the goods to the customer who offers the highest price.
This problem of revenue maximization can be formulated as
max
µ∈PΩ
n
EΩ
µ
h
max
k=1,...,n
hk(ω)
i
: EΩ
µ[ω] ≤ γ1, EΩ
µ[ω2] ≤ γ2
o
, (25)
where hk is the offer price from the k-th customer, k = 1 , . . . , n, and Ω = [0 , R], R > 0. The
support can be expressed equivalently as a spectrahedron Ω = {v ∈ R : F0 + vF1 ⪰ 0} with
F0 =

0 0
0 R

and F1 =

1 0
0 −1

.
The standard minimization form is
min (26) := min
µ∈PΩ
n
EΩ
µ[g(ω)] : EΩ
µ[ω] ≤ γ1, EΩ
µ[ω2] ≤ γ2
o
, (26)
where g is given byg(v) = − max
k=1,...,n
hk(v). The maximum expected revenue is equal to [− min (26)].
19

<!-- page 20 -->
Denoting fk(v) = −hk(v), one has g(v) = min
k=1,...,n
fk(v). Here, for each k = 1, . . . , n, we describe
the offer price by
fk(v) =
(
αk(v − bk)2 + βk(v − bk)4 + ck, if 0 ≤ v ≤ bk,
ck, if v > b k, (27)
where αk, βk, bk ≥ 0 and ck ≤ −αkb2
k − βkb4
k. Note that g is a piecewise SOS-convex function with
the following representation:
g(v) = min
k=1,...,2n
max
ℓ=1,2
gk
ℓ (v), (28)
where gk
1(v) = gk
2(v) = αk(v −bk)2+βk(v −bk)4+ck, gn+k
1 (v) = −(αkbk +βkb3
k)v+(αkb2
k +βkb4
k +ck),
and gn+k
2 (v) = ck, for k = 1, . . . , n. For details of this representation, see Appendix A.
In fact, (−fk) where fk is given as in Eqn. (27), or hk in problem (25), corresponds to a combined
quadratic-quartic utility function that is concave, non-decreasing, and smooth [9]. The coefficients
(−αk) and ( −βk) capture the rate at which the customer increases the offer price based on the
supply; bk is the maximum quantity of goods the k-th customer is willing to purchase, beyond
which the customer will no longer be willing to increase the price; and ck which excludes negative
offer prices is the maximum price.
The function in Eqn. (27) reduces to the quadratic model in [13] by setting βk = 0, k = 1, . . . , n,
where an approximating upper bound for the maximum revenue is calculated through a convex
program. Our model covers more diverse purchasing behaviours described by the quartic feature.
Suppose that there are three customers with parameters α1 = 1, α2 = 1, α3 = 1
10, β1 = 1, β2 = 1
16,
β3 = 1
100, b1 = 1, b2 = 2, b3 = 4, and c1 = −5, c2 = −7, c3 = −7.5. The supply quantity ranges
between [0 , R] = [0 , 4]. As shown in Figure 3, the piecewise SOS-convex function g is neither
convex nor concave and is not smooth.
Figure 3: Left: negative offer prices (solid lines) and negative revenue (dotted lines). Right:
piecewise SOS-convex function in Eqn. (28).
20

<!-- page 21 -->
Associate (26) with the following SOS optimization problem with n = 3 customers,
max
λ∈R2
+×R,δk∈[0,1]
Zk,Z′
k∈S2
+,k=1,...,n
− λ1γ1 − λ2γ2 − λ3
s.t. [αkb2
k + βkb4
k + ck + λ3 − tr(ZkF0)] + [−2αkbk − 4βkb3
k + λ1 − tr(ZkF1)]v+
[αk + 6βkb2
k + λ2]v2 − 4βkbkv3 + βkv4 ∈ Σ2
4(v), k = 1, . . . , n,
[δkαkb2
k + δkβkb4
k + ck + λ3 − tr(Z ′
kF0)] + [−δkαkbk − δkβkb3
k + λ1 − tr(Z ′
kF1)]v+
λ2v2 ∈ Σ2
4(v), k = 1, . . . , n,
and its equivalent SDP representation
max (29) := max
λ∈R2
+×R,δk∈[0,1]
Zk,Z′
k∈S2
+,Qk,Q′
k∈S3
+,k=1,...,n
− λ1γ1 − λ2γ2 − λ3 (29)
s.t. for each k = 1, . . . , n,
tr(QkB0) = αkb2
k + βkb4
k + ck + λ3 − tr(ZkF0),
tr(QkB1) = −2αkbk − 4βkb3
k + λ1 − tr(ZkF1),
tr(QkB2) = αk + 6βkb2
k + λ2,
tr(QkB3) = −4βkbk, tr(QkB4) = βk,
tr(Q′
kB0) = δkαkb2
k + δkβkb4
k + ck + λ3 − tr(Z ′
kF0),
tr(Q′
kB1) = −δkαkbk − δkβkb3
k + λ1 − tr(Z ′
kF1),
tr(Q′
kB2) = λ2, tr(Q′
kB3) = tr(Q′
kB4) = 0,
where B0, . . . , B4 are given in Eqn. (23).
Proposition 4.5. Assume that (26) admits a minimizer and (γ1, γ2, 1) ∈ int{{(EΩ
µ[ω], EΩ
µ[ω2], EΩ
µ[1]) :
µ ∈ X, µ ⪰B 0} + R2
+ × {0}}. Then, min (26) = max (29).
Proof. The proof is similar to Proposition 4.1.
Similar to the Newsvendor application, higher-order moments can be incorporated, and the re-
sulting SDP reformulation can be obtained.
The maximum revenue is equal to [ − max (29)]. When the first and second-order moments of the
supply quantity are at most 2, i.e., γ1 = γ2 = 2, the maximum expected revenue is $6.6495.
21

<!-- page 22 -->
We recover an optimal measure via the SDP below,
min
yk∈R5,zk∈R
k=1,...,2n
2nX
k=1
zk (30)
s.t.
2nX
k=1
yk
1 ≤ γ1,
2nX
k=1
yk
2 ≤ γ2,
2nX
k=1
yk
0 = 1,
[αkb2
k + βkb4
k + ck]yk
0 − [2αkbk + 4βkb3
k]yk
1 + [αk + 6βkb2
k]yk
2 −
4βkbkyk
3 + βkyk
4 ≤ zk, k = 1, . . . , n,
[αkb2
k + βkb4
k + ck]yn+k
0 − [αkbk + βkb3
k]yn+k
1 ≤ zn+k, k = 1, . . . , n,
ckyn+k
0 ≤ zn+k, k = 1, . . . , n,
yk
0 F0 + yk
1 F1 ⪰ 0,


yk
0 yk
1 yk
2
yk
1 yk
2 yk
3
yk
2 yk
3 yk
4

 ⪰ 0, k = 1, . . . ,2n.
Solving (30) for n = 3 gives ¯y1 = (¯y1
0, ¯y1
1, ¯y1
2, ¯y1
3, ¯y1
4) = (0 , 0, 0, 0, 0), ¯y2 = (1 , 1.4142, 2, 2.8284, 4),
¯y3 = (0, 0, 0, 0, 0), ¯y4 = (0, 0, 0, 0, 1.5563), ¯y5 = (0, 0, 0, 0, 1.5871), ¯y6 = (0, 0, 0, 0, 1.4684), ¯zk = 0,
k = 1 , 3, 4, 5, 6, and ¯z2 = −6.6495. Clearly, ¯yk
0 ≥ 0 and K = {k ∈ { 1, . . . ,6} : ¯yk
0 > 0} = {2}.
Moreover, ¯zk ≥ 0 for all k /∈ K. For h1(v) = v, it is satisfied that P4
α=0 ¯yk
α(h1)α = ¯yk
1 = 0 for
all k /∈ K, and for h2(v) = v2, P4
α=0 ¯yk
α(h2)α = ¯yk
2 = 0 for all k /∈ K. By Theorem 3.6, an
optimal probability measure for problem (26), and thus for (25) is 1 1.4142, which is the discrete
distribution P(ω = 1 .4142) = 1. This result can be interpreted by the fact that the maximum
expected revenue $6.6495 can be achieved when the merchant supplies 1 .4142 units of goods.
We want to remark that the assumption ¯yk
0 > 0, k = 1, . . . , r, in Corollary 3.7 is sufficient but not
necessary for an optimal solution to be Pr
k=1 ¯yk
0 1 buk. Indeed, with r = 2n = 6, the above revenue
problem offers an example that ¯yk
0 > 0, k = 1, . . . , r, is not a necessary condition.
Figure 4 (left) shows how the maximum expected revenue changes with respect to the mean bound
γ1. The bound for the second-order moment is set to be γ2 = γ2
1. The curve exhibits an increasing
trend. For a large mean bound γ1, the company can target the customer offering a higher price,
and thus the revenue is higher. Figure 4 (right), on the contrary, fixesγ1 = 2 and increases γ2. The
curve exhibits an upward shape, which can be justified by the fact that the larger the second-order
moment (and thus the variance), the wider the customers the company can sell products to, and
thus the higher the revenue.
22

<!-- page 23 -->
Figure 4: Maximum expected revenue by varying γ1 (left) and γ2 (right).
5 Conclusion and Outlook
We showed how to derive Sum-of-Squares (SOS) reformulations for important classes of infinite-
dimensional moment optimization problems involving piecewise SOS-convex functions. We also
showed how to recover an optimal probability measure from the associated Semi-Definite Program
(SDP) reformulation.
It is worth noting that the class of piecewise SOS-convex functions as given in Definition 2.2
facilitates tractable representations of the corresponding moment problem (Theorem 3.4) and it
is rich enough to cover a broad class of functions encountered in optimization. Our definition 2.2
of piecewise functions may be extended to more general settings where one piece of SOS-convex
polynomial is given at one partition of Rm, in line with the piecewise linear-quadratic function
studied in [31, Definition 10.20]. It would be intriguing to explore the applicability and tractability
of these generalized classes of piecewise functions.
Further, our method in Theorem 2.4 suggests that the numerically tractable SOS representation
for non-negative piecewise SOS-convex functions over a projected spectrahedron can be extended
to more general settings. For instance, it would be of interest to examine similar representations
for piecewise SOS-convex functions over any non-convex sets whose convex hulls are semi-definite-
representable, leading to SDP reformulations or relaxations for broad classes of optimization prob-
lems.
As applications, we presented numerical results for a class of generalized Newsvendor and revenue
maximization problems with higher-order moments by solving their equivalent SDP reformula-
tions. Our approach opens new avenues for further research, such as conic program reformula-
tions for distributionally robust optimization problems [33, 36] involving piecewise SOS-convex
functions, and applications to practical models such as the lot-sizing and product management
problems in the face of uncertain conditions [2]. These problems will be examined in our forth-
coming studies.
Acknowledgment
The authors are grateful to the referees for their helpful comments and valuable suggestions which
have contributed to the final preparation of the paper.
23

<!-- page 24 -->
Declarations
Data availability statement. No data was used for the research described in the article.
Conflict of interest statement . The authors have no conflict of interest to declare that are
relevant to the content of this article.
Appendix A Piecewise SOS-Convex Representation
In this appendix, we present explicit piecewise SOS-convex representations for the examples men-
tioned in the Introduction (Section 1). We note that the representation is, in general, not unique.
• The truncated ℓ1-norm, pε(v) = min{1, ε|v|}, ε > 0, can be expressed aspε(v) = min
k=1,2
max
ℓ=1,2
gk
ℓ (v)
where g1
1(v) = g1
2(v) = 1, g2
1(v) = εv, g2
2(v) = −εv. Note that the truncated ℓ1-norm is non-
convex and non-smooth. We plot max
ℓ=1,2
g1
ℓ (v) versus max
ℓ=1,2
g2
ℓ (v) in Figure 5 (left) and the
truncated ℓ1-norm in Figure 5 (right).
Figure 5: Piecewise representation (left) and truncated ℓ1-norm (right).
• The piecewise quadratic function, a, b ≥ 0, c ∈ R,
f(v) =
(
a(v − b)2 + c, if 0 ≤ v ≤ b,
c, if v > b,
can be written as f(v) = min
k=1,2
max
ℓ=1,2
gk
ℓ (v) where g1
1(v) = g1
2(v) = a(v − b)2 + c, g2
1(v) =
−abv + ab2 + c, and g2
2(v) = c. We plot max
ℓ=1,2
g1
ℓ (v) versus max
ℓ=1,2
g2
ℓ (v) in Figure 6 (left) and
the full piecewise function in Figure 6 (right).
24

<!-- page 25 -->
Figure 6: Piecewise representation (left) and the piecewise quadratic function (right).
• The piecewise quadratic-quartic function from Section 4,
f(v) =
(
α(v − b)2 + β(v − b)4 + c, if 0 ≤ v ≤ b,
c, if v > b,
with α, β, b ≥ 0, c ∈ R, can be expressed as f(v) = min
k=1,2
max
ℓ=1,2
gk
ℓ (v), where g1
1(v) = g1
2(v) =
α(v − b)2 + β(v − b)4 + c, g2
1(v) = −(αb + βb3)v + (αb2 + βb4 + c), and g2
2(v) = c. See Figure
7 for illustrations.
Figure 7: Piecewise representation (left) and the piecewise quadratic-quartic function (right).
• The Huber loss function [33] with parameter ε > 0,
Hε(v) =
(
1
2 v2, if |v| ≤ ε,
ε|v| − 1
2 ε2, otherwise,
can be expressed equivalently as Hε(v) = min
k=1,2
max
ℓ=1,2,3,4
gk
ℓ (v), where g1
ℓ (v) = 1
2 v2, ℓ = 1, . . . ,4,
g2
1(v) = εv − 1
2 ε2, g2
2(v) = −εv − 1
2 ε2, g2
3(v) = 1
2 εv, and g2
4(v) = − 1
2 εv. See Figure 8 for
illustrations.
25

<!-- page 26 -->
Figure 8: Piecewise representation (left) and Huber loss function (right).
Appendix B Conic Duality in Topological Vector Spaces
This section presents duality results for infinite-dimensional conic linear programs that were used
in Section 3.
Let X and X ′ be real topological vector spaces. They are paired if a bilinear form ⟨·, ·⟩ : X ′ ×X →
R is defined. For the linear map A : X → RJ, assume that for any λ ∈ RJ, there is a unique
x′ ∈ X ′ satisfying λ⊤Ax = ⟨x′, x⟩ for all x ∈ X. See [32, Assumption A1] for a discussion of this
assumption. For the continuous linear map A : X → RJ, the adjoint mapping A∗ : RJ → X ′ is
defined as ⟨A∗λ, x⟩ = λ⊤Ax for any x ∈ X, λ ∈ RJ. See [4, 35] for more details on the convexity
of sets and functions.
We consider the problem
min
x∈X
⟨v, x⟩ (CP)
s.t. x ∈ C, −Ax + b ∈ K,
where C ⊂ X and K ⊂ RJ are convex cones that are closed in the respective topologies, A : X →
RJ is a linear map, and v ∈ X ′, b ∈ RJ.
Associate (CP) with the set D = {x ∈ X : −Ax + b ∈ K} and the dual problem
max
λ∈RJ
− λ⊤b (CD)
s.t. v + A∗λ ∈ C+, λ ∈ K+,
where C+ = {x′ ∈ X ′ : ⟨x′, x⟩ ≥ 0, ∀x ∈ C} and K+ = {w ∈ RJ : w⊤x ≥ 0, ∀x ∈ K}.
We present a strong duality theorem under a more general constraint qualification known as the
Generalized-Sharpened strong Conical Hull Intersection Property (G-S strong CHIP [5]). A pair
of sets {C, D} is said to satisfy the G-S strong CHIP at x∗ ∈ (C ∩ D) whenever (C ∩ D − x∗)+ =
(C − x∗)+ + ∪λ∈K+{−A∗λ : λ⊤(−Ax∗ + b) = 0}.
Theorem B.1 (Strong duality under G-S strong CHIP). Suppose that x∗ is a minimizer
for the problem (CP), the linear form ⟨u, ·⟩ is continuous at x∗ for each u ∈ X ′, the linear map
A : X → RJ is continuous, the convex cones C and K are closed in the respective topologies, and
that the pair {C, D} satisfies the G-S strong CHIP at x∗. Then, min (CP) = max (CD), and the
26

<!-- page 27 -->
maximum of (CD) is attained.
Proof. Firstly, we note that, by construction, weak duality holds, that is, min (CP) ≥ max (CD).
To see min (CP) ≤ max (CD), let x∗ be a minimizer of (CP). By the necessary optimality
conditions [35], 0 ∈ ∂(⟨v, ·⟩)(x∗) − (C ∩ D − x∗)+ = {v} − (C ∩ D − x∗)+. This means v ∈
(C ∩D−x∗)+, and thus v ∈ (C −x∗)++∪λ∈K+{−A∗λ : λ⊤(−Ax∗+b) = 0} by the G-S strong CHIP
assumption. Therefore, there exist u ∈ (C − x∗)+ and ¯λ ∈ K+ satisfying u = v + A∗¯λ ∈ (C − x∗)+,
¯λ⊤(−Ax∗ + b) = 0, and
⟨v + A∗¯λ, x − x∗⟩ ≥ 0, for all x ∈ C. (31)
Since C is a convex cone, we have ⟨v + A∗¯λ, x⟩ ≥ 0 for all x ∈ C, and so v + A∗¯λ ∈ C+. In
particular, ¯λ is feasible for (CD). Moreover, letting x = 0 ∈ C in Eqn. (31) gives ⟨v + A∗¯λ, x∗⟩ =
⟨v, x∗⟩ + ¯λ⊤Ax∗ ≤ 0. But ¯λ⊤Ax∗ = ¯λ⊤b, which implies ⟨v, x∗⟩ + ¯λ⊤b ≤ 0. Hence, −¯λ⊤b ≥ ⟨v, x∗⟩,
and so,
max (CD) = max{−λ⊤b : λ ∈ K+, v + A∗λ ∈ C+} ≥ − ¯λ⊤b ≥ ⟨v, x∗⟩ = min (CP).
Therefore, we see that min (CP) = max (CD) and the maximum of (CD) is attained at ¯λ.
The G-S strong CHIP is the weakest constraint qualification guaranteeing strong (Lagrangian)
duality for convex optimization. For details, see [5] and other references therein. As we see below,
strong duality under the interior point condition [32] is a consequence of Theorem B.1.
Corollary B.2 (Strong duality under interior point condition). Suppose that x∗ is a min-
imizer for the problem (CP), the linear form ⟨u, ·⟩ is continuous at x∗ for each u ∈ X ′, the linear
map A : X → RJ is continuous, the convex cones C and K are closed in the respective topologies.
If b ∈ int(A(C) + K), then, min (CP) = max (CD), and the maximum of (CD) is attained.
Proof. The conclusion will follow from Theorem B.1 if we show that b ∈ int(A(C) + K) implies
(C ∩ D − x∗)+ = (C − x∗)+ + ∪λ∈K+{−A∗λ : λ⊤(−Ax∗ + b) = 0}.
Note that the inclusion ( C − x∗)+ + ∪λ∈K+{−A∗λ : λ⊤(−Ax∗ + b) = 0 } ⊆ (C ∩ D − x∗)+ holds
by construction. Conversely, let u ∈ (C ∩ D − x∗)+. By definition, this means ⟨u, x − x∗⟩ ≥ 0 for
all x ∈ (C ∩ D). Hence, ⟨u, x⟩ ≥ ⟨ u, x∗⟩ for all x ∈ (C ∩ D), and x∗ is a minimizer of the convex
optimization problem min {⟨u, x⟩ : x ∈ C, −Ax + b ∈ K}. By [32], there exists λ ∈ K+ with
⟨u, x∗⟩ = −λ⊤b and A∗λ + u ∈ C+. Since −Ax∗ + b ∈ K and λ ∈ K+, we have λ⊤(−Ax∗ + b) ≥ 0.
On the other hand, λ⊤(−Ax∗ + b) = ⟨−A∗λ − u, x∗⟩ ≤ 0 as x∗ ∈ C and A∗λ + u ∈ C+. These
together force λ⊤(−Ax∗ + b) = 0. In addition, ⟨A∗λ + u, x − x∗⟩ = ⟨A∗λ + u, x⟩ ≥ 0 for all x ∈ C.
This is equivalent to⟨A∗λ+u, x⟩ ≥ 0 for all x ∈ (C −x∗), which further implies A∗λ+u ∈ (C −x∗)+.
Hence, u ∈ (C − x∗)+ + ∪λ∈K+{−A∗λ : λ⊤(−Ax∗ + b) = 0}, and the proof is complete.
Appendix C SDP Representation of SOS Problems
It is known that the SOS optimization problem (D) can be expressed equivalently as an SDP. This
appendix provides technical results for Section 3.
27

<!-- page 28 -->
A monomial over v ∈ Rm of degree ¯d is vα = vα1
1 vα2
2 . . . vαm
m with ¯d = Pm
i=1 αi, and α =
(α1, . . . , αm) ∈ ({0} ∪ N)m is a multi-index. The canonical basis is
y(v) := (1, v1, . . . , vm, v2
1, v1v2, . . . , v2
m, . . . , v
¯d
1, . . . , v
¯d
m)⊤,
which is of dimension s(m, ¯d) :=
 m+ ¯d
¯d

.
Let f be a real polynomial with an even degree d = 2 ¯d written as f(v) =P
α∈N (f)αvα, where (f)α
is the α-th coefficient of f and N = {α ∈ ({0} ∪ N)m :Pm
i=1 αi ≤ d}. By [24, Proposition 2.1], f
is SOS if and only if there exists Q ∈ Ss(m,d/2)
+ such that f(v) = y(v)⊤Qy(v). This expresses the
coefficients of f(v) as linear equations of the entries in Q. If we write y(v)y(v)⊤ =P
α∈N Bαvα
for appropriate matrices Bα ∈ Ss(m,d/2), α ∈ N , checking whether f is SOS amounts to finding
Q ∈ Ss(m,d/2)
+ such that tr(QBα) = (f)α for all α ∈ N .
The following proposition shows that the SOS problem (D) and the SDP program (R) share the
same optimal values.
Proposition C.1. Let gk
ℓ , ℓ = 1, . . . , L, k = 1, . . . , r, hj, j = 1, . . . , J, be defined as in Theorem
3.4. Then, max (D) = max (R).
Proof. The SOS constraints
LX
ℓ=1
δk
ℓ gk
ℓ (v) +
JX
j=1
λjhj(v) + λJ+1 − tr(ZkF0) −
mX
i=1
vitr(ZkFi) ∈ Σ2
d(v), k = 1, . . . , r,
of (D) are equivalent to the existence of Qk ∈ Ss(m,d/2)
+ , k = 1, . . . , r, such that
LX
ℓ=1
δk
ℓ (gk
ℓ )α +
JX
j=1
λj(hj)α + λJ+1(1)α − tr(ZkFα) = tr(QkBα), for all α ∈ N , k = 1, . . . , r,
where F0 := F0 ∈ Sν, Fei := Fi ∈ Sν, i = 1, . . . , m, and Fα, α ∈ N \ { 0, e1, . . . , em}, is the zero
matrix. Thus, the conclusion follows.
References
[1] A. A. Ahmadi and P. A. Parrilo. A complete characterization of the gap between convexity and SOS-convexity.
SIAM Journal on Optimization , 23(2):811–833, 2013. https://arxiv.org/pdf/1111.4587.pdf.
[2] A. Ben-Tal, L. El Ghaoui, and A. Nemirovski. Robust Optimization, volume 28. Princeton University Press,
Princeton, 2009.
[3] D. Bertsimas and I. Popescu. On the relation between option and stock prices: a convex optimization approach.
Operations Research, 50(2):358–374, 2002.
[4] J. F. Bonnans and A. Shapiro. Perturbation Analysis of Optimization Problems . Springer Science & Business
Media, New York, 2013.
[5] N. Chieu, V. Jeyakumar, G. Li, and H. Mohebi. Constraint qualifications for convex optimization without
convexity of constraints: new connections and applications to best approximation. European Journal of
Operational Research, 265(1):19–25, 2018.
28

<!-- page 29 -->
[6] E. de Klerk, D. Kuhn, and K. Postek. Distributionally robust optimization with polynomial densities: theory,
models and algorithms. Mathematical Programming, 181:265–296, 2020.
[7] E. de Klerk and M. Laurent. A survey of semidefinite programming approaches to the generalized problem
of moments and their error analysis. In World Women in Mathematics 2018 , pages 17–56, Switzerland, 2019.
Springer.
[8] E. Delage and Y. Ye. Distributionally robust optimization under moment uncertainty with application to
data-driven problems. Operations Research, 58(3):595–612, 2010.
[9] H. U. Gerber and G. Pafum. Utility functions: from risk theory to finance. North American Actuarial Journal,
2(3):74–91, 1998.
[10] J. Goh and M. Sim. Distributionally robust optimization and its tractable approximations. Operations
Research, 58(4):902–917, 2010.
[11] M. Grant and S. Boyd. CVX: Matlab software for disciplined convex programming, version 2.0 beta. http:
//cvxr.com/cvx, 2020. Accessed: 11 Aug 2023.
[12] J. Guo, S. He, B. Jiang, and Z. Wang. A unified framework for generalized moment problems: a novel
primal-dual approach. https: // arxiv. org/ pdf/ 2201. 01445. pdf, 2022.
[13] S. Han, M. Tao, U. Topcu, H. Owhadi, and R. M. Murray. Convex optimal uncertainty quantification. SIAM
Journal on Optimization , 25(3):1368–1387, 2015.
[14] G. A. Hanasusanto, V. Roitch, D. Kuhn, and W. Wiesemann. A distributionally robust perspective on
uncertainty quantification and chance-constrained programming.Mathematical Programming, 151:35–62, 2015.
[15] J. W. Helton and J. Nie. Structured semidefinite representation of some convex sets. In 2008 47th IEEE
Conference on Decision and Control , pages 4797–4800. IEEE, 2008.
[16] J. W. Helton and J. Nie. Semidefinite representation of convex sets. Mathematical Programming, 122:21–64,
2010.
[17] Q. Y. Huang and V. Jeyakumar. A distributional Farkas’ lemma and moment optimization problems with
no-gap dual semi-definite programs. Optimization Letters, pages 1–16, 2024.
[18] V. Jeyakumar, G. M. Lee, J. H. Lee, and Q. Y. Huang. Sum-of-squares relaxations in robust DC optimization
and feature selection. Journal of Optimization Theory and Applications , 200(1):308–343, 2024.
[19] V. Jeyakumar and G. Li. Exact SDP relaxations for classes of nonlinear semidefinite programming problems.
Operations Research Letters, 40(6):529–536, 2012.
[20] V. Jeyakumar, G. Li, and J. Vicente-P´ erez. Robust SOS-convex polynomial optimization problems: exact
SDP relaxations. Optimization Letters, 9:1–18, 2015.
[21] V. Jeyakumar and J. Vicente-P´ erez. Dual semidefinite programs without duality gaps for a class of convex
minimax programs. Journal of Optimization Theory and Applications , 162:735–753, 2014.
[22] J. B. Lasserre. A semidefinite programming approach to the generalized problem of moments. Mathematical
Programming, 112:65–92, 2008.
[23] J. B. Lasserre. Convexity in semialgebraic geometry and polynomial optimization. SIAM Journal on Opti-
mization, 19(4):1995–2014, 2009.
[24] J. B. Lasserre. Moments, Positive Polynomials and Their Applications , volume 1. World Scientific, London,
2009.
[25] J. B. Lasserre. An Introduction to Polynomial and Semi-Algebraic Optimization , volume 52. Cambridge
University Press, Cambridge, 2015.
29

<!-- page 30 -->
[26] H. A. Le Thi, X. T. Vo, and T. Pham Dinh. Feature selection for linear SVMs under uncertain data: robust
optimization based on difference of convex functions algorithms. Neural Networks, 59:36–50, 2014.
[27] A. W. Lo. Semi-parametric upper bounds for option prices and expected payoffs. Journal of Financial
Economics, 19(2):373–387, 1987.
[28] MOSEK ApS. The mosek optimization toolbox for matlab manual. version 10.1. https://docs.mosek.com/
10.0/toolbox/index.html#, 2024. Accessed: 23 Apr 2024.
[29] T. Netzer and D. Plaumann. Geometry of Linear Matrix Inequalities . Springer, 2023.
[30] R. T. Rockafellar and S. Uryasev. Optimization of conditional value-at-risk. Journal of Risk , 2:21–42, 2000.
[31] R. T. Rockafellar and R. J.-B. Wets. Variational Analysis, volume 317. Springer Science & Business Media,
2009.
[32] A. Shapiro. On duality theory of conic linear problems. Semi-Infinite Programming: Recent Advances, pages
135–165, 2001.
[33] W. Wiesemann, D. Kuhn, and M. Sim. Distributionally robust convex optimization. Operations Research,
62(6):1358–1376, 2014.
[34] H. Xu, Y. Liu, and H. Sun. Distributionally robust optimization with matrix moment constraints: Lagrange
duality and cutting plane methods. Mathematical Programming, 169:489–529, 2018.
[35] C. Zalinescu. Convex Analysis in General Vector Spaces . World Scientific, Singapore, 2002.
[36] J. Zhen, D. Kuhn, and W. Wiesemann. A unified theory of robust and distributionally robust optimization
via the primal-worst-equals-dual-best principle. Operations Research, 2023.
30