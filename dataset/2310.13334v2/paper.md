# Convergence analysis on the alternating direction method of multipliers for the cosparse optimization problem

**arXiv ID:** 2310.13334v2

**Authors:** Zisheng Liu, Ting Zhang

**Abstract:** From a dual perspective of the sparse representation model, Nam et al. proposed the cosparse analysis model. In this paper, we aim to investigate the convergence of the alternating direction method of multipliers (ADMM) for the cosparse optimization problem. First, we examine the variational inequality representation of the cosparse optimization problem by introducing auxiliary variables. Second, ADMM is used to solve cosparse optimization problem. Finally, by utilizing a tight frame with a uniform row norm and building upon lemmas and the strict contraction theorem, we establish a worst-case $\mathcal{O}(1/t)$ convergence rate in the ergodic sense.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:2310.13334v2  [math.OC]  22 Nov 2023
Convergence analysis on the alternating direction method o f
multipliers for the cosparse optimization problem
Zisheng Liu a and Ting Zhang b
aSchool of Statistics and Big Data, Henan University of Economics an d Law, Zhengzhou,
China;
bSchool of Mathematics and Information Science, Henan University of Economics and Law,
Zhengzhou, China.
AR TICLE HISTOR Y
Compiled November 23, 2023
ABSTRACT
From a dual perspective of the sparse representation model, Nam et al. proposed the
cosparse analysis model. In this paper, we aim to investigat e the convergence of the
alternating direction method of multipliers (ADMM) for the cosparse optimization
problem. First, we examine the variational inequality repr esentation of the cosparse
optimization problem by introducing auxiliary variables. Second, ADMM is used
to solve cosparse optimization problem. Finally, by utiliz ing a tight frame with a
uniform row norm and building upon lemmas and the strict cont raction theorem,
we establish a worst-case O(1/t) convergence rate in the ergodic sense.
KEYWORDS
sparse representation model, cosparse analysis model, alt ernating direction method
of multipliers, variational inequality, convergence anal ysis
AMS CLASSIFICATION
90C25
1. Introduction
Low-dimensional signal recovery takes advantage of the inh erent low-dimensionality of
many natural signals, despite their high ambient dimension . Utilizing prior informa-
tion about the low-dimensional space can signiﬁcantly aid i n recovering the signal of
interest. Sparsity, a widely recognized form of prior infor mation, serves as the founda-
tion for the burgeoning ﬁeld of compressive sensing (CS [1–5 ]). The recovery of sparse
inputs has found numerous applications in areas such as imag ing, speech, radar signal
processing, sub-Nyquist sampling, and more [6–9]. A typica l sparse recovery problem
is associated with the following linear system:
y = M x, (1)
where y ∈ Rm is an observed vector, M ∈ Rm×d is a measurement matrix and x ∈ Rd
is an unknown signal which would be estimated from y. According to the Nyquist-
CONTACT Zisheng Liu Email: Z. Liu Email: liuzisheng0710@16 3.com

<!-- page 2 -->
Shannon sampling theorem, if the k-space data is undersampled so much that it fails
to meet the Nyquist sampling criterion, then reconstructin g the data can be diﬃcult
or impossible without prior knowledge of x.
1.1. Sparse synthesis model
Over the past decade, the application of compressed sensing signiﬁcantly increased
the image reconstruction speed and eﬃciency because of its c apability to reconstruct
images from highly undersampled signals. Sparse prior is wi dely used in CS-based
reconstruction methods. For the sparse synthesis model, if a vector x is suﬃcient
sparse, under the incoherence assumptions on the measureme nt matrix M , x can be
robustly estimated by the problem
min
x
∥x∥τ
s.t. y = M x,
(2)
where 0 ≤ τ ≤ 1. The advanced ideas and methods have been explored by appli cations
in signals and image processing [10–13]. After years of rese arch, this model is becoming
more and more mature and stable.
1.2. Cosparse analysis model
In the recent decade, the cosparse analysis model is an alter native approach has gained
popularity [14,15,17–20]. Within this framework, a potent ially redundant analysis op-
erator D ∈ Rn×d(n ≥ d) is employed, and the analyzed vector Dx is expected to be
sparse. This implies that a signal x ∈ Rd belongs to the cosparse analysis model with
cosparsity ℓ if ℓ = n − ∥ Dx∥0. In this paper, the quantity ℓ represents the number of
rows in D that are orthogonal to the signal. Consequently, x is referred to as ℓ-cosparse
or simply cosparse. The speciﬁc deﬁnitions of cosparse and c osupport can be found in
literature [15], for ease of reference, we have listed them b elow.
Deﬁnition 1.1 (Cosparse). A signal x ∈ Rd is said to be cosparse with respect to an
analysis operator D ∈ Rn×d if the analysis representation vector Dx contains many
zero elements. Further, the number of zero elements
ℓ = n − ∥ Dx∥0
is called the cosparsity of x, we also say x is ℓ-cosparse.
Deﬁnition 1.2 (Cosupport). For a signal x ∈ Rd and a given analysis operator
D ∈ Rn×d with its rows Dj ∈ Rd(1 ≤ j ≤ n), the cosupport is deﬁned by
Λ := {j|⟨Dj , x ⟩ = 0}.
In this paper, D is a tight frame with uniform row norm. We remind the reader
that a frame is deﬁned as below.
Deﬁnition 1.3 (Frame[21,22]). Let Φ = {ϕ i}N
i=1 ⊆ Rn be a vector sequence of the
2

<!-- page 3 -->
Hilbert space with N ≥ n. If there exist constants 0 < A ≤ B < ∞ such that
∀x ∈ Rn, A ∥x∥2 ≤
N∑
i=1
|⟨x, ϕ i⟩|2 ≤ B∥x∥2, (3)
then Φ is referred to as a ﬁnite frame of Rn. The constants A and B in the above
formula are known as the lower and upper bounds of the ﬁnite fr ame Φ, respectively.
They are considered to be the optimal bounds, with A being the supremum in the
lower bound and B being the inﬁmum in the upper bound. If A = B, then the frame
Φ is called an A-tight frame. If A = B = 1, then Φ is called a Parseval frame. If there
exists a constant C such that each meta-norm ∥ϕ i∥ = C of the frame Φ, then Φ is
called an iso-norm frame. In particular, for a tight frame, i f C = 1, it is referred to as
a uniformly tight frame.
According to the deﬁnition of cosparsity, the cosparse anal ysis model focuses on
the zero elements of the analysis representation vector Dx, rather than the non-zero
elements. This perspective contrasts with the sparse synth esis model. If the cosparsity
ℓ is signiﬁcantly large, meaning that the number of zeros ℓ is close to d, we say that x
has a cosparse representation. The cosupport set is identiﬁ ed by iteratively removing
rows from D for which ⟨Dj, x ⟩ ̸ = 0 until the index set Λ remains unchanged, with
|Λ| ≥ ℓ.
If the analysis representation vector Dx is sparse, similar to the sparse model, the
estimation of x from the measurements can be achieved by
min
x
∥Dx∥0
s.t. y = M x.
(4)
The minimization problem (4) is known to be NP-hard [15], nec essitating the use of
approximation methods. Similar to the sparse model, one opt ion is to use the greedy
analysis pursuit (GAP) approach, which is inspired by the or thogonal matching pursuit
(OMP) algorithm [14–16]. Alternatively, the nonconvex ℓ0 norm can be approximated
by the convex ℓ1 norm, leading to the relaxed problem known as analysis basis pursuit
(ABP) [23]. In this case, x can be estimated by solving a modiﬁed optimization problem
min
x
∥Dx∥1
s.t. ∥y − M x∥2 ≤ ǫ,
(5)
where ∥ · ∥1 is the ℓ1 norm that sums the absolute values of a vector and ǫ is a upper
bound on the noise level ∥v∥2.
ABP is equivalent to the unconstrained optimization
min
x
∥Dx∥1 + α
2 ∥y − M x∥2
2, (6)
which we call analysis LASSO (ALASSO). It can be said that ABP and ALASSO are
equivalent in the sense that for any ǫ > 0, there exists an α such that the optimal
solutions of ABP and ALASSO are identical. For the optimizat ion problem (6), our
previous work presented the modiﬁed GAP algorithm and error analysis [18,24]. The
simulations we conducted demonstrated the advantages of th e proposed method for the
3

<!-- page 4 -->
cosparse optimization problem. These optimization proble ms can also be solved using
interior point methods [25]. However, as the problem dimens ion increases, these tech-
niques become time-consuming since they require solutions of linear systems. Other
suggested approaches include the alternating direction me thod of multipliers (ADMM)
[26–28,37] and the accelerated alternating minimization m ethod (AAM) [29]. In this
paper, we propose a new way to analyze the convergence theory of the cosparse opti-
mization problem based on the variational inequality.
1.3. Organization of the paper
Our focus in this paper is on the cosparse optimization probl em and its convergence
study based on a variational inequality. The paper is struct ured as follows: In Section
2, we introduce auxiliary variables and investigate the var iational inequality char-
acterization of the cosparse optimization problem. In Sect ion 3, we present several
lemmas that establish the strict contraction of the ADMM for the cosparse optimiza-
tion problem. Using these lemmas and the strict contraction theorem, we establish a
worst-case O(1/t ) convergence rate in the ergodic sense. Finally, Section 4 p rovides a
brief conclusion.
2. Preliminaries
To apply the ADMM for solving the cosparse optimization prob lem (6), we convert the
unconstrained optimization problem mentioned above into a constrained optimization
problem as follows
min
x,z
∥z∥1 + α
2 ∥y − M x∥2
2
s.t. Dx − z = 0,
(7)
where an auxiliary variable z ∈ Rn is introduced in (6) to transfer Dx out of the
nondiﬀerentiable term ∥ · ∥1 and α > 0 is a penalty parameter.
In this section, we summarize the variational inequality (V I) characterization of (7).
Initially, we present the optimality condition of the const rained optimization problem
(7), which forms the foundation for our subsequent converge nce analysis [30,32]. We
then proceed to express the Lagrangian function of (7) as fol lows
L(z, x, λ ) = ∥z∥1 + α
2 ∥y − M x∥2
2 − λ T (Dx − z). (8)
In (8), we assume that x ∈ X , z ∈ Z and λ ∈ Rn where X ⊂ Rd and Z ⊂ Rn are
closed convex sets, we call ( z∗, x ∗, λ ∗) ∈ Ω := Z × X × Rn to be a saddle point of
L(z, x, λ ) if the following inequalities are satisﬁed
L(z∗, x ∗, λ ) ≤ L(z∗, x ∗, λ ∗) ≤ L(z, x, λ ∗). (9)
4

<!-- page 5 -->
Obviously, a saddle point ( z∗, x ∗, λ ∗) can be characterized by the system



z∗ = arg min{L(z, x ∗, λ ∗)|z ∈ Z} ,
x∗ = arg min{L(z∗, x, λ ∗)|x ∈ X } ,
λ ∗ = arg max{L(z∗, x ∗, λ )|λ ∈ Rn},
(10)
which can be rewritten as



z∗ ∈ Z , L (z, x ∗, λ ∗) − L(z∗, x ∗, λ ∗) ≥ 0,
x∗ ∈ X , L (z∗, x, λ ∗) − L(z∗, x ∗, λ ∗) ≥ 0,
λ ∗ ∈ Rn, L (z∗, x ∗, λ ∗) − L(z∗, x ∗, λ ) ≥ 0.
(11)
Below, we present a summary of the method for expressing the o ptimality condition
of the cosparse analysis model (7) via a variational inequal ity.
Proposition 2.1. Suppose X ⊂ Rd is a closed convex set, and θ(x) : Rd → R is a
convex function. Furthermore, let f (x) be diﬀerentiable in X . We assume that the set
of solutions for the minimization problem min{θ(x) + f (x)|x ∈ X } is nonempty, then,
x∗ = arg min{θ(x) + f (x)|x ∈ X } (12)
if and only if
x, x ∗ ∈ X , θ (x) − θ(x∗) + (x − x∗)T ∇ f (x∗) ≥ 0. (13)
The proof of Proposition 2.1 is available in [33]. Let θ1(z) = ∥z∥1 and θ2(x) =
α
2 ∥y − M x∥2
2, according to the above inequality (13), a saddle point ( z∗, x ∗, λ ∗) of
the Lagrangian function (8) can be characterized by a soluti on point of the following
variational inequality
ω, ω ∗ ∈ Ω, θ (u) − θ(u∗) + (ω − ω ∗)T F (ω ∗) ≥ 0, (14)
where
θ(u) = θ1(z) + θ2(x), Ω = Z × X × Rn, (15)
and
ω =


z
x
λ

 , u =
(
z
x
)
, F (ω ) =


λ
− DT λ
Dx − z

 , (16)
Since F is an aﬃne operator, and
F (ω ) =


0 0 I
0 0 − DT
− I D 0




z
x
λ

 , (17)
According to the antisymmetry of the aﬃne matrix, it follows that
(ω − ¯ω )T [F (ω ) − F (¯ω )] ≡ 0, ∀ ω, ¯ω ∈ Ω. (18)
5

<!-- page 6 -->
Using inequality (13) and combining (8), we derive the follo wing conclusion with
(z∗, x ∗, λ ∗) ∈ Ω,



θ1(z) − θ1(z∗) + (z − z∗)T λ ∗ ≥ 0,
θ2(x) − θ2(x∗) + (x − x∗)T (− DT λ ∗) ≥ 0,
(λ − λ ∗)T (Dx∗ − z∗) ≥ 0.
(19)
After conducting the aforementioned analysis, the linear c onstrained cosparse opti-
mization problem is reformulated as a variational inequali ty. Consequently, the task is
ultimately simpliﬁed to identifying a saddle point of the La grangian function. In the
subsequent section, the convergence analysis of the ADMM me thod for addressing the
cosparse optimization problem, as denoted by equation (7), will be discussed.
3. Convergence analysis of the cosparse optimization problem
3.1. Variational inequality characterization of ADMM
The augmented Lagrangian function of the problem (7) can be f ormulated as follows
Lβ(z, x, λ ) =∥z∥1 + α
2 ∥y − M x∥2
2 − λ T (Dx − z) + β
2 ∥Dx − z∥2
2, (20)
where λ is the Lagrange multiplier and β > 0 is a penalty parameter for the lin-
ear constraints. Thus, applying directly the augmented Lag rangian function (20) and
starting with an initial iterate ( x0, λ 0) ∈ X × Rn, the ADMM generates its sequence
via following iterative scheme



zk+1 = arg min{Lβ(z, x k, λ k)|z ∈ Z} ,
xk+1 = arg min{Lβ(zk+1, x, λ k)|x ∈ X } ,
λ k+1 = λ k − β (Dxk+1 − zk+1), λ ∈ Rn
(21)
the corresponding variational inequalities of (21) can be g iven as



θ1(z) − θ1(zk+1) + (z − zk+1)T [λ k − β (Dxk − zk+1)] ≥ 0,
θ2(x) − θ2(xk+1) + (x − xk+1)T [− DT λ k + βD T (Dxk+1 − zk+1)] ≥ 0,
(λ − λ k+1)T [(Dxk+1 − zk+1) + 1
β (λ k+1 − λ k)] ≥ 0.
(22)
For some reviews on the classical ADMM, one can refer to liter atures [28,30,31,34–36].
3.2. Assertions
To establish that {ω k} is strictly contractive with respect Ω, we ﬁrst present seve ral
lemmas.
Lemma 3.1. Let the sequence {ω k} be generated by (21). Then, we have
θ(u) − θ(uk+1) + (ω − ω k+1)T F (ω )
≥ (z − zk+1)T β (Dxk − Dxk+1) + 1
β (λ − λ k+1)T (λ k − λ k+1), ∀ω ∈ Ω. (23)
6

<!-- page 7 -->
Proof. From (22) we know that
θ1(z) − θ1(zk+1) + (z − zk+1)T [λ k − β (Dxk − zk+1)] ≥ 0, ∀z ∈ Z (24)
and
θ2(x) − θ2(xk+1) + (x − xk+1)T (− DT λ k + βD T (Dxk+1 − zk+1)) ≥ 0, ∀x ∈ X . (25)
Using λ k+1 = λ k − β (Dxk+1 − zk+1) we can easily deduce
λ k = λ k+1 + β (Dxk+1 − zk+1) (26)
and
(Dxk+1 − zk+1) = 1
β (λ k − λ k+1). (27)
Putting the formulations (26) and (27) into (24) and (25), re spectively, then we have
the following inequalities
θ1(z) − θ1(zk+1) + (z − zk+1)T [λ k+1 + β (Dxk+1 − zk+1) − β (Dxk − zk+1)] ≥ 0,
(28)
θ2(x) − θ2(xk+1) + (x − xk+1)T (− DT λ k+1) ≥ 0, (29)
and
(λ − λ k+1)T (Dxk+1 − zk+1) ≥ (λ − λ k+1)T 1
β (λ k − λ k+1). (30)
Combining (28), (29) and (30) we have



θ1(z) − θ1(zk+1) + (z − zk+1)T λ k+1 ≥ (z − zk+1)T β (Dxk − Dxk+1),
θ2(x) − θ2(xk+1) + (x − xk+1)T (− DT λ k+1) ≥ 0,
(λ − λ k+1)T (Dxk+1 − zk+1) ≥ (λ − λ k+1)T 1
β (λ k − λ k+1),
(31)
which is
θ(u) − θ(uk+1) + (ω − ω k+1)T F (ω k+1)
≥ (z − zk+1)T β (Dxk − Dxk+1) + (λ − λ k+1)T 1
β (λ k − λ k+1).
Note that the matrix in the operator F is skew-symmetric, then, using (18), we have
θ(u) − θ(uk+1) + (ω − ω k+1)T F (ω )
≥ (z − zk+1)T β (Dxk − Dxk+1) + (λ − λ k+1)T 1
β (λ k − λ k+1). (32)
The Lemma 3.1 is proved.
7

<!-- page 8 -->
Lemma 3.2. Let the sequence {ω k} be generated by (21). Then, we have
β (z − zk+1)T (Dxk − Dxk+1) + 1
β (λ − λ k+1)T (λ k − λ k+1)
= − 1
2β ∥λ k − λ∥2
2 − β
2 ∥Dxk − z∥2
2 + 1
2β ∥λ k+1 − λ∥2
2 + β
2 ∥Dxk+1 − z∥2
2
+ β
2 ∥Dxk − zk+1∥2
2.
(33)
Proof. Applying the identity
(a − b)T (c − d) = 1
2 {∥a − d∥2
2 − ∥ a − c∥2
2} + 1
2 {∥c − b∥2
2 − ∥ d − b∥2
2}
to the left-hand side in (33) with
a = z, b = zk+1, c = Dxk, d = Dxk+1,
we obtain
β (z − zk+1)T (Dxk − Dxk+1)
= β
2 {∥z − Dxk+1∥2
2 − ∥ z − Dxk∥2
2} + β
2 {∥Dxk − zk+1∥2
2 − ∥ Dxk+1 − zk+1∥2
2}.
(34)
Using the identity
bT (b − a) = 1
2 (∥b∥2
2 − ∥ a∥2
2 + ∥b − a∥2
2),
and let
a = λ − λ k, b = λ − λ k+1,
we obtain
1
β (λ − λ k+1)T (λ k − λ k+1) = 1
2β {∥λ − λ k+1∥2
2 − ∥ λ − λ k∥2
2 + ∥λ k − λ k+1∥2
2}. (35)
Using
β ∥Dxk+1 − zk+1∥2
2 = 1
β ∥λ k − λ k+1∥2
2,
and combining (34) and (35), we complete the proof of this lem ma.
Lemma 3.3. Let the sequence {xk}, {zk} and {λ k} be generated by (21), then,
β ∥Dxk − zk+1∥2
2 ≥ β ∥Dxk − Dxk+1∥2
2 + 1
β ∥λ k − λ k+1∥2
2. (36)
8

<!-- page 9 -->
Proof. Based on the second inequality of inequality (31), we can der ive the following
result
{
θ2(x) − θ2(xk+1) + (x − xk+1)T (− DT λ k+1) ≥ 0,
θ2(x) − θ2(xk) + (x − xk)T (− DT λ k) ≥ 0. (37)
Let x = xk and x = xk+1 in (37), respectively, then
{
θ2(xk) − θ2(xk+1) + (xk − xk+1)T (− DT λ k+1) ≥ 0,
θ2(xk+1) − θ2(xk) + (xk+1 − xk)T (− DT λ k) ≥ 0.
From above inequalities, we have
(λ k − λ k+1)T (Dxk − Dxk+1) ≥ 0. (38)
Using
(Dxk+1 − zk+1) = 1
β (λ k − λ k+1),
then we obtain
β ∥Dxk − zk+1∥2
2
=β ∥Dxk − Dxk+1 + Dxk+1 − zk+1∥2
2
=β ∥Dxk − Dxk+1 + 1
β (λ k − λ k+1)∥2
2
≥ β ∥Dxk − Dxk+1∥2
2 + 1
β ∥λ k − λ k+1)∥2
2.
(39)
The proof of this lemma is completed.
3.3. Strict contraction
To present the main result of the paper, it is necessary to est ablish the strict contrac-
tility of the iterative sequence. The following subsection provides a proof of the strong
contractility of the iterative sequence {ω k}, which relies on Lemma 3.1, Lemma 3.2,
and Lemma 3.3.
Theorem 3.4. Assuming that the sequence {ω k} is generated by equation (21), we
can state the following
∥vk+1 − v∗∥2
H ≤ ∥ vk − v∗∥2
H − ∥ vk − vk+1∥2
H (40)
where
v =
(
λ
x
)
, H =
( 1
β Im 0
0 βI d
)
, V ∗ = {(λ ∗, x ∗)|(z∗, x ∗, λ ∗) ∈ Ω}. (41)
9

<!-- page 10 -->
Proof. We can deduce from Lemma 3.1 and Lemma 3.2 that
θ(uk+1) − θ(u) + (ω k+1 − ω )T F (ω )
≤ 1
2β ∥λ k − λ∥2
2 + β
2 ∥Dxk − z∥2
2 − 1
2β ∥λ k+1 − λ∥2
2 − β
2 ∥Dxk+1 − z∥2
2
− β
2 ∥Dxk − zk+1∥2
2.
(42)
By utilizing Lemma 3.3, we can rewrite equation (42) as follo ws
0 ≤ θ(uk+1) − θ(u∗) + (ω k+1 − ω ∗)T F (ω ∗)
≤ 1
2β ∥λ k − λ ∗∥2
2 + β
2 ∥Dxk − z∗∥2
2 − 1
2β ∥λ k+1 − λ ∗∥2
2 − β
2 ∥Dxk+1 − z∗∥2
2
− 1
2β ∥λ k − λ k+1∥2
2 − β
2 ∥Dxk − Dxk+1∥2
2.
(43)
That is
1
β ∥λ k+1 − λ ∗∥2
2 + β ∥Dxk+1 − z∗∥2
2
≤ 1
β ∥λ k − λ ∗∥2
2 + β ∥Dxk − z∗∥2
2 − ( 1
β ∥λ k − λ k+1∥2
2 + β ∥Dxk − Dxk+1∥2
2).
(44)
Let
Dx∗ = z∗, v =
(
λ
x
)
and H =
( 1
β Im 0
0 βD T D
)
,
therefore, the left-hand side of inequality (44) becomes
1
β ∥λ k+1 − λ ∗∥2
2 + β ∥Dxk+1 − z∗∥2
2
=
(
λ k+1 − λ ∗
xk+1 − x∗
) T ( 1
β Im 0
0 βD T D
) (
λ k+1 − λ ∗
xk+1 − x∗
)
=(vk+1 − v∗)T
( 1
β Im 0
0 βD T D
)
(vk+1 − v∗)
=∥vk+1 − v∗∥2
H.
(45)
Likewise, the sum of the ﬁrst two terms on the right-hand side of inequality (44) is
1
β ∥λ k − λ ∗∥2
2 + β ∥Dxk − Dx∗∥2
2
=
(
λ k − λ ∗
xk − x∗
) T ( 1
β Im 0
0 βD T D
) (
λ k − λ ∗
xk − x∗
)
=(vk − v∗)T
( 1
β Im 0
0 βD T D
)
(vk − v∗)
=∥vk − v∗∥2
H
(46)
10

<!-- page 11 -->
and the sum of the last two terms on the right-hand side of ineq uality (44) is
1
β ∥λ k − λ k+1∥2
2 + β ∥Dxk − Dxk+1∥2
2
=
(
λ k − λ k+1
xk − xk+1
) T ( 1
β Im 0
0 βD T D
) (
λ k − λ k+1
xk − xk+1
)
=(vk − vk+1)T
( 1
β Im 0
0 βD T D
)
(vk − vk+1)
=∥vk − vk+1∥2
H.
(47)
Since D ∈ Rn×d is a unit tight frame, we have that DT D = Id. By combining formulas
(45), (46), and (47), we complete the proof of the Theorem 3.4 .
According to Theorem 3.4, we know that H is a positive deﬁnite matrix, and in-
equality (40) implies that the sequence {vk} is bounded. Assuming that the initial
vector is v0 = ( λ 0, x 0)T , we can obtain the following expression by summing both
sides of inequality (40)
∞∑
k=0
∥vk − vk+1∥2
H ≤ ∥ v0 − v∗∥2
H . (48)
The above equation indicates that lim k→∞ ∥vk − vk+1∥2
H = 0. Therefore, any subse-
quence vkj of vk also has lim j→∞ ∥vkj − vkj+1∥2
H = 0. Suppose there exists a subse-
quence that converges to ¯v, then formula (23) implies that ¯v is the solution of formula
(21). This shows that any accumulation point of the sequence vk is a solution of (21).
According to formula (40), vk cannot have more than one accumulation point, and
hence vk converges to ¯v ∈ V ∗.
3.4. Convergence rate in ergodic sense
Combining with Theorem 3.4, we prove a worst-case O(1/t ) convergence rate in a
ergodic sense of the ADMM scheme (21) for cosparse signal rec onstruction problem.
Theorem 3.5. Let the sequence {ω k} be generated by (21). Then, for any positive
integer t, we have
θ(ut) − θ(u) + (ω t − ω )T F (ω )
≤ 1
2(t + 1) [ 1
β ∥λ 0 − λ∥2
2 + β ∥Dx0 − z∥2
2], ∀ω ∈ Ω (49)
where
ω t = 1
t + 1 (
t∑
k=0
ω k+1). (50)
11

<!-- page 12 -->
Proof. For any integer k, by (42) we obtain
θ(uk+1) − θ(u) + (ω k+1 − ω )T F (ω )
≤ 1
2β ∥λ k − λ∥2
2 + β
2 ∥Dxk − z∥2
2 − 1
2β ∥λ k+1 − λ∥2
2 − β
2 ∥Dxk+1 − z∥2
2.
(51)
Suppose k = 0 , 1, 2, . . . , t are non-negative integers. By summing the left and right
ends of the inequality (51), we deduce that
t∑
k=0
θ(uk+1) − (t + 1)θ(u) +
[ t∑
k=0
ω k+1 − (t + 1)ω
]T
F (ω )
≤ 1
2β ∥λ 0 − λ∥2
2 + β
2 ∥Dx0 − z∥2
2, ∀ω ∈ Ω.
(52)
The left and right ends of the inequality (52) are multiplied by 1
t+1 at the same time,
and let
ω t = 1
t + 1
t∑
k=0
ω k+1, (53)
then, the inequality (52) is equivalent to
1
t + 1
t∑
k=0
θ(uk+1) − θ(u) + (ω t − ω )T F (ω )
≤ 1
2(t + 1) [ 1
β ∥λ 0 − λ∥2
2 + β ∥Dx0 − z∥2
2].
(54)
Given that the function θ(u) is convex, let
ut = 1
t + 1
t∑
k=0
uk+1 = 1
t + 1 (u1 + u2 + · · ·+ ut),
we can derive the following expression
θ(ut) =θ
[ 1
t + 1 (u1 + u2 + · · ·+ ut)
]
≤ 1
t + 1
[
θ(u1) + θ(u1) + · · ·+ θ(ut)
]
= 1
t + 1
t∑
k=0
θ(uk+1).
(55)
By utilizing equations (54) and (55), we complete the proof o f the Theorem 3.5.
12

<!-- page 13 -->
After t-th iterations, then ω t deﬁned by (53) satisﬁes
˜ω ∈ Ω and sup
ω∈D ˜ω
{θ(˜u) − θ(u) + (˜ω − ω )T F (ω )} ≤ d
2t = O( 1
t ),
where
D˜ω = {ω ∈ Ω|∥ω − ˜ω ∥ ≤ 1}, d := sup{ 1
β ∥λ 0 − λ∥2
2 + β ∥Dx0 − z∥2
2|ω ∈ D ˜ω},
and v0 = ( λ 0, x 0) is the initial iteration point. That means ω t is an O( 1
t ) solution of
(22).
4. Conclusions
This paper presents a novel approach to analyze the converge nce of cosparse optimiza-
tion problem. In order to complete the proof of the main theor em, we ﬁrst give the
overall framework of the ADMM to solve the cosparse optimiza tion problem, secondly,
through three lemmas, this paper gives the basic inequaliti es required for the theorem
proof, and ﬁnally, our analysis establishes a worst-case co nvergence rate of O(1/t ),
which demonstrates the eﬀectiveness of our approach.
Researchers currently rely on a range of methods to solve sep arable convex opti-
mization problems. Two popular approaches are the generali zed symmetric ADMM
and parameterizable proximal point algorithms [37–39]. Th ese methods have demon-
strated their eﬀectiveness and superiority in various exper iments. In our future work,
we plan to explore the potential of combining these methods t o solve the cosparse
signal reconstruction problem.
F unding
The authors were supported by the National Natural Science F oundation of China
Mathematics Tian Yuan Fund under grant No. 12226323 and 1222 6315, the National
Natural Science Foundation of China under grant No. 6210313 6, the Henan Province
Undergraduate College Youth Backbone Teacher Training Pro gram.
Acknowledgments
The authors wish to thank Professor Zheng-Hai Huang for prov iding his valuable
comments which have signiﬁcantly improved the quality of th is paper.
References
[1] D. L. Donoho, Compressed Sensing. IEEE Trans. Inform. Theory, 2006, 52(4), 1289-1306.
[2] Y. Tsaig, David L. Donoho, Extensions of compressed sensing . Signal Processing, 2006,
86(3), 549-571.
[3] E. J. Cand´ es, T. Tao,Decoding by linear programming. IEEE Trans. Inform. Theory, 2004,
51(12), 4203-4215.
13

<!-- page 14 -->
[4] S. Foucart, H. Rauhut, A Mathematical Introduction to Compressive Sensing , Springer,
New York. 2013.
[5] I. Mishra and S. Jain, Soft computing based compressive sensing techniques in sig nal
processing: A comprehensive review , J. Intell. Syst., 2021, 30, 312-326.
[6] A. M. Bruckstein, D. L. Donoho, M. Elad, From sparse solutions of systems of equations
to sparse modeling of signals and images , SIAM Rev. 2009, 51 34-81.
[7] M. Elad, Sparse and Redundant Representations: From Theory to Appli cations in Signal
and Image Processing , Springer press, 2010.
[8] S. Mallat, A Wavelet Tour of Signal Processing, Third Edition: The Spar se Way, 3rd
edition, Academic Press, 2008.
[9] J. L. Starck, F. Murtagh, M. J. Fadili, Sparse Image and Signal Processing-Wavelets,
Curvelets, Morphological Diversity , Cambridge University Press, 2010.
[10] M. Elad, M. Aharon, Image denoising via sparse and redundant representations o ver
learned dictionaries, IEEE Trans. Image Process, 2006, 15 (12), 3736-3745.
[11] M. Elad, J. L. Starck, Querre, P, Donoho, D.L, Simultaneous cartoon and texture im-
age inpainting using morphological component analysis (MC A), Appl. Comput. Harmon.
Anal. 2005, 19 (3), 340-358.
[12] A. L. Casanovas, G. Monaci, P. Vandergheynst, R. Gribonval, Blind audiovisual source
separation based on sparse representations , IEEE Trans. Multimedia. 2010, 12 (5), 358-
371.
[13] M. D. Plumbley, T. Blumensath, L. Daudet, R. Gribonval, M. E. Da vies, Sparse repre-
sentations in audio and music: from coding to source separat ion, Proc. IEEE., 2010, 98
(6), 995-1005.
[14] S. Nam, M. Davies, M. Elad, R. Gribonval, Cosparse analysis modeling-uniqueness and
algorithms, in: IEEE International Conference on Acoustics. Speech and Sig nal Processing,
ICASSP 2011, Prague, Czech Republic, May, 2011, 5804-5807.
[15] S. Nam, M. Davies, M. Elad, R. Gribonval, The cosparse analysis model and algorithms ,
Appl. Comput. Harmon. Anal., 2013, 34 (1), 30-56.
[16] T. Tirer and R. Giryes, Generalizing CoSaMP to signals from a union of low dimension al
linear subspaces, Appl. Comput. Harmon. Anal., 2020, 49(1), 99-122.
[17] R. Giryes, S. Nam, M. Elad, R. Gribonval, M. Davies, Greedy-like algorithms for the
cosparse analysis model , Linear Algebra Appl., 2014, 441 (15), 22-60.
[18] Z. Liu, J. Li, W. Li and P. Dai, A modiﬁed greedy analysis pursuit algorithm for the
cosparse analysis model , Numer. Algor., 2017, 74, 867-887.
[19] H. Song, X. Ren, Y. Lai, H. Meng, Sparse analysis recovery via iterative cosupport detec-
tion estimation , IEEE Access, 2021, 9, 38386-38395.
[20] M. S. Kotzagiannidis and M. E. Davies. Analysis vs synthesis with structure - An inves-
tigation of union of subspace models on graphs , Appl. Comput. Harmon Anal., 2022, 60,
293-332.
[21] V. K. Goyal, J. Kovacevic, J. A. Kelner. Quantized frame expansions with erasures , Appl.
Comput. Harmon. Anal., 2001, 10(3), 203-233.
[22] D. Li, M. Xue. Bases and frames in Banach spaces (in Chinese) , Science Press, Beijing,
2007.
[23] T. Zhao, C. E. Yonina, B Amir, and N. Arye, Smoothing and decomposition for analysis
sparse recovery, IEEE Trans. Signal Process., 2014, 62(7),1762-1774.
[24] J. Li, Z. Liu, W. Li, The reweighed greedy analysis pursuit algorithm for the cos parse
analysis model, IEEE. The 2015 11th International Conference on Natural Com putation,
2015, 1018-1022.
[25] S. Boyd and L. Vandenberghe, Convex Optimization . Cambridge, U.K, Cambridge Univ.
Press, 2004.
[26] S. Boyd, N. Parikh, E. Chu, B. Peleato and J. Eckstein, Distributed optimization and
statistical learning via alternating direction method of m ultipliers, in Found. Trends Mach.
Learn., 2010, 3(1), 1-122.
[27] M. V. Afonso, J. M. Bioucas-Dias and M. A. T. Figueiredo, Fast image recovery using
14

<!-- page 15 -->
variable splitting and constrained optimization , IEEE Trans. Image Process., 2010, 19(9),
2345-2356.
[28] D-R. Han, A survey on some Recent developments of alternating Directi on method of
multipliers, J. Oper. Res. Soc. China, 2022, 10, 1-52.
[29] J. X. Xie, A. P. Liao and Y. Lei, A new accelerated alternating minimization method for
analysis sparse recovery , Signal Process., 2018, 145, 167-174.
[30] B. S. He and X. M. Yuan, On the O(1/n ) convergence rate of the Douglas-Rachford
alternating direction method , SIAM J. Numer. Anal, 2012, 130 (3), 702-709.
[31] S-L. Hu and Z-H. Huang, Alternating direction method for bi-quadratic programmin g,
Journal of Global Optimization, 2011, 51(3), 429-446.
[32] Z-H. Huang, G. H. Lin and N. H. Xiu, Several developments of variational inequalities
and complementarity problems, bilevel programming and MPE C, Operations Research
Transactions, 2014, 18(1), 113-133.
[33] B. S. He, F. Ma and X. Yuan, Convergence study on the symmetric version of ADMM
with large step sizes , SIAM J. Imaging sciences, 2016, 9 (3), 1467-1501.
[34] S. Boyd, N. Parikh, E. Chu, B. Peleato and J. Eckstein, Distributed optimization and sta-
tistical learning via the alternating direction method of m ultipliers, Found. Trends Mach.
Learn., 2011, 3 (1), 1-122.
[35] J. Eckstein and W. Yao, Understanding the convergence of the alternating directio n
method of multipliers: Theoretical and computational pers pectives, Paciﬁc J. Optim., 2015,
11, 619-644.
[36] T. F. Chan and R. Glowinski, Finite element approximation and iterative solution of a
class of mildly nonlinear elliptic equations , Technical Report STAN-CS-78-674, Computer
Science Department, Stanford University, Stanford, CA, 1978.
[37] Z. H. Jia, G. Xue, X. J. Cai and D-R. Han, The convergence rate analysis of the symmetric
ADMM for the nonconvex separable optimization problems , J. Ind. Manag. Optim., 2021,
17(4), 1943-1971.
[38] J. Bai, X. Chang, J. Li and F. Xu, Convergence revisit on generalized symmetric ADMM ,
Optim., 2021, 70, 149-168.
[39] M. M, Alves, J. Eckstein, M. Geremia and J. G. Melo, Relative-error inertial-relaxed
inexact versions of Douglas-Rachford and ADMM splitting al gorithms, Comput. Optim.
Appl., 2020, 75, 389-422.
15