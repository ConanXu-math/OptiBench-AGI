# The exact worst-case convergence rate of the alternating direction method of multipliers

**arXiv ID:** 2206.09865v4

**Authors:** Moslem Zamani, Hadi Abbaszadehpeivasti, Etienne de Klerk

**Abstract:** Recently, semidefinite programming performance estimation has been employed as a strong tool for the worst-case performance analysis of first order methods. In this paper, we derive new non-ergodic convergence rates for the alternating direction method of multipliers (ADMM) by using performance estimation. We give some examples which show the exactness of the given bounds. We also study the linear and R-linear convergence of ADMM. We establish that ADMM enjoys a global linear convergence rate if and only if the dual objective satisfies the Polyak-Lojasiewicz (PL)inequality in the presence of strong convexity. In addition, we give an explicit formula for the linear convergence rate factor. Moreover, we study the R-linear convergence of ADMM under two new scenarios.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:2206.09865v4  [math.OC]  24 May 2023
Noname manuscript No.
(will be inserted by the editor)
The exact worst-case convergence rate of the
alternating direction method of multipliers
Moslem Zamani · Hadi
Abbaszadehpeivasti · Etienne de Klerk
Received: date / Accepted: date
Abstract Recently, semideﬁnite programming performance estimation has
been employed as a strong tool for the worst-case performance analysis of
ﬁrst order methods. In this paper, we derive new non-ergodic con vergence
rates for the alternating direction method of multipliers (ADMM) by u sing
performance estimation. We give some examples which show the exac tness
of the given bounds. We also study the linear and R-linear convergen ce of
ADMM. We establish that ADMM enjoys a global linear convergence ra te if
and only if the dual objective satisﬁes the Polyak-/suppress Lojasiewicz (P/suppressL) inequality
in the presence of strong convexity. In addition, we give an explicit f ormula for
the linear convergence rate factor. Moreover, we study the R-linear convergence
of ADMM under two new scenarios.
Keywords Alternating direction method of multipliers (ADMM) · Perfor-
mance estimation · Convergence rate · P/suppress L inequality
Mathematics Subject Classiﬁcation (2010) 90C22 · 90C25 · 65K15
This work was supported by the Dutch Scientiﬁc Council (NWO) grant
OCENW.GROOT.2019.015, Optimization for and with Machine Learning (OPTIMAL) .
M. Zamani
Tilburg University, Department of Econometrics and Operat ions Research, Tilburg, The
Netherlands
E-mail: m.zamani
1@tilburguniversity.edu
H. Abbaszadehpeivasti
Tilburg University, Department of Econometrics and Operat ions Research, Tilburg, The
Netherlands
E-mail: h.abbaszadehpeivasti@tilburguniversity.edu
E. de Klerk
Tilburg University, Department of Econometrics and Operat ions Research, Tilburg, The
Netherlands
E-mail: e.deklerk@tilburguniversity.edu

<!-- page 2 -->
2 Moslem Zamani et al.
1 Introduction
We consider the optimization problem
min
(x,z )∈ Rn× Rm
f (x) + g(z), (1)
s. t. Ax + Bz = b,
where f : Rn → R ∪ {∞} and g : Rm → R ∪ {∞} are closed proper convex
functions, 0 ̸= A ∈ Rr× n, 0 ̸= B ∈ Rr× m and b ∈ Rr. Moreover, we assume
that ( x⋆ , z ⋆ ) is an optimal solution of problem (1) and λ ⋆ is its corresponding
Lagrange multipliers. Moreover, we denote the value of f and g at x⋆ and z⋆
with f ⋆ and g⋆ , respectively.
Problem (1) appears naturally (or after variable splitting) in many ap pli-
cations in statistics, machine learning and image processing to name b ut a
few [6, 20, 25, 37]. The most common method for solving problem (1) is the
alternating direction method of multipliers (ADMM). ADMM is dual base d
approach that exploits separable structure and it may be describe d as follows.
Algorithm 1 ADMM
Set N and t > 0 (step length), pick λ0, z0.
For k = 1, 2, . . . , N perform the following step:
1. xk ∈ argmin f (x) + ⟨λk−1, Ax⟩ + t
2 ∥Ax + Bz k−1 − b∥2
2. zk ∈ argmin g(z) + ⟨λk−1, Bz⟩ + t
2 ∥Axk + Bz − b∥2
3. λk = λk−1 + t(Axk + Bz k − b).
ADMM was ﬁrst proposed in [11,13] for solving nonlinear variational p rob-
lems. We refer the interested reader to [14] for a historical review of ADMM.
The popularity of ADMM is due to its capability to be implemented parallelly
and hence can handle large-scale problems [6,19,30,40]. For exam ple, it is used
for solving inverse problems governed by partial diﬀerential equat ion forward
models [28], and distributed energy resource coordinations [26], to mention
but a few.
The convergence of ADMM has been investigated extensively in the lit er-
ature and there exist many convergence results. However, diﬀer ent perfor-
mance measures have been used for the computation of converge nce rate;
see [10, 15, 16, 21, 24, 25, 31, 39]. In this paper, we consider the d ual objective
value as a performance measure.
Throughout the paper, we assume that each subproblem in steps 1 and 2
of Algorithm 1 attains its minimum. The Lagrangian function of problem (1)
may be written as
L(x, z, λ ) = f (x) + g(z) + ⟨λ, Ax + Bz − b⟩, (2)

<!-- page 3 -->
Convergence rate of ADMM 3
and the dual objective of problem (1) is also deﬁned as
D(λ) = min
(x,z )∈ Rn× Rm
f (x) + g(z) + ⟨λ, Ax + Bz − b⟩.
We assume throughout the paper that strong duality holds for pro blem (1),
that is
max
λ ∈ Rr
D(λ) = min
Ax+Bz=b
f (x) + g(z).
Note that we have strong duality when both functions f and g are real-valued.
For extended convex functions, strong duality holds under some m ild condi-
tions; see e.g. [3, Chapter 15].
Some common performance measures for the analysis of ADMM are a s
follows,
– Objective value:
⏐
⏐f (xN ) + g(zN ) − f ⋆ − g⋆ ⏐
⏐;
– Primal and dual feasibility:

AxN + Bz N − b

 and

AT B(zN − zN − 1)

;
– Dual objective value: D(λ ⋆ ) − D(λ N );
– Distance between ( xN , z N , λ N ) and a saddle points of problem (2).
Note that the mathematical expressions are written in a non-ergo dic sense for
convenience. Each measure is useful in monitoring the progress an d conver-
gence of ADMM. The objective value is the most commonly used perfo rmance
measure for the analysis of algorithms in convex optimization [3,4,33 ]. As men-
tioned earlier, ADMM is a dual based method and it may be interpreted as a
proximal method applied to the dual problem; see [4,25] for furth er discussions
and insights. Thus, a natural performance measure for ADMM wou ld be dual
objective value. In this study, we investigate the convergence ra te of ADMM
in terms of dual objective value and feasibility. It worth noting that most
performance measures may be analyzed through the framework d eveloped in
Section 2.
Regarding dual objective value, the following convergence rate is k nown in
the literature. This theorem holds for strongly convex functions f and g; recall
that f is called strongly convex with modulus µ ≥ 0 if the function f − µ
2 ∥ · ∥2
is convex.
Theorem 1 [16, Theorem 1] Let f and g be strongly convex with moduli
µ 1 > 0 and µ 2 > 0, respectively. If t ≤ 3
√
µ 1µ 2
2
λ max(AT A)λ 2
max(BT B) , then
D(λ ⋆ ) − D(λ N ) ≤ ∥λ 1 − λ ⋆ ∥22t(N − 1) . (3)
In this study we establish that Algorithm 1 has the convergence rat e of
O( 1
N ) in terms of dual objective value without assuming the strong conv exity
of g. Under this setting, we also prove that Algorithm 1 has the converg ence
rate of O( 1
N ) in terms of primal and dual residuals. Moreover, we show that
the given bounds are exact. Furthermore, we study the linear and R-linear
convergence.

<!-- page 4 -->
4 Moslem Zamani et al.
Outline of our paper
Our paper is structured as follows. We present the semideﬁnite pro gramming
(SDP) performance estimation method in Section 2, and we develop t he perfor-
mance estimation to handle dual based methods including ADMM. In Se ction
3, we derive some new non-asymptotic convergence rates by using performance
estimation for ADMM in terms of dual function, primal and dual resid uals.
Furthermore, we show that the given bounds are tight by providing some ex-
amples. In Section 4 we proceed with the study of the linear converg ence of
ADMM. We establish that ADMM enjoys a linear convergence if and only
if the dual function satisﬁes the P/suppress L inequality when the objective function
is strongly convex. Furthermore, we investigate the relation betw een the P/suppress L
inequality and common conditions used by scholars to prove the linear con-
vergence. Section 5 is devoted to the R-linear convergence. We pr ove that
ADMM is R-linear convergent under two new scenarios which are weak er than
the existing ones in the literature.
Terminology and notation
In this subsection we review some deﬁnitions and concepts from con vex anal-
ysis. The interested reader is referred to the classical text by Ro ckafellar [36]
for more information. The n-dimensional Euclidean space is denoted by Rn.
We use ⟨·, ·⟩and ∥ · ∥to denote the Euclidean inner product and norm, respec-
tively. The column vector ei represents the i-th standard unit vector and I
stands for the identity matrix. For a matrix A, Ai,j denotes its ( i, j )-th entry,
and AT represents the transpose of A. The notation A ⪰ 0 means the matrix
A is symmetric positive semideﬁnite. We use λ max(A) and λ min(A) to denote
the largest and the smallest eigenvalue of symmetric matrix A, respectively.
Moreover, the seminorm ∥ · ∥A is deﬁned as ∥x∥A = ∥Ax∥ for any A ∈ Rm× n.
Suppose that f : Rn → (−∞ , ∞ ] is an extended convex function. The
function f is called closed if its epi-graph is closed, that is {(x, r ) : f (x) ≤ r}
is a closed subset of Rn+1. The function f is said to be proper if there exists
x ∈ Rn with f (x) < ∞ . We denote the set of proper and closed convex
functions on Rn by F0(Rn). The subgradients of f at x is denoted and deﬁned
as
∂f (x) = {ξ : f (y) ≥ f (x) + ⟨ξ, y − x⟩, ∀y ∈ Rn}.
We call a diﬀerentiable function f L -smooth if for any x1, x 2 ∈ Rn,
∥∇ f (x1) − ∇ f (x2)∥ ≤ L∥x1 − x2∥ ∀ x1, x 2 ∈ Rn.
Deﬁnition 1 Let f : Rn → (−∞ , ∞ ] be a closed proper function and let
A ∈ Rm× n. We say f is c-strongly convex relative to ∥. ∥A if the function
f − c
2 ∥. ∥2
A is convex.
In the rest of the section, we assume that A ∈ Rm× n. It is seen that
any µ -strongly convex function is µ
λ max(AT A) -strongly convex relative to ∥. ∥A.

<!-- page 5 -->
Convergence rate of ADMM 5
However, its converse does not necessarily hold unless A has full column rank.
Hence, the assumption of strong convexity relative to ∥. ∥A for a given matrix A
is weaker compared to the assumption of strong convexity. For fu rther details
on the strong convexity in relation to a given function, we refer the reader
to [29]. We denote the set of c-strongly convex functions relative to ∥. ∥A on
Rn by F A
c (Rn). We denote the distance function to the set X by dX (x) :=
inf y∈ X ∥y − x∥.
In the following sections we derive some new convergence rates for ADMM
by using performance estimation. The main idea of performance est imation
is based on interpolablity. Let I be an index set and let {(xi; gi; f i)}i∈I ⊆
Rn × Rn × R. A set {(xi; ξi; f i)}i∈I is called F A
c -interpolable if there exists
f ∈ F A
c (Rn) with
f (xi) = f i, ξ i ∈ ∂f (xi) i ∈ I .
The next theorem gives necessary and suﬃcient conditions for F A
c -interpolablity.
Theorem 2 Let c ∈ [0, ∞ ) and let I be an index set. The set {(xi; ξi; f i)}i∈I ⊆
Rn × Rn × R is F A
c -interpolable if and only if for any i, j ∈ I , we have
c
2

xi − xj
2
A ≤ f i − f j −
⟨
ξj, x i − xj ⟩
. (4)
Moreover, F0-interpolable and L-smooth if and only if for any i, j ∈ I , we
have
1
2L

gi − gj
2
≤ f i − f j −
⟨
gj, x i − xj ⟩
. (5)
Proof. The argument is analogous to that of [41, Theorem 4]. The triple
{(xi; ξi; f i)}i∈I is F A
c -interpolable if and only if the triple {(xi; ξi− cAT Axi; f i−
c
2 ∥xi∥2
A)}i∈I is F0-interpolable. By [41, Theorem 1], {(xi; ξi − cAT Axi; f i −
c
2 ∥xi∥2
A)}i∈I is F0-interpolable if and only if
f i − c
2

xi
2
A ≥ f j − c
2

xj 
2
A −
⟨
ξj − cAT Axj, x i − xj⟩
which implies inequality (4). The second part follows directly from [41, Theo-
rem 4].
Note that any convex function is 0-strongly convex relative to A. Let f ∈
F0(Rn). The conjugate function f ∗ : Rn → (−∞ , ∞ ] is deﬁned as f ∗ (y) =
supx∈ Rn⟨y, x ⟩ − f (x). We have the following identity
ξ ∈ ∂f (x) ⇔ x ∈ ∂f ∗ (ξ). (6)
Let f ∈ F 0(Rn) be µ -strongly convex. The function f is µ -strongly convex if
and only if f ∗ is 1
µ -smooth. Moreover, ( f ∗ )∗ = f .
By using conjugate functions, the dual of problem (1) may be writt en as
D(λ) = min
(x,z )∈ Rn× Rm
f (x) + g(z) + ⟨λ, Ax + Bz − b⟩
= −⟨ λ, b ⟩ − f ∗ (− AT λ) − g∗ (− BT λ). (7)

<!-- page 6 -->
6 Moslem Zamani et al.
By the optimality conditions for the dual problem, we get
b − Ax⋆ − Bz ⋆ = 0, (8)
for some x⋆ ∈ ∂f ∗ (− AT λ ⋆ ) and z⋆ ∈ ∂g ∗ (− BT λ ⋆ ). Equation (8) with (6)
imply that ( x⋆ , z ⋆) is an optimal solution to problem (1).
The optimality conditions for the subproblems of Algorithm 1 may be
written as
0 ∈ ∂f (xk) + AT λ k− 1 + tAT (
Axk + Bz k− 1 − b
)
,
0 ∈ ∂g (zk) + BT λ k− 1 + tBT (
Axk + Bz k − b
)
. (9)
As λ k = λ k− 1 + t(Axk + Bz k − b), we get
0 ∈ ∂f (xk) + AT λ k + tAT B
(
zk− 1 − zk)
, 0 ∈ ∂g (zk) + BT λ k. (10)
So, (xk, z k) is optimal for dual objective at λ k if and only if AT B
(
zk− 1 − zk)
=
0. We call AT B
(
zk− 1 − zk)
dual residual.
2 Performance estimation
In this section, we develop the performance estimation for ADMM. T he per-
formance estimation method introduced by Drori and Teboulle [9] is a n SDP-
based method for the analysis of ﬁrst order methods. Since then, many scholars
employed this strong tool to derive the worst case convergence r ate of diﬀer-
ent iterative methods; see [2, 23, 38, 41] and the references th erein. Moreover,
Gu and Yang [17] employed performance estimation to study the ext ension of
the dual step length for ADMM. Note that while there are some similar ities
between our work and [17] in using performance estimation, the fo rmulations
and results are diﬀerent.
The worst-case convergence rate of Algorithm 1 with respect to d ual ob-
jective value may be cast as the following abstract optimization prob lem,
max D(λ ⋆ ) − D(λ N )
s. t. {xk, z k, λ k}N
1 is generated by Algorithm 1 w.r.t. f, g, A, B, b, λ 0, z 0, t
(x⋆ , z ⋆ ) is an optimal solution with Lagrangian multipliers λ ⋆
∥λ 0 − λ ⋆ ∥2 + t2 
z0 − z⋆ 
2
B = ∆
f ∈ F A
c1(Rn), g ∈ F B
c2(Rm) (11)
λ 0 ∈ Rr, z 0 ∈ Rm, A ∈ Rr× n, B ∈ Rr× m, b ∈ Rr,
where f, g, A, B, b, z 0, λ 0, x ⋆ , z ⋆, λ ⋆ are decision variables and N, t, c 1, c 2, ∆ are
the given parameters.

<!-- page 7 -->
Convergence rate of ADMM 7
By using Theorem 2 and the optimality conditions (9), problem (11) ma y
be reformulated as the ﬁnite dimensional optimization problem,
max D(λ ⋆ ) − D(λ N )
s. t. {(xk; ξk; f k)}N
1 ∪ { (x⋆ ; ξ⋆ ; f ⋆ )} satisfy interpolation constraints (4)
{(zk; ηk; gk)}N
0 ∪ { (z⋆; η⋆ ; g⋆)} satisfy interpolation constraints (4)
(x⋆ , z ⋆ ) is an optimal solution with Lagrangian multipliers λ ⋆
∥λ 0 − λ ⋆ ∥2 + t2 
z0 − z⋆ 
2
B = ∆ (12)
ξk = tAT b − tAT Axk − tAT Bz k− 1 − AT λ k− 1, k ∈ { 1, ..., N }
ηk = tBT b − tBT Axk − tBT Bz k − BT λ k− 1, k ∈ { 1, ..., N }
λ k = λ k− 1 + t(Axk + Bz k − b), k ∈ { 1, ..., N }
λ 0 ∈ Rr, z 0 ∈ Rm, A ∈ Rr× n, B ∈ Rr× m, b ∈ Rr.
To handle problem (12), without loss of generality, we assume that t he
matrix
(
A B
)
has full row rank. Note this assumption does not employed
in our arguments in the following sections. In addition, we introduce s ome
new variables. As problem (1) is invariant under translation of ( x, z ), we may
assume without loss of generality that b = 0 and ( x⋆ , z ⋆ ) = (0 , 0). In addition,
due to the full row rank of the matrix
(A B )
, we may assume that λ 0 =
(
A B
) (x†
z†
)
and λ ⋆ =
(
A B
) (¯x
¯z
)
for some ¯x, x †, ¯z, z †. So,
ξ⋆ = − AT A¯x − AT B ¯z ∈ ∂f (0), η ⋆ = − BT A¯x − BT B ¯z ∈ ∂g (0),
and D(λ ⋆ ) = f ⋆ + g⋆.
By using equality constraints of problem (12) and the newly introduc ed
variables, we have for k ∈ { 1, ..., N }
λk = (Ax† + Bz †) +
k∑
i=1
t(Axi + Bz i), (13)
− (AT Ax† + AT Bz †) −
k−1∑
i=1
t(AT Axi + AT Bz i) − tAT Axk − tAT Bz k−1 ∈ ∂f (xk),
− (BT Ax† + BT Bz †) −
k∑
i=1
t(BT Axi + BT Bz i) ∈ ∂g(zk).
Note that (˜x, ˜z) ∈ argmin f (x) + g(z) + ⟨λ N , Ax + Bz − b⟩ if and only if
0 ∈ ∂f (˜x) + AT λ N , 0 ∈ ∂g (˜z) + BT λ N . (14)
It is worth noting that a point ˜ x satisfying these conditions exists, as function
f is strongly convex relative to A. In addition, one may consider ˜ z = zN by
virtue of (10). For the sake of notation convenience, we introduc e xN +1 = ˜x.
The reader should bear in mind that xN +1 is not generated by Algorithm 1.

<!-- page 8 -->
8 Moslem Zamani et al.
Therefore, D(λ N ) = f (xN +1) + g(zN ) +
⟨
λ N , Ax N +1 + Bz N ⟩
for some xN +1
with − AT λ N ∈ ∂f (xN +1). Hence, problem (12) may be written as
max f ⋆ + g⋆ − f N+1 − gN −
⟨
Ax† + Bz † +
N∑
i=1
t(Axi + Bz i), AxN+1 + Bz N
⟩
s. t. c1
2


xk − xj



2
A
≤
⟨
Ax† + Bz † +
k−1∑
i=1
t(Axi + Bz i) + tAxk + tBz k−1, A(xj − xk)
⟩
+
f j − f k, k ∈ { 1, . . . , N}, j ∈ { 1, . . . , N + 1},
c1
2


xN+1 − xj



2
A
≤
⟨
Ax† + Bz † +
N∑
i=1
t(Axi + Bz i), A
(
xj − xN+1
) ⟩
+
f j − f N+1, j ∈ { 1, . . . , N},
c2
2


zk − zj



2
B
≤
⟨
Ax† + Bz † +
k∑
i=1
t(Axi + Bz i), B
(
zj − zk
) ⟩
+
gj − gk, j, k ∈ { 1, . . . , N}, (15)
c1
2


xk



2
A
≤ f k − f ⋆ +
⣨
A¯x + B ¯z, Axk
⟩
, k ∈ { 1, . . . , N + 1},
c1
2


xk



2
A
≤ −
⟨
Ax† + Bz † +
k−1∑
i=1
t(Axi + Bz i) + tAxk + tBz k−1, Axk
⟩
+
f ⋆ − f k, k ∈ { 1, . . . , N},
c1
2


xN+1



2
A
≤ f ⋆ − f N+1 −
⟨
Ax† + Bz † +
N∑
i=1
t(Axi + Bz i), AxN+1
⟩
,
c2
2


zk



2
B
≤ gk − g⋆ +
⣨
A¯x + B ¯z, Bz k
⟩
, k ∈ { 1, . . . , N},
c2
2


zk



2
B
≤ g⋆ − gk −
⟨
Ax† + Bz † +
k∑
i=1
t(Axi + Bz i), Bz k
⟩
, k ∈ { 1, . . . , N},


Ax† + Bz † − (A¯x + B ¯z)



2
+ t2 
z0
2
B = ∆,
x† ∈ Rn, z0, z† ∈ Rm, A ∈ Rr×n, B ∈ Rr×m.
In problem (15), A, B, {xk, f k}N +1
1 , {zk, g k}N
1 , x †, z †, ¯x, f ⋆, ¯z, g ⋆ , z 0 are de-
cision variables. By using the Gram matrix method, problem (15) may b e
relaxed as a semideﬁnite program as follows. Let
U =
(
x† x1 . . . x N +1 ¯x
)
, V =
(
z† z0 . . . z N ¯z
)
.
By introducing matrix variable
Y =
(AU BV ) T (AU BV )
,

<!-- page 9 -->
Convergence rate of ADMM 9
problem (15) may be relaxed as the following SDP,
max f ⋆ + g⋆ − f N +1 − gN − tr(LoY )
s. t. tr(Lf
i,j Y ) ≤ f i − f j, i, j ∈ { 1, ..., N + 1, ⋆ }
tr(Lg
i,j Y ) ≤ gi − gj, i, j ∈ { 1, ..., N, ⋆ }
tr(L0Y ) = ∆ (16)
Y ⪰ 0,
where the constant matrices Lf
i,j , L g
i,j , L o, L 0 are determined according to the
constraints of problem (15). In the following sections, we present some new
convergence results that are derived by solving this kind of formula tion.
3 Worst-case convergence rate
In this section, we provide new convergence rates for ADMM with re spect to
some performance measures. Before we get to the theorems we n eed to present
some lemmas.
Lemma 1 Let N ≥ 4 and t, c ∈ R. Let E(t, c ) be (N + 1)× (N + 1) symmetric
matrix given by
E(t, c ) =














2c 0 0 0 . . . 0 0 . . . 0 t − c
0 α 2 β2 0 . . . 0 0 . . . 0 − t
0 β2 α 3 β3 . . . 0 0 . . . 0 t
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
0 0 0 0 . . . α k βk . . . 0 t
.
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
.
0 0 0 0 . . . 0 0 . . . α N βN
t − c − t t t . . . t t . . . β N α N +1














,
where
αk =









6c − 5t, k = 2
2
(
2k2 − 3k + 1
)
c − (4k − 1) t, 3 ≤ k ≤ N − 1
2N (N − 1)c − (2N + 1)t, k = N
2N c − (N + 1)t, k = N + 1,
βk =
{
2kt − (2k2 − k − 1)c, 2 ≤ k ≤ N − 1
3t − 2(N − 1)c, k = N ,
and k denotes row number. If c > 0 is given, then
[0, c ] ⊆ { t : E(t, c ) ⪰ 0}.

<!-- page 10 -->
10 Moslem Zamani et al.
Proof. As {t : E(t, c ) ⪰ 0} is a convex set, it suﬃces to prove the positive
semideﬁniteness of E(0, c ) and E(c, c ). Since E(0, c ) is diagonally dominant,
it is positive semideﬁnite. Now, we establish that the matrix K = E(1, 1)
is positive deﬁnite. To this end, we show that all leading principal minor s
of K are positive. To compute the leading principal minors, we perform th e
following elementary row operations on K:
i) Add the second row to the third row;
ii) Add the second row to the last row;
iii) Add the third row to the forth row;
iv) For i = 4 : N − 1
– Add i − th row to ( i + 1) − th row;
– Add 3− i
2i2− 3i− 1 times of i − th row to the last row;
v) Add N − 1
3N − 5 times of N − th row to ( N + 1) − th row.
It is seen that Kk− 1,k + Kk,k = − Kk+1,k for 2 ≤ k ≤ N − 1. Hence, by per-
forming these operations, we get an upper triangular matrix J with diagonal
Jk,k =









2, k = 1
2k2 − 3k − 1, 2 ≤ k ≤ N − 1
3N − 5, k = N
N − 2 − (N − 1)2
3N − 5 − ∑N − 1
i=4
(i− 3)2
2i2− 3i− 1 , k = N + 1.
It is seen all ﬁrst N diagonal elements of J are positive. We show that JN +1,N +1
is also positive. For i ≥ 4 we have
(i− 3)2
2i2− 3i− 1 ≤ (i− 1)2+4
2(i− 1)2 ≤ 1
2 + 2
(i− 1)(i− 2) . (17)
So,
2N 2− 9N +9
3N − 5 −
N − 1∑
i=4
(i− 3)2
2i2− 3i− 1 ≥ (N − 2)(N 2− 5N +10)
2N (3N − 5) > 0,
which implies JN +1,N +1 > 0. Since we add a factor of i − th row to j − th
row with i < j , all leading principal minors of matrices K and J are the
same. Hence K is positive deﬁnite. As E(c, c ) = cK, one can infer the positive
deﬁniteness of E(c, c ) and the proof is complete.
In the upcoming lemma, we establish a valid inequality for ADMM that
will be utilized in all the subsequent results presented in this section.
Lemma 2 Let f ∈ F A
c1(Rn), g ∈ F 0(Rm) and x⋆ = 0 , z⋆ = 0 . Suppose that
ADMM with the starting points λ 0 and z0 generates {(xk; zk; λ k)}. If N ≥ 4

<!-- page 11 -->
Convergence rate of ADMM 11
and v ∈ Rr, then
N ⟨λN , AxN + Bz N ⟩ − ⟨ λN + tAxN + tBz N−1, AxN − v⟩ + ⟨λ0 + tAx1 + tBz 0, Ax1 − v⟩+
1
2t

λ0 − λ⋆ 
2 − 1
2t


λN − λ⋆



2
+ t
2

z0
2
B − t
⣨
Ax1 − Ax2 + (N + 1)AxN + Bz N , v
⟩
−
t
N∑
k=3
⟨Axk, v⟩ + t(N−1)
2 ∥v∥2 − c1
2

x1
2
A +
N∑
k=2
α k
2


xk



2
A
+
N−1∑
k=2
βk⟨Axk, Axk+1⟩+
tN⟨Bz N−1, AxN − v⟩ + t⟨AxN , BzN ⟩ − t(N−1)2
2


zN − zN−1



2
B
− tN 2
2


AxN + Bz N



2
−
t

x2
2
A + f (x1) − f (xN ) + N
(
f (xN ) − f ⋆ + g(xN ) − g⋆
)
≥ 0, (18)
where
αk =
{
(4k − 1) t − 2
(
2k2 − 3k + 1
)
c1, 2 ≤ k ≤ N − 1,
(4N + 1) t − (2N 2 − 5N + 3) c1, k = N ,
βk = (2k2 − k − 1) c1 − 2kt.
Proof. To establish the desired inequality, we demonstrate its validity by sum -
ming a series of valid inequalities. To simplify the notation, let f k = f (xk)
and gk = g(zk) for k ∈ { 1, . . . , N }. Note that b = 0 because x⋆ = 0 , z ⋆ = 0.
By (4) and (9), we get the following inequality
N−1∑
k=1
(k2 − 1)
(
f k+1 − f k +
⣨
λk−1 + tAxk + tBz k−1, A(xk+1 − xk)
⟩
− c1
2


xk+1 − xk



2
A
)
+
N−1∑
k=1
(k2 − k)
(
f k − f k+1 +
⣨
λk + tAxk+1 + tBz k, A(xk − xk+1)
⟩
− c1
2


xk+1 − xk



2
A
)
+
N∑
k=1
(
f k − f ⋆ +
⣨
λ⋆ , Axk
⟩
− c1
2


xk



2
A
)
+
N−1∑
k=1
k2
(
gk − gk+1 +
⣨
λk+1, B(zk − zk+1)
⟩)
+
N−1∑
k=1
(k2 + k)
(
gk+1 − gk +
⣨
λk, B(zk+1 − zk)
⟩)
+
N∑
k=1
(
gk − g⋆ +
⣨
λ⋆ , Bzk
⟩)
+ t
2

Ax1 + Bz 0 − v

2 ≥ 0.

<!-- page 12 -->
12 Moslem Zamani et al.
As λ k = λ k− 1 + tAxk + tBz k, the inequality can be expressed as
N−1∑
k=1
(k2 − 1)
(⣨
tAxk + tBz k−1, A(xk+1 − xk)
⟩
− c1
2


xk+1 − xk



2
A
)
+
N−1∑
k=1
(k2 − 1)
(⣨
λk, Axk+1
⟩
−
⣨
λk−1, Axk
⟩
−
⣨
tAxk + tBz k, Axk+1
⟩)
+
N−1∑
k=1
(k2 − k)
(⣨
tAxk+1 + tBz k, A(xk − xk+1)
⟩
− c1
2


xk+1 − xk



2
A
)
+
N−1∑
k=1
(k2 − k)
(⣨
λk−1, Axk
⟩
−
⣨
λk, Axk+1
⟩
+
⣨
tAxk + tBz k, Axk
⟩)
+
N−1∑
k=1
(k2 + k)
(⣨
λk, Bz k+1
⟩
−
⣨
λk−1, Bz k
⟩
−
⣨
tAxk + tBz k, Bzk
⟩)
+
N−1∑
k=1
k2
( ⣨
λk−1, Bzk
⟩
−
⣨
λk, Bzk+1
⟩
+
⣨
tAxk + tBz k + tAxk+1 + tBz k+1, Bzk
⟩
−
⣨
tAxk+1 + tBz k+1, Bzk+1
⟩ )
+
N∑
k=1
(⣨
λ⋆ , Axk + Bz k
⟩
− c1
2


xk



2
A
)
+ t
2

Bz 0
2 +
t
2

Ax1 − v

2 + t ⟨Ax1 − v, Bz 0⟩+ f 1 − f N + N (f N − f ⋆ + gN − g⋆ ) ≥ 0.
After performing some algebraic manipulations, we obtain
N ⟨λN−1, AxN + Bz N ⟩ − ⟨ λN−1, AxN ⟩ + ⟨λ0, Ax1⟩ −
N−1∑
k=0
⟨λk − λ⋆ , Axk+1 + Bz k+1⟩+
t
2

Ax1 − v

2 + t
2

Bz 0
2 + t
⟨
Ax1 − v, Bz 0⟩
− t(N 2 − 3N + 1)⟨AxN , Bz N−1⟩−
t
N−1∑
k=1
(
(k − 1)2∥Axk∥2 − (k2 − k)⟨Axk, Axk+1⟩ − (k2 − 1)⟨Axk+1, Bz k−1⟩
)
−
t
N−1∑
k=1
(
(k2 − k + 1)∥Bz k∥2 + (− k2 + k + 1)⟨Axk, Bzk⟩ − k2⟨Bz k, Bzk+1⟩
)
−
t
N−1∑
k=2
(
(2k2 − 3k)⟨Axk, Bzk−1⟩
)
− t(N − 1)2∥Bz N ∥2 − t(N 2 − 3N + 2)∥AxN ∥2−
t(N − 1)2⟨AxN , BzN ⟩ −
N−1∑
k=1
(
(2k2 − k − 1) c1
2


xk+1 − xk



2
A
+ c1
2


xk+1



2
A
)
−
c1
2

x1
2
A + f 1 − f N + N (f N − f ⋆ + gN − g⋆ ) ≥ 0.
By using λ N − 1 = λ N − tAxN − tBz N and
2⟨λ k− λ ⋆ , Ax k+1+Bz k+1⟩ = 1
t ∥λ k+1− λ ⋆ ∥2− 1
t ∥λ k− λ ⋆ ∥2− t∥Axk+1+Bz k+1∥2,

<!-- page 13 -->
Convergence rate of ADMM 13
we get
N ⟨λN , AxN + Bz N ⟩ − ⟨ λN + tAxN + tBz N−1, AxN − v⟩ + ⟨λ0 + tAx1 + tBz 0, Ax1 − v⟩
+ 1
2t

λ0 − λ⋆ 
2 − 1
2t


λN − λ⋆



2
+ t
2

z0
2
B − t
⣨
Ax1 − Ax2 + (N + 1)AxN + Bz N , v
⟩
− t
N∑
k=3
⣨
Axk, v
⟩
− t
2
N−1∑
k=2


(k − 1)Bz k−1 − (k − 1)Bz k + kAxk − (k + 1)Axk+1 + v



2
+ t(N−1)
2 ∥v∥2 − c1
2

x1
2
A − 2t

x2
2
A + 1
2
N−1∑
k=2
(
(4k − 1) t − 2
(
2k2 − 3k + 1
)
c1
) 

xk



2
A
+
N−1∑
k=2
((
2k2 − k − 1
)
c1 − 2kt
)
⟨Axk, Axk+1⟩ +
((
2N + 1
2
)
t −
(
N 2 − 5
2 N + 3
2
)
c1
) 

xN



2
A
+ tN
⣨
Bz N−1, AxN − v
⟩
+ t
⣨
AxN , Bz N
⟩
− t(N−1)2
2


zN − zN−1



2
B
− tN 2
2


AxN + Bz N



2
+ f 1 − f N + N
(
f N − f ⋆ + gN − g⋆
)
≥ 0,
which implies the desired inequality.
We may now prove the main result of this section.
Theorem 3 Let f ∈ F A
c1(Rn) and g ∈ F 0(Rm) with c1 > 0. If t ≤ c1 and
N ≥ 4, then
D(λ ⋆ ) − D(λ N ) ≤ ∥λ 0 − λ ⋆ ∥2 + t2 
z0 − z⋆ 
2
B
4N t . (19)
Proof. As discussed in Section 2, we may assume that x⋆ = 0 and z⋆ = 0. By
(14), we have D(λ N ) = f (ˆxN ) + g(zN ) +
⟨
λ N , A ˆxN + Bz N ⟩
for some ˆxN with
− AT λ N ∈ ∂f (ˆxN ). By employing (4) and (9), we obtain
N
(
g(xN ) − g⋆ + ⟨λ⋆ , Bz N ⟩
)
+ (N − 1)
(
f (xN ) − f ⋆ + ⟨λ⋆ , AxN ⟩ − c1
2


xN



2
A
)
+
(
f (ˆxN ) − f (x1) +
⣨
λ0 + tAx1 + tBz 0, AˆxN − Ax1
⟩
− c1
2


ˆxN − x1



2
A
)
+
(2N − 2)
(
f (ˆxN ) − f (xN ) +
⣨
λN − tBz N + tBz N−1, AˆxN − AxN
⟩
−
c1
2


ˆxN − xN



2
A
)
+
(
f (ˆxN ) − f ⋆ + ⟨λ⋆ , AˆxN ⟩ − c1
2


ˆxN



2
A
)
≥ 0. (20)
By substituting v with AˆxN in inequality (18) and summing it with (20), we
get the following inequality after performing some algebraic manipulat ions
2N
(
f (ˆxN ) + g(xN ) +
⣨
λN , AˆxN + Bz N
⟩
− f ⋆ − g⋆
)
+ 1
2t

λ0 − λ⋆ 
2 + t
2

z0
2
B −
1
2t


λN − λ⋆ + t(N − 1)AxN + tAˆxN + tN BzN



2
−
t
2


(N − 1)(Bz N−1 − Bz N ) + tAxN − tAˆxN



2
− (21)
1
2 tr
(
E(t, c1)
(
Ax1 . . . A ˆxN ) T (
Ax1 . . . A ˆxN ) )
≥ 0,

<!-- page 14 -->
14 Moslem Zamani et al.
where the positive semideﬁnite matrix E(t, c 1) is given in Lemma 1. As the
inner product of positive semideﬁnite matrices is non-negative, ineq uality (21)
implies that
2N
(
D(λ ⋆ ) − D(λ N )
)
≤ 1
2t

λ 0 − λ ⋆ 
2
+ t
2

z0
2
B ,
and the proof is complete.
In comparison with Theorem 1, we could get a new convergence rate when
only f is strongly convex, i.e. g does not need to be strongly convex. Also, the
constant does not depend on λ 1. One important question concerning bound
(19) is its tightness, that is, if there is an optimization problem which a ttains
the given convergence rate. It turns out that the bound (19) is e xact. The
following example demonstrates this point.
Example 1 Suppose that c1 > 0, N ≥ 4 and t ∈ (0, c 1]. Let f, g : R → R be
given as follows,
f (x) = 1
2 |x|+ c1
2 x2, g (z) = 1
2 max{ N − 1
N (z − 1
2N t) − 1
2N t, − z}.
Consider the optimization problem
min
(x,z )∈ R× R
f (x) + g(z),
s. t. x + z = 0,
It is seen that A = B = I in this problem. Note that ( x⋆ , z ⋆) = (0 , 0) with
Lagrangian multiplier λ ⋆ = 1
2 is an optimal solution and the optimal value is
zero. One can check that Algorithm 1 with initial point λ 0 = − 1
2 and z0 = 0
generates the following points,
xk = 0 k ∈ { 1, . . . , N }
zk = 1
2N t k ∈ { 1, . . . , N }
λ k = − 1
2 + k
2N k ∈ { 1, . . . , N }.
At λ N , we have D(λ N ) = − 1
4N t = −
∥λ 0− λ ⋆ ∥2+t2∥z0− z⋆ ∥
2
B
4N t , which shows the
tightness of bound (19).
One important factor concerning dual-based methods that deter mines the
eﬃciency of an algorithm is primal and dual feasibility (residual) conve rgence
rates. In what follows, we study this subject under the setting of Theorem 3.
The next theorem gives a convergence rate in terms of primal resid ual under
the setting of Theorem 3.
Theorem 4 Let f ∈ F A
c1(Rn) and g ∈ F 0(Rm) with c1 > 0. If t ≤ c1 and
N ≥ 4, then

AxN + Bz N − b

 ≤
√
∥λ 0 − λ ⋆ ∥2 + t2 ∥z0 − z⋆ ∥2
B
tN . (22)

<!-- page 15 -->
Convergence rate of ADMM 15
Proof. The argument is similar to that used in the proof of Theorem 3. By
setting v = AxN in (18), one can infer the following inequality
N
⣨
λN , AxN + Bz N
⟩
+
⣨
λ0 + tAx1 + tBz 0, Ax1 − AxN
⟩
+ 1
2t

λ0 − λ⋆ 
2 + t
2

z0
2
B −
t
⣨
Ax1 − Ax2, AxN
⟩
+ t(N−1)
2


AxN



2
− t
N∑
k=3
⣨
Axk, AxN
⟩
− c1
2

x1
2
A − t

x2
2
A +
N−1∑
k=2
((
2k − 1
2
)
t −
(
2k2 − 3k + 1
)
c1
) 

xk



2
A
+
((3
2 N − 3
2
)
t −
(
N 2 − 5
2 N + 3
2
)
c1
) 

xN



2
A
+
N−1∑
k=2
((2k2 − k − 1) c1 − 2kt) ⟨Axk, Axk+1⟩ − t(N−1)2
2


zN − zN−1



2
B
−
tN 2
2


AxN + Bz N



2
+ f (x1) − f (xN ) + N
(
f (xN ) − f ⋆ + g(xN ) − g⋆
)
≥ 0. (23)
By employing (4) and (9), we have
N
(
f ⋆ − f (xN ) − ⟨ λN + Bz N−1 − Bz N , AxN ⟩ − c1
2


xN



2
A
)
+
(
f (xN ) − f 1 +
⣨
λ0 + tAx1 + tBz 0, AxN − Ax1
⟩
− c1
2


xN − x1



2
A
)
+ (24)
N
(
g⋆ − g(xN ) − ⟨ λN , Bz N ⟩
)
≥ 0.
By summing (23) and (24), we obtain
1
2t

λ0 − λ⋆ 
2 + t
2

z0
2
B − t(N−1)2
2


zN−1 − zN + N
(N−1)2 xN



2
B
−
tN 2
2


AxN + Bz N



2
− 1
2 tr
(
D(t, c1) (Ax1 . . . Ax N ) T (Ax1 . . . Ax N ) )
≥ 0, (25)
where the matrix D(t, c 1) is as follows,
D(t, c 1) =














2c1 0 0 0 . . . 0 0 . . . 0 t − c1
0 α 2 β2 0 . . . 0 0 . . . 0 − t
0 β2 α 3 β3 . . . 0 0 . . . 0 t
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
0 0 0 0 . . . α k βk . . . 0 t
.
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
.
0 0 0 0 . . . 0 0 . . . α N − 1 βN − 1
t − c1 − t t t . . . t t . . . β N − 1 α N














,
and
α k =







6c1 − 5t, k = 2
2
(
2k2 − 3k + 1
)
c1 − (4k − 1) t, 3 ≤ k ≤ N − 1,(
2N 2 − 4N + 4
)
c1 −
(
3N − 5 + N 2
(N − 1)2
)
t, k = N ,
βk = 2kt −
(
2k2 − k − 1
)
c1, 2 ≤ k ≤ N − 1

<!-- page 16 -->
16 Moslem Zamani et al.
As the matrix D(t, c 1) is positive semideﬁnite, see Appendix A, inequality (25)
implies that
tN 2
2

AxN + Bz N 
2
≤ 1
2t

λ 0 − λ ⋆ 
2
+ t
2

z0
2
B ,
and the proof is complete.
The following example shows the exactness of bound (22).
Example 2 Let c1 > 0, N ≥ 4 and t ∈ (0, c 1]. Consider functions f, g : R → R
given by the formulae follows,
f (x) = 1
2 |x|+ c1
2 x2,
g(z) = max {
(1
2 − 1
N
) (
z − 1
N t
)
, 1
2
( 1
N t − z
)
}.
We formulate the following optimization problem,
min
(x,z )∈ R× R
f (x) + g(z),
s. t. Ax + Bz = 0,
where A = B = I. One can verify that ( x⋆ , z ⋆) = (0 , 0) with Lagrangian
multiplier λ ⋆ = 1
2 is an optimal solution. Algorithm 1 with initial point λ 0 =
− 1
2 and z0 = 0 generates the following points,
xk = 0 k ∈ { 1, . . . , N }
zk = 1
N t k ∈ { 1, . . . , N }
λ k = 2k− N
2N k ∈ { 1, . . . , N }.
At iteration N , we have ∥AxN + Bz N ∥ = 1
tN =
√
∥λ 0− λ ⋆ ∥2+t2∥z0− z⋆ ∥2
B
tN , which
shows the tightness of bound (22).
In what follows, we study the convergence rate of ADMM in terms re sidual
dual. To this end, we investigate the convergence rate of {B
(
zk− 1 − zk)
}
as

AT B
(
zk− 1 − zk) 
 ≤ ∥ A∥

zk− 1 − zk

B. The next theorem provides a
convergence rate for the aforementioned sequence.
Theorem 5 Let f ∈ F A
c1(Rn) and g ∈ F 0(Rm) with c1 > 0. If t ≤ c1 and
N ≥ 4, then

zN − zN − 1

B ≤
√
∥λ 0 − λ ⋆ ∥2 + t2 ∥z0 − z⋆∥2
B
(N − 1)t . (26)

<!-- page 17 -->
Convergence rate of ADMM 17
Proof. Similar to the proof of Theorem 3, by setting v = AxN in (18) for N − 1
iterations, one can infer the following inequality
(N − 1)⟨λ N − 1, Ax N − 1 + Bz N − 1⟩ + 1
2t ∥λ 0 − λ ⋆ ∥2 − 1
2t ∥λ N − 1 − λ ⋆ ∥2+
t
2

z0
2
B − ⟨ λ N − 1 + tAxN − 1 + tBz N − 2, Ax N − 1 − AxN ⟩ + t(N − 2)
2 ∥xN ∥2
A+
⟨λ 0 + tAx1 + tBz 0, Ax 1 − AxN ⟩ − t
⟨
Ax1 − Ax2 + N AxN − 1 + Bz N − 1, Ax N ⟩
+ 1
2
N − 2∑
k=2
(
(4k − 1) t − 2
(
2k2 − 3k + 1
)
c1
) 
xk
2
A + t⟨AxN − 1, Bz N − 1⟩+
N − 2∑
k=2
((
2k2 − k − 1
)
c1 − 2kt
)
⟨Axk, Ax k+1⟩ + t(N − 1)⟨Bz N − 2, Ax N − 1 − AxN ⟩
+ 1
2
(
(4N − 3) t −
(
2N 2 − 9N + 10
)
c1
) 
xN − 1
2
A − t

x2
2
A − c1
2

x1
2
A −
t(N − 2)2
2

zN − 1 − zN − 2
2
B − t(N − 1)2
2 ∥AxN − 1 + Bz N − 1∥2 − t
N − 1∑
k=3
⟨Axk, Ax N ⟩+
f (x1) − f (xN − 1) + (N − 1)(f (xN − 1) − f ⋆ + g(xN − 1) − g⋆ ) ≥ 0. (27)
By using (4) and (9), we have
(N 2 − 3N + 2)
(
f (xN−1) − f (xN ) +
⣨
λN−1 + tAxN + tBz N−1, A
(
xN−1 − xN
)⟩
−
c1
2


xN − xN−1



2
A
)
+
(
f (xN ) − f (x1) +
⟨
λ0 + tAx1 + tBz 0, A
(
xN − x1
) ⟩
−
c1
2 ∥xN − x1∥2
A
)
+ N (N − 1)
(
g(zN ) − g(zN−1) +
⣨
λN−1, B
(
zN − zN−1
)⟩)
+ (28)
(N 2 − 3N + 1)
(
f (xN ) − f (xN−1) +
⟨
λN−1 − tBz N−1 + tBz N−2, A
(
xN − xN−1
) ⟩
− c1
2 ∥xN − xN−1∥2
A
)
+ (N − 1)
(
g⋆ − g(zN ) −
⣨
λN−1 + tAxN + tBz N , Bz N
⟩)
+
(N − 1)
(
f ⋆ − f (xN−1) −
⟨
λN−1 − tBz N−1 + tBz N−2, AxN−1
⟩
− c1
2 ∥xN−1∥2
A
)
+
(N − 1)2
(
g(zN−1) − g(zN ) +
⣨
λN−1 + tAxN + Bz N , B
(
zN−1 − zN
)⟩)
≥ 0.
By summing (27) and (28), we obtain
1
2t

λ0 − λ⋆ 
2 + t
2

z0
2
B − (N 2−1)t
2


 N
N+1 AxN + Bz N



2
− t(N − 1)2
2


zN − zN−1



2
B
−
(N−2)2t
2



Bz N−2 − Bz N−1 + N−1
N−2 AxN−1 −
(
1 − 1
(N−2)2
)
AxN




2
−
1
2 tr
(
F (t, c1)
(
Ax1 . . . Ax N ) T (
Ax1 . . . Ax N ) )
≥ 0,

<!-- page 18 -->
18 Moslem Zamani et al.
where the matrix F (t, c 1) is as follows,
F (t, c 1) =














2c1 0 0 0 . . . 0 0 . . . 0 t − c1
0 α 2 β2 0 . . . 0 0 . . . 0 − t
0 β2 α 3 β3 . . . 0 0 . . . 0 t
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
0 0 0 0 . . . α k βk . . . 0 t
.
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
. .
.
.
0 0 0 0 . . . 0 0 . . . α N − 1 βN − 1
t − c1 − t t t . . . t t . . . β N − 1 α N














,
and
α k =







6c1 − 5t, k = 2
2
(
2k2 − 3k + 1
)
c1 − (4k − 1) t, 3 ≤ k ≤ N − 1,(
2N 2 − 6N + 4
)
c1 − 2
(
N + 1
(N − 2)2 − 2
N +1 − 3
)
t, k = N ,
βk =
{
2kt −
(
2k2 − k − 1
)
c1, 2 ≤ k ≤ N − 2,
(N + 1
2− N − 1)t − (2N 2 − 6N + 3)c1, k = N − 1,
The rest of the proof proceeds analogously to the proof of Theor em 4.
The following example shows the tightness of this bound.
Example 3 Assume that c1 > 0, N ≥ 4 and t ∈ (0, c 1] are given, and f, g :
R → R are deﬁned by,
f (x) = 1
2 max
{
− N +1
N − 1 x, x
}
+ c1
2 x2,
g(z) = 1
2 max
{
1
t(N − 1) − z, N − 3
N − 1
(
z − 1
t(N − 1)
)}
.
Consider the optimization problem
min
(x,z )∈ R× R
f (x) + g(z),
s. t. Ax + Bz = 0.
where A = B = I. The point ( x⋆ , z ⋆ ) = (0 , 0) with Lagrangian multiplier
λ ⋆ = 1
2 is an optimal solution. After performing N iterations of Algorithm 1
with setting λ 0 = − 1
2 and z0 = 0, we have
xk = 0, k ∈ { 1, . . . , N },
zk =
{
1
t(N − 1) , k ∈ { 1, . . . , N − 1},
0, k = N,
λ k =
{
2k+1− N
2(N − 1) , k ∈ { 1, . . . , N − 1},
1
2 , k = N.

<!-- page 19 -->
Convergence rate of ADMM 19
It can be seen that

AT B
(
zN − zN − 1) 
 = 1
(N − 1)t =
√
∥λ 0− λ ⋆ ∥2+t2∥z0− z⋆ ∥2
B
(N − 1)t ,
which shows that the bound is tight.
Theorem 3 and 4 address the case that f is strongly convex relative to ∥. ∥A
and g is convex. Based on numerical results by solving performance estim ation
problems including (15) we conjecture, under the assumptions of T heorem 3,
if g is c2-strongly convex relative to ∥. ∥B, Algorithm 1 enjoys the following
convergence rates
D(λ ⋆ ) − D(λ N ) ≤ ∥λ 0 − λ ⋆ ∥2 + t2∥z0 − z⋆ ∥2
B
4N t + 2c1c2
c1+c2
,

AxN + Bz N − b

 ≤
√
∥λ 0 − λ ⋆ ∥2 + t2∥z0 − z⋆ ∥2
B
N t + c1c2
c1+c2
.
We have veriﬁed these conjectures numerically for many speciﬁc va lues of the
parameters.
4 Linear convergence of ADMM
In this section we study the linear convergence of ADMM. The linear c on-
vergence of ADMM has been addressed by some authors and some c onditions
for linear convergence have been proposed, see [8, 18, 19, 22, 2 7, 34, 42]. Two
common types of assumptions employed for proving the linear conve rgence of
ADMM are error bound property and L-smoothness. To the best knowledge
of authors, most scholars investigated the linear convergence of the sequence
{(xk, z k, λ k)} to a saddle point and there is no result in terms of dual objec-
tive value for ADMM. In line with the previous section, we study the line ar
convergence in terms of dual objective value and we derive some fo rmulas for
linear convergence rate by using performance estimation. It is not eworthy to
mention that the term ”Q-linear convergence” is also employed to de scribe the
linear convergence in the literature.
As mentioned earlier, error bound property is used by scholars for estab-
lishing the linear convergence; see e.g. [18, 22, 27, 35, 42]. Let
Da(λ) := min f (x) + g(z) + ⟨λ, Ax + Bz − b⟩ + a
2 ∥Ax + Bz − b∥2, (29)
stands for augmented dual objective for the given a > 0 and Λ ⋆ denotes the
optimal solution set of the dual problem. Note that function Da is an 1
a -smooth
function on its domain without assuming strong convexity; see [22, Lemma
2.2].
Deﬁnition 2 The function Da satisﬁes the error bound if we have
dΛ ⋆ (λ) ≤ τ∥∇ Da(λ)∥, λ ∈ Rr, (30)
for some τ > 0.

<!-- page 20 -->
20 Moslem Zamani et al.
Hong et al. [22] established the linear convergence by employing erro r
bound property (30).
Recently, some scholars established the linear convergence of gradient meth-
ods for L-smooth convex functions by replacing strong convexity with some
mild conditions, see [1, 5, 32] and references therein. Inspired by t hese results,
we prove the linear convergence of ADMM by using the so-called P/suppress L inequality.
Concerning diﬀerentiability of dual objective, by (7), we have
b − A∂f ∗ (− AT λ) − B∂g ∗ (− BT λ) ⊆ ∂ (− D(λ)) . (31)
Note that inclusion (31) holds as an equality under some mild conditions , see
e.g. [3, Chapter 3].
Deﬁnition 3 The function D is said to satisfy the P/suppress L inequality if there exists
an Lp > 0 such that for any λ ∈ Rr we have
D(λ ⋆ ) − D(λ) ≤ 1
2Lp
∥ξ∥2, ξ ∈ b − A∂f ∗ (− AT λ) − B∂g ∗(− BT λ). (32)
Note that if f and g are strongly convex, then − D is an L-smooth convex
function with L ≤ λ max(AT A)
µ 1
+ λ max(BT B)
µ 2
. Under this setting, we have Lp ≤
λ max(AT A)
µ 1
+ λ max(BT B)
µ 2
. This follows from the duality between smoothness and
strong convexity and
∥∇ D(λ) − ∇ D(ν)∥ ≤


∇ f ∗(− AT λ) − ∇ f ∗(− AT ν)



A
+


∇ g∗(− BT λ) − ∇ g∗(− BT ν)



B
≤ 1
µ 1


AT λ − AT ν



A
+ 1
µ 2


BT λ − BT ν



B
≤
(
λ max(AT A)
µ 1
+ λ max(BT B)
µ 2
)
∥λ − ν∥ .
In the next proposition, we show that deﬁnitions (30) and (32) are equiv-
alent.
Proposition 1 Let La = 1
a denote the Lipschitz constant of ∇ Da, where Da
is given in (29).
i) If Da satisﬁes the error bound (30), then D satisﬁes the P/suppress L inequality
with Lp = 1
Laτ 2 .
ii) If D satisﬁes the P/suppress L inequality, thenDa satisﬁes the error bound (30)
with τ = Lp
1+aLp
.
Proof. First we prove i). Suppose λ ∈ Rr and ξ ∈ b − A∂f ∗(− AT λ) −
B∂g ∗ (− BT λ). By identity (6), we have ξ = b − A¯x − B ¯z for some (¯x, ¯z) ∈
argmin f (x) + g(z) + ⟨λ, Ax + Bz − b⟩. Due to the smoothness of Da and (30),
we get
Da(λ ⋆ ) − Da(ν) ≤ Laτ 2
2 ∥∇ Da(ν)∥2, ν ∈ Rr. (33)
where λ ⋆ ∈ Λ ⋆ with dΛ ⋆ = ∥ν − λ ⋆ ∥. Suppose that ¯ν = λ − a(A¯x + B ¯z − b).
As we assume strong duality, we have Da(λ ⋆ ) = D(λ ⋆ ). By the deﬁnitions of
¯x, ¯y, we get
(¯x, ¯z) ∈ argmin f (x) + g(z) + ⟨¯ν, Ax + Bz − b⟩ + a
2 ∥Ax + Bz − b∥2.

<!-- page 21 -->
Convergence rate of ADMM 21
By [22, Lemma 2.1], we have ∇ Da(¯ν) = A¯x + B ¯z − b. This equality with (33)
imply
D(λ ⋆ ) − D(λ) ≤ Da(λ ⋆ ) − Da(¯ν) ≤ Laτ 2
2 ∥A¯x + B ¯z − b∥2,
and the proof of i) is complete.
Now we establish ii). Let λ be in the domain of ∇ Da. By [22, Lemma 2.1], we
have ∇ Da(λ) = A¯x + B ¯z − b for some (¯x, ¯z) ∈ argmin f (x) + g(z) + ⟨λ, Ax +
Bz − b⟩ + a
2 ∥Ax + Bz − b∥2, which implies that
0 ∈ ∂f (¯x) + AT (λ + a(A¯x + B ¯z − b)) , 0 ∈ ∂g (¯z) + BT (λ + a(A¯x + B ¯z − b)) .
(34)
Supposing ν = λ + a(A¯x + B ¯z − b). By (34), one can infer that D(ν) =
f (¯x) + g(¯z) + ⟨ν, A ¯x + B ¯z − b⟩. In addition, (6) implies that b − A¯x − B ¯z ∈
b − A∂f ∗ (− AT ν) − B∂g ∗ (− BT ν). By the P/suppress L inequality, we have
1
2Lp
∥A¯x + B ¯z − b∥2 ≥ D(λ ⋆ ) − D(ν) = Da(λ ⋆ ) − Da(λ) − a
2 ∥A¯x + B ¯z − b∥2 ,
where the equality follows from D(ν) = Da(λ)+ a
2 ∥A¯x + B ¯z − b∥2 and Da(λ ⋆ ) =
D(λ ⋆ ). Hence,
Da(λ ⋆ ) − Da(λ) ≤
(
1
2Lp
+ a
2
)
∥∇ Da(λ)∥2.
This inequality says that Da satisﬁes the P/suppress L inequality. On the other hand,
the P/suppress L inequality implies the error bound with the same constant, see [5], and
the proof is complete.
In what follows, we employ performance estimation to derive a linear c on-
vergence rate for ADMM in terms of dual objective when the P/suppress L inequality
holds. To this end, we compare the value of dual problem in two conse cutive
iterations, that is, D(λ ⋆ )− D(λ 2)
D(λ ⋆ )− D(λ 1) . The following optimization problem gives the
worst-case convergence rate,
max D(λ ⋆ )− D(λ 2)
D(λ ⋆ )− D(λ 1)
s. t. {x2, z 2, λ 2} is generated by Algorithm 1 w.r.t. f, g, A, B, b, λ 1, z 1 (35)
(x⋆ , z ⋆) is an optimal solution and its Lagrangian multipliers is λ ⋆
D satisﬁes the P/suppress L inequality
f ∈ F A
c1(Rn), g ∈ F B
c2(Rn)
λ 1 ∈ Rr, z 1 ∈ Rm, A ∈ Rr× n, B ∈ Rr× m, b ∈ Rr.
Analogous to our discussion in Section 2, we may assume without loss o f
generality b = 0, λ 1 =
(A B ) (x†
z†
)
and λ ⋆ =
(A B ) (¯x
¯z
)
for some ¯x, x †, ¯z, z †.
In addition, we assume that ˆx1 ∈ argmin f (x)+⟨λ 1, Ax ⟩ and ˆx2 ∈ argmin f (x)+
⟨λ 2, Ax ⟩. Hence,
D(λ 1) = f (ˆx1)+g(z1)+⟨λ 1, A ˆx1+Bz 1⟩, D (λ 2) = f (ˆx2)+g(z2)+⟨λ 2, A ˆx2+Bz 2⟩,

<!-- page 22 -->
22 Moslem Zamani et al.
and
− AT λ 1 ∈ ∂f (ˆx1), − BT λ 1 ∈ ∂g (z1), (36)
− AT λ 2 ∈ ∂f (ˆx2), − BT λ 2 ∈ ∂g (z2).
Moreover, by (36) and (31), we get
− Aˆx1 − Bz 1 ∈ ∂
(
− D(λ 1)
)
, − Aˆx2 − Bz 2 ∈ ∂
(
− D(λ 2)
)
.
On the other hand, λ 2 = λ 1 + tAx2 + tBz 2. Therefore, by using Theorem 2,
problem (35) may be relaxed as follows,
max f ⋆ + g⋆ − ˆf 2 − g2 − ⟨ Ax† + Bz † + tAx2 + tBz 2, Aˆx2 + Bz 2⟩
f ⋆ + g⋆ − ˆf 1 − g1 − ⟨ Ax† + Bz †, Aˆx1 + Bz 1⟩
s. t.
{ (
ˆx1, − AT Ax† − AT Bz †, ˆf 1
)
,
(
x2, − AT Ax† − AT Bz † − tAT Ax2 − tAT Bz 1, f 2
)
,
(
ˆx2, − AT Ax† − AT Bz † − tAT Ax2 − tAT Bz 2, ˆf 2
)
,
(
0, − AT A¯x − AT B ¯z, f ⋆
) }
satisfy interpolation constraints (4)
{ (
z1, − BT Ax† − BT Bz †, g1
)
,
(
z2, − BT Ax† − BT Bz † − tBT Ax2 − tBT Bz 2, g2
)
,
(
0, − BT A¯z − BT B ¯z, g⋆
) }
satisfy interpolation constraints (4)
f ∗ + g∗ − ˆf 1 − g1 −
⣨
Ax† + Bz †, Aˆx1 + Bz 1
⟩
≤ 1
2Lp

Aˆx1 + Bz 1
2 (37)
f ∗ + g∗ − ˆf 2 − g2 −
⣨
Ax† + Bz † + tAx2 + tBz 2, Aˆx2 + Bz 2
⟩
≤ 1
2Lp

Aˆx2 + Bz 2
2
A ∈ Rr×n, B ∈ Rr×m.
By deriving an upper bound for the optimal value of problem (37) in th e next
theorem, we establish the linear convergence of ADMM in the presen ce of the
P/suppress L inequality.
Theorem 6 Let f ∈ F A
c1(Rn) and g ∈ F B
c2(Rm) with c1, c 2 > 0, and let D
satisﬁes the P/suppress L inequality withLp. Suppose that t ≤ √
c1c2.
(i) If c1 ≥ c2, then
D(λ ⋆ ) − D(λ 2)
D(λ ⋆ ) − D(λ 1) ≤ 2c1c2 − t2
2c1c2 − t2 + Lpt (4c1c2 − c2t − 2t2) , (38)
in particular, if t = √ c1c2,
D(λ ⋆ ) − D(λ 2)
D(λ ⋆ ) − D(λ 1) ≤ 1
1 + Lp
(
2√ c1c2 − c2
) .
(ii) If c1 < c 2, then
D(λ ⋆ ) − D(λ 2)
D(λ ⋆ ) − D(λ 1) ≤ (39)
4c2
2 − 2c2
√
c1c2 − t2
4c2
2 − 2c2
√
c1c2 − t2 + Lpt
(
8c2
2 + 5c2t − 2√
c1c2
(
1 + t
c1
)
(2c2 + t)
) .

<!-- page 23 -->
Convergence rate of ADMM 23
Proof. The argument is based on weak duality. Indeed, by introducing suita ble
Lagrangian multipliers, we establish that the given convergence rat es are upper
bounds for problem (37). First, we prove ( i). Assume that α denotes the right
hand side of inequality (38). As 2 c1c2 − t2 > 0 and 4 c1c2 − c2t − 2t2 > 0, we
have 0 < α < 1. With some algebra, one can show that
f ⋆ + g⋆ − ˆf 2 − g2 − ⟨ Ax† + Bz † + tAx2 + tBz 2, Aˆx2 + Bz 2⟩−
α
(
f ⋆ + g⋆ − ˆf 1 − g1 − ⟨ Ax† + Bz †, Aˆx1 + Bz 1⟩
)
+
α
(
ˆf 2 − ˆf 1 + ⟨Ax† + Bz †, Aˆx2 − Aˆx1⟩ − c1
2

ˆx2 − ˆx1
2
A
)
+
α
(
f 2 − ˆf 2 + ⟨Ax† + Bz † + tAx2 + tBz 2, Ax2 − Aˆx2⟩ − c1
2

x2 − ˆx2
2
A
)
+
α
(
ˆf 2 − f 2 + ⟨Ax† + Bz † + tAx2 + tBz 1, Aˆx2 − Ax2⟩ − c1
2

ˆx2 − x2
2
A
)
+
α
(
g2 − g1 + ⟨Ax† + Bz †, Bz2 − Bz 1⟩ − c2
2

z2 − z1
2
B
)
+
(1 − α)
(
− f ⋆ − g⋆ +
⣨
Ax† + Bz † + tAx2 + tBz 2, Aˆx2 + Bz 2
⟩
+ ˆf 2 + g2+
1
2Lp

Aˆx2 + Bz 2
2 )
= −c1α
2

ˆx1 − ˆx2
2
A − c2α
2


Bz 1 − Bz 2 + t
c2
Ax2 − t
c2
Aˆx2



2
−
α(c1 − t2
2c2
)


Ax2 + tc2
2c1c2−t2 Bz 2 − tc2−2c1c2+t2
t2−2c1c2
Aˆx2



2
.
Hence, we get
f ⋆ + g⋆ − ˆf 2 − g2 − ⟨ Ax† + Bz † + tAx2 + tBz 2, A ˆx2 + Bz 2⟩ ≤
α
(
f ⋆ + g⋆ − ˆf 1 − g1 − ⟨ Ax† + Bz †, A ˆx1 + Bz 1⟩
)
for any feasible point of problem (35) and the proof of the ﬁrst par t is complete.
For (ii), we proceed analogously to the proof of ( i), but with diﬀerent Lagrange
multipliers. Let β denote the right hand side of inequality (39), i.e.
β = 4c2
2 − 2c2
√
c1c2 − t2
4c2
2 − 2c2
√
c1c2 − t2 + Lpt
(
8c2
2 + 5c2t − 2√
c1c2
(
1 + t
c1
)
(2c2 + t)
) .
It is seen that 0 < β < 1. By doing some calculus, we have

<!-- page 24 -->
24 Moslem Zamani et al.
f ⋆ + g⋆ − ˆf 2 − g2 − ⟨ Ax† + Bz † + tAx2 + tBz 2, Aˆx2 + Bz 2⟩−
β
(
f ⋆ + g⋆ − ˆf 1 − g1 − ⟨ Ax† + Bz †, Aˆx1 + Bz 1⟩
)
+
β
(
ˆf 2 − ˆf 1 + ⟨Ax† + Bz †, Aˆx2 − Aˆx1⟩ − c1
2

ˆx2 − ˆx1
2
A
)
+
√
c2
c1
β
(
f 2 − ˆf 2 + ⟨Ax† + Bz † + tAx2 + tBz 2, Ax2 − Aˆx2⟩ − c1
2

x2 − ˆx2
2
A
)
+
√
c2
c1
β
(
ˆf 2 − f 2 + ⟨Ax† + Bz † + tAx2 + tBz 1, Aˆx2 − Ax2⟩ − c1
2

ˆx2 − x2
2
A
)
+
√
c2
c1
β
(
g2 − g1 + ⟨Ax† + Bz †, Bz2 − Bz 1⟩ − c2
2

z2 − z1
2
B
)
+
(√
c2
c1
− 1
)
β
(
g1 − g2 + ⟨Ax† + Bz † + tAx2 + tBz 2, Bz1 − Bz 2⟩−
c2
2

z1 − z2
2
B
)
+ (1 − β)
(
− f ⋆ − g⋆ +
⣨
Ax† + Bz † + tAx2 + tBz 2, Aˆx2 + Bz 2
⟩
+
ˆf 2 + g2 + 1
2Lp

Aˆx2 + Bz 2
2 )
= − c1β
2

ˆx1 − ˆx2
2
A − (√
c1c2β)



Ax2 −
(
1 − t
2√ c1c2
)
Aˆx2 + t
2√ c1c2
Bz 1




2
−
(β − 1
2Lp
+ βt
(
1 − t
4√ c1c2
)) 


Aˆx2 −
(
βLp
(
− 2c2
√
c1c2 + 4c2
2 − t2)
− βLpt2 + 2√ c1c2(2βLpt + β − 1)
) 1
2
Bz 1+
(2 (2βc2Lp (t + c2) + √ c1c2 (β − βLpc2 − 1))
− βLpt2 + 2√ c1c2(2βLpt + β − 1)
) 1
2
Bz 2




2
.
The rest of the proof is similar to that of the former case.
We computed the bounds in Theorem 6 by selecting suitable Lagrangia n
multipliers and solving the semideﬁnite formulation of problem (37) by h and.
The semideﬁnite formulation is formed analogous to problem (16). No te that
the optimal value of problem (37) may be smaller than the bounds intr oduced
in Theorem 6. Indeed, our aim was to provide a concrete mathematic al proof
for the linear convergence rate. However, the linear convergenc e rate factor
is not necessarily tight. Needless to say that the optimal value of pr oblem
(37) also does not necessarily give the tight convergence factor a s it is just a
relaxation of problem (35).
Recently the authors showed that the P/suppress L inequality is necessary and suf-
ﬁcient conditions for the linear convergence of the gradient metho d with con-
stant step lengths for L-smooth function; see [1, Theorem 5]. In what follows,
we establish that the P/suppress L inequality is a necessary condition for the linear con-
vergence of ADMM. Firstly, we present a lemma that is very useful f or our
proof.
Lemma 3 Let f ∈ F A
c1(Rn) and g ∈ F B
c2(Rm). Consider Algorithm 1. If
(ˆx1, z 1) ∈ argmin f (x) + g(z) + ⟨λ 1, Ax + Bz − b⟩, then
⟨Aˆx1 + Bz 1 − b, Ax 2 + Bz 2 − b⟩ ≤

Aˆx1 + Bz 1 − b

2
. (40)

<!-- page 25 -->
Convergence rate of ADMM 25
Proof. Without loss of generality we assume that c1 = c2 = 0. By optimality
conditions, we have
f (ˆx1) − ⟨ λ 1, Ax 2 − Aˆx1⟩ ≤ f (x2), g (z1) − ⟨ λ 1, Bz 2 − Bz 1⟩ ≤ g(z2),
f (x2) − ⟨ λ 1 + t(Ax2 + Bz 1 − b), A ˆx1 − Ax2⟩ ≤ f (ˆx1),
g(z2) − ⟨ λ 1 + t(Ax2 + Bz 2 − b), Bz 1 − Bz 2⟩ ≤ g(z1).
By using these inequities, we get
0 ≤ 1
t
(
f (x2) − f (ˆx1) +
⟨
λ1, Ax2 − Aˆx1⟩)
+ 1
t
(
g(z2) − g(z1) +
⟨
λ1, Bz2 − Bz 1⟩)
+
1
t
(
f (ˆx1) − f (x2) +
⟨
λ1 + t(Ax2 + Bz 1 − b), Aˆx1 − Ax2⟩)
+
1
t
(
g(z1) − g(z2) +
⟨
λ1 + t(Ax2 + Bz 2 − b), Bz 1 − Bz 2⟩)
=

Aˆx1 + Bz 1 − b

2 − ⟨Aˆx1 + Bz 1 − b, Ax2 + Bz 2 − b⟩− 3
4

B (z1 − z2) 
2 −

A
(
ˆx1 − x2)
+ 1
2 B
(
z1 − z2) 
2
.
Hence, we have
⟨Aˆx1 + Bz 1 − b, Ax 2 + Bz 2 − b⟩
∥Aˆx1 + Bz 1 − b∥2 ≤ 1,
which completes the proof.
The next theorem establishes that the P/suppress L inequality is a necessarycondi-
tion for the linear convergence of ADMM.
Theorem 7 Let f ∈ F A
c1(Rn) and g ∈ F B
c2(Rm). If Algorithm 1 is linearly
convergent with respect to the dual objective value, then D satisﬁes the P/suppress L
inequality.
Proof. Consider λ 1 ∈ Rr and ξ ∈ b − A∂f ∗ (− AT λ 1) − B∂g ∗(− BT λ 1). Hence,
ξ = b − Aˆx1 − Bz 1 for some (ˆx1, z 1) ∈ argmin f (x) +g(z) +⟨λ, Ax + Bz − b⟩. If
one sets z0 = z1 and λ 0 = λ 1 − t(Aˆx1 + Bz 1 − b) in Algorithm 1, the algorithm
may generate λ 1. As Algorithm 1 is linearly convergent, there exist γ ∈ [0, 1)
with
D(λ ⋆ ) − D(λ 2) ≤ γ
(
D(λ ⋆ ) − D(λ 1)
)
.
So, we have
(1 − γ)
(
D(λ ⋆ ) − D(λ 1)
)
≤ D(λ 2) − D(λ 1) ≤
⟨
Aˆx1 + Bz 1 − b, λ 2 − λ 1⟩
,
where the last inequality follows from the concavity of the function D. Since
λ 2 − λ 1 = t(Ax2 + Bz 2 − b), Lemma 4 implies that
D(λ ⋆ ) − D(λ 1) ≤ t
1− γ ∥ξ∥2,
so D satisﬁes the P/suppress L inequality.

<!-- page 26 -->
26 Moslem Zamani et al.
Another assumption used for establishing linear convergence is L-smoothness;
see for example [7,8,12,34]. Deng et al. [8] show that the sequence {(xk, z k, λ k)}
is convergent linearly to a saddle point under Scenario 1 and 2 given in T able
1.
Table 1: Scenarios leading to linear convergence rates
Scenario Strong convexity Lipschitz continuity Full row ra nk
1 f, g ∇ f A
2 f, g ∇ f, ∇ g -
3 f ∇ f, ∇ g B T
It is worth mentioning that Scenario 1 or Scenario 2 implies strong con vex-
ity of the dual objective function and therefore the P/suppress L inequalityis resulted,
see [1]. Hence, Theorem 6 implies the linear convergence in terms of du al value
under Scenario 1 or Scenario 2. Deng et al. [8] studied the linear con vergence
under Scenario 3, but they just proved the linear convergence of the sequence
{(xk, Bz k, λ k)}. In the next section, we investigate the R-linear convergence
without assuming L-smoothness of f . Indeed, we establish the R-linear con-
vergence when f is strongly convex, g is L-smooth and B has full row rank.
Note that the P/suppress L inequality does not imply necessarily Scenario 1 or Sce-
nario 2. Indeed, consider the following optimization problem,
min f (x) + g(z),
s. t. x + z = 0,
x, z ∈ Rn,
where f (x) = 1
2 ∥x∥2 + ∥x∥1 and g(z) = 1
2 ∥z∥2 + ∥z∥1. With some algebra, one
may show that D(λ) = ∑n
i=1 h(λ i) with
h(s) =





− (s − 1)2, s > 1
0, |s| ≤ 1
− (s + 1)2, s < − 1.
Hence, the P/suppress L inequality holds forLp = 1
2 while neither f nor g is L-smooth.
As mentioned earlier the performance estimation problem including th e
P/suppress L inequality at ﬁnite set of points is a relaxation for computing the worst-
case convergence rate. Contrary to Theorem 6, we could not man age to prove
the linear convergence of primal and dual residuals under the assu mptions of
Theorem 6 by employing performance estimation.
5 R-linear convergence of ADMM
In this section, we study the R-linear convergence of ADMM. Recall that
ADMM enjoys R-linear convergent in terms of dual objective value if
D(λ ⋆ ) − D(λ N ) ≤ ργ N ,

<!-- page 27 -->
Convergence rate of ADMM 27
for some ρ ≥ 0 and γ ∈ [0, 1).
We investigate the R-linear convergence under the following scenar ios:
– (S1): f ∈ F A
c1(Rn) is L-smooth with c1 > 0 and A has full row rank;
– (S2): f ∈ F A
c1(Rn) with c1 > 0, g is L-smooth and B has full row rank.
Our technique for proving the R-linear convergence is based on est ablishing
the linear convergence of the sequence {V k} given by
V k = ∥λ k − λ ⋆ ∥2 + t2 
zk − z⋆ 
2
B . (41)
Note that V k is called Lyapunov function for ADMM; see [6].
First we consider the case that f is L-smooth and c1-strongly convex rela-
tive to A. The following proposition establishes the linear convergence of {V k}.
Proposition 2 Let f ∈ F A
c1(Rn) be L-smooth with c1 > 0, g ∈ F 0(Rm) and
let A has full row rank. If t <
√
c1L
λ min(AAT ) , then
V k+1 ≤
(
1 − 2c1t
c1d+2c1t+t2
)
V k, (42)
where d = L
λ min(AAT ) .
Proof. We may assume without loss of generality that x⋆ , z ⋆ and b are zero;
see our discussion in Section 2. By optimality conditions, we have
∇ f (xk+1) = − AT (
λ k + tAxk+1 + tBz k)
, η k = − BT λ k+1,
∇ f (x⋆ ) = − AT λ ⋆ , η ⋆ = − BT λ ⋆ ,
for some ηk ∈ ∂g (zk+1) and η⋆ ∈ ∂g (z⋆). Let α = 2t
c2
1d2+2c1dt2− 4c2
1t2+t4 . By
Theorem 2, we get
α
(
t2 + c1d
) 2
(
f (xk+1) − f ⋆ +
⣨
λ⋆ , Axk+1
⟩
− 1
2L


AT
(
λk + tAxk+1 + tBz k − λ⋆
) 


2)
+
2αt2(c1d + t2) (
f ⋆ − f (xk+1) − c1
2


xk+1



2
A
−
⣨
λk + tAxk+1 + tBz k, Axk+1
⟩)
+
2t
(
g(zk+1) − g⋆ +
⣨
λ⋆ , Bzk+1
⟩)
+ 2t
(
g⋆ − g(zk+1) −
⣨
λk+1, Bzk+1
⟩)
+
α (c2
1d2 − t4)
(
f ⋆ − f (xk+1) −
⣨
λk + tAxk+1 + tBz k, Axk+1
⟩
−
1
2L


AT
(
λk + tAxk+1 + tBz k − λ⋆
) 


2
)
≥ 0.
As ∥AT λ∥2 ≥ L
d ∥λ∥2 and λ k+1 = λ k + tAxk+1 + tBz k+1, we obtain the
following inequality after performing some algebraic manipulations
(
1 − 2ct
cd+2ct+t2
) (

λk − λ⋆



2
+ t2


Bz k



2)
−
(

λk+1 − λ⋆



2
+ t2


Bz k+1



2)
−
2αc2
1t


λk − λ⋆ + t2+2c1t+c1d
2c1
Axk+1 + t2+c1d
2c1
Bz k



2
≥ 0.

<!-- page 28 -->
28 Moslem Zamani et al.
The above inequality implies that
V k+1 ≤
(
1 − 2c1t
c1d+2c1t+t2
)
V k,
and the proof is complete.
Note that one can improve bound (42) under the assumptions of Pr opo-
sition 2 and the µ -strong convexity of f by employing the following known
inequality
1
2(1− µ
L )
(
1
L ∥∇ f (x) − ∇ f (y)∥2 + µ ∥x − y∥2 − 2µ
L ⟨∇ f (x) − ∇ f (y), x − y⟩
)
≤ f (y) − f (x) − ⟨∇ f (x), y − x⟩ .
Indeed, we employed the given inequality but we could not manage to o btain
a closed form formula for the convergence rate. The next theore m establishes
the R-linear convergence of ADMM in terms of dual objective value u nder the
assumptions of Proposition 2.
Theorem 8 Let N ≥ 4 and let A has full row rank. Suppose that f ∈ F c1(Rn)
is L-smooth with c1 > 0 and g ∈ F 0(Rm). If t < min{c1,
√
c1L
λ min(AAT ) }, then
D(λ ⋆ ) − D(λ N ) ≤ ρ
(
1 − 2c1t
c1d+2c1t+t2
) N
,
where d = L
λ min(AAT ) and ρ = V 0
16t
(
1 − 2c1t
c1d+2c1t+t2
) − 4
.
Proof. By Theorem 3 and Proposition 2, one can infer the following inequali-
ties,
D(λ ⋆ ) − D(λ N ) ≤ V N −4
16t
≤ V 0
16t
(
1 − 2c1t
c1d+2c1t+t2
) N − 4
,
which shows the desired inequality.
Nishihara et al. [34] showed the R-linear convergence of ADMM in term s
of {xk, z k, λ k} under the following conditions:
i) The function f is L-smooth and µ -strong with µ > 0;
ii) The matrix A is invertible and that B has full column rank.
In Theorem 8, we obtain the R-linear convergence under weaker as sumptions.
Indeed, we replace condition ii) with the matrix A having full row rank.
In the sequel, we investigate the R-linear convergence under the h ypotheses
of scenario (S2). The next proposition shows the linear convergen ce of {V k}.
Proposition 3 Let f ∈ F A
c1(Rn) with c1 > 0 and let g ∈ F 0(Rm) be L-smooth.
Suppose that B has full row rank and k ≥ 1. If t ≤ min{ c1
2 , L
2λ min(BB T ) }, then
V k+1 ≤
(
L
L+tλ min(BB T )
) 2
V k. (43)

<!-- page 29 -->
Convergence rate of ADMM 29
Proof. Analogous to the proof of Proposition 2, we assume that x⋆ = 0, z⋆ = 0
and b = 0. Due to the optimality conditions, we have
ξk+1 = − AT (
λ k + tAxk+1 + tBz k)
, ξ ⋆ = − AT λ ⋆ ,
∇ g(zk) = − BT λ k, ∇ g(zk+1) = − BT λ k+1, ∇ g(z⋆) = − BT λ ⋆ ,
for some ξk+1 ∈ ∂f (xk+1) and ξ⋆ ∈ ∂f (x⋆ ). Suppose that d = L
λ min(BB T ) and
α = 2dt
d+t . By Theorem 2, we obtain
α
(
d2 + t2)
d2 − t2
(
f ⋆ − f (xk+1) −
⣨
λk + tAxk+1 + tBz k, Axk+1
⟩
− c1
2


xk+1



2
A
)
+
α
(
d2 + t2)
d2 − t2
(
f (xk+1) − f (x⋆ ) +
⣨
λ⋆ , Axk+1
⟩
− c1
2


xk+1



2
A
)
+
α
(
g(zk+1) − g⋆ +
⣨
λ⋆ , Bzk+1
⟩
− 1
2L


BT
(
λ⋆ − λk+1
) 


2)
+
α
(
g⋆ − g(zk+1) −
⣨
λk+1, Bzk+1
⟩
− 1
2L


BT
(
λ⋆ − λk+1
) 


2)
+
α
(
g(zk) − g(zk+1) +
⣨
λk+1, Bz k − Bz k+1
⟩
− 1
2L


BT
(
λk+1 − λk
) 


2)
+
α
(
g(zk+1) − g(zk) +
⣨
λk, Bz k+1 − Bz k
⟩
− 1
2L


BT
(
λk+1 − λk
) 


2)
≥ 0.
By employing ∥BT λ∥2 ≥ L
d ∥λ∥2 and λ k+1 = λ k + tAxk+1 + tBz k+1, the
aforementioned inequality can be expressed as follows after some a lgebraic
manipulation,
−α 2
4




( 2t2
d2 − dt
)
Axk+1 + Bz k − (1 + t
d
) Bz k+1




2
− 2t
(
d2 + t2) (
cd2 − dt(c + t) − t3)
(d2 − t2)2


Axk+1



2
− α 2
4d2



λk − λ⋆ +
(2d2 − (d − t)2
d − t
)
Axk+1 + (d + t) Bz k+1




2
+
(
d
d+t
) 2 (

λk − λ⋆



2
+ t2


Bz k



2)
−
(

λk+1 − λ⋆



2
+ t2


Bz k+1



2)
≥ 0.
Hence, we have
V k+1 ≤
(
d
d+t
) 2
V k,
and the proof is complete.
As the sequence {V k} is not increasing [6, Convergence Proof], we have
V 1 ≤ V 0. Thus, by using Theorem 3 and Proposition 3, one can infer the
following theorem.
Theorem 9 Let f ∈ F A
c1(Rn) with c1 > 0 and let g ∈ F 0(Rm) be L-smooth.
Assume that N ≥ 5 and B has full row rank. If t < min{ c1
2 , L
2λ min(BB T ) }, then
D(λ ⋆ ) − D(λ N ) ≤ ρ
(
L
L+tλ min(BB T )
) 2N
, (44)
where ρ = V 0
16t
(
L
L+tλ min(BB T )
) − 10
.

<!-- page 30 -->
30 Moslem Zamani et al.
In the same line, one can infer the R-linear convergence in terms of p rimal
and dual residuals under the assumptions of Theorem 8 and Theore m 9. In this
section, we proved the linear convergence of {V k} under two scenarios (S1) and
(S2). By (7), it is readily seen that function − D is strongly convex under the
hypotheses of both scenarios (S1) and (S2). Therefore, both s cenarios imply
the P/suppress L inequality. One may wonder that if the P/suppress L inequality and the strong
convexity of f imply the linear of {V k}. By using performance estimation, we
could not establish such an implication.
As mentioned above, function − D under both scenarios are µ -strongly
convex. Hence, the optimal solution set of the dual problem is uniqu e and one
can infer the R-linear convergence of λ N by using Theorem 8 (Theorem 9) and
the known inequality,
µ
2

λ N − λ ⋆ 
2
≤ D(λ ⋆ ) − D(λ N ).
Concluding remarks
In this paper we developed performance estimation framework to h andle dual-
based methods. Thanks to this framework, we could obtain some tig ht con-
vergence rates for ADMM. This framework may be exploited for the analysis
of other variants of ADMM in the ergodic and non-ergodic sense. Mo reover,
similarly to [23], one can apply this framework for introducing and analy zing
new accelerated ADMM variants. Moreover, most results hold for a ny arbi-
trary positive step length, t, but we managed to get closed form formulas for
some interval of positive numbers.
References
1. Abbaszadehpeivasti, H., de Klerk, E., Zamani, M.: Condit ions for linear convergence of
the gradient method for non-convex optimization. Optimiza tion Letters 17(5), 1105–
1125 (2023)
2. Abbaszadehpeivasti, H., de Klerk, E., Zamani, M.: On the r ate of convergence of the
diﬀerence-of-convex algorithm (DCA). Journal of Optimiza tion Theory and Applica-
tions pp. 1–22 (2023)
3. Beck, A.: First-order methods in optimization. SIAM (201 7)
4. Bertsekas, D.: Convex optimization algorithms. Athena S cientiﬁc (2015)
5. Bolte, J., Nguyen, T.P., Peypouquet, J., Suter, B.W.: Fro m error bounds to the com-
plexity of ﬁrst-order descent methods for convex functions . Mathematical Programming
165(2), 471–507 (2017)
6. Boyd, S., Parikh, N., Chu, E., Peleato, B., Eckstein, J., e t al.: Distributed optimization
and statistical learning via the alternating direction met hod of multipliers. Foundations
and Trends® in Machine learning 3(1), 1–122 (2011)
7. Davis, D., Yin, W.: Faster convergence rates of relaxed Pe aceman-Rachford and ADMM
under regularity assumptions. Mathematics of Operations R esearch 42(3), 783–805
(2017)
8. Deng, W., Yin, W.: On the global and linear convergence of t he generalized alternating
direction method of multipliers. Journal of Scientiﬁc Comp uting 66(3), 889–916 (2016)
9. Drori, Y., Teboulle, M.: Performance of ﬁrst-order metho ds for smooth convex mini-
mization: a novel approach. Mathematical Programming 145(1), 451–482 (2014)

<!-- page 31 -->
Convergence rate of ADMM 31
10. Franca, G., Robinson, D., Vidal, R.: ADMM and accelerate d ADMM as continuous
dynamical systems. In: International Conference on Machin e Learning, pp. 1559–1567.
PMLR (2018)
11. Gabay, D., Mercier, B.: A dual algorithm for the solution of nonlinear variational prob-
lems via ﬁnite element approximation. Computers & mathemat ics with applications
2(1), 17–40 (1976)
12. Giselsson, P., Boyd, S.: Linear convergence and metric s election for Douglas-Rachford
splitting and ADMM. IEEE Transactions on Automatic Control 62(2), 532–544 (2016)
13. Glowinski, R., Marroco, A.: Sur l’approximation, par ´ e l´ ements ﬁnis d’ordre un, et la
r´ esolution, par p´ enalisation-dualit´ e d’une classe de probl` emes de dirichlet non lin´ eaires.
ESAIM: Mathematical Modelling and Numerical Analysis-Mod ´ elisation Math´ ematique
et Analyse Num´ erique9(R2), 41–76 (1975)
14. Glowinski, R., Osher, S.J., Yin, W.: Splitting methods i n communication, imaging,
science, and engineering. Springer (2017)
15. Goldfarb, D., Ma, S., Scheinberg, K.: Fast alternating l inearization methods for mini-
mizing the sum of two convex functions. Mathematical Progra mming 141(1), 349–382
(2013)
16. Goldstein, T., O’Donoghue, B., Setzer, S., Baraniuk, R. : Fast alternating direction op-
timization methods. SIAM Journal on Imaging Sciences 7(3), 1588–1623 (2014)
17. Gu, G., Yang, J.: On the dual step length of the alternatin g direction method of mul-
tipliers. arXiv preprint arXiv:2006.08309 (2020)
18. Han, D., Sun, D., Zhang, L.: Linear rate convergence of th e alternating direction method
of multipliers for convex composite programming. Mathemat ics of Operations Research
43(2), 622–637 (2018)
19. Han, D.R.: A survey on some recent developments of altern ating direction method of
multipliers. Journal of the Operations Research Society of China 10(1), 1–52 (2022)
20. Hastie, T., Tibshirani, R., W ainwright, M.: Statistica l learning with sparsity. Mono-
graphs on statistics and applied probability 143, 143 (2015)
21. He, B., Yuan, X.: On the O(1/n) convergence rate of the Dou glas-Rachford alternating
direction method. SIAM Journal on Numerical Analysis 50(2), 700–709 (2012)
22. Hong, M., Luo, Z.Q.: On the linear convergence of the alte rnating direction method of
multipliers. Mathematical Programming 162(1), 165–199 (2017)
23. Kim, D., Fessler, J.A.: Optimized ﬁrst-order methods fo r smooth convex minimization.
Mathematical programming 159(1), 81–107 (2016)
24. Li, H., Lin, Z.: Accelerated alternating direction meth od of multipliers: An optimal
O(1/k) nonergodic analysis. Journal of Scientiﬁc Computin g 79(2), 671–699 (2019)
25. Lin, Z.: Alternating Direction Method of Multipliers fo r Machine Learning. Springer
Nature
26. Liu, H., Shi, Y., W ang, Z., Ran, L., L¨ u, Q., Li, H.: A distr ibuted algorithm based on
relaxed ADMM for energy resources coordination. Internati onal Journal of Electrical
Power & Energy Systems 135, 107482 (2022)
27. Liu, Y., Yuan, X., Zeng, S., Zhang, J.: Partial error boun d conditions and the linear
convergence rate of the alternating direction method of mul tipliers. SIAM Journal on
Numerical Analysis 56(4), 2095–2123 (2018)
28. Lozenski, L., Villa, U.: Consensus ADMM for inverse prob lems governed by multiple
PDE models. arXiv preprint arXiv:2104.13899 (2021)
29. Lu, H., Freund, R.M., Nesterov, Y.: Relatively smooth co nvex optimization by ﬁrst-order
methods, and applications. SIAM Journal on Optimization 28(1), 333–354 (2018)
30. Madani, R., Kalbat, A., Lavaei, J.: ADMM for sparse semid eﬁnite programming with
applications to optimal power ﬂow problem. In: 2015 54th IEE E Conference on Decision
and Control (CDC), pp. 5932–5939. IEEE (2015)
31. Monteiro, R.D., Svaiter, B.F.: Iteration-complexity o f block-decomposition algorithms
and the alternating direction method of multipliers. SIAM J ournal on Optimization
23(1), 475–507 (2013)
32. Necoara, I., Nesterov, Y., Glineur, F.: Linear converge nce of ﬁrst order methods for
non-strongly convex optimization. Mathematical Programm ing 175(1), 69–107 (2019)
33. Nesterov, Y.: Introductory lectures on convex optimiza tion: A basic course, vol. 87.
Springer Science & Business Media (2003)

<!-- page 32 -->
32 Moslem Zamani et al.
34. Nishihara, R., Lessard, L., Recht, B., Packard, A., Jord an, M.: A general analysis of the
convergence of ADMM. In: International Conference on Machi ne Learning, pp. 343–352.
PMLR (2015)
35. Pe˜ na, J., Vera, J.C., Zuluaga, L.F.: Linear convergenc e of the Douglas-Rachford algo-
rithm via a generic error bound condition. arXiv preprint ar Xiv:2111.06071 (2021)
36. Rockafellar, R.T.: Convex analysis. In: Convex analysi s. Princeton university press
(1970)
37. Rudin, L.I., Osher, S., Fatemi, E.: Nonlinear total vari ation based noise removal algo-
rithms. Physica D: nonlinear phenomena 60(1-4), 259–268 (1992)
38. Ryu, E.K., Taylor, A.B., Bergeling, C., Giselsson, P.: O perator splitting performance
estimation: Tight contraction factors and optimal paramet er selection. SIAM Journal
on Optimization 30(3), 2251–2271 (2020)
39. Sabach, S., Teboulle, M.: Faster Lagrangian-based meth ods in convex optimization.
SIAM Journal on Optimization 32(1), 204–227 (2022)
40. Stellato, B., Banjac, G., Goulart, P., Bemporad, A., Boy d, S.: Osqp: An operator split-
ting solver for quadratic programs. Mathematical Programm ing Computation 12(4),
637–672 (2020)
41. Taylor, A.B., Hendrickx, J.M., Glineur, F.: Smooth stro ngly convex interpolation and
exact worst-case performance of ﬁrst-order methods. Mathe matical Programming
161(1-2), 307–345 (2017)
42. Yuan, X., Zeng, S., Zhang, J.: Discerning the linear conv ergence of admm for structured
convex optimization through the lens of variational analys is. J. Mach. Learn. Res. 21,
83–1 (2020)
A Appendix
Lemma 4 Let N ≥ 4 and t, c1 ∈ R. Let D(t, c1) be N × N symmetric matrix given in
Theorem 4. If c1 > 0 is given, then
[0, c1] ⊆ { t : D(t, c1) ⪰ 0}.
Proof. The argument proceeds in the same manner as in Lemma 1. Due to t he convexity
of {t : D(t, c1) ⪰ 0}, is suﬃcient to establish the positive semideﬁniteness of D(0, c1) and
D(c1, c1). As D(0, c1) is diagonally dominant, it is positive semideﬁnite. Next, we proceed
to demonstrate the positive deﬁniteness of the matrix K = D(1, 1) by computing its leading
principal minors. One can show that the claim holds for N = 4. So we investigate N ≥ 5.
To accomplish this, we perform the following elementary row operations on matrix D:
i) Add the second row to the third row;
ii) Add the second row to the last row;
iii) Add the third row to the forth row;
iv) For i = 4 : N − 2
– Add i − th row to ( i + 1) − th row;
– Add 3−i
2i2−3i−1 times of i − th row to the last row;
v) Add 2N 2 −8N+9
2N 2 −7N+4 times of ( N − 1) − th row to N − th row.
By executing these operations, we transform K into an upper triangular matirx J with
diagonal
Jk,k =







2, k = 1
2k2 − 3k − 1, 2 ≤ k ≤ N − 1
2N 2 − 7N + 8 − N 2
(N−1)2 − (2N 2−8N+9)2
2N 2−7N+4 − ∑N−2
i=4
(i−3)2
2i2−3i−1 , k = N .
It is seen all ﬁrst ( N − 1) diagonal elements of J are positive. W e show that JN,N is also
positive. By using inequality (17), we get
2N 2 − 7N + 8 − N 2
(N − 1)2 − (2N 2 − 8N + 9)2
2N 2 − 7N + 4 −
N−2∑
i=4
(i − 3)2
2i2 − 3i − 1 ≥
2N 2 − 7N + 8 − 25
16 − (2N 2 − 8N + 9) − N−5
2 − 1 + 2
N−3 ≥ N
2 − 17
16 > 0,

<!-- page 33 -->
Convergence rate of ADMM 33
for N ≥ 5, which implies JN,N > 0. Hence, D(c1, c1) ⪰ 0 and the proof is complete.