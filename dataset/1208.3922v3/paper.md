# On the Linear Convergence of the Alternating Direction Method of Multipliers

**arXiv ID:** 1208.3922v3

**Authors:** Mingyi Hong, Zhi-Quan Luo

**Abstract:** We analyze the convergence rate of the alternating direction method of multipliers (ADMM) for minimizing the sum of two or more nonsmooth convex separable functions subject to linear constraints. Previous analysis of the ADMM typically assumes that the objective function is the sum of only two convex functions defined on two separable blocks of variables even though the algorithm works well in numerical experiments for three or more blocks. Moreover, there has been no rate of convergence analysis for the ADMM without strong convexity in the objective function. In this paper we establish the global linear convergence of the ADMM for minimizing the sum of any number of convex separable functions. This result settles a key question regarding the convergence of the ADMM when the number of blocks is more than two or if the strong convexity is absent. It also implies the linear convergence of the ADMM for several contemporary applications including LASSO, Group LASSO and Sparse Group LASSO without any strong convexity assumption. Our proof is based on estimating the distance from a dual feasible solution to the optimal dual solution set by the norm of a certain proximal residual, and by requiring the dual stepsize to be sufficiently small.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:1208.3922v3  [math.OC]  26 Mar 2013
On the Linear Convergence of the Alternating
Direction Method of Multipliers ∗
Mingyi Hong and Zhi-Quan Luo †
Dedicated to the fond memories of a close friend and collabor ator, Paul Y. Tseng
August 13, 2012; Revised March 20, 2013
Abstract
We analyze the convergence rate of the alternating direction meth od of multipliers (ADMM)
for minimizing the sum of two or more nonsmooth convex separable fu nctions subject to linear
constraints. Previous analysis of the ADMM typically assumes that t he objective function
is the sum of only two convex functions deﬁned on two separable blocks of variables even
though the algorithm works well in numerical experiments for three or more blocks. Moreover,
there has been no rate of convergence analysis for the ADMM witho ut strong convexity in the
objective function. In this paper we establish the global linear conv ergence of the ADMM for
minimizing the sum of any number of convex separable functions. This result settles a key
question regarding the convergence of the ADMM when the number of blocks is more than two
or if the strong convexity is absent. It also implies the linear converg ence of the ADMM for
several contemporary applications including LASSO, Group LASSO a nd Sparse Group LASSO
without any strong convexity assumption. Our proof is based on es timating the distance from a
dual feasible solution to the optimal dual solution set by the norm of a certain proximal residual,
and by requiring the dual stepsize to be suﬃciently small.
KEY WORDS: Linear convergence, alternating directions of multiplier s, error bound, dual ascent.
AMS(MOS) Subject Classiﬁcations: 49, 90.
∗ The research is supported by the National Science Foundatio n, grant number DMS-1015346.
† Department of Electrical and Computer Engineering, Univer sity of Minnesota, Minneapolis, MN 55455, USA.
Email: {mhong, luozq }@umn.edu
1

<!-- page 2 -->
1 Introduction
Consider the problem of minimizing a separable nonsmooth co nvex function subject to linear equal-
ity constraints:
minimize f (x) = f1(x1) + f2(x2) + · · ·+ fK(xK)
subject to Ex = E1x1 + E2x2 + · · ·+ EKxK = q
xk ∈ Xk, k = 1, 2, ..., K,
(1.1)
where each fk is a nonsmooth convex function (possibly with extended valu es), x = (xT
1 , ..., x T
K )T ∈
ℜ n is a partition of the optimization variable x, X = ∏K
k=1 Xk is the feasible set for x, and
E = (E1, E 2, ..., E K ) ∈ ℜ m× n is an appropriate partition of matrix E (consistent with the partition
of x) and q ∈ ℜ m is a vector. Notice that the model (1.1) can easily accommoda te general linear
inequality constraints Ex ≥ q by adding one extra block. In particular, we can introduce a s lack
variable xK+1 ≥ 0 and rewrite the inequality constraint as Ex − xk+1 = q. The constraint xK+1 ≥ 0
can be enforced by adding a new convex component function fK+1(xK+1) = iℜ m
+ (xK+1) to the
objective function f (x), where iℜ m
+ (xK+1) is the indicator function for the nonnegative orthant ℜ m
+
iℜ m
+ (xK+1) =
{
0, if xK+1 ≥ 0 (entry wise),
∞ , otherwise.
In this way, the inequality constrained problem with K blocks is reformulated as an equivalent
equality constrained convex minimization problem with K + 1 blocks.
Optimization problems of the form (1.1) arise in many emergi ng applications involving struc-
tured convex optimization. For instance, in compressive se nsing applications, we are given an
observation matrix A and a noisy observation vector b ≈ Ax. The goal is to estimate the sparse
vector x by solving the following ℓ1 regularized linear least squares problem:
minimize ∥y∥2 + λ∥x∥1
subject to Ax + y = b,
where λ > 0 is a penalty parameter. Clearly, this is a structured conve x optimization problem
of the form (1.1) with K = 2. If the variable x is further constrained to be nonnegative, then
the corresponding compressive sensing problem can be formu lated as a three block ( K = 3) convex
separable optimization problem (1.1) by introducing a slac k variable. Similarly, in the stable version
of robust principal component analysis (PCA) [59], we are gi ven an observation matrix M ∈ ℜ m× n
which is a noise-corrupted sum of a low rank matrix L and a sparse matrix S. The goal is recover
L and S by solving the following nonsmooth convex optimization pro blem
minimize ∥L∥∗ + ρ∥S∥1 + λ∥Z∥2
F
subject to L + S + Z = M
1

<!-- page 3 -->
where ∥·∥∗ denotes the matrix nuclear norm (deﬁned as the sum of the matr ix singular eigenvalues),
while ∥ · ∥1 and ∥ · ∥F denote, respectively, the ℓ1 and the Frobenius norm of a matrix (equal to the
standard ℓ1 and ℓ2 vector norms when the matrix is viewed as a vector). In the abo ve formulation,
Z denotes the noise matrix, and ρ, λ are some ﬁxed penalty parameters. It is easily seen that the
stable robust PCA problem corresponds to the three block cas e K = 3 in the problem (1.1) with
x = (L, S, Z ) and
f1(L) = ∥L∥∗, f 2(S) = ∥S∥1, f 3(Z) = ∥Z∥2
F , (1.2)
while the coupling linear constraint is given L + S + Z = M . In image processing applications
where the low rank matrix L is additionally constrained to be nonnegative, then the abo ve problem
can be reformulated as
minimize ∥L∥∗ + ρ∥S∥1 + λ∥Z∥2
F + iℜ mn
+ (C)
subject to L + S + Z = M, L − C = 0,
where C is a slack matrix variable of the same size as L, and iℜ mn
+ (·) is the indicator function for
the nonnegative orthant ℜ mn
+ . In this case, the stable robust PCA problem is again in the fo rm
of (1.1). In particular, it has 4 block variables ( L, S, Z, C ) and the ﬁrst three convex functions are
the same as in (1.2), while the fourth convex function is give n by f4(C) = iℜ mn
+ (C). The coupling
linear constraints are L + S + Z = M, L − C = 0. Other applications of the form (1.1) include the
latent variable Gaussian graphical model selection proble m, see [9].
A popular approach to solving the separable convex optimiza tion problem (1.1) is to attach a
Lagrange multiplier vector y to the linear constraints Ex = q and add a quadratic penalty, thus
obtaining an augmented Lagrangian function of the form
L(x; y) = f (x) + ⟨y, q − Ex⟩ + ρ
2 ∥q − Ex∥2, (1.3)
where ρ ≥ 0 is a constant. The augmented dual function is given by
d(y) = min
x
f (x) + ⟨y, q − Ex⟩ + ρ
2 ∥q − Ex∥2 (1.4)
and the dual problem (equivalent to (1.1) under mild conditi ons) is
max
y
d(y). (1.5)
Moreover, if ρ > 0, then Ex is constant over the set of minimizers of (1.4) (see Lemma 2.1 in
Section 2). This implies that the dual function d(y) is diﬀerentiable with
∇ d(y) = q − Ex(y)
where x(y) is a minimizer of (1.4). Given the diﬀerentiability of d(y), it is natural to consider the
following dual ascent method to solve the primal problem (1. 1)
y := y + α ∇ d(y) = y + α (q − Ex(y)), (1.6)
2

<!-- page 4 -->
where α > 0 is a suitably chosen stepsize. Such a dual ascent strategy i s well suited for structured
convex optimization problems that are amenable to decompos ition. For example, if the objective
function f is separable (i.e., of the form given in (1.1)) and if we selec t ρ = 0, then the minimization
in (1.4) decomposes into K independent minimizations whose solutions frequently can be obtained
in a simple form. In addition, the iterations can be implemen ted in a manner that exploits the
sparsity structure of the problem and, in certain network ca ses, achieve a high degree of parallelism.
Popular choices for the ascent methods include (single) coo rdinate ascent (see [4, 8, 10, 34, 40, 42,
51, 52, 58]), gradient ascent (see [34, 42, 53]) and gradient projection [23, 32]. (See [5, 34, 49] for
additional references.)
For large scale optimization problems, it is numerically ad vantageous to select ρ > 0. Unfor-
tunately, this also introduces variable coupling in the aug mented Lagrangian (1.3), which makes
the exact minimization step in (1.4) no longer decomposable across variable blocks even if f has
a separable structure. In this case, it is more economical to minimize (1.4) inexactly by updating
the components of x cyclically via the coordinate descent method. In particula r, we can apply
the Gauss-Seidel strategy to inexactly minimize (1.4), and then update the multiplier y using an
approximate optimal solution of (1.4) in a manner similar to (1.6). The resulting algorithm is
called the Alternating Direction Method of Multipliers (AD MM) and is summarized as follows
(see [17–20]). In the general context of sums of monotone ope rators, the work of [16] describes a
large family of splitting methods for K ≥ 3 blocks which, when applied to the dual, result in similar
but not identical methods to the ADMM algorithm (1.7) below.
Alternating Direction Method of Multipliers (ADMM)
At each iteration r ≥ 1, we ﬁrst update the primal variable blocks in the Gauss-
Seidel fashion and then update the dual multiplier using the updated primal
variables:









xr+1
k = arg min
xk∈ Xk
L(xr+1
1 , ..., x r+1
k− 1, x k, x r
k+1, ..., x r
K ; yr), k = 1, 2, ..., K,
yr+1 = yr + α (q − Exr+1) = yr + α
(
q −
K∑
k=1
Ekxr+1
k
)
,
(1.7)
where α > 0 is the step size for the dual update.
Notice that if there is only one block ( K = 1), then the ADMM reduces to the standard
augmented Lagrangian method of multipliers for which the gl obal convergence is well understood
(see e.g., [1]). In particular, it is known that, under mild a ssumptions on the problem, this type of
dual gradient ascent methods generate a sequence of iterate s whose limit points must be optimal
solutions of the original problem (see [8, 49, 51]). For the s pecial case of ordinary network ﬂow
3

<!-- page 5 -->
problems, it is further known that an associated sequence of dual iterates converges to an optimal
solution of the dual (see [4]). The rate of convergence of dua l ascent methods has been studied
in the reference [37] which showed that, under mild assumpti ons on the problem, the distance
to the optimal dual solution set from any y ∈ ℜ m near the set is bounded above by the dual
optimality ‘residual’ ∥∇ d(y)∥. By using this bound, it can be shown that a number of ascent
methods, including coordinate ascent methods and a gradien t projection method, converge at least
linearly when applied to solve the dual problem (see [35, 36] ; also see [2, 11, 29] for related analysis).
(Throughout this paper, by ‘linear convergence’ we mean roo t–linear convergence (denoted by
R-linear convergence) in the sense of Ortega and Rheinboldt [41].)
When there are two blocks ( K = 2), the convergence of the ADMM was studied in the context of
Douglas-Rachford splitting method [13–15] for ﬁnding a zer o of the sum of two maximal monotone
operators. It is known that in this case every limit point of t he iterates is an optimal solution of the
problem. The recent work of [21, 22, 26] have shown that, unde r some additional assumptions, the
objective values generated by the ADMM algorithm and its acc elerated version (which performs
some additional line search steps for the dual update) conve rge at a rate of O(1/r ) and O(1/r 2)
respectively. Moreover, if the objective function f (x) is strongly convex and the constraint matrix
E is row independent, then the ADMM is known to converge linear ly to the unique minimizer of
(1.1) [33]. [One notable exception to the strong convexity r equirement is in the special case of
linear programming for which the ADMM is linearly convergen t [14].] More recent convergence
rate analysis of the ADMM still requires at least one of the co mponent functions ( f1 or f2) to
be strongly convex and have a Lipschitz continuous gradient . Under these and additional rank
conditions on the constraint matrix E, some linear convergence rate results can be obtained for a
subset of primal and dual variables in the ADMM algorithm (or its variant); see [6,12,24]. However,
when there are more than two blocks involved ( K ≥ 3), the convergence (or the rate of convergence)
of the ADMM method is unknown, and this has been a key open ques tion for several decades. The
recent work [38] describes a list of novel applications of th e ADMM with K ≥ 3 and motivates
strongly for the need to analyze the convergence of the ADMM i n the multi-block case. The recent
monograph [7] contains more details of the history, converg ence analysis and applications of the
ADMM and related methods.
A main contribution of this paper is to establish the global ( linear) convergence of the ADMM
method for a class of convex objective functions involving a ny number of blocks ( K is arbitrary).
The key requirement for the global (linear) convergence is t he satisfaction of a certain error bound
condition that is similar to that used in the analysis of [37] . This error bound estimates the
distance from an iterate to the optimal solution set in terms of a certain proximity residual. The
class of objective functions that are known to satisfy this e rror bound condition include many of
the compressive sensing applications, such as LASSO [48], G roup LASSO [56] or Sparse Group
LASSO [57].
4

<!-- page 6 -->
2 Technical Preliminaries
Let f be a closed proper convex function in ℜ n, let E be an m × n matrix, let q be a vector in ℜ m.
Let dom f denote the eﬀective domain of f and let int(dom f ) denote the interior of dom f . We
make the following standing assumptions regarding f :
Assumption A.
(a) The global minimum of (1.1) is attained and so is its dual o ptimal value. The intersection
X ∩ int(dom f ) ∩ { x |Ex = q} is nonempty.
(b) f = f1(x1) + f2(x2) + · · ·+ fK(xK ), with each fk further decomposable as
fk(xk) = gk(Akxk) + hk(xk)
where gk and hk are both convex and continuous over their domains, and Ak’s are some given
matrices (not necessarily full column rank, and can be zero) .
(c) Each gk is strictly convex and continuously diﬀerentiable on int(do m gk) with a uniform
Lipschitz continuous gradient
∥∇ AT
k gk(Axk) − AT
k ∇ gk(Ax′
k)∥ ≤ L∥xk − x′
k∥, ∀ xk, x ′
k ∈ Xk
where L > 0 is a constant.
(d) Each hk satisﬁes either one of the following conditions
1. The epigraph of hk(xk) is a polyhedral set.
2. hk(xk) = λ k∥xk∥1 + ∑
J wJ ∥xk,J ∥2, where xk = (· · ·, x k,J , · · ·) is a partition of xk with
J being the partition index.
3. Each hk(xk) is the sum of the functions described in the previous two ite ms.
(e) For any ﬁxed and ﬁnite y and ξ, ∑
k hk(xk) is ﬁnite for all x ∈ { x : L(x; y) ≤ ξ} ∩ X.
(f) Each submatrix Ek has full column rank.
(g) The feasible sets Xk, k = 1, · · ·, K are compact polyhedral sets.
We have the following remarks regarding to the assumptions m ade.
1. Each fk may only contain convex function hk. That is, the strongly convex part gk can be
absent. Also, since the matrices Ak’s are not required to have full column rank, the overall
objective function f (x) is not necessarily strongly convex. In fact, under Assumpt ion A, the
optimization problem (1.1) can still have multiple primal o r dual optimal solutions. This
makes the convergence (and rate of convergence) analysis of ADMM diﬃcult.
5

<!-- page 7 -->
2. Assumption ( e) does allow hk(·) to be an indicator function, as in this case the set {xk |∑
k hk(xk) = ∞} ∩ X is not a subset of {x : L(x; y) ≤ ξ} ∩ X for any given y and ξ.
3. Linear term of the form ⟨bk, x k⟩ is already included in hk, as its ephigraph is polyhedral.
Moreover, from the assumption that Xk is polyhedral, the feasibility constraint xk ∈ Xk can
be absorbed into hk by adding to it an indicator function iXk (xk). To simplify notations, we
will not explicitly write xk ∈ Xk in the ADMM update (1.7) from now on.
4. Assumption (f) is made to ensure that the subproblems for e ach xk is strongly convex. This
assumption will be relaxed later when the subproblems are so lved inexactly; see Section 4.1.
5. Assumption (g) requires the feasible set of the variables to be compact, which is needed to
ensure that certain error bounds of the primal and dual probl ems of (1.1) hold. This assump-
tion is usually satisﬁed in practical applications (e.g. th e consensus problems) whenever a
priori knowledge on the variable domain is available. This assumpt ion can be further relaxed;
see the discussion at the end of Section 3.
Under Assumption A, both the primal optimum and the dual opti mum values of (1.1) are attained
and are equal (i.e., the strong duality holds for (1.1)) so th at
d∗ = max
y
L(x; y) = max
y
(
f (x) + ⟨y, q − Ex⟩ + ρ
2 ∥Ex − q∥2
)
= min
Ex=q
f (x),
where d∗ is the optimal value of the dual of (1.1).
Under Assumption A, there may still be multiple optimal solu tions for both the primal problem
(1.1) and its dual problem. We ﬁrst claim that the dual functi onal
d(y) = min
x
L(x; y) = min
x
f (x) + ⟨y, q − Ex⟩ + ρ
2 ∥q − Ex∥2, (2.1)
is diﬀerentiable everywhere. Let X(y) denote the set of optimal solutions for (2.1).
Lemma 2.1 For any y ∈ ℜ m, both Ex and Akxk, k = 1 , 2, ..., K , are constant over X(y). More-
over, the dual function d(y) is diﬀerentiable everywhere and
∇ d(y) = q − Ex(y),
where x(y) ∈ X(y) is any minimizer of (2.1).
Proof. Fix y ∈ ℜ m. We ﬁrst show that Ex is invariant over X(y). Suppose the contrary, so that
there exist two optimal solutions x and x′ from X(y) with the property that Ex ̸= Ex′. Then, we
have
d(y) = L(x; y) = L(x′; y).
6

<!-- page 8 -->
Due to the convexity of L(x; y) with respect to the variable x, the solution set X(y) must be convex,
implying ¯x = (x + x′)/ 2 ∈ X(y). By the convexity of f (x), we have
1
2
[
(f (x) + ⟨y, q − Ex⟩) + (f (x′) + ⟨y, q − Ex′⟩)
]
≥ f (¯x) + ⟨y, q − E ¯x⟩.
Moreover, by the strict convexity of ∥ · ∥2 and the assumption Ex ̸= Ex′, we have
1
2
(
∥Ex − q∥2 + ∥Ex′ − q∥2)
> ∥E ¯x − q∥2.
Multiplying this inequality by ρ/ 2 and adding it to the previous inequality yields
1
2
[
L(x; y) + L(x′; y)
]
> L(¯x; y),
which further implies
d(y) > L(¯x; y).
This contradicts the deﬁnition d(y) = min x L(x; y). Thus, Ex is invariant over X(y). Notice that
d(y) is a concave function and its subdiﬀerential is given by [1]
∂d(y) = Closure of the convex hull { q − Ex(y) |x(y) ∈ X(y) }.
Since Ex(y) is invariant over X(y), the subdiﬀerential ∂d(y) is a singleton. By Danskin’s Theorem,
this implies that d(y) is diﬀerentiable and the gradient is given by ∇ d(y) = q − Ex(y), for any
x(y) ∈ X(y).
A similar argument (and using the strict convexity of gk) shows that Akxk is also invariant over
X(y). The proof is complete. Q.E.D.
By using Lemma 2.1, we show below a Lipschitz continuity prop erty of ∇ d(y), for y over any
level set of d.
Lemma 2.2 Fix any scalar η ≤ f ∗ and let U = { y ∈ ℜ m | d(y) ≥ η }. Then there holds
∥∇ d(y′) − ∇ d(y)∥ ≤ 1
ρ ∥y′ − y∥, ∀ y′ ∈ U , y ∈ U .
Proof. Fix any y and y′ in U . Let x = x(y) and x′ = x(y′) be two minimizers of L(x; y) and
L(x; y′) respectively. By convexity, we have
z − ET y + ρET (Ex − q) = 0 and z′ − ET y′ + ρET (Ex′ − q) = 0 ,
7

<!-- page 9 -->
where z and z′ are some subgradient vectors in the subdiﬀerential ∂f (x) and ∂f (x′) respectively.
Thus, we have
⟨z − ET y + ρET (Ex − q), x ′ − x⟩ = 0
and
⟨z′ − ET y′ + ρET (Ex′ − q), x − x′⟩ = 0.
Adding the above two equalities yields
⟨z − z′ + ET (y′ − y) − ρET E(x′ − x), x ′ − x⟩ = 0.
Upon rearranging terms and using the convexity property
⟨z′ − z, x ′ − x⟩ ≥ 0,
we get
⟨y′ − y, E (x′ − x)⟩ = ⟨z′ − z, x ′ − x⟩ + ρ∥E(x′ − x)∥2 ≥ ρ∥E(x′ − x)∥2.
Thus, ρ∥E(x′ − x)∥ ≤ ∥ y′ − y∥ which together with ∇ d(y′) − ∇ d(y) = E(x − x′) (cf. Lemma 2.1)
yields
∥∇ d(y′) − ∇ d(y)∥ = ∥E(x′ − x)∥ ≤ 1
ρ ∥y − y′∥.
The proof is complete. Q.E.D.
To show the linear convergence of the ADMM method, we need cer tain local error bounds around
the optimal solution set X(y) as well as around the dual optimal solution set Y ∗. To describe these
local error bounds, we ﬁrst deﬁne the notion of a proximity op erator. Let h : dom ( h) ↦→ ℜ be
a (possibly nonsmooth) convex function. For every x ∈ dom (h), the proximity operator of h is
deﬁned as [44]
proxh(x) = argmin
u∈ℜ n
h(u) + 1
2 ∥x − u∥2.
Notice that if h(x) is the indicator function of a closed convex set X, then
proxh(x) = proj X (x),
so the proximity operator is a generalization of the project ion operator. In particular, it is known
that the proximity operator satisﬁes the nonexpansiveness property:
∥proxh(x) − proxh(x′)∥ ≤ ∥ x − x′∥, ∀ x, x ′. (2.2)
The proximity operator can be used to characterize the optim ality condition for a nonsmooth
convex optimization problem. Suppose a convex function f is decomposed as f (x) = g(Ax) + h(x)
8

<!-- page 10 -->
where g is strongly convex and diﬀerentiable, h is a convex (possibly nonsmooth) function, then
we can deﬁne the proximal gradient of f with respect to h as
˜∇ f (x) := x − proxh(x − ∇ (f (x) − h(x))) = x − proxh(x − AT ∇ g(Ax)).
If h ≡ 0, then the proximal gradient ˜∇ f (x) = ∇ f (x). In general, ˜∇ f (x) can be used as the
(standard) gradient of f for the nonsmooth minimization min x∈ X f (x). For example, ˜∇ f (x∗) = 0
iﬀ x∗ is a global minimizer.
For the Lagrangian minimization problem (2.1) and under Ass umption A, the work of [37,50,57]
suggests that the size of the proximal gradient
˜∇ xL(x; y) := x − proxh (x − ∇ x(L(x; y) − h(x)))
= x − proxh
(
x − AT ∇ g(Ax) + ET y − ρET (Ex − q)
)
(2.3)
can be used to upper bound the distance to the optimal solutio n set X(y) of (2.1). Here
h(x) :=
K∑
k=1
hk(xk), g (Ax) :=
K∑
k=1
gk(Akxk)
represent the nonsmooth and the smooth parts of f (x) respectively.
In our analysis of ADMM, we will also need an error bound for th e dual function d(y). Notice
that a y ∈ ℜ m solves (1.5) if and only if y satisﬁes the system of nonlinear equations
∇ d(y) = 0 .
This suggests that the norm of the ‘residual’ ∥∇ d(y)∥ may be a good estimate of how close y is
from solving (1.5). The next lemma says if the nonsmooth part of fk takes certain forms, then the
distance to the primal and dual optimal solution sets can ind eed be bounded.
Lemma 2.3 Suppose assumptions A(a)—A(e) hold.
1. If in addition X is a polyhedral set, then there exists a positive scalar τ and δ such that the
following error bound holds
dist (x, X (y)) ≤ τ∥ ˜∇ xL(x; y)∥, (2.4)
for all (x, y ) such that ∥ ˜∇ xL(x; y)∥ ≤ δ, where the proximal gradient ˜∇ xL(x; y) is given by
(2.3). Furthermore, if X is also a compact set, then there exists some τ > 0 such that the
error bound (2.4) holds for all x ∈ X ∩ dom(h).
9

<!-- page 11 -->
2. Similarly, if assumption A-(g) also holds, then for any sc alar ζ, there exist positive scalars δ
and τ such that
dist (y, Y ∗) = ∥y − y∗∥ ≤ τ∥∇ d(y)∥, whenever d(y) ≥ ζ and ∥∇ d(y)∥ ≤ δ. (2.5)
Moreover, in both cases the constant τ is independent of the choice of y and x.
For any ﬁxed y, the proof for the ﬁrst part of Lemma 2.3 is identical to those of [37, 50, 57],
each of which shows the error bound with diﬀerent assumptions on the objective function f . In
particular, it was shown that (2.4) holds for all x with ∥ ˜∇ xL(x; y)∥ ≤ δ (i.e., suﬃciently close to
X(y)). An important new ingredient is the claim that the error bo und holds over the compact
set X ∩ dom(h). This can be seen in two steps as follows: (1) for all x ∈ X ∩ dom(h) such that
∥ ˜∇ xL(x; y)∥ ≤ δ, the error bound (2.4) is already known to hold; (2) for all x ∈ X ∩ dom(h) such
that ∥ ˜∇ xL(x; y)∥ ≥ δ, the ratio
dist (x, X (y))
∥ ˜∇ xL(x; y)∥
is a continuous function and well deﬁned over the compact set X ∩ dom(h)∩
{
x | ∥˜∇ xL(x; y)∥ ≥ δ
}
.
Thus, the above ratio must be bounded from above by a constant τ′ (independent of y). Combining
(1) and (2) yields the desired error bound over the set X ∩ dom(h).
Another new ingredient in Lemma 2.3 is the additional claim t hat the constants δ, τ are both
independent of the choice of y. This property follows directly from a similar property of H oﬀman’s
error bound [28] (on which the error bounds of [37, 50, 57] are based) for a feasible linear system
P := {x |Ax ≤ b}:
dist (x, P ) ≤ τ∥[Ax − b]+∥, ∀ x ∈ ℜ n,
where τ is independent of b. In fact, a careful checking of the proofs of [37, 50, 57] show s that the
corresponding error constants δ and τ for the augmented Lagrangian function L(x; y) can be indeed
made independent of y. We omit the proof of the ﬁrst part of Lemma 2.3 for space consi deration.
Dual error bounds like the one stated in the second part of the lemma have been studied
previously by Pang [43] and by Mangasarian and Shiau [39], th ough in diﬀerent contexts. The
above error bound is ‘local’ in that it holds only for those y that are bounded or near Y ∗ (i.e., when
∥∇ d(y)∥ ≤ δ as opposed to a ‘global’ error bound which would hold for all y in ℜ m). However if in
addition y also lies in some compact set Y , then the dual error bound hold true for all y ∈ Y (using
the same argument as in the preceding paragraph). In the appe ndix, we include a proof showing
that the dual error bound holds true, for the case where the ep igraph of hk is polyhedral (which
includes ℓ1 norm and indicator function for polyhedral sets). We note th at from this proof it is
clear that indeed the value of τ in the dual error bound does not depend on the choice of either x
or y.
10

<!-- page 12 -->
Under Assumption A(f), the augmented Lagrangian function L(x; y) (cf. (1.3)) is strongly con-
vex with respect to each subvector xk. As a result, each alternating minimization iteration of
ADMM (1.7)
xr+1
k = argmin
xk
L(xr+1
1 , ..., x r+1
k− 1, x k, x r
k+1, ..., x r
K ; yr), k = 1, ..., K.
has a unique optimal solution. Thus the sequence of iterates {xr} of the ADMM are well deﬁned.
The following lemma shows that the alternating minimizatio n of the Lagrangian function gives a
suﬃcient descent of the Lagrangian function value.
Lemma 2.4 Suppose Assumptions A(b) and A(f ) hold. Then ﬁx any index r, we have
L(xr; yr) − L(xr+1; yr) ≥ γ∥xr − xr+1∥2, (2.6)
where the constant γ > 0 is independent of r and yr.
Proof. By assumptions A(b) and A(f) , the augmented Lagrangian func tion
L(x; y) =
K∑
k=1
(fk(xk) + ⟨yk, q k − Ekxk⟩) + ρ
2





K∑
k=1
Ekxk − q





2
is strongly convex in each variable xk and has a uniform modulus ρλ min(ET
k Ek) > 0. Here, the
notation λ min(·) denotes the smallest eigenvalue of a symmetric matrix. Thi s implies that, for each
k, that
L(x; y) − L(x1, .., x k− 1, ¯xk, x k+1, ..., x K ; y) ≥ ρλ min(ET
k Ek)∥xk − ¯xk∥2, (2.7)
for all x, where ¯xk is the minimizer of min xk L(x; y) (when all other variables {xj}j̸=k are ﬁxed).
Fix any index r. For each k ∈ { 1, ..., K }, by ADMM (1.7), xr+1
k is the minimizer of
L(xr+1
1 , ..., x r+1
k− 1, x k, x r
k+1, x r
k+2, ..., x r
K ; yr). It follows from (2.7)
L(xr+1
1 , ..., x r+1
k− 1, x r
k, ..., x r
K ; yr) − L(xr+1
1 , ..., x r+1
k , x r
k+1, ..., x r
K ; yr) ≥ γ∥xr
k − xr+1
k ∥2, ∀ k, (2.8)
where
γ = ρ min
k
λ min(ET
k Ek)
is independent of r and yr. Summing this over k, we obtain the
suﬃcient decrease condition
L(xr; yr) − L(xr+1; yr) ≥ γ∥xr − xr+1∥2.
This completes the proof of Lemma 2.4.
Q.E.D.
To prove the linear convergence of the ADMM algorithm, we als o need the following lemma
which bounds the size of the proximal gradient ˜∇ L(xr; yr) at an iterate xr.
11

<!-- page 13 -->
Lemma 2.5 Suppose assumptions A(b)—A(c) hold. Let {xr} be generated by the ADMM algo-
rithm (1.7). Then there exists some constant σ > 0 ( independent of yr) such that
∥ ˜∇ L(xr; yr)∥ ≤ σ ∥xr+1 − xr∥ (2.9)
for all r ≥ 1.
Proof. Fix any r ≥ 1 and any 1 ≤ k ≤ K. According to the ADMM procedure (1.7), the variable
xk is updated as follows
xr+1
k = argmin
xk

 hk(xk) + gk(Akxk) − ⟨ yr, E kxk⟩ + ρ
2






Ekxk +
∑
j<k
Ejxr+1
j +
∑
j>k
Ejxr
j − q






2
 .
The corresponding optimality condition can be written as
xr+1
k = proxhk

 xr+1
k − AT
k ∇ xkgk(Akxr+1
k ) + ET
k yr − ρET
k

 ∑
j≤ k
Ejxr+1
j +
∑
j>k
Ejxr
j − q



 . (2.10)
Therefore, we have

xr+1
k − proxhk
(
xr
k − AT
k ∇ xk gk(Akxr
k) + ET
k yr − ρET
k (Exr − q)
) 
 =




proxhk

 xr+1
k − AT
k ∇ xk gk(Akxr+1
k ) + ET
k yr + ρET
k

 ∑
j≤ k
Ejxr+1
j +
∑
j>k
Ejxr
j − q




− proxhk
(
xr
k − AT
k ∇ xk gk(Akxr
k) + ET
k yr + ρET
k (Exr − q)
)





≤


(xr+1
k − xr
k) − AT
k (∇ xk gk(Akxr+1
k ) − ∇ xkgk(Akxr
k)) + ρET
k
∑
j≤ k
Ej(xr+1
j − xr
j )



≤ ∥ xr+1
k − xr
k∥ + L∥AT
k ∥∥Ak∥∥xr+1
k − xr
k∥ + ρ∥ET
k ∥
∑
j≤ k
∥Ej∥∥xr+1
j − xr
j ∥
≤ c∥xr+1 − xr∥, for some c > 0 independent of yr, (2.11)
where the ﬁrst inequality follows from the nonexpansive pro perty of the prox operator (2.2), and the
second inequality is due to the Lipschitz property of the gra dient vector ∇ gk (cf. Assumption A-(c)).
Using this relation and the deﬁnition of the proximal gradie nt ˜∇ L(xr; yr), we have
∥ ˜∇ xk L(xr; yr)∥ =

xr
k − proxhk
(
xr
k − AT
k ∇ xk gk(Akxr
k) + ET
k yr − ρET
k (Exr − q)
) 

≤ ∥ xr
k − xr+1
k ∥ +

xr+1
k − proxhk
(
xr
k − AT
k ∇ xk gk(Akxr
k) + ET
k yr − ρET
k (Exr − q)
) 

≤ (c + 1)∥xr+1 − xr∥, ∀ k = 1, 2, ..., K.
12

<!-- page 14 -->
This further implies that the entire proximal gradient vect or can be bounded by ∥xr+1 − xr∥:
∥ ˜∇ L(xr; yr)∥ ≤ (c + 1)
√
K∥xr+1 − xr∥.
Setting σ = (c + 1)
√
K (which is independent of yr) completes the proof. Q.E.D.
3 Linear Convergence of ADMM
Let d∗ denote the dual optimal value and {xr, y r} be the sequence generated by the ADMM method
(1.7). Due to assumption A(a), d∗ also equals to the primal optimal value. Further we denote
∆r
d = d∗ − d(yr) (3.1)
which represents the gap from dual optimality at the r-th iteration. The primal gap to optimality
at iteration r is deﬁned as
∆r
p = L(xr+1; yr) − d(yr), r ≥ 1. (3.2)
Clearly, we have both ∆ r
d ≥ 0 and ∆ r
p ≥ 0 for all r. To establish the linear convergence of ADMM,
we need several lemmas to estimate the sizes of the primal and dual optimality gaps as well as their
respective decrease.
Let X(yr) denote the set of optimal solutions for the following optim ization problem
min
x
L(x; yr) = min
x
f (x) + ⟨yr, q − Ex⟩ + ρ
2 ∥Ex − q∥2.
We denote
¯xr = argmin
¯x∈ X(yr)
∥¯x − xr∥.
We ﬁrst bound the sizes of the dual and primal optimality gaps .
Lemma 3.1 Suppose assumptions A(a)—A(e) and A(g) hold. Then for any sca lar δ > 0, there
exists a positive scalar τ′ such that
∆r
d ≤ τ′∥∇ d(yr)∥2 = τ′∥Ex(yr) − q∥2, (3.3)
for any yr ∈ ℜ m with ∥∇ d(yr)∥ ≤ δ. Moreover, there exist positive scalars ζ and ζ
′
(independent
of yr) such that
∆r
p ≤ ζ∥xr+1 − xr∥2 + ζ
′
∥xr − ¯xr∥2, for all r ≥ 1. (3.4)
13

<!-- page 15 -->
Proof. Fix any yr, and let y∗ be the optimal dual solution closest to yr. Then it follows from the
mean value theorem that there exists some ˜ y in the line segment joining yr and y∗ such that
∆r
d = d(y∗) − d(yr)
= ⟨∇ d(˜y), y ∗ − yr⟩
= ⟨∇ d(˜y) − ∇ d(y∗), y ∗ − yr⟩
≤ ∥∇ d(˜y) − ∇ d(y∗)∥∥y∗ − yr∥
≤ 1
ρ ∥˜y − y∗∥∥y∗ − yr∥
≤ 1
ρ ∥yr − y∗∥∥y∗ − yr∥
= 1
ρ ∥y∗ − yr∥2
where the second inequality follows from Lemma 2.2. Recall f rom the second part in Lemma 2.3
that there exists some τ such that
dist (yr, Y ∗) = ∥yr − y∗∥ ≤ τ∥∇ d(yr)∥.
Combining the above two inequalities yields
∆r
d = d(y∗) − d(yr) ≤ τ′∥∇ d(yr)∥2,
where τ′ = τ2/ρ is a constant. This establishes the bound on the size of dual g ap (3.3).
It remains to prove the bound on the primal gap (3.4). For nota tional simplicity, let us separate
the smooth and nonsmooth part of the augmented Lagrangian as follows
L(x; y) = g(x) + h(x) + ⟨y, q − Ex⟩ + ρ
2 ∥q − Ex∥2 := ¯L(x; y) + h(x).
Let xr+1
k denote the k-th subvector of the primal vector xr+1. From the way that the variables
are updated (2.10), we have
xr+1
k = proxhk
[
xr+1
k − ∇ xk
¯L
(
{xr+1
j≤ k}, {xr
j }j>k ; yr
)]
= proxhk
[
xr
k − ∇ xk
¯L(xr; yr) − xr
k + xr+1
k − ∇ xk
¯L
(
{xr+1
j≤ k}, {xr
j }j>k ; yr
)
+ ∇ xk
¯L(xr; yr)
]
:= proxhk
[
xr
k − ∇ xk
¯L(xr; yr) − er
k
]
(3.5)
where the gradient vector ∇ xk
¯L
(
{xr+1
j≤ k}, {xr
j }j>k ; yr
)
can be explicitly expressed as
∇ xk
¯L
(
{xr+1
j≤ k}, {xr
j }j>k ; yr
)
= AT
k ∇ xk g(Akxr+1
k ) − ET
k yr + ρET
k

 ∑
j≤ k
Ejxr+1
j +
∑
j>k
Ejxr
j − q


14

<!-- page 16 -->
and the error vector er
k is deﬁned by
er
k := xr
k − xr+1
k + ∇ xk
¯L
(
{xr+1
j≤ k}, {xr
j }j>k ; yr
)
− ∇ xk
¯L(xr; yr). (3.6)
Note that we can bound the norm of er
k as follows
∥er
k∥ ≤ ∥ xr
k − xr+1
k ∥ + ∥∇ xk
¯L
(
{xr+1
j≤ k}, {xr
j }j>k ; yr
)
− ∇ xk
¯L(xr; yr)∥
≤ ∥ xr
k − xr+1
k ∥ +






AT
k
(
∇ xkg(Akxr+1
k ) − ∇ xk g(Akxr
k)
)
+ ρET
k

 ∑
j≤ k
Ej(xr+1
j − xr
k)








≤ c∥xr − xr+1∥, (3.7)
where the constant c > 0 is independent of yr, and can take the same value as in (2.11).
Using (3.5), and by the deﬁnition of the proximity operator, we have the following
hk(xr+1
k ) + ⟨xr+1
k − xr
k, ∇ xk
¯L(xr; yr) + er
k⟩ + 1
2 ∥xr+1
k − xr
k∥2
≤ hk(¯xr
k) + ⟨¯xr
k − xr
k, ∇ xk
¯L(xr; yr) + er
k⟩ + 1
2 ∥¯xr
k − xr
k∥2. (3.8)
Summing over all k = 1, · · ·, K , we obtain
h(xr+1) + ⟨xr+1 − xr, ∇ x ¯L(xr; yr) + er⟩ + 1
2 ∥xr+1 − xr∥2
≤ h(¯xr) + ⟨¯xr − xr, ∇ x ¯L(xr; yr) + er⟩ + 1
2 ∥¯xr − xr∥2.
Upon rearranging terms, we obtain
h(xr+1) − h(¯xr) + ⟨xr+1 − ¯xr, ∇ x ¯L(xr; yr)⟩ ≤ 1
2 ∥¯xr − xr∥2 − ⟨ xr+1 − ¯xr, e r⟩. (3.9)
Also, we have from the mean value theorem that there exists so me ˜x in the line segment joining
xr+1 and ¯xr such that
¯L(xr+1; yr) − ¯L(¯xr; yr) = ⟨∇ x ¯L(˜x; yr), x r+1 − ¯xr⟩.
15

<!-- page 17 -->
Using the above results, we can bound ∆ r
p by
∆r
p = L(xr+1; yr) − L(¯xr; yr)
= ¯L(xr+1; yr) − ¯L(¯xr; yr) + h(xr+1) − h(¯xr)
= ⟨∇ x ¯L(˜x; yr), x r+1 − ¯xr⟩ + h(xr+1) − h(¯xr)
= ⟨∇ x ¯L(˜x; yr) − ∇ x ¯L(xr; yr), x r+1 − ¯x⟩ + ⟨∇ x ¯L(xr; yr), x r+1 − ¯x⟩ + h(xr+1) − h(¯xr)
≤ ⟨∇ x ¯L(˜x; yr) − ∇ x ¯L(xr; yr), x r+1 − ¯x⟩ + 1
2 ∥¯xr − xr∥2 + c
√
K∥xr+1 − xr∥∥xr+1 − ¯xr∥
≤
( K∑
k=1
L∥Ak∥T ∥Ak∥ + ρ∥ET E∥
)
∥˜x − xr∥∥xr+1 − ¯xr∥
+ 1
2 ∥¯xr − xr∥2 + c
√
K∥xr+1 − xr∥∥xr+1 − ¯xr∥
≤
( K∑
k=1
L∥Ak∥T ∥Ak∥ + ρ∥ET E∥
)
(
∥xr+1 − xr∥ + ∥¯xr − xr∥
) 2
+ 1
2 ∥¯xr − xr∥2 + c
√
K∥xr+1 − xr∥
(
∥xr+1 − xr∥ + ∥¯xr − xr∥
)
≤ ζ∥xr+1 − xr∥2 + ζ
′
∥¯xr − xr∥2, for some ζ, ζ ′ > 0,
where the ﬁrst inequality follows from (3.9) and (3.7), the s econd inequality is due to the Cauchy-
Schwartz inequality and the Lipschitz continuity of ∇ ¯Lx(x; yr), while the third inequality follows
from the fact that ˜x lies in the line segment joining xr+1 and ¯xr so that ∥˜x − xr∥ ≤ ∥ xr+1 − xr∥ +
∥¯xr − xr∥. This completes the proof. Q.E.D.
We then bound the decrease of the dual optimality gap.
Lemma 3.2 For each r ≥ 1, there holds
∆r
d − ∆r− 1
d ≤ − α (Exr − q)T (E ¯xr − q). (3.10)
Proof. The reduction of the optimality gap in the dual space can be bo unded as follows:
∆r
d − ∆r− 1
d = [ d∗ − d(yr)] − [d∗ − d(yr− 1)]
= d(yr− 1) − d(yr)
= L(¯xr− 1; yr− 1) − L(¯xr; yr)
= [ L(¯xr; yr− 1) − L(¯xr; yr)] + [L(¯xr− 1; yr− 1) − L(¯xr; yr− 1)]
= ( yr− 1 − yr)T (q − E ¯xr) + [L(¯xr− 1; yr− 1) − L(¯xr; yr− 1)]
= − α (Exr − q)T (E ¯xr − q) + [L(¯xr− 1; yr− 1) − L(¯xr; yr− 1)]
≤ − α (Exr − q)T (E ¯xr − q), ∀ r ≥ 1,
16

<!-- page 18 -->
where the last equality follows from the update of the dual va riable yr− 1, and the fact that ¯ xr− 1
minimizes L(·, y r− 1). Q.E.D.
Lemma 3.2 implies that if q − Exr is close to the true dual gradient ∇ d(yr) = q − E ¯xr, then
the dual optimal gap is reduced after each ADMM iteration. Ho wever, since ADMM updates the
primal variable by only one Gauss-Seidel sweep, the primal i terate xr is not necessarily close the
minimizer ¯xr of L(x; yr). Thus, unlike the method of multipliers (for which xr = ¯xr for all r), there
is no guarantee that the dual optimality gap ∆ r
d is indeed reduced after each iteration of ADMM.
Next we proceed to bound the decrease in the primal gap ∆ r
p.
Lemma 3.3 Suppose assumptions A(b) and A(f ) hold. Then for each r ≥ 1, we have
∆r
p − ∆r− 1
p ≤ α ∥Exr − q∥2 − γ∥xr+1 − xr∥2 − α (Exr − q)T (E ¯xr − q) (3.11)
for some γ independent of yr.
Proof. Fix any r ≥ 1, we have
L(xr; yr− 1) = f (xr) + ⟨yr− 1, q − Exr⟩ + ρ
2 ∥Exr − q∥2
and
L(xr+1; yr) = f (xr+1) + ⟨yr, q − Exr+1⟩ + ρ
2 ∥Exr+1 − q∥2.
By the update rule of yr (cf. (1.7)), we have
L(xr; yr) = f (xr) + ⟨yr− 1, q − Exr⟩ + ρ
2 ∥Exr − q∥2 + α ∥Exr − q∥2.
This implies
L(xr; yr) = L(xr; yr− 1) + α ∥Exr − q∥2.
Recall from Lemma 2.4 that the alternating minimization of t he Lagrangian function gives a suﬃ-
cient descent. In particular, we have
L(xr+1; yr) − L(xr; yr) ≤ − γ∥xr+1 − xr∥2,
for some γ > 0 that is independent of r and yr. Therefore, we have
L(xr+1; yr) − L(xr; yr− 1) ≤ α ∥Exr − q∥2 − γ∥xr+1 − xr∥2, ∀ r ≥ 1.
Hence, we have the following bound on the reduction of primal optimality gap
∆r
p − ∆r− 1
p = [ L(xr+1; yr) − d(yr)] − [L(xr; yr− 1) − d(yr− 1)]
= [ L(xr+1; yr) − L(xr; yr− 1)] − [d(yr) − d(yr− 1)]
≤ α ∥Exr − q∥2 − γ∥xr+1 − xr∥2 − α (Exr − q)T (E ¯xr − q), ∀ r ≥ 1,
17

<!-- page 19 -->
where the last step is due to Lemma 3.2. Q.E.D.
Notice that when α = 0 (i.e., no dual update in the ADMM algorithm), Lemma 3.3 red uces to
the suﬃcient decrease estimate (2.6) in Lemma 2.4. When α > 0, the primal optimality gap is not
necessarily reduced after each ADMM iteration due to the pos itive term α ∥Exr − q∥2 in (3.11).
Thus, in general, we cannot guarantee a consistent decrease of either the dual optimality gap ∆ r
d
or the primal optimality gap ∆ r
p. However, somewhat surprisingly, the sum of the primal and d ual
optimality gaps decreases for all r, as long as the dual step size α is suﬃciently small. This is used
to establish the linear convergence of ADMM method.
Theorem 3.1 Suppose the conditions in Assumption A hold. Then the sequenc e of iterates
{(xr, y r)} generated by the ADMM algorithm (1.7) converges linearly to an optimal primal-dual
solution for (1.1), provided the stepsize α is suﬃciently small. Moreover, the sequence of feasibility
violation {∥Exr − q∥} also converges linearly.
Proof. We show by induction that the sum of optimality gaps ∆ r
d + ∆ r
p is reduced after each
ADMM iteration, as long as the stepsize α is chosen suﬃciently small. For any r ≥ 1, we denote
¯xr = argmin
¯x∈ X(yr)
∥¯x − xr∥. (3.12)
By induction, suppose ∆ r− 1
d + ∆r− 1
p ≤ ∆0
d + ∆0
p for some r ≥ 1. Recall that each Xk is compact
and that the indicator function iXk (xk) is included in hk(xk) (see the discussion after Assumption
A), it follows that xr ∈ X, implying the boundedness of xr. Thus, we obtain from Lemma 2.3 that
∥xr − ¯xr∥ ≤ τ∥ ˜∇ L(xr; yr)∥ (3.13)
for some τ > 0 (independent of yr). To prove Theorem 3.1, we combine the two estimates (3.10)
and (3.11) to obtain
[∆r
p + ∆r
d] − [∆r− 1
p + ∆r− 1
d ] = [∆ r
p − ∆r− 1
p ] + [∆r
d − ∆r− 1
d ]
≤ α ∥Exr − q∥2 − γ∥xr+1 − xr∥2 − 2α (Exr − q)T (E ¯xr − q)
= α ∥Exr − E ¯xr∥2 − α ∥E ¯xr − q∥2 − γ∥xr+1 − xr∥2. (3.14)
Now we invoke (3.13) and Lemma 2.5 to lower bound ∥xr+1 − xr∥:
∥xr − ¯xr∥ ≤ τ∥ ˜∇ L(xr; yr)∥ ≤ τ σ∥xr+1 − xr∥. (3.15)
Substituting this bound into (3.14) yields
[∆r
p + ∆r
d] − [∆r− 1
p + ∆r− 1
d ] ≤ (α ∥E∥2τ2σ 2 − γ)∥xr+1 − xr∥2 − α ∥E ¯xr − q∥2. (3.16)
18

<!-- page 20 -->
Thus, if we choose the stepsize α suﬃciently small so that
0 < α < γτ − 2σ − 2∥E∥− 2, (3.17)
then the above estimate shows that
[∆r
p + ∆r
d] ≤ [∆r− 1
p + ∆r− 1
d ], (3.18)
which completes the induction. Moreover, the induction arg ument shows that if the stepsize α
satisﬁes the condition (3.17), then the descent condition ( 3.16) holds for all r ≥ 1.
By the descent estimate (3.16), we have
∥xr+1 − xr∥ → 0, ∥∇ d(yr)∥ = ∥E ¯xr − q∥ → 0. (3.19)
We now show that the sum of optimality gaps ∆ r
d + ∆r
p in fact contracts geometrically after
a ﬁnite number of ADMM iterations. By (3.19), for any δ > 0, there must exist a ﬁnite integer
¯r > 0 such that for all r ≥ ¯r, ∥∇ d(yr)∥ ≤ δ. Since ∆ r
d, ∆ r
p are nonnegative and bounded from
above (see (3.18)), it follows that d(yr) is bounded from below by a constant ζ independent of
r. Applying the second part of Lemma 2.3, we have that for all r ≥ ¯r, the dual error bound
dist(yr, Y ∗) ≤ τ∥∇ d(yr)∥ holds true.
Therefore, it follows from Lemma 3.1 that we have the followi ng cost-to-go estimate
∆r
d = d∗ − d(yr) ≤ τ′∥∇ d(yr)∥2 = τ′∥E ¯xr − q∥2, (3.20)
for some τ′ > 0 and for all r ≥ ¯r.
Moreover, we can use Lemma 3.1 to bound ∥xr+1 − xr∥2 from below by ∆ r
p. In particular, we
have from (3.15) and Lemma 3.1 that
∆r
p ≤ ζ∥xr+1 − xr∥2 + ζ
′
∥¯xr − xr∥2
≤ ζ∥xr+1 − xr∥2 + ζ
′
τ2σ 2∥xr+1 − xr∥2
=
(
ζ + ζ
′
τ2σ 2
)
∥xr+1 − xr∥2.
Substituting this bound and (3.20) into (3.16), and assumin g that α > 0 satisﬁes (3.17), we obtain
[∆r
p + ∆r
d] − [∆r− 1
p + ∆r− 1
d ] ≤ (α ∥E∥2τ2σ 2 − γ)∥xr+1 − xr∥2 − α ∥E ¯xr − q∥2
≤ − (γ − α ∥E∥2τ2σ 2)
ζ + ζ
′
τ2σ 2 ∆r
p − α (τ′)− 1∆r
d
≤ − min
{ (γ − α ∥E∥2τ2σ 2)
ζ + ζ
′
τ2σ 2 , α (τ′)− 1
}
[∆r
p + ∆r
d].
19

<!-- page 21 -->
Since α > 0 is chosen small enough such that (3.17) holds, we have
λ := min
{ γ − α ∥E∥2τ2σ 2
ζ + ζ
′
τ2σ 2 , α (τ′)− 1
}
> 0.
Consequently, we have
[∆r
p + ∆r
d] − [∆r− 1
p + ∆r− 1
d ] ≤ − λ[∆r
p + ∆r
d]
which further implies
0 ≤ [∆r
p + ∆r
d] ≤ 1
1 + λ [∆r− 1
p + ∆r− 1
d ].
This shows that the sequence {∆r
p + ∆r
d}r≥ ¯r converges to zero Q-linearly 1. As a result, we conclude
that {∆r
p + ∆r
d} and hence both ∆ r
p and ∆ r
d globally converge to zero R-linearly 2.
We next show that the dual sequence {yr} is also R-linearly convergent. To this end, notice
that the inequalities (3.15) and (3.16) imply
[∆r
p + ∆r
d] − [∆r− 1
p + ∆r− 1
d ] ≤ (α ∥E∥2 − γτ − 2σ − 2)∥xr − ¯xr∥2 − α ∥E ¯xr − q∥2. (3.21)
Then by (3.21), we see that both ∥xr − ¯xr∥ → 0 and ∥E ¯xr − q∥ → 0 R-linearly. This implies that
Exr − q → 0 R-linearly and ∇ d(yr) → 0 R-linearly. Using the fact that dist( yr, Y ∗) ≤ τ∥∇ d(yr)∥,
we conclude that yr converges R-linearly to an optimal dual solution.
We now argue that the primal iterates {xr} converge to an optimal solution of (1.1). By the
inequality (3.16), we can further conclude that
∥xr+1 − xr∥2 → 0, ∥E ¯xr − q∥ → 0
R-linearly. Notice that the R-linear convergence of ∥xr+1 − xr∥2 → 0 implies that ∥xr+1 − xr∥ → 0
R-linearly. This further shows that xr → x∞ R-linearly for some x∞ . Denote the limit of dual
sequence {yr} by y∞ . By the preceding argument, we know y∞ is a dual optimal solution of (1.1).
To show that x∞ is a primal optimal solution of (1.1), it suﬃces to prove that x∞ ∈ X(y∞ ). Using
(3.15), and the fact that ∥xr − ¯xr∥ → 0, we have
∥x∞ − ¯xr∥ ≤ ∥ xr − x∞ ∥ + ∥xr − ¯xr∥ → 0.
Since ¯xr ∈ X(yr), we have L(¯xr, y r) ≤ L(x, y r) for all x ∈ X. Passing limit, we obtain L(x∞ , y ∞ ) ≤
L(x, y ∞ ) for all x ∈ X, that is, x∞ ∈ X(y∞ ). It then follows that the sequence {xr} converges
R-linearly to a primal optimal solution. Q.E.D.
1A sequence {xr} is said to converge Q-linearly to some ¯x if ∥xr+1 − ¯x∥/∥xr − ¯x∥ ≤ µ for all r, where µ ∈ (0, 1)
is some constant. A sequence {xr} is said to converge to ¯ x R -linearly if ∥xr − ¯x∥ ≤ cµr for all r and for some c > 0.
2To see that such R-linear convergence is in fact global, note that ¯r > 0 is ﬁnite, and ∆ r
p + ∆ r
d is Q-linearly
convergent for r ≥ ¯r. Then one can always ﬁnd an appropriate constant c such that ∆ r
p + ∆r
d ≤ c(1 + λ)− r for all
r = 1, 2, . . . .
20

<!-- page 22 -->
The following corollary relaxes the compactness assumptio n of the feasible set X (Assumption
A-(g) in Theorem 3.1), and replaces it by the boundedness of t he primal-dual iterates.
Corollary 3.1 Suppose assumptions A(a)—A(f ) hold, and that X is a polyhedral set. Further
assume that either one of the following two assumptions hold s true
1. The sequence of dual iterates {yr} lies in a compact set, and that the set {x |L(x, y ) ≤ ζ} is
compact for any ﬁnite y and ζ.
2. The sequence of primal-dual iterates {(xr, y r)} lies in a compact set.
Then if α is chosen suﬃciently small (cf. (3.17)), all the conclusions stated in Theorem 3.1 still
hold true. Moreover, the sequence of function values {f (xr)} also converges linearly.
Proof. Suppose the ﬁrst assumption is true, so that all the dual iter ates {yr} lie in a compact set
˜Y . Let us deﬁne
δ := sup
{
∥x∥ | [L(x, y ) − d(y)] + [d∗ − d(y)] ≤ ∆0
p + ∆0
p, ∀ y ∈ ˜Y
}
.
Clearly we have [ d∗ − d(y)] ≥ 0 and [ L(x, y ) − d(y)] ≥ 0 for all y ∈ ˜Y . Moreover δ < ∞ due to
the compactness of ˜Y as well as the compactness of the set {x |L(x, y ) ≤ ζ} for any ﬁnite y.
Again we show by induction that the sum of optimality gaps ∆ r
d + ∆ r
p is reduced after each
ADMM iteration, as long as the stepsize α is chosen suﬃciently small. For any r ≥ 1, by induction
we suppose ∆ r− 1
d + ∆r− 1
p ≤ ∆0
d + ∆0
p. Then ∥xr∥ ≤ δ < ∞ , and by the ﬁrst part of Lemma 2.3,
the primal error bound holds. Then we can carry out exactly th e same analysis as the proof of
Theorem 3.1 to arrive at the same conclusion.
Additionally, since
f (xr+1) − d∗ = [ f (xr+1) − d(yr)] + [d(yr) − d∗]
= [ L(xr+1; yr) − d(yr)] − [d∗ − d(yr)] − ⟨ yr, q − Exr+1⟩ − ρ
2 ∥Exr+1 − q∥2
= ∆ r
p − ∆r
d − ⟨ yr, q − Exr+1⟩ − ρ
2 ∥Exr+1 − q∥2
and
∆r
p → 0, ∆r
d → 0, ∥q − Exr∥ → 0,
linearly, it follows that f (xr+1) − d∗ → 0 R-linearly (recall that yr is now in a compact set). The
proof is complete.
21

<!-- page 23 -->
Similarly, when the second assumption is true, then the erro r bound condition in Lemma 2.3
again holds, and by using the same argument above we arrive at the desired conclusion. Q.E.D.
As a remark, we point out that the proof of Corollary 3.1 also s hows that the same linear
convergence of f (xr) → d∗ also holds under the assumptions in Theorem 3.1.
We close this section by providing a few examples that satisf y the assumptions in Theorem 3.1
and Corollary 3.1. First consider the following ℓ1 minimization problem
min
x
∥x∥1, s. t. Ex = b, a ≤ xk ≤ b, k = 1, · · ·, K (3.22)
which can be equivalently written as a K-block problem
min
{xk}
K∑
k=1
|xk|, s. t.
K∑
k=1
ekxk = b, a ≤ xk ≤ b, k = 1, · · ·, K, (3.23)
where ek is the k-th column of E, a and b are some scalars. It is easy to verify that this problem
meets all the conditions listed in Assumption A, hence the li near convergence result in Theorem 3.1
applies. The same is true for the following mixed ℓ1/ℓ 2 minimization problem
min
{xk}
K∑
k=1
∥xk∥, s. t.
K∑
k=1
Ekxk = b, a ≤ xk ≤ b, k = 1, · · ·, K, (3.24)
where xk, a and b are now n-dimensional vectors, and Ek ∈ ℜ m× n is some matrix with full column
rank.
Furthermore, the boundedness assumption in Corollary 3.1 i s satisﬁed by many examples of the
two-block ADMM described in [7]. This is because when K = 2, α/β ∈ (0, 1
2 (1 +
√
5)) and Ek’s are
full column rank, it is known that both the primal and dual ite rates generated by the two-block
ADMM algorithm indeed lie in a bounded set [19]. Therefore th e second assumption made in the
Corollary 3.1 holds true, hence we only require assumptions A(a)–A(f), which are in fact quite
mild. For example, they are met by the following instance of t he consensus problem (see [7, Section
7] for introduction of the consensus problem)
min
{xk},z
∑
k=1,···,K
∥Axk − b∥2 + w∥xk∥1 s. t. x k − z = 0, k = 1, · · ·, K, (3.25)
where w > 0 is some constant. Thus, the two block ADMM algorithm conver ges linearly for (3.25)
regardless of the rank of A. Note that when A is full column rank, the objective is strongly convex.
Consequently the error bound condition in Lemma 2.3 holds tr ue globally, and the coeﬃcient τ can
be at least greater than λ min(AT A).
22

<!-- page 24 -->
4 Variants of ADMM
The convergence analysis of Section 3 can be extended to some variants of the ADMM. We brieﬂy
describe two of them below.
4.1 Proximal ADMM
In the original ADMM (1.7), each block xk is updated by solving a convex optimization subproblem
exactly. For large scale problems, this subproblem may not be easy to solve unless the matrix Ek
is unitary (i.e., ET
k Ek = I) in which case the variables in xk can be further decoupled (assuming fk
is separable). If the matrix Ek is not unitary, we can still employ a simple proximal gradien t step
to inexactly minimize L(xr+1
1 , ..., x r+1
k− 1, x k, x r
k+1, ..., x r
K ). More speciﬁcally, we update each block of
xk according to the following procedure
xr+1
k = argmin
xk
{
hk(xk) + ⟨yr, q − Ekxk⟩ +
⟨
AT
k ∇ gk(Akxr
k), x k − xr
k
⟩
+ β
2 ∥xk − xr
k∥2
+
⣨
ρET
k
(∑
j<k
Ejxr+1
j +
∑
j≥ k
Ejxr
j − q
)
, x k − xr
k
⟩}
(4.1)
in which the smooth part of the objective function in the k-th subproblem, namely,
gk(Akxk) + ⟨yr, q − Ekxk⟩ + ρ
2


Ekxk +
∑
j<k
Ejxr+1
j +
∑
j>k
Ejxr
j − q



2
is linearized locally at xr
k, and a proximal term β
2 ∥xk − xr
k∥2 is added. Here, β > 0 is a positive
constant. With this change, updating xk is easy when hk (the nonsmooth part of fk) is separable.
For example, this is the case for compressive sensing applic ations where hk(xk) = ∥xk∥1, and the
resulting subproblem admits a closed form solution given by the component-wise soft thresholding
(also known as the shrinkage operator). We note that the prox imal ADMM algorithm described
here is slightly more general than the proximal ADMM algorit hm seen in the literature, in which
only the penalization term ρ
2


Ekxk + ∑
j<k Ejxr+1
j + ∑
j>k Ejxr
j − q



2
is linearized locally at xr
k;
see e.g., [54, 55].
We claim that Theorem 3.1 holds for the proximal ADMM algorit hm without requiring assump-
tion A-(f ) (the full rankness of Ek’s). Indeed, to establish the (linear) convergence of the pr oximal
ADMM (4.1), we can follow the same proof steps as that for Theo rem 3.1, with the only changes
being in the proof of Lemmas 2.4-2.5 and Lemma 3.1. We ﬁrst sho w that Lemma 2.4 holds without
assumption A(f). Clearly subproblem (4.1) is now strongly convex without the full column rank
assumption of Ek’s made in A(f). In the following, we will show that as long as β is large enough,
23

<!-- page 25 -->
there is a suﬃcient descent:
L(xr+1; yr) − L(xr; yr) ≤ − γ∥xr+1 − xr∥2, for some γ > 0 independet of yr. (4.2)
This property can be seen by bounding the smooth part of L(xr+1
1 , ..., x r+1
k− 1, x k, x r
k+1, ..., x r
K ), which
is given by
¯Lk(xk) := gk(Akxk) + ⟨yr, q − Ekxk⟩ + ρ
2



∑
j<k
Ejxr+1
j +
∑
j>k
Ejxr
j + Ekxk − q



2
,
with the Taylor expansion at xr
k:
¯Lk(xr+1
k ) ≤ ¯Lk(xr
k) + ⟨∇ ¯Lk(xr
k), x r+1
k − xr
k⟩ + ν
2 ∥xr+1
k − xr
k∥2 (4.3)
where
ν := L∥Ak∥∥AT
k ∥ + ρ∥ET
k Ek∥
is the Lipschitz constant of ¯Lk(·) and L is the Lipschitz constant of ∇ gk(·). Making the above
inequality more explicit yields
L(xr+1
1 , ..., x r+1
k− 1, x r+1
k , x r
k+1, ..., x r
K ; yr) − L(xr+1
1 , ..., x r+1
k− 1, x r
k, x r
k+1, ..., x r
K ; yr)
≤ hk(xr+1
k ) − hk(xr
k) + ⟨yr, E k(xr
k − xr+1
k )⟩ +
⟨
AT
k ∇ gk(Akxr
k), x r+1
k − xr
k
⟩
+
⟨
ρET
k

 ∑
j<k
Ejxr+1
j +
∑
j≥ k
Ejxr
j − q

 , x r+1
k − xr
k
⟩
+ ν
2 ∥xr+1
k − xr
k∥2
≤ − β
2 ∥xr+1
k − xr
k∥2 + ν
2 ∥xr+1
k − xr
k∥2
= − γ∥xr+1
k − xr
k∥2, ∀ k, (4.4)
provided the regularization parameter β satisﬁes
γ := 1
2 (β − ν) > 0.
In the above derivation of (4.4), the ﬁrst step is due to (4.3) , while the second inequality follows
from the deﬁnition of xr+1
k (cf. (4.1)). Summing (4.4) over all k yields the desired estimate of
suﬃcient descent (4.2).
To verify that Lemma 2.5 still holds for the proximal ADMM alg orithm, we note from the
corresponding optimality condition for (4.1)
xr+1
k = proxhk

 xr+1
k − AT
k ∇ xk gk(Akxr
k) + ET
k yr − ρET
k

 ∑
j<k
Ejxr+1
j +
∑
j≥ k
Ejxr
j − q

 − β (xr+1
k − xr
k)

 .
24

<!-- page 26 -->
Using this relation in place of (2.10) and following the same proof steps, we can easily prove that
the bound (2.9) in Lemma 2.5 can be extended to the proximal AD MM algorithm. Thus, the
convergence results in Theorem 3.1 remain true for the proxi mal ADMM algorithm (4.1).
It remains to verify that Lemma 3.1 still holds true. In fact t he ﬁrst part of Lemma 3.1 can
be shown to be independent of the iterates, thus it trivially holds true for the proximal ADMM
algorithm. To show that the second part of Lemma 3.1 is true, n ote that the optimality condition
of the proximal ADMM algorithm implies that
xr+1
k = proxhk
[
xr+1
k − ∇ xk
¯L
(
{xr+1
j<k }, {xr
j }j≥ k; yr
)
− β (xr+1
k − xr
k)
]
:= proxhk
[
xr
k − ∇ xk
¯L(xr; yr) − er
k
]
where in this case er
k is given as
er
k := xr
k − xr+1
k + ∇ xk
¯L
(
{xr+1
j<k }, {xr
j }j≥ k; yr
)
− ∇ xk
¯L(xr; yr) + β (xr+1
k − xr
k).
It is then straightforward to show that the norm of er
k can be bounded by c
′
∥xr − xr+1∥ for some
constant c
′
> 0. The rest of the proof follows the same steps as in Lemma 3.1.
4.2 Jacobi Update
Another popular variant of the ADMM algorithm is to use a Jaco bi iteration (instead of a Gauss-
Seidel iteration) to update the primal variable blocks {xk}. In particular, the ADMM iteration
(1.7) is modiﬁed as follows:
xr+1
k = argmin
xk

 hk(xk) + gk(Akxk) − ⟨ yr, E kxk⟩ + ρ
2






Ekxk +
∑
j̸=k
Ejxr
j − q






2
 , ∀ k. (4.5)
The convergence for this direct Jacobi scheme is unclear, as the augmented Lagrangian func-
tion may not decrease after each Jacobi update. In the follow ing, we consider a modiﬁed Jacobi
scheme with an explicit stepsize control. Speciﬁcally, let us introduce an intermediate variable
w = (wT
1 , · · ·, w T
K )T ∈ ℜ n. The modiﬁed Jacobi update is given as follows:
wr+1
k = argmin
xk

 hk(xk) + gk(Akxk) − ⟨ yr, E kxk⟩ + ρ
2






Ekxk +
∑
j̸=k
Ejxr
j − q






2
 , ∀ k, (4.6)
xr+1
k = xr
k + 1
K
(
wr+1
k − xr
k
)
, ∀ k. (4.7)
where a stepsize of 1 /K is used in the update of each variable block.
25

<!-- page 27 -->
With this modiﬁcation, we claim that Lemmas 2.4-2.5 and Lemm a 3.1 still hold. In particular,
Lemma 2.4 can be argued as follows. The strong convexity of L(x; y) with respect to the variable
block xk implies that
L
(
xr
1, · · ·, x r
k− 1, x r
k, x r
k+1 · · ·, x r
K ; yr)
− L
(
xr
1, · · ·, x r
k− 1, w r
k, x r
k+1 · · ·, x r
K ; yr)
≥ γ∥wr+1
k − xr
k∥2, ∀ k.
Using this inequality we obtain
L(xr; yr) − L(xr+1; yr)
= L(xr; yr) − L
(K − 1
K xr + 1
K wr+1; yr
)
= L(xr; yr) − L
(
1
K
K∑
k=1
(xr
1, · · ·, x r
k− 1, w r+1
k , x r
k+1 · · ·, x r
K); yr
)
≥ L(xr; yr) − 1
K
K∑
k=1
L
(
xr
1, · · ·, x r
k− 1, w r+1
k , x r
k+1 · · ·, x r
K ; yr)
= 1
K
K∑
k=1
(
L(xr; yr) − L
(
xr
1, · · ·, x r
k− 1, w r+1
k , x r
k+1 · · ·, x r
K ; yr))
≥ γ
K
K∑
k=1
∥wr+1
k − xr
k∥2
= γ
K ∥wr+1 − xr∥2.
where the ﬁrst inequality comes from the convexity of the aug mented Lagrangian function.
From the update rule (4.7) we have K(xr+1
k − xr
k) = ( wr+1
k − xr
k), which combined with the
previous inequality yields
L(xr; yr) − L(xr+1; yr) ≥ γK ∥xr+1 − xr∥2.
The proof of Lemma 2.5 also requires only minor modiﬁcations . In particular, we have the
following optimality condition for (4.5)
wr+1
k = proxhk

 wr+1
k − AT
k ∇ xkgk(Akwr+1
k ) + ET
k yr − ρET
k

 ∑
j̸=k
Ejxr
j + Ekwr+1
k − q




Similar to the proof of Lemma 2.5, we have

wr+1
k − proxhk
[
xr
k − AT
k ∇ xk gk(Akxr
k) + ET
k yr − ρET
k (Exr − q)
]
 ≤ c∥wr+1 − xr∥.
26

<!-- page 28 -->
Utilizing the relationship K(xr+1
k − xr
k) = ( wr+1
k − xr
k), we can establish Lemma 2.5 by following
similar proof steps (which we omit due to space reason).
Lemma 3.1 can be shown as follows. We ﬁrst express wr+1
k as
wr+1
k = proxhk
[
wr+1
k − ∇ xk
¯L
(
{xr
j̸=k}, w r+1
k ; yr)]
= proxhk
[
xr
k − ∇ xk
¯L (xr; yr) − er
k
]
where we have deﬁned
er
k := ∇ xk
¯L
(
{xr
j̸=k}, w r+1
k ; yr)
− ∇ xk
¯L (xr; yr) + xr
k − wr+1
k .
Again by using the relationship K(xr+1
k − xr
k) = ( wr+1
k − xr
k), we can bound the norm of er
k by
c
′
∥xr+1 − xr∥, for some c
′
> 0. The remaining proof steps are similar to those in Lemma 3.1 .
Since Lemmas 2.4-2.5 and Lemma 3.1 hold for the Jacobi versio n of the ADMM algorithm with
a step size control, we conclude that the convergence result s of Theorem 3.1 remain true in this
case.
5 Concluding Remarks
In this paper we have established the convergence and the rat e of convergence of the classical
ADMM algorithm when the number of variable blocks are more th an two and in the absence of
strong convexity. Our analysis is a departure of the convent ional analysis of ADMM algorithm
which relies on a contraction argument involving a weighted (semi-)norm of ( xr − x∗, y r − y∗),
see [17–20, 25, 26, 30, 47]. In our analysis, we require neith er the strong convexity of the objective
function nor the row independence assumption of the constra ined matrix E. Instead, we use a
local error bound to show that when the stepsize of dual updat e is made suﬃciently small, the
sum of the primal and the dual optimality gaps decreases afte r each ADMM iteration, although
separately they may individually increase. An interesting issue for further research is to identify
good practical stepsize rules for dual update. While (3.17) does suggest a dual stepsize rule in
terms of error bound constants, it may be too conservative an d cumbersome to compute unless the
objective function is strongly convex. One possibility may be to use an adaptive dual stepsize rule
to guarantee the decrease of the sum of the primal and dual opt imality gaps.
Acknowledgement: The authors are grateful to Xiangfeng Wang and Dr. Min Tao of N anjing
University for their constructive comments.
27

<!-- page 29 -->
6 Appendix
6.1 Proof of Dual Error Bound (2.5)
The augmented Lagrangian dual function can be expressed as
d(y) = min
x∈ X
⟨y, q − Ex⟩ + ρ
2 ∥q − Ex∥2 + g(Ax) + h(x). (6.1)
For convenience, deﬁne p(Ex) := ρ
2 ∥q − Ex∥2, and let ℓ(x) := p(Ex) +g(Ax) +h(x). For simplicity,
in this proof we further restrict ourselves to the case where the nonsmooth part has polyhedral level
sets, i.e., {x : h(x) ≤ ξ} is polyhedral for each ξ. More general cases can be shown along similar
lines, but the arguments become more involved.
Let us deﬁne
x(y) ∈ arg min
x∈ X
ℓ(x) + ⟨y, q − Ex⟩.
Let ( x∗, y ∗) denote a primal and dual optimal solution pair. Let X ∗ and Y ∗ denote the primal
and dual optimal solution set. The he following properties w ill be useful in our subsequent analysis.
(a) There exist positive scalars σg, Lg such that ∀ x(y), x (y′) ∈ X
a-1) ⟨AT ∇ g(Ax(y′)) − AT ∇ g(Ax(y)), x (y′) − x(y)⟩ ≥ σg∥Ax(y′) − Ax(y)∥2.
a-2) g(Ax(y′)) − g(Ax(y)) − ⟨ AT ∇ g(Ax(y)), x (y′) − x(y)⟩ ≥ σg
2 ∥Ax(y′) − Ax(y)∥2.
a-3) ∥AT ∇ g(Ax(y′)) − AT ∇ g(Ax(y))∥ ≤ Lg∥Ax(y′) − Ax(y)∥.
b) All a-1)–a-3) are true for p(·) as well, with some constants σp and Lp.
c) ∇ d(y) = q − Ex(y), and ∥∇ d(y′) − ∇ d(y)∥ ≤ 1
ρ ∥y′ − y∥.
Part (a) is true due to the assumed Lipchitz continuity and st rong convexity of the function
g(·). Part (b) is from the Lipchitz continuity and strong convex ity of the quadratic penalization
p(·). Part (c) has been shown in Lemma 2.1 and Lemma 2.2.
To proceed, let us rewrite the primal problem equivalently a s
d(y) = min
(x,s):x∈ X,h(x)≤ s
⟨y, q − Ex⟩ + p(Ex) + g(Ax) + s. (6.2)
Let us write the polyhedral set {(x, s ) : x ∈ X, h (x) ≤ s} compactly as Cxx + Css ≥ c for some
matrices Cx ∈ Rj× n, Cs ∈ Rj× 1 and c ∈ Rj× 1, where j is some integer. For any ﬁxed y, let
(x(y), s (y)) denote one optimal solution for (6.2), note we must have h(x(y)) = s(y). Due to
equivalence, if y∗ ∈ Y ∗, we must also have x(y∗) ∈ X ∗.
28

<!-- page 30 -->
Deﬁne a set-valued function M that assigns the vector ( d, e ) ∈ Rn × Rm to the set of vectors
(x, s, y, λ ) ∈ Rn × R × Rm × Rj that satisfy the following system of equations
ET y + C T
x λ = d,
C T
s λ = 1,
q − Ex = e,
λ ≥ 0, (Cxx + Css) ≥ c, ⟨Cxx + Css − c, λ ⟩ = 0.
It is easy to verify by using the optimality condition for pro blem (6.2) that
(x, s, y, λ ) ∈ M (ET ∇ p(Ex) + AT ∇ g(Ax), e ) for some λ
if and only if x = x(y), e = ∇ d(y). (6.3)
We can take e = 0, and use the fact that x(y∗) ∈ X ∗, we see that ( x, s, y, λ ) ∈ M (ET ∇ p(Ex) +
AT ∇ g(Ax), 0) if and only if x ∈ X ∗ and y ∈ Y ∗.
The following result states a well-known local upper Lipsch itzian continuity property for the
polyhedral multifunction M; see [28, 36, 37].
Proposition 6.1 There exists a positive scalar θ that depends on A, E, C x, C s only, such that for
each ( ¯d, ¯e) there is a positive scalar δ′ satisfying
M(d, e ) ⊆ M ( ¯d, ¯e) + θ∥(d, e ) − ( ¯d, ¯e)∥B, (6.4)
whenever ∥(d, e ) − ( ¯d, ¯e)∥ ≤ δ′. (6.5)
where B denotes the unit Euclidean ball in Rn × Rm × R × Rj.
The following is the main result for this appendix. Note that the scalar τ in the claim is
independent the choice of y, x, s, and is independent on the coeﬃcients of the linear term s.
Claim 6.1 Suppose all the assumptions in Assumption A are satisﬁed. The n there exits positive
scalars δ, τ such that for all y ∈ U and ∥∇ d(y)∥ ≤ δ, there holds dist(y, Y ∗) ≤ τ∥∇ d(y)∥.
Proof. By the previous claim, M is locally Lipschitzian with modulus θ at ( ∇ ℓ(x∗), 0) =
(ET ∇ p(Ex∗) + AT ∇ g(Ax∗), 0).
Let δ ≤ δ′/ 2. We ﬁrst show that if ∥∇ d(y)∥ ≤ δ, then we must have ∥∇ ℓ(x(y))−∇ ℓ(x∗)∥ ≤ δ′/ 2.
To this end, take a sequence y1, y 2, · · ·, such that er := ∇ d(yr) → 0. By assumption A(g) {x(yr)}
lies in a compact set. Due to the fact that s(yr) = h(x(yr)), so the sequence {s(yr)} also lies in
29

<!-- page 31 -->
a compact set (cf. Assumption A(e)). By passing to a subseque nce if necessary, let ( x∞ , s ∞ ) be
a cluster point of {x(yr), s (yr)}. In light of the continuity of ∇ ℓ(·), we have ( ∇ ℓ(x(yr)), e r) →
(∇ ℓ(x∞ ), 0). Now for all r, {(x(yr), s (yr), ∇ ℓ(x(yr)), e r)} lies in the set
{(x, s, d, e ) |(x, s, y, λ ) ∈ M (d, e ) for some ( y, λ )}
which is polyhedral and thus is closed. Then we can pass limit to it and conclude (cf. Proposition
6.1)
(x∞ , s ∞ , y ∞ , λ ∞ ) ∈ M (∇ ℓ(x∞ ), 0)
for some ( y∞ , λ ∞ ) ∈ Rm × Rj. Thus by (6.3) and the discussions that follow, we have x∞ ∈ X ∗
and y∞ ∈ Y ∗. By Lemma 2.1, we have ∇ ℓ(x∗) = ∇ ℓ(x∞ ), which further implies that ∇ ℓ(x(yr)) →
∇ ℓ(x∗). This shows that the desired δ exists.
Then we let e = ∇ d(y), and suppose ∥e∥ ≤ δ. From the previous argument we have
∥∇ ℓ(x(y)) − ∇ ℓ(x∗)∥ + ∥e∥ ≤ δ′/ 2 + δ′/ 2 = δ′.
Using the results in Proposition 6.1, we have that there exis ts ( x∗, s ∗, y ∗, λ ∗) ∈ M (∇ ℓ(x∗), 0)
satisfying
∥(x(y), s, y, λ ) − (x∗, s ∗, y ∗, λ ∗)∥ ≤ θ (∥∇ ℓ(x∗) − ∇ ℓ(x(y))∥ + ∥e∥) .
Since ( x(y), s, y, λ ) ∈ M (∇ ℓ(x(y)), e ), it follows from the deﬁnition of M that
ET y + C T
x λ = ∇ ℓ(x(y)), (6.6)
C T
s λ = 1, (6.7)
q − Ex(y) = e, (6.8)
λ ≥ 0, (Cxx(y) + Css(y)) ≥ c, ⟨Cxx(y) + Css(y) − c, λ ⟩ = 0. (6.9)
Since ( x∗, s ∗, y ∗, λ ∗) ∈ M (∇ ℓ(x∗), 0), we have from the deﬁnition of M
ET y∗ + C T
x λ ∗ = ∇ ℓ(x∗), (6.10)
C T
s λ ∗ = 1, (6.11)
q − Ex∗ = 0, (6.12)
λ ∗ ≥ 0, (Cxx∗ + Css∗) ≥ c, ⟨Cxx∗ + Css∗ − c, λ ∗⟩ = 0. (6.13)
Moreover, we have
σg∥A(x(y) − x∗)∥2 + σp∥E(x(y) − x∗)∥2
≤ ⟨ AT ∇ g(Ax(y)) − AT ∇ g(Ax(y∗)), x (y) − x(y∗)⟩ + ⟨ET ∇ p(Ex(y)) − ET ∇ p(Ex(y∗)), x (y) − x(y∗)⟩
= ⟨∇ ℓ(x(y)) − ∇ ℓ(x(y∗)), x (y) − x(y∗)⟩
= ⟨λ − λ ∗, C xx(y) − Cxx∗⟩ + ⟨y − y∗, Ex (y) − Ex∗⟩
30

<!-- page 32 -->
where the ﬁrst inequality comes from the strong convexity of g(·) and p(·); the last equality is from
(6.6) and (6.10). Moreover, we have
⟨λ − λ ∗, C xx(y) − Cxx∗⟩
= ⟨λ − λ ∗, C xx(y) − Cxx∗⟩ + ⟨λ − λ ∗, C ss − Css∗⟩
= ⟨λ − λ ∗, (Cxx(y) + Css) − (Cxx∗ + Css∗)⟩
= −⟨ λ ∗, C xx(y) + Css − c⟩ − ⟨ λ, C xx∗ + Css∗ − c⟩ ≤ 0 (6.14)
where in the ﬁrst equality we have used the fact that C T
s λ − C T
s λ ∗ = 0; see (6.7) (6.11); in the third
equality and in the last inequality we have used the compleme ntary conditions (6.13) and (6.9). As
a result, we have
σg∥A(x(y) − x∗)∥2 + σp∥E(x(y) − x∗)∥2
≤ ⟨ y − y∗, (Ex(y) − q) − (Ex∗ − q)⟩ ≤ ∥ y − y∗∥∥e∥, (6.15)
where the last step is due to ∇ d(y) = Ex(y) − q and ∇ d(y∗) = Ex∗ − q = 0. Finally we have from
Proposition 6.1
∥(x(y), s, y, λ ) − (x∗, s ∗, y ∗, λ ∗)∥2
≤ θ2 (∥∇ ℓ(x∗) − ∇ ℓ(x(y))∥ + ∥e∥)2
≤ θ2 (
2∥∇ ℓ(x∗) − ∇ ℓ(x(y))∥2 + 2∥e∥2)
≤ 2θ2 (
2∥∇ g(x∗) − ∇ g(x(y))∥2 + 2∥∇ p(x∗) − ∇ p(x(y))∥2 + ∥e∥2)
≤ 2θ2 (
L2
g∥AT (x(y) − x∗)∥2 + L2
p∥ET (x(y) − x∗)∥2 + ∥e∥2)
≤ 2θ2 max
(
2L2
g
σg
, 2L2
p
σp
, 1
)
(
σg∥AT (x(y) − x∗)∥2 + σp∥ET (x(y) − x∗)∥2 + ∥e∥2)
≤ 2θ2 max
(
2L2
g
σg
, 2L2
p
σp
, 1
)
(
∥e∥∥y − y∗∥ + ∥e∥2)
≤ 2θ2 max
(
2L2
g
σg
, 2L2
p
σp
, 1
)
(
∥e∥∥(x(y), s, y, λ ) − (x∗, s ∗, y ∗, λ ∗)∥ + ∥e∥2)
,
where the second inequality is due to ∇ ℓ(x) = ∇ g(x) + ∇ p(x) and the fourth step follows from
properties a-3) and b).
We see that the above inequality is quadratic in ∥(x(y), s, y, λ )− (x∗ , s ∗, y ∗, λ ∗)∥/ ∥e∥, so we have
∥(x(y), s, y, λ ) − (x∗, s ∗, y ∗, λ ∗)∥/ ∥e∥ ≤ τ
for some scalar τ depending on θ, Lg, Lp, σg, σp. It is worth noting that τ does not depend on the
choice of the coeﬃcients of the linear term s. We conclude dist( y, Y ∗) ≤ τ∥∇ d(y)∥. Q.E.D.
31

<!-- page 33 -->
References
[1] Bertsekas, D.P.: Nonlinear Programming. Athena Scient iﬁc (2010).
[2] Bertsekas, D. P. and Gafni, E.: Projection methods for va riational inequalities with ap-
plication to the traﬃc assignment problem. Math. Prog. Stud y. 17, 139–159 (1982)
[3] Bertsekas, D. P. and Gallager, R.: Data Networks. Prenti ce–Hall, Englewood Cliﬀs, New
Jersey. (1987)
[4] Bertsekas, D. P., Hosein, P. A., and Tseng, P.: Relaxatio n methods for network ﬂow
problems with convex arc costs. SIAM J. Control Optim. 25, 12 19–1243 (1987)
[5] Bertsekas, D. P. and Tsitsiklis, J. N.: Parallel and Dist ributed Computation: Numerical
Methods. Prentice–Hall, Englewood Cliﬀs, New Jersey (1989)
[6] Boley, D.: Linear convergence of ADMM on a model problem. TR 12-009, Department of
Computer Science and Engineering, University of Minnesota (2012)
[7] Boyd, S., Parikh, N., Chu, E., Peleato, B. and Eckstein, J .: Distributed Optimization and
Statistical Learning via the Alternating Direction Method of Multipliers. Foundations and
Trends in Machine Learning, Michael Jordan, Editor in Chief , 3, 1–122 (2011)
[8] Bregman, L. M.: The relaxation method of ﬁnding the commo n point of convex sets and
its application to the solution of problems in convex progra mming. USSR Comp. Math.
and Math. Physics 7, 200–217 (1967)
[9] Chandrasekaran, V., Parrilo, P., and Willsky, A.: Laten t variable graphical model via
convex optimization. Preprint. (2010)
[10] Cottle, R. W. Duvall, S. G., and Zikan, K.: A Lagrangian r elaxation algorithm for the
constrained matrix problem. Naval Res. Logistics Quarterl y. 33, 55–76 (1986)
[11] De Pierro, A. R. and Iusem, A. N.: On the convergence prop erties of Hildreth’s quadratic
programming algorithm. Math. Prog. 47, 37–51 (1990)
[12] Deng, W. and Yin, W.: On the global linear convergence of alternating direction methods.
Preprint. (2012)
[13] Douglas, J. and Rachford, H.H.: On the numerical soluti on of the heat conduction problem
in 2 and 3 space variables. Transactions of the American Math ematical Society. 82, 421-439
(1956)
32

<!-- page 34 -->
[14] Eckstein, J.: Splitting methods for monotone operator s with applications to parallel op-
timization. Ph.D Thesis, Operations Research Center, MIT. (1989)
[15] Eckstein, J. and Bertsekas, D.P.: On the Douglas-Rachf ord splitting method and the prox-
imal point algorithm for maximal monotone operators. Mathe matical Programming.55,
293-318 (1992)
[16] Eckstein, J. and Svaiter, B.F.: General projective spl itting methods for sums of maximal
monotone operators. SIAM Journal on Control and Optimizati on. 48, 787–811 (2010)
[17] Gabay, D.: Application of the method of multipliers to v aruational inequalities, In:
Fortin, M., Glowinski, R., eds., Augmented Lagrangian meth ods: Application to the
numerical solution of Boundary-Value Problem, North-Holl and, Amsterdam, The Nether-
lands. 299-331 (1983)
[18] Gabay, D. and Mercier, B.: A dual algorithm for the solut ion of nonlinear variational
problems via ﬁnite-element approximations. Computer and M athematics with Applica-
tions. 2, 17-40 (1976)
[19] Glowinski, R.: Numerical methods for nonlinear variat ional problems. Springer-Verlag,
New York, Berlin, Heidelberg, Tokyo (1984)
[20] Glowinski, R. and Le Tallec, P.: Augmented Lagrangian a nd operator splitting methods
in nonlinear mechanics. SIAM Studies in Applied Mathematic s, Philadelphia, PA. (1989)
[21] Goldfarb, D. and Ma, S.: Fast multiple splitting algori thms for convex optimization. Arxiv
preprint arXiv:0912.4570 (2009)
[22] Goldfarb, D., Ma, S. and Scheinberg, K.: Fast alternati ng linearization methods for min-
imizing the sum of two convex functions. Arxiv preprint arXi v:0912.4571 (2009)
[23] Goldstein, A. A.: Convex programming in hilbert space. Bull. Am. Math. Soc. 70, 709–710
(1964)
[24] Goldstein, T., O’Donoghue, B. and Setzer, S.: Fast alte rnating direction optimization
methods. CAM report 12-35, UCLA, (2012).
[25] He, B. S., Tao, M. and Yuan, X.M.: Alternating direction method with gaussian back
substitution for separable convex programming. SIAM J. Opt im. 22, 313-340 (2012)
[26] He, B. S. and Yuan, X. M.: On the O(1/n ) convergence rate of the Douglas-Rachford
alternating direction method. SIAM J. Numer. Anal. 50, 700- 709 (2012)
33

<!-- page 35 -->
[27] Herman, G. T.: Image Reconstruction from Projection: T he Fundamentals of Computer-
ized Tomography. Academic Press, New York. (1980)
[28] Hoﬀman, A. J.: On Approximate Solutions of Systems of Lin ear Inequalities. J. Res. Nat.
Bur. Standards. 49, 263-265 (1952)
[29] Iusem, A. N.: On dual convergence and the rate of primal c onvergence of bregman’s
convex programming method. SIAM J. Control Optim. 1, 401–42 3 (1991)
[30] Kontogiorgis, S. and Meyer, R. R.: A variable-penalty a lternating directions method for
convex optimization, Mathematical Programming. 83, 29-53 (1998)
[31] Lamond, B., and Stewart, N. F.: Bregman’s balancing met hod. Transportation Res. 15B,
239–248 (1981)
[32] Levitin, E. S. and Poljak, B. T.: Constrained minimizat ion methods. Z. Vycisl. Mat. i
Mat. Fiz. 6, 787-823 (1965). English translation in USSR Com put. Math. Phys. 6, 1–50
(1965)
[33] Lions, P. L., and Mercier., B.: Splitting algorithms fo r the sum of two nonlinear operators.
SIAM Journal on Numerical Analysis, 16, 964–979, (1979)
[34] Lin, Y. Y., and Pang, J.-S.: Iterative methods for large convex quadratic programs: a
survey. SIAM J. Control Optim. 18, 383-411 (1987)
[35] Luo, Z.-Q., and Tseng, P.: On the convergence of the coor dinate descent method for
convex diﬀerentiable minimization. J. Optim. Theory and App l. 72, 7–35 (1992)
[36] Luo, Z.-Q., and Tseng, P.: On the Linear convergence of d escent methods for convex
essentially smooth minimization. SIAM J. Control Optim. 30 , 408-425 (1992)
[37] Luo, Z.-Q. and Tseng, P.: On the convergence rate of dual ascent methods for strictly
convex minimization. Mathematics of Operations Research 1 8, 846–867 (1993)
[38] Ma, S.: Alternating proximal gradient method for conve x minimization. Preprint (2012)
[39] Mangasarian, O. L., and Shiau, T.-H.: Lipschitz Contin uity of solutions of linear inequal-
ities, programs and complementarity problems. SIAM J. Cont rol Optim. 25 583–595
(1987)
[40] Ohuchi, A., and Kaji, I.: Lagrangian dual coordinatewi se maximization algorithm for
network transportation problems with quadratic costs. Net works 14, 515–530 (1984)
[41] Ortega, J. M., and Rheinboldt, W. C.: Iterative Solutio n of Nonlinear Equations in
Several Variables. Academic Press, New York, New York. (197 0)
34

<!-- page 36 -->
[42] Pang, J.-S.: On the convergence of dual ascent methods f or large scale linearly constrained
optimization problems. The University of Texas, School of M anagement, Dallas, Texas.
(1984)
[43] Pang, J.-S.: A posteriori error bounds for the linearly –constrained variational inequality
problem. Math. Oper. Res. 12, 474–484 (1987)
[44] Rockafellar, R. T.: Convex Analysis, Princeton Univer sity Press, Princeton, New Jersey.
(1970)
[45] Schneider, M. H., and Zenios, S. A.: A comparative study of algorithms for matrix bal-
ancing. Oper. Res. 38, 439–455 (1990)
[46] Stephenson, D.: Pipeﬂow Analysis, Elsevier, Amsterda m. (1984)
[47] Tao, M. and Yuan, X.M.: Recovering low-rank and sparse c omponents of matrices from
incomplete and noisy observations. SIAM J. Optim. 21, 57-81 (2011)
[48] Tibshirani, R.: Regression shrinkage and selection vi a the Lasso. Jourbal of the Royal
Statistical Society Ser. B. 58, 267-288 (1996)
[49] Tseng, P.: Dual ascent methods for problems with strict ly convex costs and linear con-
straints: a uniﬁed approach. SIAM J. Control Optim. 28, 214- 242 (1990)
[50] Tseng, P.: Approximation accuracy, gradient methods, and error bound for structured
convex optimization. Technical report. (2009)
[51] Tseng, P., and Bertsekas D. P.: Relaxation methods for p roblems with strictly convex
separable costs and linear constraints. Math. Prog. 38, 303 –321 (1987)
[52] Tseng, P., and Bertsekas D. P.: Relaxation methods for p roblems with strictly convex
costs and linear constraints. Math. Oper. Res. 16, 462–481 ( 1991)
[53] Ventura, J. A., and Hearn, D. W.: Computational develop ment of a Lagrangian dual ap-
proach for quadratic networks. University of Florida, Indu strial and Systems Engineering
Department Report 88–14 (1988)
[54] Wang, X. F., and Yuan, X. M.: The linearized alternating direction method of multipliers
for dantzig selector. SIAM Journal on Scientiﬁc Computing, 34, 2792–2811, (2012).
[55] Yang, J. F., and Zhang, Y.: Alternating direction algor ithms for l1-problems in compres-
sive sensing. SIAM Journal on Scientiﬁc Computing, 33, 250– 278 (2011).
35

<!-- page 37 -->
[56] Yuan, M. and Lin, Y.: Model selection and estimation in r egression with grouped vari-
ables. Journal of the Royal Statistical Society: Series B (S tatistical Methodology). 68,
49-67 (2006)
[57] Zhang, H., Jiang, J.J. and Luo, Z.-Q.: On the linear conv ergence of a proximal gradient
method for a class of nonsmooth convex minimization problem s. Preprint. (2012)
[58] Zenios, S. A., and Mulvey, J. M.: Relaxation techniques for strictly convex network
problems. Ann. Oper. Res. 5, 517–538 (1986)
[59] Zhou, Z., Li, X., Wright, J., Candes, E.J., and Ma, Y.: St able principal component pur-
suit. Proceedings of 2010 IEEE International Symposium on I nformation Theory (2010)
36