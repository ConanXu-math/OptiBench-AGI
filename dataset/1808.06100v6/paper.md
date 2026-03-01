# On the solution existence and stability of polynomial optimization problems

**arXiv ID:** 1808.06100v6

**Authors:** Vu Trung Hieu

**Abstract:** This paper introduces and investigates a regularity condition in the asymptotic sense for optimization problems whose objective functions are polynomial. Under this regularity condition, the normalization argument in asymptotic analysis enables us to see the solution existence as well as the solution stability of these problems. We prove a Frank-Wolfe type theorem for regular optimization problems and an Eaves type theorem for non-regular pseudoconvex optimization problems. Moreover, we show results on the stability such as upper semicontinuity and local upper-Hölder stability of the solution map of polynomial optimization problems. At the end of the paper, we discuss the genericity of the regularity condition.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:1808.06100v6  [math.OC]  9 Aug 2021
Noname manuscript No.
(will be inserted by the editor)
On the solution existence and stability of polynomial
optimization problems
V u T rung Hieu
Received: date / Accepted: date
Abstract In this paper, we introduce and investigate a new regularity condi-
tion in the asymptotic sense for optimization problems whose object ive func-
tions are polynomial. The normalization argument in asymptotic analys is en-
ables us to study the existence as well as the stability of solutions of these
problems. We prove a Frank-Wolfe type theorem for regular optimization prob-
lems and an Eaves type theorem for non-regular pseudoconvex op timization
problems. Moreover, under the regularity condition, we show resu lts on the
stability such as upper semicontinuity and local upper-H¨ older stab ility of the
solution map of polynomial optimization problems. At the end of the pa per,
we discuss the genericity of the regularity condition.
Keywords Polynomial optimization · Regularity condition · Asymptotic
cone ·Frank–Wolfe type theorem ·Eaves type theorem ·Upper semicontinuity ·
Local upper-H¨ older stability· Genericity
Mathematics Subject Classiﬁcation (2000) 90C30 · 14P10
1 Introduction
We consider the following optimization problem
minimize f (x) subject to x ∈ K,
where K is a nonempty, closed subset of Rn and f : Rn → R is a polynomial
in n variables of degree d ≥ 2. The problem and its solution set are denoted
by OP( K, f ) and Sol( K, f ) respectively. Let fd be the homogeneous compo-
nent of degree d of f , and let K∞ be the asymptotic cone of K that will be
Vu Trung Hieu
Sorbonne Universit´ e,CNRS, LIP6, F-75005, Paris, France
Division of Mathematics, Phuong Dong University, 171 Trung Kinh Street, Hanoi, Vietnam
E-mail: trung-hieu.vu@lip6.fr

<!-- page 2 -->
2 Vu Trung Hieu
introduced in Section 2. We say that OP( K, f ) is regular if the solution set
of the asymptotic problem OP( K∞ , f d) is bounded, and the problem is non-
regular otherwise. The regularity condition has appeared in studies about the
solution existence and stability in quadratic programming (see, e.g., [1 , 2] and
the references therein).
Asymptotic cones and functions play an important role in optimization
and variational inequalities [3]. The normalization argument in asympto tic
analysis enables us to study the existence and stability of solutions n ot only
for quadratic programming, linear complementarity problems, and a ﬃne vari-
ational inequalities (see, e.g., [1, 4]), but also for polynomial compleme ntarity
problems and polynomial variational inequalities that have unbounde d con-
straint sets (see, e.g., [5, 6]). In this paper, the normalization arg ument is used
as the main technique to investigate the existence as well as the sta bility of
solutions to polynomial optimization problems.
In 1956, Frank and Wolfe [7] proved that if K is polyhedral and f is
quadratic and bounded from below on K, then Sol( K, f ) is nonempty. Sev-
eral versions of the Frank-Wolfe theorem for quadratic, cubic, a nd polynomial
optimization problems have been shown in [1, 2, 8, 9, 10, 11, 12]. Belou sov and
Klatte [9], and Obuchowska [10] have proved Frank-Wolfe type the orems for
convex and quasiconvex polynomial optimization problems. Recently , by us-
ing a technique from semi-algebraic geometry, Dinh, Ha and Pham [11] have
shown a Frank-Wolfe type theorem for nondegenerate problems. The present
paper gives another Frank-Wolfe type theorem, which says that if OP(K, f )
is regular and f is bounded from below on K, then the problem has a solu-
tion. Besides, the Eaves theorem [13] provides us with another cr iterion for
the existence of solutions to quadratic optimization problems. Exte nsions of
this theorem for quadratically constrained quadratic problems hav e been in-
vestigated in [1, 2, 14, 15]. This paper introduces an Eaves type the orem for
non-regular pseudoconvex optimization problems, where the cons traint sets
are convex.
Under the assumption that the constraint set K is compact and semi-
algebraic, some stability and genericity results for polynomial optimiz ation
problems have been shown by Lee and Pham [16]. If K is compact, then its
asymptotic cone is trivial, i.e., K∞ = {0}; Hence that OP( K, f ) satisﬁes the
regularity condition obviously. In the present paper, K may be unbounded.
Under the regularity condition, we prove several local properties of the solution
map of polynomial optimization problems such as local boundedness a nd upper
semicontinuity. Furthermore, based on an error bound for a polyn omial system
in [20], we prove the local upper-H¨ older stability of the solution map.
We denote by Rd[x] the space of all polynomials of degree at most d and
by Rd the set of all polynomials g of degree d such that OP( K, g ) is regular.
The set Rd is an open cone in Rd[x]. At the end of this work, K is deﬁned by
convex polynomials, we prove that Rd is generic in Rd[x].
The organization of the paper is as follows. Section 2 gives a brief intr oduc-
tion to asymptotic cones, polynomials, and the regularity condition. Section 3
proves two criteria of the solution existence. Section 4 investigate s properties

<!-- page 3 -->
On the solution existence and stability of polynomial optim ization problems 3
of the solution map. The last section discusses the genericity of the regularity
condition.
2 Preliminaries
Recall that the asymptotic cone [3] of a nonempty closed subset S in Rn is
deﬁned and denoted by
S∞ =
{
v ∈ Rn : ∃tk → +∞ , ∃xk ∈ S with lim
k→∞
xk
tk
= v
}
.
Clearly, the cone S∞ is closed and contains 0. The set S is bounded if and
only if S∞ is trivial. Furthermore, if S is convex then S∞ is a closed convex
cone and S∞ = 0 +S, where 0 +S is the recession cone of S, that consists of all
vectors v ∈ Rn such that x + tv ∈ S for any x ∈ S and t ≥ 0. Thus, one has
S = S + S∞ when S is convex.
Let d ≥ 2 be given. The dimension of the space Rd[x] is ﬁnite; its dimension
is denoted by ρ. Let X(x) be the vector consisting of ρ monomials of degree
at most d which is listed by lexicographic ordering
X(x) := (1 , x 1, x 2, . . . , x n, x 2
1, x 1x2, . . . , x 1xn, . . . , x d
1, x d− 1
1 x2, . . . , x d
n)T .
For every g ∈ Rd[x], there exists a unique vector a = ( a1, . . . , a ρ ) ∈ Rρ such
that g(x) = aT X(x). We denote by ∥g∥ the ℓ2–norm of the polynomial g,
namely
∥g∥ := ∥a∥ =
√
a2
1 + · · ·+ a2
ρ .
The Cauchy–Schwarz inequality yields |g(x)| ≤ ∥ X(x)∥∥g∥. Furthermore, if
{gk} is a convergent sequence in Rd[x] with gk → g, then gk
d → gd.
Throughout the paper, we assume that the constraint set K ⊂ Rn is
nonempty and closed, and the objective function f : Rn → R is a polyno-
mial of degree d ≥ 2.
We say that OP(K, f ) is a polynomial optimization problem if K is given by
polynomials. With the given set K and the given integer d, the solution map
of polynomial optimization problems OP( K, g ), where g ∈ Rd[x], is deﬁned by
SolK(·) : Rd[x] ⇒ Rn, g ↦→ Sol(K, g ).
Assume that g ∈ Rd[x] with deg g = d and g = gd + · · ·+ g1 + g0, where
gl is a homogeneous polynomial of degree l, i.e., gl(tx) = tlgl(x) for all t ≥ 0
and x ∈ Rn, l ∈ [d] := {1, . . . , d }, and g0 ∈ R. Then, gd is the leading term (or
the recession polynomial) of the polynomial g (of degree d). Clearly, one has
gd(x) = lim
λ → +∞
g(λx)
λ d , ∀x ∈ Rn .
For the pair (K, f ), the asymptotic pair ( K∞ , f d) is unique. The asymptotic
optimization problem OP( K∞ , f d) plays a vital role in the investigation of
behavior of OP( K, f ) at inﬁnity. The following remarks point out (without
proof) the basic properties of the asymptotic problem.

<!-- page 4 -->
4 Vu Trung Hieu
Remark 2.1 Since fd is a homogeneous polynomial and K∞ is a closed cone,
the asymptotic optimization problem OP( K∞ , f d) has a solution if and only
if fd is non-negative on K∞ .
Remark 2.2 Assume that Sol( K∞ , f d) is nonempty. Then, this set is a closed
cone with 0 ∈ Sol(K∞ , f d). In addition, Sol( K∞ , f d) coincides with the zero
set of fd in K∞ , i.e.,
Sol(K∞ , f d) = {x ∈ K∞ : fd(x) = 0 }.
Now, we introduce the regularity notion concerning the boundedne ss of the
solution set of OP( K∞ , f d).
Deﬁnition 2.1 The problem OP( K, f ) is said to be regular if Sol( K∞ , f d) is
bounded and non-regular otherwise.
Denote by E d (resp., Od, U d) the set of all polynomials g of degree d such
that Sol(K∞ , f d) is the empty set (resp., the trivial cone, an unbounded cone).
Clearly, Rd = E d ∪ O d, and one has the following disjoint union:
Rd[x] = Rd− 1[x] ∪ E d ∪ O d ∪ U d . (2.1)
Remark 2.3 The boundedness of Sol( K∞ , f d) implies that of Sol( K, f ). In-
deed, assume to the contrary that Sol( K, f ) is unbounded. There exists an
unbounded sequence {xk} ⊂ Sol(K, f ). Without loss of generality, we can as-
sume xk is nonzero for all k, ∥xk∥ → +∞ , and ∥xk∥− 1xk → ¯x for some ¯x ∈ Rn
with ∥¯x∥ = 1. Note that f (xk) = f ∗ , where f ∗ ∈ R is the minimum of f over
K, for all k. By dividing the last equation by ∥xk∥d and letting k → +∞ ,
we obtain fd(¯x) = 0. It follows that ¯ x ∈ Sol(K∞ , f d). Since ¯x ̸= 0, the cone
Sol(K∞ , f d) is unbounded, which contradicts our assumption. Thus, the claim
is proved.
Remark 2.4 We observe that the set Rd is nonempty. If K is bounded then
Rd coincides with the set of all g ∈ Rd[x] such that deg g = d. Hence, we
can suppose that K is unbounded. Clearly, the cone K∞ also is unbounded.
Let ¯x ∈ K∞ be nonzero. There exists l ∈ [n] such that ¯xl ̸= 0. Let us deﬁne
a homogeneous polynomial of degree d as f (x) := − (¯xlxl)d. For any t > 0,
one has t¯x ∈ K∞ , and f (t¯x) = − (¯x2
l )dtd → −∞ as t → +∞ . Then, f is not
bounded from below on K∞ . This yields Sol( K∞ , f ) = ∅ and f ∈ R d.
Example 2.1 Consider the case that n = 1, K = R and d = 2. One has R2[x] =
{a2x2+a1x+a0 : (a2, a 1, a 0) ∈ R3}. Since K∞ = R, an easy computation shows
that E 2 = {a2x2 + a1x + a0 : a2 < 0, a 1 ∈ R, a 0 ∈ R}, O2 = {a2x2 + a1x + a0 :
a2 > 0, a 1 ∈ R, a 0 ∈ R}, and R2 = {a2x2 + a1x + a0 : a2 ̸= 0, a 1 ∈ R, a 0 ∈ R}.
3 Two criteria for the solution existence
We now introduce two criteria for the solution existence of OP( K, f ). In the
proofs, the normalization argument in asymptotic analysis plays a vit al role;
meanwhile, the semi-algebraicity of K is not required.

<!-- page 5 -->
On the solution existence and stability of polynomial optim ization problems 5
3.1 A Frank-Wolfe type theorem for regular problems
The following theorem provides us a criterion for the solution existen ce of
regular optimization problems.
Theorem 3.1 (F rank-Wolfe type theorem) If OP(K, f ) is regular and f
is bounded from below on K, then its solution set is nonempty and compact.
Proof Suppose that f ∈ R d, i.e. f ∈ E d ∪ O d, and there exists γ ∈ R such
that γ ≤ f (x), for all x ∈ K. For any given v ∈ K∞ , there are two sequences
{tk} ⊂ R+ and {xk} ⊂ K such that tk → +∞ and t− 1
k xk → v as k → +∞ .
For any k, one has γ ≤ f (xk). Dividing both sides of the last inequality by td
k
and letting k → +∞ , we obtain 0 ≤ fd(v). Thus, fd is non-negative over K∞ .
It follows from Remark 2.1 that f does not belong to E d; hence, we conclude
that f must be in Od.
Let ¯x ∈ K be given, and M := {x ∈ K : f (x) ≤ f (¯x)}. It is easy to check
that Sol( M, f ) = Sol( K, f ). Hence, we need only to prove that Sol( M, f ) is
nonempty and compact.
Clearly, M is closed. We claim that M is bounded. On the contrary, we
suppose that there exists an unbounded sequence {xk} ⊂ M such that xk is
nonzero for all k, ∥xk∥ → +∞ , and ∥xk∥− 1xk → v for some v ∈ Rn with
∥v∥ = 1. One has
γ ≤ f (xk) ≤ f (¯x), (3.1)
for all k. Dividing the values in (3.1) by ∥xk∥d and letting k → +∞ , we get
fd(v) = 0. This yields v ∈ Sol(K∞ , f d). Because of v ̸= 0, Sol( K∞ , f d) is
unbounded. This contradicts our assumption and, thus, the claim is proved.
The compactness of M and Bolzano-Weierstrass’ Theorem allow us to con-
clude that Sol( M, f ) is nonempty and compact. ⊓ ⊔
Remark 3.1 From the proof of Theorem 3.1, we see that if f is bounded from
below on K then Sol( K∞ , f d) is nonempty, i.e. f ∈ O d ∪ U d. Hence, if f ∈ E d
then OP( K, f ) has no solution.
Corollary 3.1 Assume that f = α 1xd
1 +· · ·+α nxd
n +p where d is even, α ℓ > 0
for all ℓ ∈ [n], p is a polynomial with deg p < d . Then, Sol(K, f ) is nonempty
and compact.
Proof Clearly, fd = α 1xd
1 + · · ·+ α nxd
n is non-negative over Rn. It follows that
fd is also non-negative over K∞ . From Remarks 2.1 and 2.2, it is clear that
Sol(K∞ , f d) is nonempty and
Sol(K∞ , f d) = {x ∈ K∞ : α 1xd
1 + · · ·+ α nxd
n = 0} = {0}.
This means that f ∈ O d. Clearly, f is bounded from below on K, and the
condition of Theorem 3.1 holds. Therefore, Sol( K, f ) is nonempty and com-
pact. ⊓ ⊔
The following example illustrates Theorem 3.1, in which the constraint s et
is neither convex nor semi-algebraic.

<!-- page 6 -->
6 Vu Trung Hieu
Example 3.1 Consider the optimization problem OP( K, f ), where the polyno-
mial f is given by f (x1, x 2) = x3
2 − x1x2 and the constraint set K is given
by
K = {(x1, x 2) ∈ R2 : x1 ≥ 0, x 2 − x1 ≥ 0, e x1 − x2 ≥ 0}.
Since f3(x1, x 2) = x3
2 and K∞ = {(x1, x 2) ∈ R2 : x1 ≥ 0, x 2 − x1 ≥ 0}, one
has Sol( K∞ , f 3) = {(0, 0)}. According to Theorem 3.1, Sol( K, f ) is nonempty
and compact.
3.2 An Eaves type theorem for non-regular problems
In this subsection, we investigate the solution existence of non-re gular opti-
mization problems, where the objective functions are pseudoconv ex on the
constraint sets.
Assume that U is an open subset of Rn. One says the polynomial f is
pseudoconvex on U if, for any x, y ∈ U such that ⟨∇ f (x), y − x⟩ ≥ 0, here ∇ f
is the gradient of f , we have f (y) ≥ f (x). Recall that f is pseudoconvex on
U if and only if ∇ f is pseudomonotone on U [17, Theorem 3.1], i.e. if, for any
x, y ∈ U such that ⟨∇ f (x), y − x⟩ ≥ 0, we have ⟨∇ f (y), y − x⟩ ≥ 0.
Lemma 3.1 Assume that K is convex and f is pseudoconvex on an open set
U containing K. If x0 ∈ Sol(K, f ), then ⟨∇ f (x), x − x0⟩ ≥ 0 for all x ∈ K.
Proof Since f is pseudoconvex on the set U , the gradient ∇ f is pseudomono-
tone on U . Suppose that x0 ∈ Sol(K, f ), one has ⟨∇ f (x0), x − x0⟩ ≥ 0 for all
x ∈ K (see, e.g., [18, Proposition 5.2]). The pseudomonotonicity of the gra di-
ent implies that ⟨∇ f (x), x − x0⟩ ≥ 0 for all x ∈ K. The lemma is proved. ⊓ ⊔
Theorem 3.2 (Eaves type theorem) Assume that K is convex and f is
pseudoconvex on an open set containing K. If OP(K, f ) is non-regular, then
the following statements are equivalent:
(a) If v ∈ Sol(K∞ , f d) \ {0}, then there exists x ∈ K such that ⟨∇ f (x), v ⟩ > 0;
(b) Sol( K, f ) is nonempty and compact.
Proof Suppose that OP( K, f ) is non-regular. We prove (a) ⇒ (b). Assume
that (a) holds. For each k ∈ N, we denote
Kk = {x ∈ Rn : x ∈ K, ∥x∥ ≤ k}.
Clearly, Kk is compact and convex. Without loss of generality, we can assume
that Kk is nonempty. According to Bolzano-Weierstrass’ Theorem, OP( Kk, f )
has a solution, denoted by xk.
We assert that the sequence {xk} is bounded. Indeed, suppose on the con-
trary that {xk} is unbounded, xk ̸= 0 for all k, ∥xk∥ → +∞ , and ∥xk∥− 1xk →
v, where v ∈ K∞ and ∥v∥ = 1. For each k, one has
f (xk) ≤ f (x), ∀x ∈ Kk. (3.2)

<!-- page 7 -->
On the solution existence and stability of polynomial optim ization problems 7
Let y ∈ K be given. For k large enough, y ∈ Kk and f (xk) ≤ f (y). By
dividing two sides of the last inequality by ∥xk∥d and letting k → +∞ , we
obtain fd(v) ≤ 0. This leads to v ∈ Sol(K∞ , f d) \ {0}. Furthermore, since f is
pseudoconvex on Kk, from Lemma 3.1 we have
⟨∇ f (y), y − xk⟩ ≥ 0. (3.3)
Dividing both sides of the inequality in (3.3) by ∥xk∥ and letting k → +∞ , we
obtain ⟨∇ f (y), v ⟩ ≤ 0. The conclusion holds for any x ∈ K, i.e., ⟨∇ f (x), v ⟩ ≤ 0
for all x ∈ K. This contradicts to our assumption. Hence, {xk} is bounded.
We can assume that xk → ¯x. From (3.2), by the continuity of f , it easy to
check that ¯x solves OP( K, f ), so Sol( K, f ) is nonempty.
To prove the compactness of the solution set, we can repeat the p revious
argument by supposing that there is an unbounded solution sequen ce {xk},
and can show that there exists v ∈ Sol(K∞ , f d) \ {0} such that ⟨∇ f (x), v ⟩ ≤ 0
for all x ∈ K. This contradicts to (b).
(b) ⇒ (a) Since K is convex, one has K∞ = 0 +K and K = K + K∞ .
Suppose that Sol( K, f ) is nonempty and compact, but (b) is wrong, i.e. there
exists v ∈ Sol(K∞ , f d) \ {0} such that ⟨∇ f (x), v ⟩ ≤ 0 for all x ∈ K. Let ¯x be a
solution of OP(K, f ). For any t ≥ 0, one has ¯x+tv ∈ K and ⟨∇ f (¯x+tv), v ⟩ ≤ 0.
Thus, we have
⟨∇ f (¯x + tv), ¯x − (¯x + tv)⟩ ≥ 0.
The pseudoconvexity of f yields f (¯x) ≥ f (x0 + tv). Hence, x0 + tv belongs
to Sol( K, f ), for any t ≥ 0. This shows that Sol( K, f ) is unbounded which
contradicts to our assumption. Thus (a) holds, and the proof is co mplete. ⊓ ⊔
Corollary 3.2 Assume that K is convex and f is convex on an open convex
set containing K. If OP(K, f ) is non-regular, then the following statements
are equivalent:
(a) If v ∈ Sol(K∞ , f d) \ {0}, then there exists x ∈ K such that ⟨∇ f (x), v ⟩ > 0;
(b) Sol( K, f ) is nonempty and compact.
Proof Since the convexity implies the pseudoconvexity, by applying Theore m
3.2 for the convex polynomial f , we have the assertion. ⊓ ⊔
The following example illustrates Corollary 3.2.
Example 3.2 Consider the polynomial optimization problem OP( K, f ) with
the objective function f (x1, x 2) = 1
6 x3
2 + 1
2 x2
1 − x1x2 and the constraint set
K = {(x1, x 2) ∈ R2 : x1x2 ≥ 1, x 2 ≥ 2}. The gradient and the Hessian matrix
of f , respectively, are given by
∇ f =
[ x1 − x2
− x1 + 1
2 x2
2
]
, H =
[ 1 − 1
− 1 x2
]
.
It is easy to check that K is convex and H is positive semideﬁnite on the open
set U = {(x1, x 2) ∈ R2 : x1x2 > 0, x 2 > 1} ⊃ K; hence f is convex on K. One
has K∞ = R2
+ and f3(x1, x 2) = 1
6 x3
2. This yields
Sol(K∞ , f 3) = {(x1, x 2) ∈ R2 : x1 ≥ 0, x 2 = 0}.

<!-- page 8 -->
8 Vu Trung Hieu
For every v = (α, 0) in Sol( K∞ , f 3)\ {0}, one has α > 0. By choosing the point
x = (3, 2) in the constraint set, we get ⟨∇ f (x), v ⟩ = α > 0. Finally, according
to Corollary 3.2, the solution set of OP( K, f ) is nonempty and compact.
4 Stability of the solution map
We investigate the local boundedness, the upper semicontinuity, a nd the local
upper-H¨ older stability of the solution map under the regularity con dition.
4.1 Upper semicontinuity of the solution map
To prove the local boundedness and the upper semicontinuity of th e solution
map, we need the following lemma.
Lemma 4.1 The set Rd is open in Rd[x].
Proof To prove the openness of Rd, we only need to show that the complement
Rd[x]\ Rd is closed. Clearly, Rd[x]\ Rd = Rd− 1[x]∪U d. Let {gk} be a sequence
in Rd[x] \ Rd such that gk → g. From the deﬁnition of Rd, if deg g < d , i.e.,
g ∈ Rd− 1[x], then g ∈ Rd[x] \ Rd. Thus, we can suppose that deg g = d. One
has gk
d → gd, here gk
d is the component of degree d of gk.
We now prove that g belongs to U d. For each k, Sol(K∞ , g k
d ) is unbounded.
There exists an unbounded sequence {xk} such that xk ∈ Sol(K∞ , g k
d ), ∥xk∥ →
+∞ , ∥xk∥− 1xk → ¯x with ∥¯x∥ = 1. Let v ∈ K∞ be given. One has ∥xk∥v ∈ K∞
and gk
d (∥xk∥v) ≥ gk
d (xk), for any k. Dividing the last inequality by ∥xk∥d and
letting k → +∞ , one has gk
d (v) ≥ gk
d (¯x). This conclusion holds for every
v ∈ K∞ . This yields ¯x ∈ Sol(K∞ , g d). As ∥¯x∥ = 1, we have ¯x ̸= 0. It follows
that g belongs to U d. The closedness of Rd[x] \ Rd is proved. ⊓ ⊔
Recall that a set-valued map Ψ : Rm ⇒ Rn is locally bounded at ¯u if
there exists an open neighborhood U of ¯u such that ∪ u∈ U Ψ (u) is bounded [19,
Deﬁnition 5.14]. The map Ψ is upper semicontinuous at ¯u ∈ T if for any open
set V ⊂ Rn such that Ψ (¯u) ⊂ V there exists a neighborhood U of ¯u such that
Ψ (u) ⊂ V for all u ∈ U . Recall that if Ψ is closed, namely, the graph
gph(Ψ ) :=
{
(u, v ) ∈ Rm × Rn : v ∈ Ψ (u)
}
is closed in Rm × Rn, and locally bounded at u, then Ψ is upper semicontinuous
at u [19, Theorem 5.19].
Proposition 4.1 Assume that K is convex. If OP(K, f ) is regular, then the
following statements hold:
(a) The solution map SolK(·) is locally bounded at f , i.e., there exists ε > 0
such that the set
Oε :=
⋃
g∈B (ε,d )
Sol(K, f + g), (4.1)
where B(ε, d ) is the open ball in Rd[x] with center 0 and radius ε, is
bounded.

<!-- page 9 -->
On the solution existence and stability of polynomial optim ization problems 9
(b) The solution map SolK(·) is upper semicontinuous at f .
Proof (a) According to Lemma 4.1, Rd is open in Rd[x]. There is a closed ball
B(ε, d ) such that
f + B(ε, d ) ⊂ R d . (4.2)
Assume to the contrary that Oε is unbounded. Then, there exists an un-
bounded sequence {xk} and a sequence {gk} ⊂ B (ε, d ) such that xk solves
OP(K, f + gk) with xk ̸= 0 for every k, ∥xk∥ → +∞ , and ∥xk∥− 1xk → ¯x with
∥¯x∥ = 1. By the compactness of B(ε, d ), without loss of generality, we can
assume that gk → g with g ∈ B(ε, d ).
From assumptions, for every k, one has
(f + gk)(y) ≥ (f + gk)(xk), (4.3)
for any y ∈ K. Let y ∈ K be ﬁxed and assume that v ∈ K∞ . By the convexity
of K, one has y + ∥xk∥v ∈ K for any k. From (4.3), we conclude that
(f + gk)(y + ∥xk∥v) ≥ (f + gk)(xk).
Dividing this inequality by ∥xk∥d and taking k → +∞ , we obtain
(f + g)d(v) ≥ (f + g)d(¯x).
The conclusion hold for any v ∈ K∞ . It follows that ¯ x ∈ Sol (K∞ , (f + g)d).
Because of (4.2), Sol ( K∞ , (f + g)d) is contained in {0}, which contradicts to
∥¯x∥ = 1. Hence, Oε must be bounded.
(b) It is not diﬃcult to prove that the graph
gph(Sol) :=
{
(g, x ) ∈ Rd[x] × Rn : x ∈ Sol(K, g )
}
is closed in Rd[x] × Rn. Since Sol K(·) is locally bounded on Rd, according to
[19, Theorem 5.19], Sol K(·) is upper semicontinuous at f . ⊓ ⊔
4.2 Local upper-H¨ older stability of the solution map
When the constraint set K is convex and given by polynomials, we can inves-
tigate the local upper-H¨ older stability of the solution map under th e regular
condition. To prove the stability, we need the following lemma.
Lemma 4.2 ([20]) Let U be a semi-algebraic subset in Rn, represented by
U = {x ∈ Rn : ui(x) = 0 , i ∈ [l], v j(x) ≤ 0, j ∈ [m]} ,
where ui(x), i ∈ [l], and vj(x), j ∈ [m], are polynomials. For any compact set
V ⊂ Rn, there are constants c > 0 and H > 0 such that
d(x, U ) ≤ c
( l∑
i=1
|ui(x)|+
m∑
j=1
[vj(x)]+
) H
,
for all x ∈ V , here [r]+ := max{0, r } and d(x, U ) the usual distance from x to
the set U .

<!-- page 10 -->
10 Vu Trung Hieu
Theorem 4.1 Assume that OP(K, f ) is regular and K is a convex set given
by
K = {x ∈ Rn : pi(x) = 0 , i ∈ [l], q j(x) ≤ 0, j ∈ [m]} ,
where all pi, q j are polynomials. If Sol(K, f ) is nonempty, then the map SolK(·)
is locally upper-H¨ older stable at f , i.e., there exist ℓ > 0, H > 0 and ε > 0
such that
Sol(K, g ) ⊂ Sol(K, f ) + ℓ∥g − f ∥H B, (4.4)
for all g ∈ Rd[x] satisfying ∥g − f ∥ < ε, where B is the closed unit ball in Rn.
Proof Suppose Sol(K, f ) is nonempty and its optimal value is f ∗ . Since OP(K, f )
is regular and K is convex, according to Proposition 4.1, there exists ε > 0
such that Sol( K, f ) ⊂ Oε, deﬁned by (4.1), is bounded. Let V be the closure
of Oε. It follows that V is a nonempty compact set. By the assumptions, we
see that
Sol(K, f ) = {x ∈ Rn : f (x) − f ∗ = 0, p i(x) = 0 , i ∈ [l], q j(x) ≤ 0, j ∈ [m]}.
From this equality, by applying Lemma 4.2 for U = Sol(K, f ) and the compact
set V , there are constants c0 > 0 and H > 0 such that
d(x, Sol(K, f )) ≤ c0A(x)H ∀x ∈ V, (4.5)
where
A(x) := |f (x) − f ∗ |+
l∑
i=1
|pi(x)|+
m∑
j=1
[qj (x)]+.
Let g ∈ Rd[x] be arbitrary given such that ∥g − f ∥ < ε. From the deﬁnition
of V , Sol( K, f ) and Sol( K, g ) are subsets of V . Here, Sol( K, g ) may be empty.
By the compactness of V , we deﬁne the constant L := max{∥X(x)∥ : x ∈ V }.
Hence, one has
|g(x) − f (x)| ≤ ∥ X(x)∥∥g − f ∥ ≤ L∥g − f ∥ ∀ x ∈ V. (4.6)
If Sol(K, g ) is empty, then (4.4) is obvious. Thus, we consider the case that
Sol(K, g ) ̸= ∅. Since both Sol( K, f ) and Sol( K, g ) are nonempty and compact,
for any xg ∈ Sol(K, g ), there is xf ∈ Sol(K, f ) such that
∥xg − xf ∥ = d(xg, Sol(K, f )). (4.7)
Because of pi(xg) = 0 for i ∈ [l] and qj(xg) ≤ 0 for j ∈ [m], from the deﬁnition
of A(x), one has A(xg) = |f (xg) − f ∗ |. By (4.7) and (4.5), we see that
∥xg − xf ∥ ≤ c0A(xg)H = c0|f (xg) − f ∗ |H .
Since xf ∈ Sol(K, f ), we have f (xf ) = f ∗ ≤ f (xg). Therefore, we obtain
∥xg − xf ∥ ≤ c0|f (xg) − f ∗|H = c0(f (xg) − f (xf ))H . (4.8)

<!-- page 11 -->
On the solution existence and stability of polynomial optim ization problems 11
It follows from xg ∈ Sol(K, g ) that g(xg) − g(xf ) ≤ 0. Since xg, x f ∈ V , we
conclude from (4.6) that
f (xg) − f (xf ) = ( f (xg) − g(xg)) + (g(xg) − g(xf )) + (g(xf ) − f (xf ))
≤ (f (xg) − g(xg)) + (g(xf ) − f (xf ))
≤ 2L∥g − f ∥.
The inequality (4.8) and the last result lead to
∥xg − xf ∥ ≤ c0(2L)H ∥g − f ∥H,
consequently,
d(xg, Sol(K, f )) = ∥xg − xf ∥ ≤ ℓ∥g − f ∥H,
where ℓ = c0(2L)H.
The conclusion holds for any xg in Sol( K, g ). Hence, the inclusion (4.4) of
the theorem is proved. ⊓ ⊔
5 Genericity of the regularity condition
In this section, we discuss the genericity of the regularity condition of polyno-
mial optimization problems.
A subset A is called generic in Rm if A contains a countable intersection
of dense and open sets in Rm. If A is generic in Rm and A ⊂ B then B also is
generic in Rm. Let T be a topological space. It is known that if h : Rm → T
is a homeomorphism and A is generic in Rm, then h(A) is generic in T .
Let U ⊂ Rm be a semi-algebraic set. Then, there exists a decomposition
of U into a disjoint union [21, Theorem 2.3.6], U = ∪ s
i=1Ui, where each Ui
is semi-algebraically homeomorphic to (0 , 1)di. Here, let (0 , 1)0 be a point,
(0, 1)k ⊂ Rk be the set of points x = ( x1, . . . , x k) such that xj ∈ (0, 1) for all
j ∈ [k]. The dimension of U is deﬁned by dim( U ) := max {d1, . . . , d s}. The
dimension is well-deﬁned and does not depends on the decomposition o f U .
Recall that if U is nonempty and dim( U ) is zero, then U has ﬁnitely many
points. Furthermore, if dim( Rm \U ) < m , then U is generic in Rm (see, e.g.
[22, Lemma 2.3]).
The space generated by all monomials of degree d listed by lexicographic
ordering {xd
1, x d− 1
1 x2, x d− 1
1 x3, . . . , x d
n} is denoted by Hd. One has the direct
sum Rd[x] = Hd ⊕ Rd− 1[x]. The dimension of Hd is denoted by η. For every
homogeneous polynomial h ∈ H d, one has a unique b ∈ Rη , such that h(x) =
bT Xd(x), where
X T
d (x) = ( xd
1, x d− 1
1 x2, x d− 1
1 x3, . . . , x d
n).
Here, ∇ (bT Xd(x)) is the gradient vector of bT Xd(x) and Db[∇ (bT Xd(x))] is
the Jacobian matrix of bT Xd(x) with respect to b.
Lemma 5.1 One has rank(Db[∇ (bT Xd(x))]) = n for all x ∈ Rn \{0}.

<!-- page 12 -->
12 Vu Trung Hieu
Proof In the proof, we are only interested in the monomials xd− 1
i xj, where
i, j ∈ [n]. Hence, for convenience, we rewrite X T
d (x) and bT respectively as
follows:
(xd
1, x d− 1
1 x2, . . . , x d− 1
1 xn; . . . ; xd− 1
n x1, x d− 1
n x2, . . . , x d
n; . . . )
and
(b11, b 12, . . . , b 1n; . . . ; bn1, b n2, . . . , b nn; . . . ).
Then, we have
bT Xd(x) =
∑
j∈ [n]
b1jxd− 1
1 xj + · · ·+
∑
j∈ [n]
bnjxd− 1
n xj + Q, (5.1)
where Q is a homogeneous polynomial of degree d.
From (5.1), an easy computation shows that
∂(bT Xd(x))
∂x i
= dbiixd− 1
i + (d − 1)
∑
j̸=i
bijxd− 2
i xj +
∑
j̸=i
bjixd− 1
j + ∂Q
∂x i
,
and the n × η-matrix Db[∇ (bT Xd(x))] can be described as follows
Db[∇ (bT Xd(x))] =
[
M1, M 2, · · ·, M n, · · ·
]
,
where the submatrix Mi, for i ∈ [n], is deﬁned by
Mi =



xd− 1
i Ii− 1 Oi× 1 Oi× (n− i)
L1× (i− 1) dxd− 1
i R1× (n− i)
Oi× (n− i) O(n− i)× 1 xd− 1
i In− i



with Ik being the unit k × k-matrix, Ok× s being the zero k × s-matrix,
L1× (i− 1) =
(
(d − 1)xd− 2
i x1, . . . , (d − 1)xd− 2
i xi− 1
)
,
and
R1× (n− i) =
(
(d − 1)xd− 2
i xi+1, . . . , (d − 1)xd− 2
i xn
)
.
We observe that det( Mi) = dxd(d− 1)
i , for all i ∈ [n]. Since x ̸= 0, there exists
l ∈ [n] such that xl ̸= 0. This implies that rank( Ml) = n. Hence, the rank of
Db[∇ (bT Xd(x))] is n, for any x ̸= 0. ⊓ ⊔
Suppose that C is a polyhedral cone given by
C = {x ∈ Rn : Ax ≥ 0} , (5.2)
where A = ( aij) ∈ Rp× n. Let KKT( C, g ), where g ∈ Rd[x], be the set of the
Karush-Kuhn-Tucker points of OP( C, g ), i.e., x ∈ KKT(C, g ) if and only if
there exists λ ∈ Rp such that
{ ∇ g(x) − AT λ = 0,
λ T (Ax) = 0 , λ ≥ 0, Ax ≥ 0. (5.3)

<!-- page 13 -->
On the solution existence and stability of polynomial optim ization problems 13
From the Karush-Kuhn-Tucker conditions, we see that Sol( C, g ) ⊂ KKT(C, g )
for all g ∈ Rd[x].
For each index set α ⊂ [p], we associate the pseudo-face Cα of C, which is
denoted and deﬁned by
Cα :=
{
x ∈ Rn :
n∑
j=1
aijxj = 0 ∀i ∈ α,
n∑
j=1
aijxj > 0 ∀i ∈ [p] \α
}
,
where aij is the element in the i-th row and the j-th column of A. The number
of pseudo-faces of C is ﬁnite. These pseudo-faces establish a disjoint decom-
position of C. So, we obtain
KKT(C, g ) =
⋃
α ⊂ [p]
(KKT(C, g ) ∩ Cα ) , (5.4)
The following proposition shows that the Karush-Kuhn-Tucker set -valued
map of homogeneous polynomial optimization problems
KKTC : Rη ⇒ Rn, b ↦→ KKTC (b) = KKT( C, b T Xd(x)),
is ﬁnite-valued, i.e., the cardinal # KKT C (b) is ﬁnite, on a generic semi-
algebraic set of Rη .
Proposition 5.1 Assume that C is a polyhedral cone given by (5.2) and the
matrix A is full rank. Then, there exists a generic semi-algebraic se t S ⊂ Rη
such that # KKTC (b) < ∞ for any b ∈ S.
Proof Let Cα be a nonempty pseudo-face of C and 0 /∈ Cα . This implies that
Xd(x) is nonzero on this pseudo-face. We consider the function
Φ α : Rη × Cα × R|α |
+ → Rn+|α |,
which is deﬁned by
Φ α (b, x, λ α ) =
(
∇ (bT Xd(x)) +
∑
i∈ α
λ iAi, A α x
)
,
where Aα x = ( Ai1 x, . . . , A i|α|x), i j ∈ α. Clearly, Cα is smooth and Φ α is a
semi-algebraic function of class C∞ . The Jacobian matrix of Φ α is determined
as follows
DΦ α =
[ Db[∇ (bT Xd(x))] ∗ AT
α
O|α |× η Aα O|α |×|α |
]
.
From Lemma 5.1, for all x ∈ Cα , the rank of Db[∇ (bT Xd(x))] is n. Since the
rank of Aα is |α |, we conclude that the rank of the matrix DΦ α is n + |α | for
all x ∈ Cα . Therefore, 0 ∈ Rn+|α |+|J| is a regular value of Φ α . According to
the Sard Theorem with parameter [22, Theorem 2.4], there exists a g eneric
semi-algebraic set Sα ⊂ Rη such that if b ∈ Sα then 0 is a regular value of the
map
Φ α,b : Cα × R|α | → Rn+|α |, Φ α,b (x, λ α ) = Φ α (b, x, λ α ).

<!-- page 14 -->
14 Vu Trung Hieu
We see that Ω (α, b ) := Φ − 1
α,b (0) is a semi-algebraic set. From the Regular Level
Set Theorem [23, Theorem 9.9], we can claim that if Ω (α, b ) is nonempty
then the set is a 0 − dimensional semi-algebraic set. It follows that Ω (α, b ) is a
ﬁnite set. Moreover, from (5.3), one has KKT C (b) ∩ Cα = π (Ω (α, b )), where
π is the projection Rn+|α | → Rn which is deﬁned by π (x, λ α ) = x. Hence,
KKTC (b) ∩ Cα is a ﬁnite set.
We consider the case that 0 ∈ Cα and deﬁne U := Cα \ {0}. As is clear, U
is semi-algebraic. From (5.3), we see that 0 ∈ KKTC (b). Hence,
KKTC (b) ∩ Cα = {0} ∪ (KKTC (b) ∩ U ).
From the previous argument, KKTC (b)∩ U is a ﬁnite set. By the decomposition
(5.4), KKT C (b) is a ﬁnite set.
Take S = ∩ α ⊂ [p] Sα , it follows that S is generic in Rη and KKT C (b) has
ﬁnite points for any b ∈ S. Hence, # KKT C (b) < ∞ for all b in S. The proof
is complete. ⊓ ⊔
Corollary 5.1 Assume that C is a polyhedral cone given by (5.2) and the
matrix A is full rank. Then there exists a generic set Gd in Hd such that
# Sol(C, g ) < ∞ for any g ∈ G d.
Proof Since Rη and Hd are homeomorphic, with the isomorphism Π : Rη →
Hd deﬁned by Π (b) = bT Xd(x). According to Proposition 5.1, there exists a
generic set S ⊂ Rη such that the Karush-Kuhn-Tucker set KKT( C, b ) is ﬁnite,
for any b ∈ S. Clearly, Gd := Π (S) is generic in Hd. Since Sol( C, b T Xd(x)) ⊂
KKT(C, b ), one has # Sol( C, g ) < ∞ , for any g ∈ G d. ⊓ ⊔
Remark 5.1 If the constraint K is represented by
K = {x ∈ Rn : q1(x) ≤ 0, . . . , q m(x) ≤ 0} , (5.5)
where q1, . . . , q m are convex polynomials, then the recession cone of K is a
nonempty polyhedral cone. We denote
K j = {x ∈ Rn : qj(x) ≤ 0} , j ∈ [m].
For each j ∈ [m], K j is closed convex set, and K j
∞ is polyhedral (see [9, p.39]).
Since K = K 1 ∩ · · · ∩K m, according to [3, Proposition 2.1.9], one has
K∞ =
⋂
j∈ [m]
K j
∞ .
If follows that K∞ is a nonempty polyhedral cone. Hence, there exists a matrix
A ∈ Rp× n such that
K∞ = {x ∈ Rn : Ax ≥ 0}. (5.6)
Theorem 5.1 Assume that K be represented by (5.5) and the cone K∞ rep-
resented by (5.6), where A is full rank. Then, the set Rd is generic in Rd[x].

<!-- page 15 -->
On the solution existence and stability of polynomial optim ization problems 15
Proof From Remark 5.1, the recession cone K∞ is a nonempty polyhedral cone,
where K∞ = {x ∈ Rn : Ax ≥ 0}. According to Corollary 5.1, there exists a
generic set Gd in Hd such that # Sol( K∞ , g ) < ∞ for any g ∈ G d. Because of
the direct sum Rd[x] = Hd ⊕ Rd− 1[x], the set Gd ⊕ Rd− 1[x] is generic in Rd[x].
It is easy to check that Gd ⊕ Rd− 1[x] ⊂ R d. Hence, Rd is generic in Rd[x]. ⊓ ⊔
Example 5.1 Consider the problem OP( K, f ) given in Example 2.1, we see
that
R2 = {a2x2 + a1x + a0 : a2 ̸= 0, a 1 ∈ R, a 0 ∈ R}
is open and dense in R2[x].
Perspectives
The regularity condition enables us to investigate the stability of the optimal
value function of polynomial optimization problems. Furthermore, t he regu-
larity condition is useful to study the connectedness of the solutio n sets of
convex polynomial vector optimization problems.
Acknowledgements The author would like to thank the anonymous referees for the ir
corrections and comments. This work has been supported by Eu ropean Union’s Horizon
2020 research and innovation programme under the Marie Sk/suppress lodowska-Curie Actions, grant
agreement 813211 (POEMA).
References
1. Lee, G.M., Tam, N.N., Yen, N.D.: Quadratic Programming an d Aﬃne Variational In-
equalities: A Qualitative Study. Springer Verlag, New York (2005)
2. Tam, N.N., Nghi, T.V.: On the solution existence and stabi lity of quadratically con-
strained nonconvex quadratic programs. Optim. Lett. 12, 1045–1063 (2018)
3. Auslender, A., Teboulle, M.: Asymptotic Cones and Functi ons in Optimization and
Variational Inequalities. Springer-Verlag, New York (200 3)
4. Cottle, R.W., Pang, J.-S., Stone, R.E.: The Linear Comple mentarity Problem. Aca-
demic, Boston (1992)
5. Gowda, M.S.: Polynomial complementarity problems. Pac. J. Optim. 13, 227–241 (2017)
6. Hieu, V.T.: Solution maps of polynomial variational ineq ualities, J. Global Optim. 77,
807–824(2020)
7. Frank, M., W olfe, P.: An algorithm for quadratic programm ing. Naval Res. Logist. Q.
3, 95–110 (1956)
8. Luo, Z.-Q., Zhang, S.: On extensions of the Frank-W olfe th eorems. Comput. Optim.
Appl. 13, 87–110 (1999)
9. Belousov, E.G., Klatte, D.: A Frank-W olfe type theorem fo r convex polynomial pro-
grams. Comput. Optim. Appl. 22, 37–48 (2002)
10. Obuchowska, W.T.: On generalizations of the Frank–W olf e theorem to convex and quasi-
convex programmes. Comput. Optim. Appl. 33, 349–364 (2006)
11. Dinh, S.T., Ha, H.V., Pham, T.S.: A Frank-W olfe type theo rem for nondegenerate
polynomial programs. Math. Program. 147, 519–538 (2014)
12. Klatte, D.: On a Frank-W olfe type theorem in cubic optimi zation. Optimization 68,
539–547 (2019)
13. Eaves, B.C.: On quadratic programming. Manag. Sci. 17, 698–711 (1971)

<!-- page 16 -->
16 Vu Trung Hieu
14. Kim, D.S., Tam, N.N., Yen, N.D.: Solution existence and s tability of quadratically
constrained convex quadratic programs. Optim. Lett. 6, 363–373 (2012)
15. Nguyen, H.Q., Nguyen, V.B., Sheu, R.L.: Extension of Eav es Theorem for determining
the boundedness of convex quadratic programming problems. Taiwanese J. Math., 24,
1551–1563 (2020)
16. Lee, G.M., Pham, T.S.: Stability and genericity for semi -algebraic compact programs.
J. Optim. Theory Appl. 169, 473–495 (2016)
17. Karamardian S.: Complementarity problems over cones wi th monotone and pseu-
domonotone maps. J. Optim. Theory Appl. 18, 445–454 (1976)
18. Ansari, Q.H,, Lalitha, C.S,, Mehta, M.: Generalized Con vexity, Nonsmooth Variational
Inequalities, and Nonsmooth Optimization. Chapman and Hal l/CRC (2014)
19. Rockafellar, R.T., W ets, R.J.-B.: Variational Analysi s. Springer-Verlag, Berlin (2009)
20. Li, C., Mordukhovich, B.S., Pham, T.S.: New fractional e rror bounds for polynomial
systems with applications to Holderian stability in optimi zation and spectral theory of
tensors. Math. Program. Ser. A 153, 333–362 (2014)
21. Bochnak, R., Coste, M., Roy, M.F.: Real Algebraic Geomet ry. Springer-Verlag, Berlin
(1998)
22. Dang, V.D., Ha, H.V., Pham, T.S.: W ell-posedness in unco nstrained polynomial opti-
mization problems. SIAM J. Optim. 26, 1411–1428 (2016)
23. Tu, L.W.: An Introduction to Manifolds. Springer-Verla g, New York (2010)