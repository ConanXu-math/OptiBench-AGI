# A projected primal-dual splitting for solving constrained monotone inclusions

**arXiv ID:** 1805.11687v1

**Authors:** Luis Briceño-Arias, Sergio López Rivera

**Abstract:** In this paper we provide an algorithm for solving constrained composite primal-dual monotone inclusions, i.e., monotone inclusions in which a priori information on primal-dual solutions is represented via closed convex sets. The proposed algorithm incorporates a projection step onto the a priori information sets and generalizes the method proposed in [Vũ, B.C.: A splitting algorithm for dual monotone inclusions involving cocoercive operators. Adv. Comput. Math. 38, 667-681 (2013)]. Moreover, under the presence of strong monotonicity, we derive an accelerated scheme inspired on [Chambolle, A.; Pock, T.: A first-order primal-dual algorithm for convex problems with applications to imaging. J. Math. Imaging Vis. 40, 120-145 (2011)] applied to the more general context of constrained monotone inclusions. In the particular case of convex optimization, our algorithm generalizes the methods proposed in [Condat, L.: A primal-dual splitting method for convex optimization involving lipschitzian, proximable and linear composite terms. J. Optim. Theory Appl. 158,460-479 (2013)] allowing a priori information on solutions and we provide an accelerated scheme under strong convexity. An application of our approach with a priori information is constrained convex optimization problems, in which available primal-dual methods impose constraints via Lagrange multiplier updates, usually leading to slow algorithms with unfeasible primal iterates. The proposed modification forces primal iterates to satisfy a selection of constraints onto which we can project, obtaining a faster method as numerical examples exhibit. The obtained results extend and improve several results in the literature.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
A projected primal-dual splitting for solving constrained monotone
inclusions∗
Luis M. Brice˜ no-Arias and Sergio L´ opez Rivera1
1Universidad T´ ecnica Federico Santa Mar´ ıa
Departamento de Matem´ atica
March 3, 2022
Abstract
In this paper we provide an algorithm for solving constrained composite primal-dual mono-
tone inclusions, i.e., monotone inclusions in which a priori information on primal-dual solutions
is represented via closed convex sets. The proposed algorithm incorporates a projection step onto
the a priori information sets and generalizes the method proposed in [2]. Moreover, under the
presence of strong monotonicity, we derive an accelerated scheme inspired on [3] applied to the
more general context of constrained monotone inclusions. In the particular case of convex opti-
mization, our algorithm generalizes the methods proposed in [1, 3] allowing a priori information
on solutions and we provide an accelerated scheme under strong convexity. An application of our
approach with a priori information is constrained convex optimization problems, in which avail-
able primal-dual methods impose constraints via Lagrange multiplier updates, usually leading
to slow algorithms with unfeasible primal iterates. The proposed modiﬁcation forces primal iter-
ates to satisfy a selection of constraints onto which we can project, obtaining a faster method as
numerical examples exhibit. The obtained results extend and improve several results in [1, 2, 3].
Keywords: accelerated schemes, constrained convex optimization, monotone operator theory,
proximity operator, splitting algorithms
∗Contact author: L. M. Brice˜ no-Arias, luis.briceno@usm.cl,
1
arXiv:1805.11687v1  [math.OC]  29 May 2018

<!-- page 2 -->
1 Introduction
This paper is devoted to the numerical resolution of composite primal-dual monotone inclusions in
which a priori information on solutions is known. The relevance of monotone inclusions and convex
optimization is justiﬁed via the increasing number of applications in several ﬁelds of engineering
and applied mathematics as image processing, evolution inclusions, variational inequalities, learning,
partial diﬀerential equations, Mean Field Games, among others (see, e.g., [4, 5, 6, 7, 8, 9, 10] and
references therein). The a priori information on primal-dual solutions is represented via closed
convex sets in primal and dual spaces, following some ideas developed in [11, 12]. We force primal-
dual iterates to belong to these information sets by adding additional projections on primal-dual
iterates in each iteration of our proposed method.
An important instance, in which the advantage of our formulation arises, is composite convex op-
timization with aﬃne linear equality constraints. In this context, the primal-dual methods proposed
in [1, 2, 3, 13, 14, 15, 16] impose feasibility through Lagrange multiplier updates. A disadvantage of
this approach is that such algorithms are usually slow and their primal iterates do not necessarily
satisfy any of the constraints (see, e.g., [7]), leading to unfeasible approximate primal solutions.
By projecting onto the aﬃne subspace generated by the constraints, previous problem is solved.
However, in several applications this projection is not easy to compute because of singularity or bad
conditioning on the linear system (see, e.g. [17]). In this context, the a priori information on primal
solutions can be set as any selection of the aﬃne linear constraints. Indeed, since any solution is
feasible, we know it must satisfy any selection of the constraints. Even if in the previous context
the formulation with a priori information may be seen as artiﬁcial, in a practical point of view this
formulation allows us to propose a method with an additional projection onto an arbitrary selection
of the constraints, which improves the eﬃciency of the method (see Section 4.2). This method forces
primal iterates to satisfy the selection of the constraints, which can be chosen in order to compute
the projection easily.
In this paper we provide a new projected primal-dual splitting method for solving constrained
monotone inclusions, i.e., inclusions in which we count on a priori information on primal-dual solu-
tions. We also provide an accelerated scheme of our method in the presence of strong monotonicity
and we derive linear convergence in the fully strongly monotone case. In the case without a pri-
ori information, our results give an accelerated scheme of the method proposed in [2] for strongly
monotone inclusions. In the context of convex optimization, our method generalize the algorithms
proposed in [1, 3] and [16] without inertia, by incorporating a projection onto an the a priori primal-
dual information set. This method is applied in the context of convex optimization with equality
constraints, when the a priori information set is chosen as a selection of the aﬃne linear constraints
in which it is easy to project. The advantages of this approach with respect to classical primal-dual
approaches are justiﬁed via numerical examples. Our acceleration scheme in the convex optimiza-
tion context is obtained as a generalization of [3], complementing the ergodic rates obtained in the
case without projection in [15] and, as far as we know, have not been developed in the literature
and are interesting in their own right.
The paper is organized as follows. In Section 2 we set our notation and we give a brief background.
In section 3, we set the constrained primal-dual monotone inclusion and we propose our algorithm
together with the main results. We also provide connections with existing methods in the literature.
In Section 4, we apply previous results to convex optimization problems with equality aﬃne linear
constraints, together with numerical experiences illustrating the improvement in the eﬃciency of
2

<!-- page 3 -->
the algorithm with the additional projection. We ﬁnish with some conclusions in Section 5.
2 Notation and preliminaries
LetH andG be real Hilbert spaces. We denote the scalar products of H andG by⟨·|·⟩ and the
associated norms by ∥·∥ . The projector operator onto a nonempty closed convex set C⊂H is
denoted by PC and, for a set-valued operator M :H→ 2H we use ran(M) for the range of M,
gra(M) for its graph, M−1 for its inverse, JM = (Id +M)−1 for its resolvent, and □ stands for the
parallel sum as in [18]. Moreover, M isρ-strongly monotone if, for every (x,u ) and (y,v ) in gra(M),
⟨x−y,u−v⟩≥ ρ∥x−y∥2, it isρ−cocoercive ifM−1 isρ−strongly monotone,M is monotone if it is
ρ-strongly monotone withρ = 0, and it is maximally monotone if its graph is maximal in the sens of
inclusions inH×H , among the graphs of monotone operators. The class of all lower semicontinuous
convex functions f :H→ (−∞, +∞] such that dom(f) ={x∈H| f(x)< +∞}̸ = ∅ is denoted by
Γ0(H) and, for every f∈ Γ0(H), the Fenchel conjugate of f is denoted by f∗, its subdiﬀerential by
∂f , and its proximity operator by proxf, as in [18]. We recall that ( ∂f )−1 =∂f∗ and J∂f = proxf.
In addition, when C⊂H is a convex closed subset, we have that J∂ιC = proxιC = PC, where ιC
is the indicator function of C, which is 0 in C and +∞ otherwise. Given α∈]0, 1[, an operator
T :H→H satisfying FixT̸= ∅ isα−averaged quasi-nonexpansive if, for everyx∈H andy∈ FixT
we have∥Tx−y∥2≤∥ x−y∥2− ( 1−α
α )∥x−Tx∥2. We refer the reader to [18] for deﬁnitions and
further results in monotone operator theory and convex optimization.
3 Problem and main results
We consider the following problem.
Problem 3.1 LetT :H→H be an α−averaged quasi-nonexpansive operator with α∈]0, 1[, let V
be a closed vector subspace of G, let L:H→G be a nonzero linear bounded operator satisfying
ranL⊂ V , let A :H→ 2H and D :G→ 2G be maximally monotone operators which are ρ and
δ−strongly monotone, respectively, and let B :G→ 2G and C :H→H be χ and β−cocoercive,
respectively, for (ρ,χ )∈ [0, +∞[2 and (δ,β )∈ [0, +∞]2. The problem is to solve the primal and
dual inclusions
ﬁnd ˆ x∈ FixT such that 0 ∈Aˆx +L∗(B □D)(Lˆx) +Cˆx (P)
ﬁnd ˆ u∈V such that ( ∃ ˆx∈ FixT )
{
−L∗ˆu∈Aˆx +Cˆx
ˆu∈ (B □D)(Lˆx), (D)
under the assumption that solutions exist.
When A = ∂f , B = ∂g, C =∇h, and D = ∂ℓ, where f ∈ Γ0(H), g∈ Γ0(G), h:H→ R is a
diﬀerentiable convex function with β−1−Lipschitz gradient, and ℓ∈ Γ0(G) is δ−strongly convex,
Problem 3.1 reduces to
ﬁnd ˆx∈ FixT∩ argminx∈HF (x) :=f(x) + (g □ℓ)(Lx) +h(x) ( P0)
3

<!-- page 4 -->
together with the dual problem
ﬁnd ˆu∈V∩ argminu∈Gg∗(u) + (f∗ □h∗)(−L∗u) +ℓ∗(u), (D0)
assuming that some qualiﬁcation condition holds. Note that, when T = PX, any solution to (P0)
is a solution to min x∈XF (x), but the converse is not true. The set X in this case represents an a
priori information on the primal solution. As you can see in the next section, an application of this
formulation is constrained convex optimization, in which X may represent a selection of the aﬃne
linear constraints. Even if, in this case, the formulation can be set without considering the set X,
its artiﬁcial appearance has a practical relevance: the method obtained include a projection onto
X which helps to the performance of the method as stated in Section 4.
When ρ = χ = 0, V =G and T = Id , (P0)-(D0) can be solved by using [20, Theorem 4.2] or
[16, Theorem 5]. In the last method, inertial terms are also included. In the case when ℓ∗ = 0 the
algorithm in [1] can be used and if ℓ∗ =h = 0, (P0)-(D0) can be solved by [3, 19] or a version of [3]
with linesearch proposed in [23]. In [3], the strong convexity is exploited via acceleration schemes.
Moreover, whenT =PX,X⊂H is nonempty, closed and convex,V =G andℓ∗ =h = 0, (P0)-(D0)
is solved in [7, Theorem 3.1]. When ρ >0 or χ >0, ergodic convergence rates are derived in [15]
whenV =G andT = Id . In its whole generality, as far as we know, (P0)-(D0) has not been solved
and strong convexity has not been exploited.
In Problem 3.1 set T = Id ,V =G =G1⊕···⊕ Gm,L: x↦→ (L1x,...,L mx),B : (u1,...,u m)↦→
ω1B1u1×···× ωmBmum and D : (u1,...,u m)↦→ ω1D1u1×···× ωmDmum, where for every i∈
{1,...,m }, Li :H→ Gi is linear and bounded, Bi : Gi↦→ 2Gi and Di : Gi↦→ 2Gi are maximally
monotone operators such that Di is strongly monotone, and ωi > 0 satisﬁes ∑m
i=1ωi = 1. Then,
Problem 3.1 reduces to [20, Problem 1.1] (see also [2, Problem 1.1]). We prefer to set m = 1
for simplicity. In [20], previous problem is solved when C is monotone and Lipschitz by applying
the method in [21] to the product primal-dual space. Accelerated versions of previous algorithm
under strong monotonicity are proposed in [22]. The cocoercivity of C is exploited in [2], where
an algorithm is proposed for solving Problem 3.1 when ρ = χ = 0, V =G and T = Id . In the
following theorem we provide an algorithm for solving Problem 3.1 in its whole generality with
weak convergence to a solution when the stepsizes are ﬁxed. Moreover, when A orB−1 are strongly
monotone (ρ >0 or χ >0), we provide an accelerated version inspired on (and generalizing) [3,
Section 5.1]. Finally, we generalize [3, Section 5.2] for obtaining linear convergence when ρ> 0 and
χ> 0.
Theorem 3.2 Letγ0∈]0, 2δ[ and τ0∈]0, 2β[ be such that
∥L∥2≤
( 1
τ0
− 1
2β
) ( 1
γ0
− 1
2δ
)
(3.1)
and let (x0, ¯x0,u 0)∈H×H×G such that ¯x0 =x0. Let (θk)k∈N, (γk)k∈N and (τk)k∈N be sequences
in ]0, 1], ]0, 2δ[ and ]0, 2β[, respectively, and consider
(∀k∈ N)

ηk+1 =JγkB−1(uk +γk(L¯xk−D−1uk))
uk+1 =PV ηk+1
pk+1 =JτkA(xk−τk(L∗uk+1 +Cxk))
xk+1 =Tp k+1
¯xk+1 =xk+1 +θk(pk+1−xk).
(3.2)
Then, the following hold.
4

<!-- page 5 -->
(i) For every k∈ N and for every solution (ˆx, ˆu) to Problem 3.1, we have
∥xk− ˆx∥2
τk
+∥uk− ˆu∥2
γk
≥ (2ρτk + 1)∥pk+1− ˆx∥2
τk
+
pk+1−xk

2 ( 1
τk
− 1
2β
)
+ (2χγk + 1)∥ηk+1− ˆu∥2
γk
+∥ηk+1−uk∥2
( 1
γk
− 1
2δ
)
+ 2
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
− 2θk−1
⣨
L(pk−xk−1)|ηk− ˆu
⟩
− 2θk−1∥L∥∥pk−xk−1∥∥ηk+1−uk∥. (3.3)
(ii) Suppose that ρ = 0 and χ = 0. If we set θk≡ 1, τk≡ τ, γk≡ γ and we assume that (3.1)
holds with strict inequality, we obtain xk ⇀ ˆx and uk ⇀ ˆu, for some solution (ˆx, ˆu) to
Problem 3.1.
(iii) Suppose that ρ> 0, χ = 0, and D−1 = 0. If we set
(∀k∈ N) θk = 1√1 + 2ρτk
, τ k+1 =θkτk, γ k+1 =γk/θk, (3.4)
and we assume that (3.1) holds with equality, we obtain, for every solution (ˆx, ˆu) to Prob-
lem 3.1, (∀ε> 0)(∃N0∈ N)(∀k≥N0)
xk− ˆx
2≤ 1 +ε
k2
( x0− ˆx
2
ρ2τ 2
0
+ 2β∥L∥2
ρ2(2β−τ0)
u0− ˆu
2
)
. (3.5)
(iv) Suppose that ρ> 0 and χ> 0 and deﬁne
µ = 2√ρχ
∥L∥ and α = min
{
µρ
ρ + µ
4β
, µχ
χ + µ
4δ
}
. (3.6)
If we set θk≡θ∈ ((1 +α)−1, 1], τk≡τ and γk≡γ with
τ = 2βµ
µ + 4βρ and γ = 2µδ
µ + 4δχ, (3.7)
we obtain linear convergence. That is, for every k∈ N,
(
χ(1−ω) + µ
4δ
)uk− ˆu

2
+
(
ρ + µ
4β
) xk− ˆx

2
≤ωk
((
χ + µ
4δ
)u0− ˆu
2
+
(
ρ + µ
4β
) x0− ˆx
2
)
, (3.8)
whereω = (1 +θ)/(2 +α)∈ [(1 +α)−1,θ ).
Proof. (i): Fix k∈ N and let (ˆx, ˆu) be a solution to Problem 3.1. We have ˆx∈ FixT , ˆu∈V and,
using B □D = (B−1 +D−1)−1, we deduce−(L∗ˆu +Cˆx)∈Aˆx and Lˆx−D−1ˆu∈B−1ˆu. Therefore,
since (3.2) yields {xk−pk+1
τk
−L∗uk+1−Cxk∈Apk+1
uk−ηk+1
γk
+L¯xk−D−1uk∈B−1ηk+1 (3.9)
5

<!-- page 6 -->
and A and B−1 are ρ and χ-strongly monotone, respectively, we deduce
⟨xk−pk+1
τk
−L∗(uk+1− ˆu)
⏐⏐⏐⏐pk+1− ˆx
⟩
+
⟨uk−ηk+1
γk
+L(¯xk− ˆx)
⏐⏐⏐⏐ηk+1− ˆu
⟩
−
⣨
Cxk−Cˆx|pk+1− ˆx
⟩
−
⣨
D−1uk−D−1ˆu|ηk+1− ˆu
⟩
≥ρ
pk+1− ˆx

2
+χ∥ηk+1− ˆu∥2. (3.10)
From the cocoercivity of C and D−1 we have from ab≤βa2 +b2/(4β) that
⣨
Cxk−Cˆx|pk+1− ˆx
⟩
=
⣨
Cxk−Cˆx|pk+1−xk
⟩
+
⣨
Cxk−Cˆx|xk− ˆx
⟩
≥−∥Cxk−Cˆx∥∥pk+1−xk∥ +β∥Cxk−Cˆx∥2
≥−∥pk+1−xk∥2
4β , (3.11)
and, analogously,
⟨
D−1uk−D−1ˆu|ηk+1− ˆu
⟩
≥− ∥ηk+1−uk∥2
4δ . Hence, by using [18, Lemma 2.12(i)]
in (3.10), we deduce
∥xk− ˆx∥2
τk
+∥uk− ˆu∥2
γk
≥
(
2ρ + 1
τk
) pk+1− ˆx

2
+
(
2χ + 1
γk
)
∥ηk+1− ˆu∥2
+ 2
[⣨
L(pk+1− ˆx)|uk+1− ˆu
⟩
−
⣨
L(¯xk− ˆx)|ηk+1− ˆu
⟩]
+∥ηk+1−uk∥2
(1
γk
− 1
2δ
)
+∥pk+1−xk∥2
(1
τk
− 1
2β
)
. (3.12)
Moreover, (3.2), ran(L)⊂V and uk−ηk∈V⊥, for every k∈ N, yield
⣨
L(pk+1− ˆx)|uk+1− ˆu
⟩
−
⣨
L(¯xk− ˆx)|ηk+1− ˆu
⟩
=
⣨
L(pk+1− ˆx)|uk+1− ˆu
⟩
−
⣨
L(xk− ˆx)|ηk+1− ˆu
⟩
−θk−1
⣨
L(pk−xk−1)|ηk+1− ˆu
⟩
=
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
−θk−1
⣨
L(pk−xk−1)|ηk+1− ˆu
⟩
=
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
−θk−1
⣨
L(pk−xk−1)|ηk+1−uk
⟩
−θk−1
⣨
L(pk−xk−1)|ηk− ˆu
⟩
≥
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
−θk−1∥L∥∥pk−xk−1∥∥ηk+1−uk∥
−θk−1
⣨
L(pk−xk−1)|ηk− ˆu
⟩
, (3.13)
which, together with (3.12), yield (3.3).
(ii): For every k∈ N, it follows from Theorem 3.2((i)), [24, Lemma 2.1], ρ = χ = 0, θk≡ 1,
6

<!-- page 7 -->
τk≡τ, γk≡γ, and the properties of T and PX that
∥pk− ˆx∥2
τ +∥ηk− ˆu∥2
γ ≥∥pk+1− ˆx∥2
τ +
(1−α
α
)∥xk−pk∥2
τ +∥uk−ηk∥2
γ
+∥ηk+1− ˆu∥2
γ +∥pk+1−xk∥2
(1
τ− 1
2β
)
+∥ηk+1−uk∥2
(1
γ− 1
2δ
)
+ 2
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
− 2
⣨
L(pk−xk−1)|ηk− ˆu
⟩
− 2∥L∥∥pk−xk−1∥∥ηk+1−uk∥
≥∥pk+1− ˆx∥2
τ +∥uk−ηk∥2
γ +∥ηk+1−uk∥2
(1
γ− 1
2δ− 1
ν
)
+
(1−α
α
)∥xk−pk∥2
τ +∥ηk+1− ˆu∥2
γ +∥pk+1−xk∥2
(1
τ− 1
2β
)
+ 2
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
− 2
⣨
L(pk−xk−1)|ηk− ˆu
⟩
−ν∥L∥2∥pk−xk−1∥2, (3.14)
for every ν > 0. If we let ε =
[(
1
τ− 1
2β
) (
1
γ− 1
2δ
)
−∥L∥2
] (
βτ
2β−τ
)
> 0, and we choose ν =
( 1
γ− 1
2δ−ε)−1> 0, we have ν∥L∥2 = ( 1
τ− 1
2β )−νε( 1
τ− 1
2β ). Hence, from (3.14) we have
Υk +
pk− ˆx
2
τ ≥ Υk+1 +
pk+1− ˆx
2
τ +
(1−α
α
)∥xk−pk∥2
τ +∥uk−ηk∥2
γ
+ε∥ηk+1−uk∥2 +νε
(1
τ− 1
2β
)
∥pk−xk−1∥2, (3.15)
where, for every k∈ N,
Υk =∥ηk− ˆu∥2
γ + 2
⣨
L(pk−xk−1)|ηk− ˆu
⟩
+
(1
τ− 1
2β
)
∥pk−xk−1∥2. (3.16)
Note that from (3.1) we have, for every k∈ N,
Υk≥∥ηk− ˆu∥2
γ + 2
⣨
L(pk−xk−1)|ηk− ˆu
⟩
+ ∥L∥2
(
1
γ− 1
2δ
)∥pk−xk−1∥2
≥∥ηk− ˆu∥2
γ + 2
⣨
L(pk−xk−1)|ηk− ˆu
⟩
+γ∥L∥2∥pk−xk−1∥2
≥ 1
γ∥ηk− ˆu +γL(pk−xk−1)∥2≥ 0, (3.17)
and, hence, from (3.15) we deduce that (Υk +∥pk− ˆx∥2/τ)k∈N is a F´ ejer sequence. We deduce from
[18, Lemma 5.31] that (ηk)k∈N and (pk)k∈N are bounded,
xk−pk→ 0, u k−ηk→ 0, η k+1−uk→ 0, and pk−xk−1→ 0. (3.18)
Therefore, there exist weak accumulation points ¯ x and ¯u of the sequences ( pk)k∈N and (ηk)k∈N,
respectively, say pkn ⇀ ¯x and ηkn ⇀ ¯u and, from (3.18), we have ukn ⇀ ¯u, ukn+1 ⇀ ¯u, pkn ⇀ ¯x,
7

<!-- page 8 -->
pkn+1 ⇀ ¯x, xkn−1 ⇀ ¯x and ¯xkn = xkn +pkn−xkn−1 ⇀ ¯x. Since T and PV are nonexpansive,
Id−T and Id−PV are maximally monotone [18, Example 20.29] and, therefore, they have weak-
strong closed graphs [18, Proposition 20.38]. Hence, it follows from (3.18) that (Id −T )pk→ 0 and
(Id−PV )ηk→ 0 and, hence, (¯x, ¯u)∈ FixT×V . Moreover, (3.9) can be written equivalently as
(vkn,wkn)∈ (M + Q)(pkn+1,ηkn+1), (3.19)
where M : (p,η )↦→ (Ap +L∗η)× (B−1η−Lp) is maximally monotone [19, Proposition 2.7(iii)],
Q: (p,η )↦→ (Cp,D−1η) is min{β,δ}−cocoercive, and
{
vk := xk−pk+1
τ −L∗(uk+1−ηk+1) +Cpk+1−Cxk
wk := uk−ηk+1
γ +L(xk−pk+1 +pk−xk−1) +D−1ηk+1−D−1uk.
(3.20)
It follows from [18, Corollary 25.5] that M + Q is maximally monotone and, since (3.18) and the
uniform continuity of C, D and L yields vkn→ 0 and wkn→ 0, we deduce from the weak-strong
closedness of the graph of M + Q that (¯x, ¯u) is a solution to Problem 3.1, and the result follows.
(iii): Fix k∈ N. Since ρ> 0, δ = +∞, χ = 0, from Theorem 3.2((i)) we have
∥xk− ˆx∥2
τk
+∥uk− ˆu∥2
γk
≥ (2ρτk + 1)τk+1
τk
∥pk+1− ˆx∥2
τk+1
+ γk+1
γk
∥ηk+1− ˆu∥2
γk+1
+ 2
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
− 2θk−1
⣨
L(pk−xk−1)|ηk− ˆu
⟩
+
pk+1−xk

2 ( 1
τk
− 1
2β
)
−θ2
k−1γk∥L∥2∥pk−xk−1∥2, (3.21)
where we have used 2ab≤a2/γ +γb2. Moreover, it follows from (3.4) that
(∀k∈ N) (1 + 2 ρτk)τk+1
τk
= (1 + 2ρτk)θk = 1
θk
= γk+1
γk
, (3.22)
which, combined with (3.21), yields
∥xk− ˆx∥2
τk
+∥uk− ˆu∥2
γk
≥ 1
θk
(∥pk+1− ˆx∥2
τk+1
+∥ηk+1− ˆu∥2
γk+1
)
+ 2
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
− 2θk−1
⣨
L(pk−xk−1)|ηk− ˆu
⟩
+
pk+1−xk

2 ( 1
τk
− 1
2β
)
−θ2
k−1γk∥L∥2∥pk−xk−1∥2. (3.23)
Now deﬁne
(∀k∈ N) ∆ k =
xk− ˆx
2
τk
+
uk− ˆu
2
γk
. (3.24)
Dividing (3.23) by τk and using θkτk =τk+1 we obtain from the nonexpansivity of PV and T that
∆k
τk
≥ ∆k+1
τk+1
+ 2
τk
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
− 2
τk−1
⣨
L(pk−xk−1)|ηk− ˆu
⟩
+
pk+1−xk2
τ 2
k
(
1− τk
2β
)
−γkτk∥L∥2∥pk−xk−1∥2
τ 2
k−1
. (3.25)
8

<!-- page 9 -->
In addition, (3.1) with equality reduces to
∥L∥2 =
( 1
τ0
− 1
2β
) 1
γ0
⇔ γ0τ0∥L∥2 =
(
1− τ0
2β
)
. (3.26)
Since, for every k∈ N\{ 0}, γkτk =γ0τ0 and{τk}k∈N is decreasing (see (3.4)), we have from (3.26)
that
γkτk∥L∥2 =γ0τ0∥L∥2 =
(
1− τ0
2β
)
≤
(
1− τk−1
2β
)
, (3.27)
and (3.25) yields
∆k
τk
≥ ∆k+1
τk+1
+
pk+1−xk2
τ 2
k
(
1− τk
2β
)
−∥pk−xk−1∥2
τ 2
k−1
(
1− τk−1
2β
)
+ 2
τk
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
− 2
τk−1
⣨
L(pk−xk−1)|ηk− ˆu
⟩
. (3.28)
Now ﬁx N≥ 1. By adding from k = 0 to k = N− 1 in (3.28), using that ¯x0 = x0 and deﬁning
p0 =x0 =:x−1, we obtain from uN =PVηN, and ranL⊂V that
∆0
τ0
≥ ∆N
τN
+
pN−xN−12
τ 2
N−1
(
1− τN−1
2β
)
+ 2
τN−1
⟨
L(pN−xN−1)|uN− ˆu
⟩
≥ ∆N
τN
− ∥L∥2
(1− τN−1
2β )
uN− ˆu
2
= 1
τN
(
∆N− γNτN∥L∥2
(1− τN−1
2β )
uN− ˆu
2
γN
)
≥∥xN− ˆx∥2
τ 2
N
, (3.29)
where the last inequality follows from (3.27) and (3.24). Multiplying (3.29) byτ 2
N and usingγNτN =
γ0τ0 and (3.26), we conclude that
xN− ˆx
2
≤τ 2
N
( x0− ˆx
2
τ 2
0
+ ∥L∥2
(1− τ0
2β )
u0− ˆu
2
)
. (3.30)
The result follows from limN→∞NρτN = 1 [3, Corollary 1].
(iv): Fix k∈ N. Note that (3.7) yields
(
1
τ− 1
2β
) (
1
γ− 1
2δ
)
=∥L∥2 and, from (3.12), uk+1 =
PVηk+1 and ranL⊂V we have
uk− ˆu
2
2γ +
xk− ˆx
2
2τ ≥ (2ρτ + 1)∥pk+1− ˆx∥2
2τ + (2χγ + 1)∥ηk+1− ˆu∥2
2γ
+
pk+1−xk2
2
(1
τ− 1
2β
)
+∥ηk+1−uk∥2
2
(1
γ− 1
2δ
)
+
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
. (3.31)
Hence, by deﬁning
(∀k∈ N) Ω k :=
(
χ + µ
4δ
)
∥uk− ˆu∥2 +
(
ρ + µ
4β
)
∥xk− ˆx∥2, (3.32)
9

<!-- page 10 -->
multiplying (3.31) by µ and using (3.6), (3.7), uk+1 = PV (ηk+1), ran(L)⊂ V , and the nonexpan-
sivity of T and PV we have
Ωk≥ Ωk+1 +µρ
pk+1− ˆx
2 +µχ
ηk+1− ˆu
2 +ρ
pk+1−xk2
+χ
ηk+1−uk2 +µ
⣨
L(pk+1−xk)|ηk+1− ˆu
⟩
≥ (1 +α)Ωk+1 +ρ
pk+1−xk2 +χ
uk+1−uk2
+µ
⣨
L(pk+1−xk)|uk+1− ˆu
⟩
. (3.33)
Moreover, for every ω,λ> 0 we have
µ
⣨
L(pk+1−xk)|uk+1− ˆu
⟩
=µ
⣨
L(pk+1−xk)|uk+1− ˆu
⟩
−µθ
⣨
L(pk−xk−1)|uk+1− ˆu
⟩
=µ
⣨
L(pk+1−xk)|uk+1− ˆu
⟩
−ωµ
⣨
L(pk−xk−1)|uk− ˆu
⟩
−ωµ
⣨
L(pk−xk−1)|uk+1−uk
⟩
− (θ−ω)µ
⣨
L(pk−xk−1)|uk+1− ˆu
⟩
≥µ
⣨
L(pk+1−xk)|uk+1− ˆu
⟩
−ωµ
⣨
L(pk−xk−1)|uk− ˆu
⟩
−ωµ∥L∥
(
λ
pk−xk−12
2 +
uk+1−uk2
2λ
)
− (θ−ω)µ∥L∥
(
λ
pk−xk−12
2 +
uk+1− ˆu
2
2λ
)
=µ
⣨
L(pk+1−xk)|uk+1− ˆu
⟩
−ωµ
⣨
L(pk−xk−1)|uk− ˆu
⟩
−µθλ∥L∥
pk−xk−12
2 − ωµ∥L∥
uk+1−uk2
2λ
− (θ−ω)µ∥L∥
uk+1− ˆu
2
2λ . (3.34)
By choosing λ =ω
√ρ
χ, from (3.33), (3.34) and (3.6), we obtain
Ωk≥ Ωk+1
ω +
(
1 +α− 1
ω
)
Ωk+1 +ρ
pk+1−xk

2
+µ
⣨
L(pk+1−xk)|uk+1− ˆu
⟩
−ωµ
⣨
L(pk−xk−1)|uk− ˆu
⟩
−ωθρ
pk−xk−1

2
−
(θ−ω
ω
)
χ
uk+1− ˆu

2
.
(3.35)
Since θ∈
]
(1 +α)−1, 1
]
, by setting ω = 1 +θ
2 +α∈
[
(1 +α)−1,θ
[
, we have 1 + α− 1
ω = θ−ω
ω > 0.
Hence, since (3.32) yields Ωk+1≥χ
uk+1− ˆu
2, from (3.35) and θ≤ 1 we have
Ωk≥ Ωk+1
ω +ρ
pk+1−xk

2
−ωρ
pk−xk−1

2
+µ
⣨
L(pk+1−xk)|uk+1− ˆu
⟩
−ωµ
⣨
L(pk−xk−1)|uk− ˆu
⟩
.
(3.36)
10

<!-- page 11 -->
Moreover, using p0 =x0 =:x−1, multiplying (3.36) by ω−k and adding from k = 0 to k =N− 1,
we conclude from the deﬁnition of µ that
Ω0≥ω−NΩN +ω−N +1ρ
pN−xN−12
+µω−N +1⟨
L(pN−xN−1)|uN− ˆu
⟩
≥ω−NΩN +ω−N +1ρ
pN−xN−12
−µω−N +1∥L∥
(√ρ
χ
pN−xN−12
2 +
√χ
ρ
uN− ˆu
2
2
)
=ω−NΩN−ω−N +1χ
uN− ˆu
2
, (3.37)
or, equivalently,
ωN
((
χ + µ
4δ
)u0− ˆu
2
+
(
ρ + µ
4β
) x0− ˆx
2
)
≥
(
χ(1−ω) + µ
4δ
)uN− ˆu
2
+
(
ρ + µ
4β
) xN− ˆx
2
, (3.38)
which proves the linear convergence since ω <θ≤ 1.
Remark 3.3 (i) Note that condition (3.1) is weaker than the condition needed in [2]. Indeed,
this condition in our case reads 2 ρ min{β,δ}> 1, where ρ = min
{
γ−1,τ−1}
(1−
√
τγ∥L∥2),
which implies 2 min{δ,β}> 1
ρ > max{γ,τ},
(
1− τ
2β
)
>
√
τγ∥L∥2 and
(
1− γ
2δ
)
>
√
τγ∥L∥2. (3.39)
Thus, by multiplying last expressions we obtain
(
1− τ
2β
)(
1− γ
2δ
)
> τγ∥L∥2, which implies
(3.1). Our condition is strictly weaker, as it can be seen in Figure 1, in which we plot the case
∥L∥ = 1 and δ =β =b, for b = 1, b = 1/2 and b = 1/4. That is, we compare regions
Rb =
{
(τ,γ )∈ [0, 2b]× [0, 2b]
⏐⏐ min
{ 1−√τγ
τ , 1−√τγ
γ
}
> 1
2b
}
(3.40)
Sb =
{
(τ,γ )∈ [0, 2b]× [0, 2b]
⏐⏐
(
1− τ
2b
) (
1− γ
2b
)
>τγ
}
. (3.41)
(ii) It is not diﬃcult to extend our method by replacing the averaged quasi-nonexpansive operator
T by (αk)k∈N−averaged quasi-nonexpansive operators (Tk)k∈N varying at each iteration and
satisfying supk∈Nαk < 1. Indeed, as in [24], we have to assume that xk−Tkxk → 0 and
xk ⇀ ximpliesx∈∩ k∈N FixTk, which is satisﬁed in several cases. In particular, if we set, for
everyk∈ N,Tk :=JγkM(Id−γkN), where M :H→ 2H is maximally monotone and N :H→
H is ξ-cocoercive, if γk ∈]0, 2ξ[, Tk is γk/2ξ−cocoercive and ∩k∈N FixTk = zer(M +N).
Therefore, our method using these operators leads to the common solution to zer( M +N)
and zer(A +L∗◦B◦L +C). Previous example can also be tackled by Theorem 3.2 if we use
γk≡ γ and Tk≡ T := JγM (Id−γN ) and we prefer to keep the constant operator case for
avoiding additional hypotheses and for the sake of simplicity.
(iii) The method proposed in [22] is an accelerated version of the method proposed in [20] under
the assumption that A +C is strongly monotone. Of course, this weaker assumption can also
be used in our context, but we prefer to keep the statement of Theorem 3.2 simpler.
11

<!-- page 12 -->
0.0 0.5 1.0 1.5 2.0
0.0
0.5
1.0
1.5
2.0
0.0 0.2 0.4 0.6 0.8 1.0
0.0
0.2
0.4
0.6
0.8
1.0
0.0 0.1 0.2 0.3 0.4 0.5
0.0
0.1
0.2
0.3
0.4
0.5
Figure 1: We plot regions Rb in blue and Sb in orange. Left: case b = 1, Center: case b = 1/2,
Right: case b = 1/4. Note that in the case τ =γ the regions coincide.
(iv) Theorem 3.2((iii)) generalizes the acceleration scheme proposed in [3] to monotone inclusions
with a priori information. From this results we derive an accelerated version of the methods
in [2] in the strongly monotone case when T = Id and V =G. These accelerated version, as
far as we know, have not been developed in the literature.
(v) In the context of primal-dual problem ( P0)-(D0), (3.2) reduces to
(∀k∈ N)

ηk+1 = proxγkg∗(uk +γk(L¯xk−∇ℓ∗(uk)))
uk+1 =PV ηk+1
pk+1 = proxτkf
(
xk−τk
(
L∗uk+1 +∇h(xk)
))
xk+1 =Tp k+1
¯xk+1 =xk+1 +θk(pk+1−xk),
(3.42)
and our conditions on the parameters coincide with [15, 16]. Without strong convexity of f
andg∗, we deduce from Theorem 3.2((ii)) the weak convergence of the sequences generated by
(3.42), generalizing results in [1, 7, 16]. When f or g∗ is strongly convex, Theorem 3.2((iii))
yields an accelerated and projected version of [1]. When V =G and T = Id , this result
complements the ergodic convergence rates obtained in [15] and generalizes [3]. When ℓ∗ = 0,
V =G, T = Id and f and g∗ are strongly convex Theorem 3.2((iv)) yields non-ergodic linear
convergence of [1], complementing the ergodic linear convergence in [15]. The advantage of
the algorithm (3.42) with respect to [1, 3] is that primal-dual iterates of the former are forced
to be in X×V when T = PX. This feature leads to a faster algorithm in the context of
constrained convex optimization, by choosing X to be some of the constraints. This can be
observed in the particular instance developed in [7] and in Section 4, in which we provide
some numerical simulations.
4 Application to constrained convex optimization
In this section, we explore the advantages of the proposed method in constrained convex optimiza-
tion.
12

<!-- page 13 -->
4.1 Constrained convex optimization problem
Problem 4.1 Let f∈ Γ0(RN), let R and S be m×N and n×N real matrices, respectively, and
let c∈ Rm and d∈ Rn. The problem is to
min
x∈RN
f(x) s.t. Rx =c Sx =d, (4.1)
under the assumption that solutions exist.
Note that (4.1) can be written equivalently as min x∈RNf(x) +ι{b}(Lx), where L: x↦→ (Rx,Sx )
andb = (c,d )∈ Rm+n. Assume that 0∈ sri(L(domf)−b). Note that, since proxγι∗
{b}
= Id−γb [18,
Proposition 24.8(ix)], the method proposed in [3, Algorithm 1] in this case reads: given x0 = ¯x0∈H
and u0∈G ,
(∀k∈ N)

uk+1 =uk +γ(L¯xk−b)
xk+1 = proxτf (xk−τL∗uk+1)
¯xk+1 = 2xk+1−xk,
(4.2)
where γτ∥L∥2 < 1. The constraint is imposed via the Lagrange multiplier update in the ﬁrst
step of (4.2). This implies that the primal sequence {xk}k∈N does not necessarily satisfy any of
the constraints. For ensuring feasibility, we should project onto L−1b by considering the problem
minx∈RNf(x) +ιL−1b(x). However, this is not always possible since, in several applications, the
matrices involved are singular or very bad conditioned (see discussion in [7, 17]). If it is diﬃcult
to compute PL−1b but we can project onto R−1c, we can rewrite (4.1) as the problem of ﬁnding
ˆx∈R−1c∩ argminx∈RNf(x) +ι{b}(Lx), which is (P0) when X =R−1c, h =ℓ∗ = 0, and g =ι{b}.
Next corollary follows from Theorem 3.2, (3.42) and PX : x↦→x−R∗(RR∗)−1(Rx−c).
Corollary 4.2 Letγ >0 andτ >0 be such thatγτ∥L∥2< 1 and let (x0, ¯x0,u 0)∈ RN×RN×Rm+n
be such that x0 = ¯x0. Consider the routine
(∀k∈ N)

uk+1 =uk +γ(L¯xk−b)
pk+1 = proxτf (xk−τL∗uk+1)
xk+1 =pk+1−R∗(RR∗)−1(Rpk+1−c)
¯xk+1 =xk+1 +pk+1−xk.
(4.3)
Then, there exist a solution ˆx to Problem 4.1 and an associated multiplier ˆu such that xk→ ˆx and
uk→ ˆu.
4.2 Numerical experiences
In this section we consider some particular instances of Problem 4.1. We consider the case whenf =
∥·∥1∈ Γ0(RN),N = 1000,τ = 0.99
γ∥L∥2 and the relative error in (4.3) isrk =
√
∥uk+1−uk∥
2
+∥xk+1−xk∥
2
∥uk∥
2
+∥xk∥
2 ,
for everyk∈ N. We set γ = 10−2 and (x0, ¯x0,u 0) = (0, 0, 0)∈ RN× RN× Rm+n and, in each test we
show the average execution time and the number of average iterations of both methods, obtained
by considering 20 random realizations of matrices R,S and vectors c∈ Rm andd∈ Rn. Here PCP
and CP denote the algorithms (4.3) and (4.2), respectively.
13

<!-- page 14 -->
m = 1, n = 100 e = 10−4 e = 5· 10−5 e = 10−5
iter time (s) iter time (s) iter time (s)
PCP 9265 22.28 14570 37.02 46191 116.26
CP 9732 23.04 15718 39.21 50544 125.49
%improv. 4.8 3.3 7.3 5.6 8.6 7.4
Table 1: Average time and number of iterations when m = 1 for obtaining rk <e .
T est 1. In Problem 4.1, Table 1 show the eﬃciency of CP and PCP for the case m = 1 and
n = 100. We see that both algorithms are similar in terms of the execution time and the number of
iterations, with a small advantage for the PCP algorithm. In addition, by decreasing the tolerance
e, the percentage of improvement, computed as 100 · (xCP−xPCP)/xCP, slightly increases.
T est 2. In Problem 4.1, Table 2 show the eﬃciency of CP and PCP for the case m = 10 and
n = 100. In this case, there are clear diﬀerences between both algorithms and, as before, PCP is
Table 2: Average time and number of iterations when m = 10 for obtaining rk <e .
m = 10, n = 100 e = 10−4 e = 5· 10−5 e = 10−5
iter time (s) iter time (s) iter time (s)
PCP 6865 18.65 10229 27.86 22855 65.05
CP 9280 23.72 16033 39.13 49526 129.78
%improv. 26.0 21.4 36.2 28.8 53.9 49.9
more eﬃcient as tolerance decreases. In fact, when tolerance is 10 −5, there is an improvement of
approximately 50% with respect to the CP in the execution time and the number of iterations is
less than a half.
T est 3. Finally, in Problem 4.1, Table 3 show the eﬃciency of CP and PCP for the case m = 30
andn = 100. We note that the improvement in execution times are considerably higher than in the
Table 3: Average time and number of iterations when m = 30 for obtaining rk <e .
m = 30, n = 100 e = 10−4 e = 5· 10−5 e = 10−5
iter time (s) iter time (s) iter time (s)
PCP 5146 7.68 7143 10.67 13421 19.70
CP 9941 12.93 16438 21.37 50841 64.23
%improv. 48.2 40.6 56.5 50.1 73.6 69.3
previous cases. For example, in the case of e = 10−4 the improvement increases by approximately
20% with respect to the case m = 10 and by approximately 40% in the case of m = 1. As in the
previous cases, if we decrease the tolerance to 10−5, PCP has better eﬃciency reaching almost 70%
improvement with respect to CP. Table 4 summarizes the percentage of improvements for each test.
We observe that for larger values ofm we obtain a better relative performance of PCP with respect
to CP. The larger is m, the larger is the proportion of constraints on which we project.
14

<!-- page 15 -->
Table 4: Comparison of improvement of average iterations and average times.
% improv. m = 1 m = 10 m = 30
iter time (s) iter time (s) iter time (s)
e = 10−4 4.8 3.3 26.0 21.4 48.2 40.6
e = 5· 10−5 7.3 5.6 36.2 28.8 56.5 50.1
e = 10−5 8.6 7.4 53.9 49.9 73.6 69.3
5 Conclusions
In this paper we provide a projected primal-dual method for solving composite monotone inclusions
with a priori information on solutions. We provide acceleration schemes in the presence of strong
monotonicity and we derive linear convergence in the fully strongly monotone case. The importance
of the a priori information set is illustrated via a numerical example in convex optimization with
equality constraints, in which the proposed method outperforms [3].
Acknowledgements The authors thank the “Programa de ﬁnanciamiento basal” from CMM–
Universidad de Chile and the project DGIP-UTFSM PI-M-18.14 of Universidad T´ ecnica Federico
Santa Mar´ ıa.
References
[1] Condat, L.: A primal-dual splitting method for convex optimization involving lipschitzian, prox-
imable and linear composite terms. J. Optim. Theory Appl. 158,460-479 (2013)
[2] V˜ u, B.C.: A splitting algorithm for dual monotone inclusions involving cocoercive operators.
Adv. Comput. Math. 38, 667-681 (2013)
[3] Chambolle, A., Pock, T.: A ﬁrst-order primal-dual algorithm for convex problems with applica-
tions to imaging. J. Math. Imaging Vis. 40, 120-145 (2011)
[4] Attouch, H., Brice˜ no-Arias, L.M., Combettes, P.L.: A parallel splitting method for coupled
monotone inclusions. SIAM J. Control Optim. 48, 3246-3270 (2010)
[5] Attouch, H., Bolte, J., Redont, P., Soubeyran, A.: Alternating proximal algorithms for weakly
coupled convex minimization problems–applications to dynamical games and PDE’s. J. Convex
Anal. 15, 485-506 (2008)
[6] Attouch, H., Brice˜ no-Arias, L.M., Combettes, P.L.: A strongly convergent primal-dual method
for nonoverlapping domain decomposition, Numer. Math. 133, 433-470 (2016)
[7] Brice˜ no-Arias, L.M., Kalise, D., Silva, F.J.: Proximal methods for stationary Mean Field Games
with local couplings, SIAM J. Control Optim. 56(2), 801-836 (2018)
[8] Facchinei, F., Pang, J.-S.: Finite-Dimensional Variational Inequalities and Complementarity
Problems. Springer, New York (2003)
[9] Gabay, D.: Applications of the method of multipliers to variational inequalities. In: Fortin, M.,
Glowinski, R. (eds.): Augmented Lagrangian Methods: Applications to the Numerical Solutionnof
15

<!-- page 16 -->
Boundary Value Problems, pp. 299–331. North-Holland, Amsterdam (1983)
[10] Mercier, B.: Topics in Finite Element Solution of Elliptic Problems (Lectures on Mathematics,
no. 63). Tata Institute of Fundamental Research, Bombay (1979)
[11] Tseng, P.: A modiﬁed forward-backward splitting method for maximal monotone mappings.
SIAM J. Control Optim. 38, 431-446 (2000)
[12] Brice˜ no-Arias, L.M., Davis, D.: Forward-backward-half forward algorithm for solving monotone
inclusions, SIAM J. Optim., in press (2018)
[13] Esser, E.; Zhang, X.; Chan, T.F.: A general framework for a class of ﬁrst order primal-dual
algorithms for convex optimization in imaging science. SIAM J. Imaging Sci. 3, 1015-1046 (2010)
[14] He, B., Yuan, X.: Convergence analysis of primal-dual algorithms for a saddle-point problem:
from contraction perspective. SIAM J. Imaging Sci. 5, 119-149 (2012)
[15] Chambolle, A.; Pock, T.: On the ergodic convergence rates of a ﬁrst-order primal-dual algo-
rithm. Math. Program. 159, 253–287 (2016)
[16] Lorenz, D.A., Pock, T.: An inertial forward-backward algorithm for monotone inclusions. J.
Math. Imaging Vision 51, 311-325 (2015)
[17] Brice˜ no-Arias, L.M., Kalise, D., Kobeissi, Z., Lauri` ere, M., Gonz´ alez, A.M., Silva, F.J.: On
the implementation of a primal-dual algorithm for second order time-dependent mean ﬁeld games
with local couplings (https://arxiv.org/abs/1802.07902)
[18] Bauschke, H.H., Combettes, P.L.: Convex Analysis and Monotone Operator Theory in Hilbert
Spaces, 2nd ed., Springer, New York (2017)
[19] Brice˜ no-Arias, L.M., Combettes, P.L.: A monotone+skew splitting model for composite mono-
tone inclusions in duality. SIAM J. Optim. 21, 1230-1250 (2011)
[20] Combettes, P. L.; Pesquet, J.-C.: Primal-dual splitting algorithm for solving inclusions with
mixtures of composite, Lipschitzian, and parallel-sum type monotone operators. Set-Valued Var.
Anal. 20, 307-330 (2012)
[21] Tseng, P.: A modiﬁed forward-backward splitting method for maximal monotone mappings.
SIAM J. Control Optim. 38, 431-446 (2000)
[22] Bot ¸, R. I.; Hendrich, C.: Convergence analysis for a primal-dual monotone + skew splitting
algorithm with applications to total variation minimization. J. Math. Imaging Vision 49, 551-568
(2014)
[23] Malitsky, Y.; Pock, T.: A ﬁrst-order primal-dual algorithm with linesearch. SIAM J. Optim.
28, 411-432 (2018)
[24] Combettes, P. L.: Solving monotone inclusions via compositions of nonexpansive averaged
operators. Optimization 53, 475-504 (2004)
16