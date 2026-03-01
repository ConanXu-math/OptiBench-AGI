# A Robust Version of Convex Integral Functionals

**arXiv ID:** 1305.6023v3

**Authors:** Keita Owari

**Abstract:** We study the pointwise supremum of convex integral functionals $\mathcal{I}_{f,γ}(ξ)= \sup_{Q} \left( \int_Ωf(ω,ξ(ω))Q(dω)-γ(Q)\right)$ on $L^\infty(Ω,\mathcal{F},\mathbb{P})$ where $f:Ω\times\mathbb{R}\rightarrow\overline{\mathbb{R}}$ is a proper normal convex integrand, $γ$ is a proper convex function on the set of probability measures absolutely continuous w.r.t. $\mathbb{P}$, and the supremum is taken over all such measures. We give a pair of upper and lower bounds for the conjugate of $\mathcal{I}_{f,γ}$ as direct sums of a common regular part and respective singular parts; they coincide when $\mathrm{dom}(γ)=\{\mathbb{P}\}$ as Rockafellar's result, while both inequalities can generally be strict. We then investigate when the conjugate eliminates the singular measures, which a fortiori yields the equality in bounds, and its relation to other finer regularity properties of the original functional and of the conjugate.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
A Robust V ersion of Convex Integral Functionals
Keita Owari
Graduate School of Economics, The University of Tokyo
7-3-1 Hongo, Bunkyo-ku, Tokyo 113-0033, Japan
owari@e.u-tokyo.ac.jp
First Version: 26.05.2013
Major Revision: 25.06.2014
This version: 21.05.2015
Abstract
We study the pointwise supremum of convex integral functionals
I f,γ(ξ) = sup
Q
(∫
Ω
f (ω, ξ(ω))Q(dω)− γ(Q)
)
, ξ ∈ L∞(Ω, F ,P),
where f : Ω×R→R is a proper normal convex integrand,γ is a proper convex function on the set of probability
measures absolutely continuous w.r.t. P, and the supremum is taken over all such measures. We give a pair
of upper and lower bounds for the conjugate of I f,γ as direct sums of a common regular part and respective
singular parts; they coincide when dom( γ) ={P} as Rockafellar’s classical result, while both inequalities can
generally be strict. We then investigate when the conjugate eliminates the singular measures, which a fortiori
yields the equality in bounds, and its relation to other ﬁner regularity properties of the original functional and of
the conjugate.
Key Words: convex integral functionals, duality, robust stochastic optimization, ﬁnancial mathematics
MSC (2010): 46N10, 46E30, 49N15, 52A41, 91G80
1. Introduction
Let (Ω, F ,P) be a probability space and f : Ω×R→ (−∞,∞] a measurable mapping with
f (ω,·) being proper, convex, lsc for a.e. ω. Then ξ↦→ I f (ξ) :=
∫
Ω f (ω, ξ(ω))P(dω) deﬁnes a
convex functional on L∞ := L∞(Ω, F ,P), called the convex integral functional. Among many
others, R.T. Rockafellar obtained in [22] that under mild integrability assumptions on f , the
conjugate I∗
f (ν) = supξ∈L∞(ν(ξ)− I f (ξ)) of I f is expressed as the direct sum of regular and
singular parts (which we call the Rockafellar theorem):
(1) I∗
f (ν) = I f∗(dνr/dP) + sup
ξ∈dom(I f )
νs(ξ), ∀ν∈ (L∞)∗,
where f∗(ω, y) := supx(xy− f (ω, x)) and νr (resp. νs) denotes the regular (resp. singular) part
of ν∈ (L∞)∗ in the Yosida-Hewitt decomposition. In particular, if I f is ﬁnite everywhere on
L∞, the conjugate I∗
f “eliminates” the singular elements of (L∞)∗ in that I∗
f (ν) =∞ unless ν is
σ-additive. The latter property implies the weak-compactness of all the sublevel sets of I∗
f|L1
and the continuity of I f for the Mackey topology τ(L∞, L1) and so on (e.g. [23, Th. 3K]).
This paper is concerned with a robust version of integral functionals of the form
I f,γ(ξ) := sup
Q∈Q
(∫
Ω
f (ω, ξ(ω))Q(dω)− γ(Q)
)
, ξ ∈ L∞.
where Q is the set of all probability measures Q absolutely continuous w.r.t. P, and γ is a
convex penalty function on Q. (see Section 2 for precise formulation). As the pointwise
1
arXiv:1305.6023v3  [math.FA]  21 May 2015

<!-- page 2 -->
2 K. Owari
supremum of convex functions, I f,γ is convex, lower semicontinuous if so is each I f,Q(ξ) :=∫
Ω f (ω, ξ(ω))Q(dω), and is even norm-continuous as soon as it is ﬁnite everywhere. On the
other hand, it is less obvious what the convex conjugate I∗
f,γ is, when singular measures are
eliminated, and whether the latter property yields ﬁner regularity properties ofI f,γ and I∗
f,γ.
Our main result (Theorem 3.1) is a robust version of the Rockafellar theorem which con-
sists of a pair of upper and lower bounds (instead of a single equality) for I∗
f,γ on (L∞)∗. Both
bounds are of forms analogous to (1) with a common regular part, but a di ﬀerence appears in
the singular parts. They coincide in the classical case I f = I f,δ{P}, while the subsequent Exam-
ple 3.3 shows that both inequalities can be strict and one may not hope for sharper bounds in
general. The same example shows also that the everywhere ﬁniteness ofI f,γ is not enough for
the property that I∗
f,γ eliminates the singular measures, while the lower bound in Theorem 3.1
provides us a simple suﬃcient condition in a form of “uniform integrability”. In Theorem 3.4
and its corollaries, it is shown that given dom( I f,γ) = L∞ (plus a technical assumption), the
condition is even necessary, and is equivalent to other ﬁner properties ofI f,γ and I∗
f,γ, includ-
ing the weak compactness of sublevels of I∗
f,P|L1, the Mackey continuity of I f,γ on L∞ among
others, which are guaranteed solely by the ﬁniteness in the classical case.
Certain class of robust optimization problems are formulated as (or reduced to) the mini-
mization of a robust integral functional I f,γ over a convex set C⊂ L∞. Our initial motivation
of this work lies in the Fenchel duality for this type of problems with⟨L∞, L1⟩ pairing:
(2) inf
ξ∈C
I f,γ(ξ) =− min
η∈L1
(
I∗
f,γ(η) + sup
ξ∈C
⟨−ξ, η⟩
)
,
The Rockafellar-type result provides us the precise form of I∗
f,γ|L1, hence of the dual problem,
while the classical Fenchel duality theorem tells us that a su ﬃcient condition for (2) is the
τ(L∞, L1)-continuity of I f,γ at some point ξ0∈ C∩ dom(I f,γ), or in another view, the ﬁniteness
implies (via the norm continuity) the⟨L∞, (L∞)∗⟩-duality which reduces to the⟨L∞, L1⟩-duality
if I∗
f,γ eliminates the singular elements of (L∞)∗.
Our motivating example of minimization of I f,γ is the robust utility maximization
maximize u(ξ) := inf
Q∈Q
(EQ[U(·, ξ)] + γ(Q)) over a convex cone C⊂ L∞(3)
where U : Ω×R→R is a random utility function. In this context, each Q∈ Q is considered
as a model used to evaluate the quality of control ξ via the expected utility EQ[U(·, ξ)] =∫
Ω U(ω, ξ(ω))Q(dω), but we are not sure which model is true; so we optimize the worst case
among all models Q penalized by γ(Q) according to the likelihood. Note that (3) is equivalent
to minimize I fU ,γ with fU(ω, x) =−U(ω,−x) over the cone−C. The corresponding Fenchel
duality in⟨L∞, L1⟩ constitute a half of what we call the martingale duality in mathematical
ﬁnance; given the⟨L∞, L1⟩-duality with a “good” cone C having the polar generated by so-
called local martingale measures, the theory of (semi)martingales takes care of the other half.
2. Preliminaries
We use the probabilistic notation. Let ( Ω, F ,P) be a complete probability space and L0 :=
L0(Ω, F ,P) denote the space of (equivalence classes modulo P-a.s. equality of) R-valued
random variables deﬁned on it. As usual, we do not distinguish a random variable and the class
it generates, and a constant c∈R is regarded as the random variable c1Ω. Here 1A denotes
the indicator of a set A in measure theoretic sense while δA :=∞1Ac is the one in convex
analysis. The P-expectation of ξ ∈ L0 is denoted by E[ξ] : =
∫
Ω ξ(ω)P(dω) and we write

<!-- page 3 -->
Robust Integral Functionals 3
Lp := Lp(Ω, F ,P) and∥·∥ p :=∥·∥ Lp for each 1≤ p≤∞ . Any probabilistic notation without
reference to a measure is to be understood with respect toP. Especially, “a.s.” means “ P-a.s.”,
and identiﬁcation of random variables is always made byP. For other probability measures Q
absolutely continuous with respect toP (Q≪P), we writeEQ[·] for Q-expectation, Lp(Q) :=
{ξ∈ L0 : EQ[|ξ|p] <∞} (which is a set of P-classes!) etc, explicitly indicating the measure
that involves. We write Q∼ P to mean Q and P are equivalent (Q≪ P and P≪ Q).
The norm dual of L∞ is ba := ba(Ω, F ,P), the space of allbounded ﬁnitely additive signed
measures ν respectingP-null sets, i.e., supA∈F|ν(A)| <∞, ν(A∪ B) = ν(A) + ν(B) if A∩ B =∅,
and ν(A) = 0 ifP(A) = 0 (see [12, pp. 354-357]). The bilinear form of (L∞, ba) is given by the
(Radon) integral ν(ξ) =
∫
Ω ξdν which coincides with the usual integral when ν is σ-additive.
We regard any σ-additive ν∈ ba as an element of L1 via the mapping ν↦→ dν/dP which is an
isometry from the subspace of such ν’s onto L1. In particular, the set
Q :={Q : probability measures on ( Ω, F) with Q≪P}
is regarded as{η∈ L1 : η≥ 0,E[η] = 1}. A ν∈ ba is said to be purely ﬁnitely additive if
there exists a sequence (An) in F such thatP(An)↗ 1 but|ν|(An) = 0 for all n, and we denote
by bas the totality of such ν∈ ba. Any ν∈ ba admits a unique Yosida-Hewitt decomposition
ν = νr + νs where νr is the regular (σ-additive) part, and νs is the purely ﬁnitely additive part
(e.g. [7, Th. III.7.8]), thus ba = L1⊕ bas with the above identiﬁcation. We denote byba+ the
set of positive elements of ba and bas
+ := ba+∩ bas etc.
2.1. Convex Penalty Function and Associated Orlicz Spaces
We make the following assumption on the penalty function γ:
Assumption 2.1. γ : Q→R+ is a σ(L1, L∞)-lsc proper convex function such that
inf
Q∈Q
γ(Q) = 0; (4)
∃Q0∈ Q with Q0∼P and γ(Q0) <∞; (5)
{Q∈ Q : γ(Q)≤ 1} is σ(L1, L∞)-compact(6)
(⇔ closed and uniformly integrable).
We denote Qγ :={Q∈ Q : γ(Q) <∞}, the eﬀective domain of γ.
Remark 2.2. (4) and (5) are normalizing assumptions and the latter one is equivalent to saying
that Qγ∼P, i.e., for each A∈ F,P(A) = 0 iﬀ Q(A) = 0 for every Q∈ Qγ. Given the lower
semicontinuity of γ, (6) implies that the inﬁmum in (4) is attained, and (6) is equivalent to
apparently stronger
{Q∈ Q : γ(Q)≤ c} is σ(L1, L∞)-compact,∀c > 0.
Indeed, for any c > 1, pick a Qc ∈ Qγ with γ(Qc) < 1/c (by (4)), then γ(Q) ≤ c im-
plies γ
(
1
1+cQ + c
1+cQc
)
≤ 1
1+c· c + c
1+c· 1
c = 1, thus (6) implies the weak compactness of{
1
1+cQ + c
1+cQc : γ(Q)≤ c
}
, hence of{Q∈ Qγ : γ(Q)≤ c}. ♦
In general, any lower semicontinuous proper convex function γ on Q deﬁnes a function
ργ(ξ) := sup
Q∈Qγ
(EQ[ξ]− γ(Q)) on
{
ξ∈ L0 : ξ−∈ ⋂
Q∈Qγ L1(Q)
}
⊃ L∞∪ L0
+.

<!-- page 4 -->
4 K. Owari
Regardless of Assumption 2.1, ργ restricted to L∞ is a σ(L∞, L1)-lsc ﬁnite-valued monotone
convex function with ργ(ξ + c) = ργ(ξ) + c if c is a constant, whose conjugate on L1 is
γ(η)1Q(η) +∞1L1\Q(η). In ﬁnancial mathematics, such a function is called a convex risk
measure (up to a change of sign). Then (4) reads as ργ(0) = 0, (5) as ργ(ε1A) = 0 for some
ε > 0⇒P(A) = 0, while (6) is equivalent to saying that ργ has the Lebesgue property on L∞
(see [13, 5, 17]):
sup
n
∥ξn∥∞ <∞, ξ n→ ξ a.s. ⇒ ργ(ξ) = lim
n
ργ(ξn).
To a penalty function γ satisfying Assumption 2.1, we associate the gauge norm
∥ξ∥ργ = inf
{
λ > 0 : ργ(|ξ|/λ)≤ 1
}
= sup
Q∈Qγ
EQ[|ξ|]
1 + γ(Q) , ξ ∈ L0.
In view of (5), this is indeed a norm on the Orlicz space
Lργ =
{
ξ∈ L0 :∃α > 0 with ργ(α|ξ|) <∞
}
={ξ∈ L0 :∥ξ∥ργ <∞},
which is a solid subspace (lattice ideal) of L0 (i.e., ξ∈ Lργ and|ζ|≤| ξ| a.s.⇒ ζ∈ Lργ), and
(Lργ,∥·∥ ργ) is a Banach lattice. We consider the following subspaces of Lργ too:
Mργ :=
{
ξ∈ L0 :∀α > 0, ρ γ(α|ξ|) <∞
}
,
M
ργ
u :=
{
ξ∈ L0 :∀α > 0, lim
N
ργ
(
α|ξ|1{|ξ|>N}
)
= 0
}
.
Both Mργ and M
ργ
u are solid as subspaces of L0. These Orlicz-type spaces are studied in [19]
where γ corresponds to ϕ∗
0 (more precisely ϕ∗
0 = γ1Q +∞1L1\Q) and ˆϕ = ργ in the notation of
the current paper. In general under Assumption 2.1,
L∞⊂ M
ργ
u ⊂ Mργ⊂ Lργ⊂ ⋂
Q∈Qγ L1(Q),
while all inclusions can generally be strict ([19, Examples 3.6 and 3.7]). From the last inclu-
sion, ργ is well-deﬁned on Lργ as a proper monotone convex function, and it is ﬁnite-valued
on Mργ while not on Lργ (unless Mργ = Lργ), and M
ργ
u is the maximum solid subspace of L0 to
which ργ retain the Lebesgue property (as the order-continuity; see [19, Theorem 3.5]), and is
characterized by a uniform integrability property [19, Theorem 3.8]: for ξ∈ Mργ,
(7) ξ∈ M
ργ
u ⇔{ ξdQ/dP : γ(Q)≤ c} is uniformly integrable,∀c > 0.
Actually, by the same argument as in Remark 2.2, we have only to consider the case c = 1.
2.2. Robust V ersion of Integral Functionals
In the sequel, let f : Ω×R→R be a proper normal convex integrandonR, i.e.,
f is F⊗ B(R)-measurable and
f (ω,·) is an lsc proper convex function for a.e. ω.(8)
Since F is assumedP-complete, the ﬁrst part is equivalent to the measurability of epigraphical
mapping (see [24, Ch. 14] for a general reference). As immediate consequences of (8), f (·, ξ)
is F-measurable for each ξ∈ L0, and f∗ is also a proper normal convex integrand where
f∗(ω, y) := sup
x∈R
(xy− f (ω, x)), ω ∈ Ω, y∈R,

<!-- page 5 -->
Robust Integral Functionals 5
Given such f and γ, we deﬁne a robust analogue of convex integral functional
(9) I f,γ(ξ) := sup
Q∈Qγ
(EQ[ f (·, ξ)]− γ(Q))= ργ ( f (·, ξ)) , ξ ∈ L∞.
This functional is well-deﬁned as a proper convex functional onL∞ as soon as:
∃ξ0∈ L∞ s.t. f (·, ξ0)+∈ Mργ,(10)
∃η0∈ Lργ s.t. f∗(·, η0)+∈ Lργ.(11)
Indeed, since f (·, ξ)≥ ξη0− f∗(·, η0)+, (11) implies f (·, ξ)−∈ Lργ⊂ ⋂
Q∈Qγ L1(Q), soEQ[ f (·, ξ)]
is well-deﬁned with values in (−∞,∞] for each Q∈ Qγ and ξ ∈ L∞, thus so is I f,γ(ξ) =
supQ∈Qγ
(EQ[ f (·, ξ)]− γ(Q)), while (10) shows that I f,γ(ξ0) <∞, thus I f,γ .∞. Also, I f,γ is
convex as the pointwise supremum of convex functions. Note also that
{ξ∈ L∞ : f (·, ξ)+∈ Mργ}⊂ dom(I f,γ)⊂{ ξ∈ L∞ : f (·, ξ)+∈ Lργ}.(12)
where dom(I f,γ) :={ξ∈ L∞ : I f,γ(ξ) <∞}. In particular, by (5),
(13) ξ∈ dom(I f,γ)⇒ ξ(ω)∈ dom f (ω,·) for a.e. ω.
If Mργ = Lργ, all three sets in (12) coincide, so ργ ( f (·, ξ)+) <∞ as soon as I f,γ(ξ) <∞. In
general, however, ργ ( f (·, ξ)+) =∞ may happen even if ξ∈ dom(I f,γ).
We next check that I f,γ has a nice regularity on L∞.
Lemma 2.3. Under Assumption 2.1, (10) and (11),I f,γ is σ(L∞, L1)-lower semicontinuous, or
equivalently I f,γ has the Fatou property:
(14) sup
n
∥ξn∥∞ <∞, ξ n→ ξ a.s. ⇒ I f,γ(ξ)≤ lim inf
n
I f (ξn).
Proof. Suppose a := supn∥ξn∥∞ <∞ and ξn→ ξ a.s., then automatically ξ∈ L∞. With η0 as
in (11), f (·, ξn)≥− a|η0|− f∗(·, η0)+∈ Lργ⊂ ⋂
Q∈Qγ L1(Q). Hence for each Q∈ Qγ, Fatou’s
lemma shows thatEQ[ f (·, ξ)]≤ lim infnEQ[ f (·, ξn)], and we deduce (14) as
sup
Q∈Qγ
(EQ[ f (·, ξ)]− γ(Q))≤ sup
Q∈Qγ
(
lim inf
n
EQ[ f (·, ξn)]− γ(Q)
)
≤ lim inf
n
sup
Q∈Qγ
(EQ[ f (·, ξn)]− γ(Q)).
The equivalence between (14) and the weak*-lower semicontinuity follows from a well-known
consequence of the Krein-Šmulian and Mackey-Arens theorems that (see e.g. [11]) a convex
set C⊂ L∞ is weak*-closed if and only if C∩{ ξ :∥ξ∥∞≤ a} is L0-closed for all a > 0.
2.3. Robust f ∗-divergence
We proceed to the functional that plays the role of I f∗ in the classical case. Let
˜f∗(ω, y, z) := sup
x∈dom f (ω,·)
(xy− z f(ω, x)), ω ∈ Ω, (y, z)∈R×R+.
Noting that (a−
f∗, a+
f∗)⊂ dom f⊂ [a−
f∗, a+
f∗] where a±
f∗ := limk→± f∗(·, k)/k, we have
˜f∗(ω, y, z) =

0 if y = z = 0,
y· a±
f∗(ω) if y≷ 0, z = 0,
z f∗(ω, y/z) if z > 0.
(15)

<!-- page 6 -->
6 K. Owari
Lemma 2.4. Suppose (8). Then ˜f∗ : Ω×R×R+→R is a proper normal convex integrand on
R×R+, i.e., it is F⊗ B(R)⊗ B(R+)-measurable and (y, z)↦→ ˜f∗(ω, y, z) is a lower semicon-
tinuous proper convex function for a.e. ω∈ Ω. Also, for a.e. ω∈ Ω,
(16) xy≤ z f(ω, x) + ˜f∗(ω, y, z), ∀x∈ dom f (ω,·),∀y∈R,∀z≥ 0.
Proof. Since f is normal, there exists a sequence of measurable functions ( ξn)n∈N⊂ L0 such
that{ξn(ω)}n∩ dom f (ω,·) is dense in dom f (ω,·) ([23], Proposition 2D). Modifying the se-
quence as ¯ξn := ξn1{ f (·,ξn)<∞} + ξ01{ f (·,ξn)=∞} where ξ0∈ domI f,γ, we have for a.e. ω, ¯ξn(ω)∈
dom f (ω,·) and{¯ξn(ω)}n is dense in dom f (ω,·). Thus
˜f∗(ω, y, z) = sup
n
(¯ξn(ω)y− z f(·, ¯ξn(ω))).
Consequently, ˜f∗ is a normal convex integrand as the countable supremum of aﬃne integrands
with ˜f∗(·, 0, 0) = 0. (16) is obvious from the deﬁnition.
Now we deﬁneH f∗(η|Q) :=E[ ˜f∗ (·, η, dQ/dP)] for η∈ L1, Q∈ Qγ, and
H f∗,γ(η) := inf
Q∈Qγ
(
H f∗(η|Q) + γ(Q)
)
, η ∈ L1.(17)
In view of identiﬁcation{ν∈ ba : σ-additive} = L1, we deﬁne also for anyσ-additive ν∈ ba,
H f∗(ν|Q) := H f∗ (dν/dP|Q) , H f∗,γ(ν) = H f∗,γ(dν/dP).
Since ˜f∗(·, y, 1) = f∗(·, y), we recover H f∗(η|P) = H f∗,δ{P}(η) = I f∗(η) in the classical case.
Under the above assumptions, H f∗(·|·) and H f∗,γ are well-deﬁned.
Lemma 2.5. Under Assumption 2.1, (8), (10) and (11), H f∗(·|·) and H f∗,γ are well-deﬁned as
proper convex functionals respectively on L1× Qγ and L1, and it holds
(18) E[ηξ]≤ I f,γ(ξ) + H f∗,γ(η), ∀ξ∈ L∞,∀η∈ L1.
Proof. Let ψQ := dQ/dP for each Q∈ Qγ. In view of (13), we see that
˜f∗(·, η, ψQ)≥ ξη− ψQ f (·, ξ)≥ ξη− ψQ f (·, ξ)+∈ L1,(19)
for any ξ∈ dom(I f,γ) (,∅), η∈ L1 and Q∈ Qγ, thus H f∗(η|Q) =E[ ˜f∗(·, η, ψQ)] >−∞ is well-
deﬁned, convex onL1× Qγ (since ˜f∗ is convex), and proper sinceH f∗(η0ψQ|Q) =EQ[ f∗(·, η0)]
with η0∈ Lργ as in (11) (then η0ψQ∈ L1). Taking the expectation in (19),
H(η|Q) + γ(Q)≥E[ξη]− (EQ[ f (·, ξ)]− γ(Q))≥E[ξη]− I f,γ(ξ) >−∞,
for any ξ∈ dom(I f,γ), η∈ L1 and Q∈ Qγ, so H f∗,γ,(η) = infQ∈Qγ
(
H f∗(η|Q) + γ(Q)
)
>−∞
and we have (18) (which is trivially true when I f,γ(ξ) =∞). The convexity of H f∗,γ follows
from that of H f∗(·|·) + γ(·) and of Qγ.
Remark 2.6. Under assumption (22) in the (main) Theorem 3.1 below, H f∗ + γ (resp. H f∗,γ)
is weakly lower semicontinuous on L1× Qγ (resp. L1), and for each η∈ L1, the inﬁmum
infQ∈Qγ
(
H f∗(η|Q) + γ(Q)
)
is attained (see Appendix A). Though we could give direct proofs

<!-- page 7 -->
Robust Integral Functionals 7
here, we will not use these properties, and the lower semicontinuity of H f∗,γ will be obtained
as Corollary 3.2 to the main theorem which does not internally use that property.
Note also that when f (hence f∗ too) is non-random, all the integrability assumptions are
trivialized, in which case H f∗(·|·) is called the f∗-divergence while H f∗,γ is a slight general-
ization of robust f∗-divergence (the latter is a special case with γ(Q) = δP(Q) with P⊂ Q),
and the joint lower semicontinuity etc are found e.g. in [9, Lemma 2.7]. In fact, if f is ﬁnite
(P( f (·, x) <∞,∀x) = 1⇔ lim|y|→∞ f∗(·, y)/|y| =∞ a.s.), we have from (15) that
H f∗(ν|Q) =

EQ[ f∗(·, dν/dQ)] if ν≪ Q,
+∞ otherwise. ♦
Here are some typical examples of penalty function γ and associated integral functionals.
Example 2.7 (Classical case). The “classical” integral functional I f,P(ξ) : = E[ f (·, ξ)] corre-
sponds to the penalty function γP(Q) := δ{P}(Q) which clearly satisﬁes Assumption 2.1, and
ργP(ξ) =E[ξ]. Then M
ργ
u = Mργ = Lργ = L1, and the integrability assumptions (10) and (11)
are identical to the ones in [22]:
∃ξ0∈ L∞ with f (·, ξ0)+∈ L1 and ∃η0∈ L1 with f∗(·, η0)+∈ L1,
under which (1) holds true ([22], Th. 1), and since I f∗,P(η) = H f∗(η|P), it reads as
I∗
f,P(ν) = H f∗(νr|P) + sup
ξ∈dom(I f,P)
νs(ξ), ∀ν∈ ba. ♦
We can also consider other probability P≪P and I f,P := I f,γP where γP = δ{P}, but we need
a little care when P≁P (then (5) is violated): here we consider I f,P as a functional on L∞(P)
rather than the space of P-equivalence classes of P-essentially bounded random variables.
Equivalently, I f,P = Ig,P with g(·, x) = f (·, x) dP
dP 1{dP/dP>0}, and its conjugate is
(20) I∗
f,P(ν) = H f∗(νr|P) +∞1{νr3P} + sup
ξ∈dom(I f,P)
νs(ξ).
Example 2.8 (Homogeneous case). The formulation (9) covers the following form
I f,P(ξ) := sup
Q∈P
EQ[ f (·, ξ)] = I f,δP (ξ),
where P⊂ Q is a nonempty convex set. δP clearly satisﬁes (4), and (6) (resp. (5)) is equiva-
lent to the weak compactness of P itself (resp.∃Q∈ P with Q∼P), and ρP(ξ) := ρδP (ξ) =
supQ∈PEQ[ξ] is a positively homogeneous monotone convex function, called sublinear expec-
tation or coherent risk measure (modulo change of sign). In particular, Mργ = Lργ, while
M
ργ
u ⊊ Mργ is possible (see [19, Example 3.7]). ♦
Example 2.9 (Polyhedral case). This is a special case of Example 2.8. Suppose we are given
a ﬁnite number of probability measuresP1, ..., Pn≪P which generate a polyhedral convex set
P = conv(P1, ..., Pn). This P is clearly (convex and) weakly compact in L1, (5) is equivalent
to 1
n(P1, ..., Pn)∼P, and we have MρP
u = MρP = LρP = ⋂
k≤n L1(Pk) since
1
n
∑
k≤n
EPk[|ξ|]≤∥ ξ∥ρP = max
1≤k≤n
EPk[|ξ|]≤
∑
k≤n
EPk[|ξ|].

<!-- page 8 -->
8 K. Owari
In particular, noting that I f,(λ1P1+···+λnPn) = λ1I f,P1 +··· + λI f,Pn,
I f,P(ξ) = sup
Q∈P
I f,Q(ξ) = max
1≤k≤n
I f,Pk(ξ), dom(I f,P) = ⋂
k≤n dom(I f,Pk). ♦
In [28, Cor. 2.8.11], the conjugate of pointwise maximum ofﬁnitely manyconvex functions is
obtained, which reads in our context as (compare to (25) below): if dom(I f,P) ,∅,
I∗
f,P(ν) = min
{
ϕ∗(ν) : ϕ∈ conv(I f,P1, ..., I f,Pn)
}
= min
Q∈P
I∗
f,Q(ν)
(20)
= min
Q∈P
H f∗(νr|Q) +∞1{νr3Q} + sup
ξ∈dom(I f,Q)
νs(ξ)
 .
(21)
Example 2.10 (Entropic penalty). Let γent(Q) be the relative entropy of Q w.r.t.P:
γent(Q) := Hx log x(Q|P) :=E
[dQ
dP log dQ
dP
]
.
This function satisﬁes Assumption 2.1: γent(P) = 0, hence (4) and (5), while (6) follows from
the de la Vallée-Poussin theorem since limx→∞
x log x
x = limx→∞ log x =∞. Let
LΦexp :={ξ∈ L0 :∃α > 0,E[exp(α|ξ|)] <∞} (exponential Orlicz space),
MΦexp :={ξ∈ L0 :∀α > 0,E[exp(α|ξ|)] <∞} (Morse subspace).
In this case, ρent(ξ) : = ργent(ξ) = logE [exp(ξ)] whenever ξ−∈ LΦexp. Hence, M
ργ
u = Mργ =
MΦexp ⊂ LΦexp = Lργ and MΦexp ⊊ LΦexp if ( Ω, F ,P) is atomless (e.g. exponential random
variable). The integrability assumptions (10) for ξ0 ∈ L∞ and (11) for η0 ∈ LΦexp read as
E [exp(α f (·, ξ0))] <∞ for all α > 0 andE [exp(ε f∗(·, η0))] <∞ for some ε > 0, respectively.
Moreover, the corresponding integral functional is explicitly written as
I f,γent(ξ) = logE [exp ( f (·, ξ))] , ∀ξ∈ L∞. ♦
3. Statements of Main Results
3.1. A Rockafellar-Type Theorem for the Convex Conjugate
The following is a robust analogue of the Rockafellar theorem for the conjugate of I f,γ
I∗
f,γ(ν) := sup
ξ∈L∞
(ν(ξ)− I f,γ(ξ)), ν ∈ ba = ba(Ω, F ,P).
Theorem 3.1. Suppose Assumption 2.1, (8), (11) and
∃ξ0∈ L∞ such that f (·, ξ0)+∈ M
ργ
u .(22)
Then for any ν∈ ba with the Yosida-Hewitt decomposition ν = νr + νs,
(23) H f∗,γ(νr) + sup
ξ∈D f,γ
νs(ξ)≤ I∗
f,γ(ν)≤ H f∗,γ(νr) + sup
ξ∈dom(I f,γ)
νs(ξ),
where D f,γ :={ξ∈ L∞ : f (·, ξ)+∈ M
ργ
u}⊂ dom(I f,γ) (by (12)).

<!-- page 9 -->
Robust Integral Functionals 9
A proof is given in Section 4.1. In contrast to the classical Rockafellar theorem (1), our
robust version (23) consists of two inequalities instead of a single equality. But the possible
diﬀerence appears only in the singular part, thus
Corollary 3.2 (Restriction to L1). Under the same assumptions as in Theorem 3.1,
I∗
f,γ(η) = H f∗,γ(η) = inf
Q∈Qγ
(
H f∗(η|Q) + γ(Q)
)
, ∀η∈ L1.
In particular, H f∗,γ is weakly lower semicontinuous on L1, and
(24) I f (ξ) = sup
η∈L1
(
E[ξη]− H f∗,γ(η)
)
, ξ ∈ L∞.
Proof. The ﬁrst assertion is clear from (23), by which the conjugate H f∗,γ is lower semicon-
tinuous for any topology consistent with the duality ⟨L∞, L1⟩, while (24) is a consequence of
σ(L∞, L1)-lower semicontinuity of I f,γ (Lemma 2.3) via the Fenchel-Moreau theorem.
In the classical case of Example 2.7, D f,γ = dom(I f,γ) ={ξ∈ L∞ : f (·, ξ)+∈ L1} (since
M
ργ
u = Lργ = L1), hence (23) reduces to a single equality which is exactly (1) as in the Rock-
afellar theorem [22, Theorem 1]. The original version of [22] is slightly more general, where
the integral functional is deﬁned with respect to aσ-ﬁnite(rather than probability) measure µ
forRd-valued random variables ξ∈ L∞(Ω, F , µ;Rd). There are also some extensions replac-
ing L∞(Ω, F , µ;Rd) by some decomposable spaces of measurable functions taking values in a
Banach space. See in this line [14], [4], and [23] for a general reference.
In the polyhedral case of Example 2.9 ( γ = δP and P = conv(P1, ..., Pn)), we still have
D f,γ = dom(I f,γ) = ⋂
k≤n dom(I f,Pk). Thus (23) reduces to
I∗
f,γ(ν) = min
Q∈P
H f∗(νr|Q) + sup
{
νs(ξ) : ξ∈ ⋂
k≤n dom(I f,Pk)
}
.(25)
This is slightly sharper than (21) in the sense that regular and singular parts are separated.
To the best of our knowledge, Rockafellar-type result for the robust form (9) of convex
integral functionals (including the homogeneous case of Example 2.8) is new. A possible
complaint would be the di ﬀerence between singular parts in the upper and lower bounds in
(23). In the full generality of Theorem 3.1, however, both inequalities can really be strict
and one can not hope for sharper bounds as the next example illustrates (see Appendix B for
details).
Example 3.3 (Badly Behaving Integrand). Let (Ω, F) := (N, 2N) withP given byP({n}) =
2−n, and ( Pn)n a sequence of probability measures on 2 N speciﬁed by P1({1}) = 1; Pn({1}) =
1− 1/n, Pn({n}) = 1/n. Then P = conv(Pn; n∈N) is weakly compact in L1(N, 2N,P), thus
γ = δP is a penalty function satisfying Assumption 2.1 and ργ(ξ) = supnEPn[ξ] if ξ∈ L0
+.
In this case, L∞ is regarded as the sequence space ℓ∞ with the norm∥ξ∥∞ = supn|ξ(n)|, and
ν∈ bas
+(N, 2N,P) if and only if ν vanishes on any ﬁnite set, or equivalently, for anyν∈ ba+,
(26) ν∈ bas
+ ⇔∥ ν∥· lim inf
n
ξ(n)≤ ν(ξ)≤∥ ν∥· lim sup
n
ξ(n), ∀ξ∈ ℓ∞ = L∞.
(Such ν , 0 exists, thus bas
+\{ 0} ,∅; see [1, Lemmas 16.29 and 16.30]). Now we set
(27) f (n, x) = nx+ex, n∈N = Ω, x∈R.

<!-- page 10 -->
10 K. Owari
Then I f,γ(ξ) = supn
((
1− 1
n
)
ξ(1)+eξ(1) + ξ(n)+eξ(n)+)
≤ 2∥ξ∥∞e∥ξ∥∞, so dom( I f,γ) = L∞, and
limN→∞ supnEPn[ f (·, ξ)1{ f (·,ξ)≥N}] = lim supn ξ(n)+eξ(n)+
(Lemma B.1), thus
(28) 0 ∈ D f,γ ={ξ∈ ℓ∞ : lim sup
n
ξ(n)≤ 0}⊊ dom(I f,γ) = L∞.
As for I∗
f,γ, dom( I∗
f,γ) ⊂ ba+ since I f,γ is increasing, and H f∗,γ(0) = 0 since f∗(·, 0) =
inf x f (·, x) = 0, thus (23) reads as sup ξ∈D f,γ ν(ξ) ≤ I∗
f,γ(ν) ≤ supξ∈dom(I f,γ) ν(ξ) on bas
+. On
the other hand, for ν∈ bas
+, supξ∈D f,γ ν(ξ) = 0, supξ∈domI f ν(ξ) = +∞, and (Lemma B.2):
I∗
f,γ(ν) = sup
x≥0
x(∥νs∥− ex),∀ν∈ bas
+.
In particular, I∗
f,γ(ν) = 0 if ν∈ U s
+ :={ν∈ bas
+ : ∥ν∥ = ν(N)≤ 1}, 0 < I∗
f,γ(ν) <∞ if
ν∈ bas
+\ U s
+ :={ν∈ bas
+ :∥ν∥ > 1}, and lim∥ν∥→∞,ν∈bas
+\U s
+ I∗
f,γ(ν) =∞. In summary,
– I∗
f,γ coincides with the lower bound supξ∈D f,γ ν(ξ) = 0 on U s
+, while
– on bas
+\ U s
+, I∗
f,γ is strictly between the upper and lower bounds and it runs through the
whole interval of these bounds (in this speciﬁc case, [0,∞]). ♦
3.2. Finer Properties in the Finite-V alued Case
We now consider the regularities of I f,γ and H f∗,γ in terms of the dual paring⟨L∞, L1⟩. In the
classical case of Example 2.7, the singular part of I∗
f in (1) is trivialized (i.e., δ{0}) as soon as
I f := I f,{P} is ﬁnite-valued, then I∗
f reduces entirely to I f∗. It implies that all the sublevels of
I f∗ are σ(L1, L∞)-compact (see [23, Th. 3K]), which is equivalent to the continuity of I f for
the Mackey topology τ(L∞, L1), and I f admits a σ-additive subgradient at every point (weak*
subdiﬀerentiable). Consequently, we can work entirely with the dual pair⟨L∞, L1⟩.
In the robust case, the “triviality of singular part of I∗
f,γ” should be understood as
(29) ∀ν∈ ba, I∗
f,γ(ν) <∞⇒ ν is σ-additive,
i.e. that I∗
f,γ eliminates the singular measures, which still makes sense even though I∗
f,γ itself
need not be the direct sum of regular and singular parts. This guarantees in particular that (23)
reduces to a single equality (of course). In the case of Example 3.3, D f,γ⊊ dom(I f,γ) = L∞,
and I∗
f,γ(νs) <∞ as long as νs≥ 0. Thus dom( I f,γ) = L∞ is not enough for (29), while from
(23), D f,γ = L∞ is clearly suﬃcient. In fact, given the ﬁniteness (plus a technical assumption),
D f,γ = L∞ is also necessary for (29), and equivalent to other basic⟨L∞, L1⟩-regularities of I f,γ
and I∗
f,γ which follow solely from the ﬁniteness in the classical case.
Theorem 3.4. In addition to the assumptions of Theorem 3.1, suppose dom(I f,γ) = L∞ and
(30) ∃ξ′
0∈ L∞ with f (·, ξ′
0)−∈ Mργ.
Then the following are equivalent:
(i) D f,γ = L∞, i.e., f (·, ξ)+∈ M
ργ
u for all ξ∈ L∞;
(ii)R⊂ D f,γ, i.e., f (·, x)+∈ M
ργ
u for all x∈R;
(iii) I∗
f,γ eliminates singular measures in the sense of (29);
(iv){η∈ L1 : H f∗,γ(η)≤ c} is σ(L1, L∞)-compact for all c∈R;
(v) I f,γ is continuous for the Mackey topology τ(L∞, L1);
(vi) I f,γ(ξ) = limn I f,γ(ξn) if supn∥ξn∥∞ <∞ and ξn→ ξ a.s. (the Lebesgue property);

<!-- page 11 -->
Robust Integral Functionals 11
(vii) I f,γ(ξ) = maxη∈L1
(
E[ξη]− H f∗,γ(η)
)
, i.e., the supremum is attained in (24),∀ξ∈ L∞.
Here implications (i)⇔ (ii)⇒ (iii)⇒ (iv)⇔ (v)⇒ (vi) and (v)⇒ (vii) are true without (30).
A proof will be given in Section 4.2.
Remark 3.5. The ﬁniteness ofI f,γ already implies that f (·, ξ)+∈ Lργ,∀ξ∈ L∞ (in particular f
is ﬁnite-valued), while the additional assumption (30) is made to guarantee f (·, ξ)+∈ Mργ for
all ξ∈ L∞ (see (49) and the subsequent paragraph). These coincide in the homogeneous case
(Example 2.8; including the classical case) since then Lργ = Mργ, but in general, Mργ⊊ Lργ is
possible. A suﬃcient condition for (30) is that
(31) ∃η0∈ Mργ with f∗(·, η0)+∈ Mργ,
(then f (·, ξ)− ∈ Mργ for all ξ ∈ L∞), which is identical to (11) (contained in the standing
assumptions) in the homogeneous case. ♦
Remark 3.6. The equivalence between (iv), (vi) and (vii) for convex risk measures on L∞,
i.e., for ργ|L∞ with penalty function γ as in Assumption 2.1 is known ([13] and [5]) followed
by some generalizations: [16, 17] for convex risk measures on Orlicz spaces, and [19, 20]
for ﬁnite-valuedmonotone convex functions on solid spaces of measurable functions among
others. ♦
From (ii)⇒ (v), we derive a simple criterion in terms of integrability off for the⟨L∞, L1⟩-
Fenchel duality for the minimization of robust integral functionalI f,γ. Here we recall Fenchel’s
duality theorem (see [21, Th. 1]): if⟨E, E′⟩ is a dual pair, ϕ, ψ are proper convex functions on
E, and if either ϕ or ψ is τ(E, E′)-continuous at some x∈ dom(ϕ)∩ dom(ψ), then
inf
x∈E
(ϕ(x) + ψ(x)) =− min
x′∈E′
(ϕ∗(x′) + ψ∗(−x′)).
Putting E = L∞, E′ = L1, ϕ = I f,γ and ψ = δC with C⊂ L∞ convex, (ii)⇒ (v) tells us that
Corollary 3.7 (Fenchel Duality). Let γ be a penalty function satisfying Assumption 2.1 and
f a proper normal convex integrand. If f (·, x)+∈ M
ργ
u for all x∈R, and f ∗(·, η0)+∈ Lργ for
some η0∈ Lργ, then for any convex set C⊂ L∞,
inf
ξ∈C
I f,γ(ξ) =− min
η∈L1
(
H f∗,γ(−η) + sup
ξ′∈C
E[ξ′η]
)
.
If in addition C is a convex cone, the right hand side is equal to − minη∈C◦ H f∗,γ(−η), where
C◦ ={η∈ L1 : E[ξη]≤ 1,∀ξ∈ C} (the one-sided polar of C in⟨L∞, L1⟩).
The subdiﬀerential of I f,γ at ξ∈ L∞ is the following set of ν∈ (L∞)∗ called subgradients:
∂I f,γ(ξ) :={ν∈ (L∞)∗ : ν(ξ)− I f,γ(ξ)≥ ν(ξ′)− I f,γ(ξ′),∀ξ′∈ L∞}.
We say thatI f,γ is subdiﬀerentiable at ξ if ∂I f,γ(ξ) ,∅. In view of (24), η∈ ∂I f,γ(ξ)∩ L1 (then
η is called a σ-additive subgradient of I f,γ at ξ) if and only if it maximizes η′↦→ E[ξη′]−
H f∗,γ(η′), thus (vii) is equivalent to saying that for every ξ∈ L∞, ∂I f,γ(ξ)∩ L1 ,∅. Note
also that ∂I f,γ(ξ)⊂ dom(I∗
f,γ) since I∗
f,γ(ν) = supξ′∈L∞(ν(ξ′)− I f,γ(ξ′))≤ ν(ξ)− I f,γ(ξ) <∞ if
ν∈ ∂I f,γ(ξ). Thus (iii) implies ∂I f,γ(ξ)⊂ L1. Summing up,

<!-- page 12 -->
12 K. Owari
Corollary 3.8. Under the assumptions of Theorem 3.4, (i) – (vii) are equivalent also to
(viii)∅ , ∂I f,γ(ξ)⊂ L1 for every ξ∈ L∞.
The weak compactness of the sublevels of H f∗,γ can be viewed as a generalization of the
de la Vallée-Poussin theorem which asserts that a set C ⊂ L1 is uniformly integrable if and
only if there exists a function g :R→ (−∞,∞] which is coercive: lim|y|→∞ g(y)/y =∞ and
supη∈CE[g(|η|)] = supη∈C Hg,δ{P}(|η|) < ∞ (e.g. [6, Th. II.22]). The coercivity condition is
equivalent to saying that dom(g∗) =R. Now we have as a consequence of (ii)⇔ (iv):
Corollary 3.9 (cf. [9] when g is non-random). A set C ⊂ L1 is uniformly integrable if and
only if there exists a convex penalty function γ on Q satisfying Assumption 2.1 as well as a
proper normal convex integrand g with g∗(·, x)∈ M
ργ
u ,∀x∈R, such that supη∈C Hg,γ(η) <∞.
Proof. Let f = g∗, then f∗ = g∗∗ = g (since normal). Then D f,γ = L∞ by assumption and
Hg,γ = H f∗,γ is well-deﬁned while sup η∈C Hg,γ(η) < ∞ guarantees (11) as well. Now the
suﬃciency is nothing but (ii)⇒ (iv), while the necessity is clear from the above paragraph.
3.3. Examples of “Nice” Integrands and Robust Utility Maximization
When f is non-random and ﬁnite,R ⊂ D f,γ is automatic, while f∗(y) ∈ M
ργ
u for any y ∈
dom f∗ ,∅ (since constant). Here are some ways to generate “nice” random integrands.
Example 3.10 (Random scaling). Let g :R→R be a (non-random) ﬁnite convex function
. 0, and W∈ L0 be strictly positive (i.e.,P(W > 0) = 1). Then put
f (ω, x) := g(W(ω)x), ∀(ω, x)∈ Ω×R.
In this case, f∗(ω, y) = g∗(y/W(ω)) andR⊂ D f,γ is true if
(32) ∃δ > 0, p > 1 such that g(−δW p)+∨ g(δW p)+∈ M
ργ
u .
Note that|W x| = δ
2|W(2x/δ)| ≤1
2
(
δW p + 2q
δq−1|x|q
)
where 1
p + 1
q = 1. Applying the (quasi)
convexity of g twice, (32) implies for each x∈R
g(W x)≤ g(−δW p)+∨ g(δW p)+∨ g
(
− 2q
δq−1|x|
)+
∨ g
( 2q
δq−1|x|
)+
∈ M
ργ
u .
Also, since g . 0, domg∗\{ 0} ,∅. If y∈ domg∗ and y > 0 (resp. y < 0),
0≤ W≤ 1 + W p≤ 1 + g(δW p)+ + g∗(y)
yδ ; resp. ≤ 1 + g(−δW p)+ + g∗(y)
−yδ .
In both cases, (32) implies W∈ Mργ, and consequently, ηy = yW∈ Mργ and f∗(·, ηy) = g∗(y)∈
L∞. Thus (31) (⇒ (11)) follows from (32) as well. If in addition g is monotone increasing,
g(−W p)+≤ g(0)+, thus the half of (32) is automatically true. ♦
Example 3.11 (Random parallel shift). Let f be a ﬁnite normal convex integrand satisfying
(11) andR⊂ D f , and B∈ L0. Then put
(33) fB(ω, x) = f (ω, x + B(ω)), (ω, x)∈ Ω×R.

<!-- page 13 -->
Robust Integral Functionals 13
By convexity of f , f (·, x + B)≤ ε
1+ε f
(
·, 1+ε
ε x
)
+ 1
1+ε f (·, (1 + ε)B) and f
(
·, ε
1+ε x
)
≤ ε
1+ε f (·, x +
B) + 1
1+ε f (·,−εB), thus putting Γα(x) = f (·, αx)+/α,
1 + ε
ε f
(
·, ε
1 + ε x
)
− Γε(−B)≤ fB(·, x)≤ ε
1 + ε f
(
·, 1 + ε
ε x
)
+ Γ1+ε(B),(34)
f∗(·, y)− Γ1+ε(B)≤ f∗
B(·, y)≤ f∗(·, y) + Γε(−B),(35)
where f∗
B(·, y) = f∗(·, y)− yB, and (35) follows from (34) by taking conjugates. Thus if
∃ε > 0 such that Γ1+ε(B)∈ M
ργ
u and Γε(−B)∈ Lργ,(36)
thenR⊂ D fB,γ and f∗
B satisﬁes (11) ((31) if f∗ does). Moreover, (35) implies in this case
♦
H f∗
B,γ(η) <∞⇔ H f∗,γ(η) <∞⇒ ηB∈ L1,
and H f∗
B,γ is explicitly given in terms of H f∗,γ as
(37) H f∗
B,γ(η) = H f∗,γ(η)−E[ηB], ∀η∈ dom(H f∗
B,γ) = dom(H f∗,γ).
We can combine the preceding two examples:
Example 3.12. Let g :R→R, W and B be as in Examples 3.10 and 3.11, and put
h(·, x) := g(W x+ B) = g(W(x + B/W))
This h satisﬁes (11) andR⊂ Dh,γ if (g, W) satisﬁes (32) and (36) holds with f = g. Note that
if we apply Example 3.11 to f (·, x) = g(W x) and B/W, then h(·, x) = f (·, x + B/W) and e.g.
f (·, (1 + ε)B/W) = g((1 + ε)B). ♦
Our initial motivation was a duality method for robust utility maximization of the general
form (3) with random utility function U : Ω×R→R. See [10] for the ﬁnancial background of
the problem. A motivational example of random utility is of the typeUD,B(·, x) = U(D−1x+ B)
where U :R→R is a proper concave increasing function and D, B∈ L0 (with D > 0 a.s.)
correspond respectively to the discount factor and a payoﬀ of a claim. Then the problem is to
(38) maximize uD,B,γ(ξ) := inf
Q∈Q
(
EQ[U(D−1ξ + B)] + γ(Q)
)
=−I fD,B,γ(−ξ)
over a convex cone C⊂ L∞ where fD,B(·, x) =−U(−D−1x + B) is a proper normal convex inte-
grand of the form in Example 3.12. The full detail of this problem in more concrete ﬁnancial
setup will be given in a separate future paper together with an application to a robust version
of utility indiﬀerence valuation. Here we just give a criterion in terms of “integrabilities” of
D and B for the duality without singular term as well as its explicit form when U is ﬁnite on
R. It constitutes a half of what we call the martingale duality method (see e.g. [2, 25, 3]
for the other half in the classical case and [18] 1 for a partial result in the robust case). For
the case dom(U) =R+, see [26] when D, B are constants; [27] with bounded B, and [9] for
dom(U) =R with constant D, B; see also [10] for more thorough references. The following is
an immediate consequence of Corollary 3.7 and Example 3.12.
1There an earlier version of this paper (still available as arXiv:1101.2968) was used.

<!-- page 14 -->
14 K. Owari
Corollary 3.13. In the above notation, suppose U is ﬁnite onR, and
(39) ∃δ, ε > 0 with U(−δD−(1+ε))−∈ M
ργ
u , U(−(1 + ε)B)−∈ M
ργ
u , U(εB)−∈ Lργ.
Then for any nonempty convex cone C⊂ L∞, it holds that
sup
ξ∈C
uB,D,γ(ξ) = min
η∈C◦
V
(
HV,γ(Dη) +E[DBη]
)
,(40)
where V(y) := supx∈R(U(x)− xy) and C◦
V :={η∈ C◦ : HV,γ(Dη) <∞}.
Here C can be any convex cone and the possibility of both sides being∞ is not excluded;
it does not happen iﬀ C◦
V ,∅. If in addition C◦,e
V :={η∈ C◦
V : η > 0 a.s.} ,∅, we can replace
“min η∈C◦
V” by “inf η∈C◦,e
V
” etc with a little more e ﬀort and certain regularities of U. Choosing a
“good” cone C, these conditions as well as the dual problem have clear ﬁnancial interpretations
and consequences (see [2] for a good exposition in the classical case).
A couple of features deserve attention: (i) We directly invoke Fenchel’s theorem to the
functional uD,B,γ = −I fD,B,γ(−·) by means of our main theorem, instead of interchanging
“sup ξ∈C” and “inf Q∈Q”, and invoke a classical duality (see [2, 3] and its references) under
each Q, so we do not need to mind what happens under “extreme” Q∈ Qγ; (ii) embedding the
randomness D, B to the utility function U instead of transforming the domain C to D−1C + B,
we retain the “good form” of C (which is essential for the probabilistic techniques to work
well), and obtain a criterion for the duality in terms solely of B and D (when U is ﬁnite).
Those integrability conditions are weak even in the classical case where γ = δ{P} and D≡ 1
(then (39) reads as U(−(1 + ε)B)−, U(εB)−∈ L1 for some ε > 0 complementing the result of
[3]).
4. Proofs
4.1. Proof of Theorem 3.1
The upper bound is simply a consequence of Young’s inequality (18) andν = νr + νs:
I∗
f,γ(ν) = sup
ξ∈dom(I f,γ)
(
νr(ξ)− I f,γ(ξ) + νs(ξ)
)(18)
≤ H f∗,γ(νr) + sup
ξ∈dom(I f,γ)
νs(ξ).
The lower bound is more involved. First, ﬁxν∈ ba and deﬁne
Lν(Q, ξ) := ν(ξ)−EQ[ f (·, ξ)] + γ(Q), Q∈ Qγ, ξ∈ D f,γ.
It is ﬁnite-valued onQγ× D f,γ, convex in Q∈ Qγ and concave in ξ∈ D f,γ. Moreover,
Lemma 4.1. With the notation above, it holds that
inf
Q∈Qγ
sup
ξ∈D f,γ
Lν(Q, ξ) = sup
ξ∈D f,γ
inf
Q∈Qγ
Lν(Q, ξ).(41)
Proof. We claim that for any ξ∈ D f,γ and c > 0, Λc(ξ) :={Q∈ Qγ : Lν(Q, ξ)≤ c} is weakly
compact in L1. Again let ψQ := dQ/dP for each Q∈ Qγ. For the uniform integrability of
Λc(ξ), pick a Q0∈ Qγ with γ(Q0) = 0 (see Remark 2.2) and observe that
EQ[ f (·, ξ)]≤E
[
2 f (·, ξ)+ ψQ + ψQ0
2
]
≤ ργ
(2 f (·, ξ)+)+ 1
2γ(Q) <∞,(42)

<!-- page 15 -->
Robust Integral Functionals 15
since ξ∈ D f,γ. Thus Λc(ξ)⊂
{
Q∈ Qγ : γ(Q)≤ 2
(
c− ν(ξ) + ργ(2 f (·, ξ)+)
)}
, hence Λc(ξ) is
uniformly integrable by (6) and so isΛ′
c(ξ) := { f (·, ξ)+ψQ : Lν(Q, ξ)≤ c} by (7) since ξ∈ D f,γ.
To see that Λc(ξ) is weakly closed, it su ﬃces to show that it is norm-closed since con-
vex. So let ( Qn)n be a sequence in Λc converging in L1 to some Q, then passing to a sub-
sequence, we can suppose without loss that ψQn → ψQ a.s. too. From the previous para-
graph, ( f (·, ξ)+ψQn)n is uniformly integrable, hence by the (reverse) Fatou lemma, we have
EQ[ f (·, ξ)]≥ lim supnEQn[ f (·, ξ)] and consequently,
Lν(Q, ξ) = ν(ξ)−EQ[ f (·, ξ)] + γ(Q)
≤ ν(ξ) + lim inf
n
−EQn[ f (·, ξ)] + lim inf
n
γ(Qn)
≤ lim inf
n
(ν(ξ)−EQn[ f (·, ξ)] + γ(Qn))≤ c.
We deduce that Λc(ξ) is weakly closed, hence weakly compact. Now (41) follows from a
minimax theorem (see [19, Appendix A]).
Noting that I∗
f,γ(ν)≥ supξ∈D f,γ(ν(ξ)− I f,γ(ξ)) = supξ∈D f,γ infQ∈Qγ Lν(Q, ξ), we deduce from
Lemma 4.1 that for any ν∈ ba,
(43) I∗
f,γ(ν)≥ inf
Q∈Qγ
sup
ξ∈D f,γ
Lν(Q, ξ) = inf
Q∈Qγ
sup
ξ∈D f,γ
(ν(ξ)−EQ
[ f (·, ξ)] + γ(Q)).
Lemma 4.2. Let η, ζ, ψ∈ L1 such that ψ≥ 0 a.s. and
(44) ˜f∗ (·, η, ψ) = sup
x∈dom f
(xη− ψ f (·, x)) > ζ a.s.
Then there exists a ˆξ∈ L0 such that
(45) f (·, ˆξ) <∞ a.s. and ˆξη− ψ(·, ˆξ)≥ ζ a.s.
Proof. This amounts to proving that the multifunction
S (ω) :={x∈ dom f (ω,·) : xη(ω)− ψ(ω) f (ω, x)≥ ζ(ω)}
admits a measurable selection. S is nonempty valued by (44), and measurable sinceg(ω, x) :=
ψ(ω) f (ω, x)− xη(ω) (with the convention 0·∞ = 0) is a normal convex integrand (see [24,
Prop. 14.44, Cor. 14.46]), and S (ω) = dom f (ω,·)∩{ x : g(ω, x)≤− ζ(ω)}. On {ψ > 0}, we
have simply S ={x : f (·, x)− x η
ψ≤− ζ
ψ} which is closed since f is normal. Thus
S′(ω) =

S (ω) if ω∈{ ψ > 0},
∅ if ω∈{ ψ = 0},
is a closed-valued measurable multifunction with dom S′ ={ω : S′(ω) ,∅} ={ψ > 0}, thus
the standard measurable selection theorem (see [24], Cor. 14.6) shows the existence ofξ′∈ L0
such that ξ′(ω)∈ S′(ω) = S (ω) for ω∈{ ψ > 0}.
On{ψ = 0}, the multifunction S need not be closed-valued. So we explicitly construct a
selector. First, on {ψ = η = 0}, we can take any ξ0∈ domI f,γ(,∅ by assumption), which
satisﬁesξ(ω)∈ dom f (ω,·) for a.e. ω by (13), and ξ0η− ψ f (·, ξ0) = 0 > ζ on{ψ = η = 0} since
(44) reads as ζ < 0 when ψ = η = 0.

<!-- page 16 -->
16 K. Owari
Next, we put
ξ′′ := 1
2
((
a+
f∧ ζ
η
)
+
(
a−
f∨ ζ
η
))
on{ψ = 0, η , 0},
where a±
f∗(ω) = limx→±∞ f∗(ω, x)/x as in (15). Recalling that dom f (ω,·) = {a+
f∗(ω)} if
a−
f∗(ω) = a+
f∗(ω), and otherwise ( a−
f∗(ω), a+
f∗(ω))⊂ dom f (ω,·)⊂ [a−
f∗(ω), a+
f∗(ω)], we have
ξ′′(ω)∈ dom f (ω,·) for ω∈{ ψ = 0, η , 0}. Also, (44) reads as a+
f η > ζ when ψ = 0 and η > 0,
hence
ξ′′η− ψ f (·, ξ′′) = ξ′′η = 1
2ζ + 1
2
(
a−
f∨ ζ
)
≥ ζ on{ψ = 0, η > 0}.
Similarly, (44) reads as a−
f η > ζ on{ψ = 0, η < 0}, hence
ξ′′η = 1
2
((
a+
f η∨ ζ
)
+
(
a−
f η∧ ζ
))
= 1
2
(
a+
f∨ ζ
)
+ 1
2ζ≥ ζ on{ψ = 0, η < 0}.
Now ˆξ := ξ′1{ψ>0} + ξ01{ψ=0,η=0} + 1
2ξ′′1{ψ=0,η,0} is a desired measurable selection.
Proof of Theorem 3.1. The upper bound is already established at the beginning of this section.
For the lower bound, it suﬃces to show that for any ν = νr + νs∈ ba,
(46) α < H f∗,γ(νr), β < sup
ξ∈D f,γ
νs(ξ)⇒ α + β < I∗
f,γ(ν).
In the sequel, we ﬁx an arbitraryν = νr + νs∈ ba and α, β as in (46), and we denote
η := dνr
dP and ψQ := dQ
dP ,∀Q∈ Qγ.
By the assumption on β, there exists a ξs∈ D f,γ with νs(ξs) > β, and by the singularity of
νs, there exists an increasing sequence (An)n in F withP(An)↑ 1 and|νs|(An) = 0 for all n, so
that νs(ξs1Acn) = νs(ξs) > β. On the other hand, the assumption on α implies that
α < H f∗(νr|Q) + γ(Q) =E
[ ˜f∗ (·, η, ψQ
)]
+ γ(Q), ∀Q∈ Qγ.
Then for each Q∈ Qγ, there exists a ζQ∈ L1 with
E[ζQ] > α− γ(Q) and ζQ < ˜f∗(·, η, ψQ) a.s.
(even if Φ := ˜f∗ (·, η, ψQ
)< L1: since Φ−∈ L1, choosing ε > 0 so thatE[Φ]− ε > α , we have
limNE[(Φ− ε)∧ N] > α by the monotone convergence theorem, so ( Φ− ε)∧ N0 with a big
N0 does the job.) Therefore, Lemma 4.2 implies that there exists a ξ0
Q∈ L0 such that
(47) f (·, ξ0
Q) <∞ and ξ0
Qη− ψQ f (·, ξ0
Q)≥ ζQ a.s.
Note that this does not guarantees that ξ0
Q is in D f,γ (if it was, there would be nothing to
prove anymore). So we approximate ξ0
Q by elements of D f,γ. Let Bn :={|ξ0
Q|≤ n}∩{| f (·, ξ0
Q)|≤
n} and Cn := An∩ Bn, thenP(Cn)↑ 1 since f (·, ξ0
Q) <∞ a.s. Put
ξn
Q := ξ0
Q1Cn + ξs1Ccn, ∀n∈N.

<!-- page 17 -->
Robust Integral Functionals 17
Then for each n, ξn
Q ∈ D f,γ since∥ξ0
Q∥∞ ≤ n +∥ξs∥∞ <∞ and f (·, ξn
Q)+ = f (·, ξ0
Q)+1Cn +
f (·, ξs)+1Ccn≤ n + f (·, ξs)+∈ M
ργ
u by the solidness of the space. On the other hand,
νr(ξn
Q)−EQ
[
f (·, ξn
Q)
]
=E
[
ξn
Qη− ψQ f (·, ξn
Q)
]
=E
[
1Cn
(
ξ0
Qη− ψQ f (·, ξ0
Q)
)]
+E
[
1Ccn
(ξsη− ψQ f (·, ξs))]
≥E[ζQ1Cn] +E
[
1Ccn
(ξsη− ψQ f (·, ξs))]
=E[ζQ] +E
[
1Ccn
(ξsη− ψQ f (·, ξs)− ζQ                                        
=:ΞQ
)]
.
Note that ΞQ ∈ L1 since f (·, ξs) ∈ ⋂
Q∈Qγ L1(Q) (by ξs ∈ D f,γ), thus lim nE[1CcnΞQ] = 0.
Therefore, noting that νs(ξn
Q) = νs(1Ccnξn
Q) = νs(1Ccnξs) = νs(ξs),
sup
ξ∈D f,γ
(ν(ξ)−EQ[ f (·, ξ)])≥ sup
n
(
νr(ξn
Q)−EQ[ f (·, ξn
Q)] + νs(ξn
Q)
)
≥ lim sup
n
(
E[ζQ] +E[1CcnΞQ] + νs(ξs)
)
=E[ζQ] + νs(ξs) > α− γ(Q) + β.
In view of (43), we have
I∗
f,γ(ν)
(43)
≥ inf
Q∈Qγ
sup
ξ∈D f,γ
(ν(ξ)−EQ[ f (·, ξ)] + γ(Q))≥ α + β.
Since α < H f∗,γ(νr) and β < supξ∈D f,γ νs(ξ) are arbitrary, this completes the proof.
4.2. Proof of Theorem 3.4
In the sequel, the assumptions of Theorem 3.4excepting (30) are supposed without notice. The
implication (i)⇒ (ii) is trivial sinceR⊂ L∞, and (i)⇒ (iii) is clear from (23) of Theorem 3.1,
while (ii)⇒ (i) follows from (cf. [22, 23])
(48) a≤ ξ≤ b a.s., a, b∈R⇒ f (·, ξ)≤ f (·, a)+ + f (·, b)+.
Indeed, the assumption implies the existence of a [0 , 1]-valued random variable α such that
ξ = αa+(1−α)b a.s. Since f (ω,·) is convex for a.e.ω, we see that f (ω, ξ(ω))≤ α(ω) f (ω, a)+
(1− α(ω)) f (ω, b)≤ f (ω, a)+ + f (ω, b)+ for a.e. ω∈ Ω.
Given dom(I f,γ) = L∞ and σ(L∞, L1)-lsc of I f,γ (Lemma 2.3) as well as I∗
f,γ|L1 = H f∗,γ
(Corollary 3.2), (iv)⇔ (v) is a special case of the following (with E = L∞ and E′ = L1):
Lemma 4.3 ([15], Propositions 1 and 2). Let⟨E, E′⟩ be a dual pair and ϕ a σ(E, E′)-lsc ﬁ-
nite convex function on E with the conjugate ϕ∗ on E′. Then ϕ is τ(E, E′)-continuous on E if
and only if{x′∈ E′ : ϕ∗(x′)≤ c} is σ(E′, E)-compact for each c∈R.
Proof of (iii)⇒ (iv). Put E = L∞, E′ = ba, then τ(L∞, ba) is the norm-topology while the
ﬁnite lsc convex function I f,γ on the Banach space L∞ is norm-continuous (see [8, Ch.1,
Cor. 2.5]). Thus Λc :={ν∈ ba : I f,γ(ν)≤ c} is σ(ba, L∞)-compact from Lemma 4.3, a fortiori
{η∈ L1 : H f∗,γ(η)≤ c} = Λc∩ L1 = Λc (by Corollary 3.2 and (iii)) is σ(L1, L∞)-compact since
σ(L1, L∞) = σ(ba, L∞)|L1.

<!-- page 18 -->
18 K. Owari
Proof of (v)⇒ (vi). This follows from the observation that
sup
n
∥ξn∥∞ <∞ and ξn→ ξ a.s. ⇒ ξn→ ξ for τ(L∞, L1).
Indeed, for any weakly compact (⇒ uniformly integrable) subset C⊂ L1,
sup
η∈C
E[|ξ− ξn||η|]≤ sup
η∈C
E[|ξ− ξn||η|1{|η|>N}] + sup
η∈C
E[|ξ− ξn||η|1{|η|≤N}]
≤ 2 sup
n
∥ξn∥∞ sup
η∈C
E[|η|1{|η|>N}] + NE[|ξ− ξn|].
Taking a diagonal, we see qC(ξ− ξn) := supη∈C|E[(ξ− ξn)η]|→ 0, while qC with C running
through (convex, circled) weakly compact subsets of L1 generates τ(L∞, L1).
Proof of (iv)⇒ (vii). Fix an arbitrary η0∈ dom(H f∗,γ) (,∅) and ξ∈ L∞. Observe that
E[ξ(η + η0)] =E
[
2ξ η + η0
2
]
≤ I f,γ(2ξ) + 1
2
(
H f∗,γ(η) + H f∗,γ(η0)
)
,∀η∈ L1.
HenceE[ξη]− H f∗,γ(η) ≤ I f,γ(2ξ) +∥ξ∥∞∥η0∥1 + 1
2 H f∗,γ(η0)− 1
2H f∗,γ(η). Putting Cc,ξ :=
2
(
c + I f,γ(2ξ) +∥ξ∥∞∥η0∥1
)
+ H f∗,γ(η0) which does not depend on η, we see that
E[ξη]− H f∗,γ(η)≥ c⇒ H f∗,γ(η)≤ Cc,ξ.
Consequently,{η∈ L1 : E[ξη]− H f∗,γ(η)≥ c} is weakly compact for each c > 0, since it is
contained in a weakly compact set {η∈ L1 : H f∗,γ(η)≤ Cc,ξ} and η↦→ E[ξη]− H f∗,γ(η) is
weakly upper semicontinuous. Therefore, supη∈L1
(
E[ξη]− H f∗,γ(η)
)
is attained.
From now on, we assume (30) (thus all the assumptions of Theorem 3.4) which implies
that
(49) f (·, ξ)+∈ Mργ,∀ξ∈ L∞.
Indeed, by dom(I f,γ) = L∞ and∃ξ0∈ D f,γ, we have for all ξ∈ L∞ that
ργ (λ f (·, ξ))≤ 1
2ργ ( f (·, 2λξ− (2λ− 1)ξ0)) + 1
2ργ
((2λ− 1) f (·, ξ0)+)<∞, λ > 1.
Here the ﬁrst term in the right hand side is 1
2 I f,γ(2λξ− (2λ− 1)ξ0) <∞. On the other hand,
if f (·, ξ′
0)−∈ Mργ as in (30), putting A ={ f (·, ξ)≥ 0} with an arbitrary ξ∈ L∞, we see that
f (·, ξ1A+ξ′
01Ac) = f (·, ξ)++ f (·, ξ′
0)1Ac, hence f (·, ξ)+≤ f (·, ξ1A+ξ′
01Ac)+ f (·, ξ′
0)−. Therefore,
ργ (λ f (·, ξ)+)≤ 1
2 ργ
(
2λ f (·, ξ1A + ξ′
01Ac)
)
+ 1
2ργ
(
2λ f (·, ξ′
0)−
)
,∀λ > 1.
Proof of (vi)⇒ (i). In view of (7) and (49), it suﬃces to show that for each ξ∈ L∞,
lim
N
sup
γ(Q)≤1
EQ
[
f (·, ξ)+1{ f (·,ξ)+≥N}
]
= 0, ∀c > 0.
Let AN :={ f (·, ξ)≥ N}∈ F (N∈N), thenP(AN)→ 0 since I f,γ is ﬁnite. Pick a ξ0∈ D f,γ and
put ξN := (ξ−ξ0)1AN so that ξN +ξ0 = ξ1AN +ξ01Ac
N. Then f (·, ξN +ξ0) = f (·, ξ)1AN + f (·, ξ0)1Ac
N,
while for λ > 1, f (·, ξN + ξ0)≤ 1
λ f (·, λξN + ξ0) + λ−1
λ f (·, ξ0), hence
f (·, ξ)+1AN = f (·, ξ)1AN≤ 1
λ f (·, λξN + ξ0) + 1
λ f (·, ξ0)− + f (·, ξ0)+1AN , ∀λ > 1.

<!-- page 19 -->
Robust Integral Functionals 19
Note that since ξ0∈ D f,γ and since f (·, ξ0)−∈ Lργ,
C1 := sup
γ(Q)≤c
EQ[ f (·, ξ0)−] <∞ and lim
N
sup
γ(Q)≤c
EQ[ f (·, ξ0)+1AN] = 0.
Also,EQ[ f (·, λξN + ξ0)]≤ ργ( f (·, λξN + ξ0)) + γ(Q) = I f,γ(λξN + ξ0) + γ(Q), while for each
λ > 1, we have sup N∥λξN + ξ0∥∞≤ λ∥ξ∥∞ + (λ− 1)∥ξ0∥∞ <∞ and λξN + ξ0→ ξ0 a.s., thus
I f,γ(λξN + ξ0)→ I f,γ(ξ0) by (vi), hence for some big Nλ, I f,γ(λξN + ξ0)≤ I f,γ(ξ0) + 1 =: C2
for N > Nλ. Summing up,
sup
γ(Q)≤c
EQ[ f (·, ξ)+1AN]≤ C2 + c + C1
λ + sup
γ(Q)≤c
EQ[ f (·, ξ0)+1AN],∀N > Nλ, λ > 1.
Now a diagonal argument yields that supγ(Q)≤cEQ[ f (·, ξ)+1AN]→ 0.
Proof of (vii)⇒ (iv). Under (49), lim∥η∥1→∞ H f∗,γ(η)/∥η∥1 =∞, i.e., H f∗,γ is coercive on L 1.
For, since∥η∥1 =E[sgn(η)η] where sgn(η) = 1{η≥0}− 1{η<0}∈ L∞, f (·, nsgn(η))≤ f (·,−n)+ +
f (·, n)+∈ Mργ by (48), (49), and H f∗,γ(η) = supξ∈L∞(E[ξη]− I f,γ(ξ)), we have
H f∗,γ(η)≥E[nsgn(η)η]− I f,γ(nsgn(η))
≥ n∥η∥1− ργ (2 f (·,−n)+)− ργ (2 f (·, n)+)
2 >−∞,∀n.
Then the coercive James’ theorem [16, Th. 2] applied to the coercive function H f∗,γ on the
Banach space E = L1 implies the relative weak compactness of all the sublevel sets{η∈ L1 :
H f∗,γ(η)≤ c} which are weakly closed since H f∗,γ is weakly lower semicontinuous.
Appendix
A. Lower semicontinuity of H f ∗ + γ
Lemma A.1. Under the assumptions of Theorem 3.1, H f∗ + γ is jointly weakly lower semi-
continuous on L1× Qγ, and infQ∈Qγ
(
H f∗(η|Q) + γ(Q)
)
is attained for every η∈ L1.
Proof. Let ξ0 be as in (22) (so f (·, ξ0)+ ∈ M
ργ
u ), ψQ = dQ/dP for each Q ∈ Qγ, and put
K(c, a) := 2
(
c + ργ (2 f (·, ξ0)+) + a∥ξ0∥∞
)
<∞. Then by (42) in the proof of Lemma 4.1,
∥η∥1≤ a and H f∗(ν|Q) + γ(Q)≤ c⇒ γ(Q)≤ K(c, a).(50)
Also, from (19) with ξ0 above and (7), we see that whenever (ηn)n is uniformly integrable and
supn γ(Qn) <∞, { ˜f∗(·, ηn, ψQn)−}
n is uniformly integrable.
To see that H f∗ + γ is weakly lsc, it su ﬃces by convexity that for each c ∈ R, Λc :=
{(η, ψQ)∈ L1× L1 : Q∈ Qγ, H f∗(η|Q) + γ(Q)≤ c} is norm-closed. Let ( ηn, ψQn)n⊂ Λc be
norm convergent to (η, ψQ). Passing to a subsequence, we can assume a.s. convergence too.
By the above paragraph and the norm convergence, ( ηn)n and { ˜f∗(·, ηn, ψQn)−}
n are uniformly
integrable, while γ(Q)≤ supn γ(Qn) <∞ by (50) and lsc of γ. Thus by Fatou’s lemma,
H f∗(η|Q) + γ(Q) =E
[ ˜f∗(η, ψQ)
]
+ γ(Q)
≤ lim inf
n
E
[ ˜f∗(ηn, ψQn)
]
+ lim inf
n
γ(Qn)
≤ lim inf
n
(
E
[ ˜f∗(ηn, ψQn)
]
+ γ(Qn)
)
≤ c,

<!-- page 20 -->
20 K. Owari
hence (η, ψQ)∈ Λc, obtaining the lower semicontinuity ofH f∗ + γ. In particular, H f∗(η|·)+ γ(·)
is weakly lower semicontinuous onQγ, so another application of (50) as well as (6) imply that
the inﬁmum infQ∈Qγ
(
H f∗(η|Q) + γ(Q)
)
is attained for each η∈ L1.
B. Some Details of Example 3.3
Let (N, 2N,P) and ( Pn)n be as in Example 3.3 with P = conv(Pn, n ∈ N) (closed convex
hull). The weak compactness of P follows from sup n Pn({k, k + 1, k + 2, ...}) = sup{1/n :
n≥ k} = 1/k → 0 as k → ∞. Also ργ(ξ) = supP∈PEP[ξ] = supnEPn[ξ] if ξ ≥ 0 since
P↦→ EP[ξ∧ N] is continuous for any N∈ N, so sup P∈PEP[ξ] = supN supP∈PEP[ξ∧ N] =
supN supP∈conv(Pn;n∈N)EP[ξ∧ N], while if P = α1Pn1 +··· αlPnl, thenEP[ξ∧ N] = α1EPn1
[ξ∧
N] +··· + αlEPnl
[ξ∧ N]≤ max1≤i≤l EPni
[ξ∧ N]≤ supn EPn[ξ∧ N].
Lemma B.1. Let f be given by (27) in Example 3.3. Then we have
(51) lim
N→∞
sup
n
EPn[ f (·, ξ)1{ f (·,ξ)≥N}] = lim sup
n
ξ(n)+eξ(n)+
.
Proof. Let h(x) : = x+ex and ﬁx ξ∈ ℓ∞. If ∥ξ∥∞ = 0, then both sides of (51) are 0, thus we
assume∥ξ∥∞ > 0 (⇔ h(∥ξ∥∞) > 0). Note that for any N, n∈N, we have (by deﬁnition)
EPn[ f (·, ξ)1{ f (·,ξ)≥N}] =
(
1− 1
n
)
h(ξ(1))1{h(ξ(1))≥N} + h(ξ(n))1{nh(ξ(n))≥N}.
In particular,EPn[ f (·, ξ)1{ f (·,ξ)≥N}]≤ h(ξ(1))1{h(ξ(1))≥N} + h(ξ(n))1{nh(∥ξ∥∞)≥N}, thus
lim
N→∞
sup
n
EPn[ f (·, ξ)1{ f (·,ξ)≥N}]≤ lim
N→∞
sup
n≥N/h(∥ξ∥∞)
h(ξ(n)) = lim sup
n
h(ξ(n)) =: α.
This is “≤” in (51). If α = 0, we are done; otherwise, for any ε > 0 and N∈N, there exists
some nε
N > N/(α− ε) > 0 with h(ξ(nε
N)) > α− ε. In particular, nε
Nh(ξ(nε
N)) > N, hence
sup
n
EPn[ f (·, ξ)1{ f (·,ξ)≥N}]≥ sup
n
h(ξ(n))1{nh(ξ(n))≥N}≥ h(ξ(nε
N)) > α− ε.
This proves “≥” in (51).
Since h(x) = x+ex is increasing, continuous and h(0) = 0, lim sup n h(ξ(n)) = 0⇔
lim supn ξ(n) = 0, hence (28). Recall that dom(I∗
f,γ)⊂ ba+ (since I f,γ is increasing). Finally,
Lemma B.2. The conjugate I∗
f,γ is explicitly given on bas
+ as:
I∗
f,γ(ν) = sup
x≥0
x(∥ν∥− ex), ∀ν∈ bas
+.(52)
Proof. Let ν∈ bas
+. SinceEPn[ f (·, ξ)] =
(
1− 1
n
)
h(ξ(1)) + h(ξ(n)), we have
h(ξ(n))
(∗)
≤EPn[ f (·, ξ)]
(∗∗)
≤ h(ξ(1)) + h(ξ(n)).
From (∗) and (26), ν(ξ)− I f (ξ)≤ ∥ν∥ lim supn ξ(n)− supn h(ξ(n))≤ ∥ξ+∥∞(∥ν∥− e∥ξ+∥∞)≤
supx≥0 x(∥ν∥− ex) which shows “≤” in (52). Considering ¯x0 := (0, x, x,··· )∈ ℓ∞ with x≥ 0
(then∥ ¯x0∥∞ = x and ν( ¯x0) = x∥ν∥ since ν vanishes on ﬁnite sets), we deduce from (∗∗) that
I∗
f,γ(ν)≥ sup
ξ∈ℓ∞
(ν(ξ)− h(ξ(1))− h(∥ξ∥∞))≥ sup
x≥0
(x∥ν∥− xex).

<!-- page 21 -->
Robust Integral Functionals 21
Acknowledgements
The author gratefully acknowledges the ﬁnancial support of the Center for Advanced Research
in Finance (CARF) at the Graduate School of Economics of the University of Tokyo.
References
[1] Aliprantis, C. D. and K. C. Border (2006): Inﬁnite dimensional analysis: A hitchhiker’s guide.
Springer, Berlin, 3rd ed.
[2] Bellini, F. and M. Frittelli (2002): On the existence of minimax martingale measures. Math.
Finance 12, 1–21.
[3] Biagini, S., M. Frittelli and M. Grasselli (2011): Indi ﬀerence price with general semimartingales.
Math. Finance 21, 423–446.
[4] Castaing, C. and M. Valadier (1977): Convex analysis and measurable multifunctions, Lecture
Notes in Math., vol. 580. Springer-Verlag, Berlin.
[5] Delbaen, F. (2009): Di ﬀerentiability properties of utility functions. In: Optimality and risk—
modern trends in mathematical ﬁnance, Springer, Berlin, pp. 39–48.
[6] Dellacherie, C. and P.-A. Meyer (1978): Probabilities and potential, North-Holland Mathematics
Studies, vol. 29. North-Holland Publishing Co., Amsterdam.
[7] Dunford, N. and J. T. Schwartz (1988): Linear operators. Part I. General theory. Wiley Classics
Library. John Wiley & Sons Inc., New York. With the assistance of William G. Bade and Robert
G. Bartle, Reprint of the 1958 original, A Wiley-Interscience Publication.
[8] Ekeland, I. and R. Temam (1976): Convex analysis and variational problems, Studies in Mathe-
matics and its Applications, vol. 1. North-Holland Publishing Co., Amsterdam.
[9] Föllmer, H. and A. Gundel (2006): Robust projections in the class of martingale measures.Illinois
J. Math. 50, 439–472 (electronic).
[10] Föllmer, H., A. Schied and S. Weber (2009): Robust preferences and robust portfolio choice.
In: A. Bensoussan, Q. Zhang and P. G. Ciarlet (eds.), Mathematical Modelling and Numerical
Methods in Finance, Handbook of Numerical Analysis, vol. 15, North-Holland, pp. 29–88.
[11] Grothendieck, A. (1973): Topological vector spaces. Gordon and Breach Science Publishers, New
York. Translated from the French by Orlando Chaljub, Notes on Mathematics and its Applications.
[12] Hewitt, E. and K. Stromberg (1975): Real and abstract analysis: A modern treatment of the theory
of functions of a real variable, Graduate Texts in Mathematics , vol. 25. Springer-Verlag, New
York. Third printing.
[13] Jouini, E., W. Schachermayer and N. Touzi (2006): Law invariant risk measures have the Fatou
property. Advances in Mathematical Economics 9, 49–71.
[14] Kozek, A. (1980): Convex integral functionals on Orlicz spaces. Comment. Math. Prace Mat. 21,
109–135.
[15] Moreau, J.-J. (1964): Sur la fonction polaire d’une fonction semi-continue supérieurement. C. R.
Acad. Sci. Paris 258, 1128–1130.
[16] Orihuela, J. and M. Ruiz Galán (2012): A coercive James’s weak compactness theorem and non-
linear variational problems. Nonlinear Anal. 75, 598–611.
[17] Orihuela, J. and M. Ruiz Galán (2012): Lebesgue property for convex risk measures on Orlicz
spaces. Math. Financ. Econ. 6, 15–35.
[18] Owari, K. (2012): On admissible strategies in robust utility maximization. Math. Financ. Econ.
6, 77–92.
[19] Owari, K. (2014): Maximum Lebesgue extension of monotone convex functions. J. Funct. Anal.
266, 3572–3611.

<!-- page 22 -->
22 K. Owari
[20] Owari, K. (2014): On the Lebesgue property of monotone convex functions. Math. Financ. Econ.
8, 159–167.
[21] Rockafellar, R. T. (1966): Extension of Fenchel’s duality theorem for convex functions. Duke
Math. J. 33, 81–89.
[22] Rockafellar, R. T. (1971): Integrals which are convex functionals. II.Paciﬁc J. Math.39, 439–469.
[23] Rockafellar, R. T. (1976): Integral functionals, normal integrands and measurable selections. In:
Nonlinear operators and the calculus of variations (Summer School, Univ. Libre Bruxelles, Brus-
sels, 1975), Lecture Notes in Math., vol. 543, Springer, Berlin, pp. 157–207.
[24] Rockafellar, R. T. and R. J.-B. Wets (1998): Variational analysis, Grundlehren der Mathematis-
chen Wissenschaften, vol. 317. Springer-Verlag, Berlin.
[25] Schachermayer, W. (2003): A super-martingale property of the optimal portfolio process.Finance
Stoch. 7, 433–456.
[26] Schied, A. (2007): Optimal investments for risk- and ambiguity-averse preferences: a duality
approach. Finance Stoch. 11, 107–129.
[27] Wittmüss, W. (2008): Robust optimization of consumption with random endowment. Stochastics
80, 459–475.
[28] Z ˘alinescu, C. (2002): Convex analysis in general vector spaces. World Scientiﬁc Publishing Co.,
Inc., River Edge, NJ.