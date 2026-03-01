# Convex Integration and Legendrian Approximation of Curves

**arXiv ID:** 1507.07661v2

**Authors:** Norbert Hungerbühler, Thomas Mettler, Micha Wasem

**Abstract:** Using convex integration we give a constructive proof of the well-known fact that every continuous curve in a contact $3$-manifold can be approximated by a Legendrian curve.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:1507.07661v2  [math.DG]  21 Jul 2016
CONVEX INTEGRA TION AND LEGENDRIAN
APPROXIMA TION OF CUR VES
NORBERT HUNGERB ¨UHLER, THOMAS METTLER AND MICHA W ASEM
Abstract. Using convex integration we give a constructive proof of the well-
known fact that every continuous curve in a contact 3-manifo ld can be approx-
imated by a Legendrian curve.
1. Introduction
A contact structure on a 3-manifold M is a maximally non-integrable rank 2
subbundle ξ of the tangent bundle of M . If α is a 1-form on M whose kernel is
ξ, then ξ is a contact structure if and only if α ∧ dα ̸= 0. A curve η in a contact
3-manifold ( M, ξ) is called Legendrian, whenever η∗α = 0 for some (local) 1-form
α deﬁning ξ.
The purpose of this note is to give a detailed proof of the foll owing statement
which is often used in contact geometry and Legendrian knot t heory.
Theorem 1.1. Any continuous map from a compact 1-manifold to a contact 3-
manifold can be approximated by a Legendrian curve in the C 0-Whitney topology.
Whereas this theorem is a special case of Gromov’s h-principle for Legendrian
immersions [Gro86], the curve-case can be treated by more elementary techniqu es.
Sketches of proofs of Theorem 1.1 have already appeared in the literature, see for
example [ Etn05, p.6-7], [ Gei06, p.40] or [ Gei08, p.102]. Exploiting the fact that
every contact 3-manifold is locally contactomorphic to R3 equipped with the stan-
dard contact structure deﬁned by α = dz − ydx, Etnyre and Geiges indicate that
either the front-projection (x, z) of a given curve ( x, y, z) can be approximated
by a zig-zag-curve whose slope approximates the y-component of the curve or the
Lagrangian projection (x, y) can be approximated by a curve whose area integral
approximates the z component of the curve, which can be achieved by adding
Date: April 30, 2015.
1

<!-- page 2 -->
small negatively or positively oriented loops.
Here, we give a diﬀerent and analytically rigorous proof of Th eorem 1.1 by
using convex integration. Our proof has the advantage of pro viding a constructive
approximation. In particular, in the case of a continuous cu rve in R3 equipped
with the standard contact structure, we obtain an explicit L egendrian curve given
in terms of an elementary integral. For instance, we obtain a n explicit solution
to the “parallel parking problem” in Example 3.1. Example 3.2 shows how our
technique recovers the zig-zag-curves and the small loops i n the front - respectively
Lagrangian projections.
2. Proof of the Theorem
We start by ﬁrst treating the case where the contact manifold is R3 equipped
with the standard contact structure, that is, we aim to prove the following:
Proposition 2.1. Let υ ∈ C 0([0, 2π], R3). For every ε > 0 there exists a Legen-
drian curve η ∈ C ∞([0, 2π], R3) such that ∥υ − η∥C0([0,2π]) ⩽ ε.
Remark 2.2. Here, as usual, ∥γ∥C0(I) := sup t∈I |γ(t)| and ∥γ∥C1(I) := ∥γ∥C0(I) +
∥dγ∥C0(I).
Let the curve we wish to approximate be given by ( x, y, z) ∈ C ∞([0, 2π], R3).
The regularity is no restriction due to a standard approxima tion argument us-
ing convolution. Let η = ( a, b, c) ∈ C ∞([0, 2π], R3) denote the approximating
Legendrian curve. For every choice of smooth functions ( a, c) ∈ C ∞([0, 2π], R2)
satisfying ˙a ̸= 0, we obtain a Legendrian curve by deﬁning b = ˙c/ ˙a. Therefore, if
( ˙a(t), ˙c(t)) lies in the set
Rt,ε :=
{
(u, v) ∈ R2, |v − y(t)u|⩽ ε min{|u|, |u|2}
}
,
for every t ∈ [0, 2π], then ∥b − y∥C0([0,2π]) ⩽ ε. This condition can be achieved by
deﬁning
(a(t), c(t)) := ( x(0), z(0)) +
∫ t
0
γ(u, nu) du,
with γ ∈ C ∞([0, 2π]× S1, R2) and n ∈ N, provided that γ(t, ·) ∈ R t,ε. Furthermore,
if γ additionally satisﬁes
1
2π
∮
S1
γ(t, s) ds = ( ˙x(t), ˙z(t)),
2

<!-- page 3 -->
for all t ∈ [0, 2π], then – as we will show below – ( a(t), c(t)) approaches ( x(t), z(t))
as n gets suﬃciently large.
The set Rt,ε is ample, i.e., the interior of its convex hull is all of R2. For
any given point ( ˙x(t), ˙z(t)) ∈ R2 we will thus be able to ﬁnd a loop in Rt,ε having
( ˙x(t), ˙z(t)) as its barycenter. This fact is sometimes referred to as th e fundamental
lemma of convex integration (see for instance [ Spr10, Prop. 2.11, p. 28]). In the
particular case studied here we obtain an explicit formula f or γ:
Lemma 2.3. There exists a family of loops γ ∈ C ∞([0, 2π] × S1, R2) satisfying
γ(t, ·) ∈ R t,ε and such that
(1) 1
2π
∮
S1
γ(t, s) ds = ( ˙x(t), ˙z(t)),
for all t ∈ [0, 2π].
Proof. The map γ := (γ1, γ2), where
γ1(t, s) := r cos s + ˙x(t)
and
γ2(t, s) := γ1(t, s)
(
y(t) + 2( ˙z(t) − y(t) ˙x(t))
r2 + 2 ˙x(t)2 γ1(t, s)
)
satisﬁes ( 1) for every r > 0. If r is large enough one obtains γ(t, ·) ∈ R t,ε, where
r can be chosen independently of t by compactness of [0 , 2π]. □
We now have:
Proof of Proposition 2.1. With the deﬁnitions above we obtain
b(t) := ˙c(t)
˙a(t) = y(t) + 2( ˙z(t) − y(t) ˙x(t))
r2 + 2 ˙x(t)2 γ1(t, nt).
We are left to show that |(a, c) − (x, z)| ⩽ ε provided n is large enough. This
follows from the following estimate
(2) ∥(a, c) − (x, z)∥C0([0,2π]) ⩽ 4π2
n ∥γ∥C1([0,2π]×S1).
The estimate is in fact a geometric property of the derivativ e and can be inter-
preted as follows: Since ( ˙ a, ˙c) and ( ˙x, ˙z) coincide “in average” on shorter and
shorter intervals when n gets bigger and bigger, ( a, c) and ( x, z) tend to become
close: Let
Ik :=
[ 2πk
n , 2π(k + 1)
n
]
for k = 0, . . . ,
⌊ nt
2π
⌋
− 1 and J :=
[⌊ nt
2π
⌋ 2π
n , t
]
.
3

<!-- page 4 -->
Then we can estimate D = |(a(t), c(t)) − (x(t), z(t))|:
D =
⏐
⏐
⏐
⏐
∫ t
0
γ(u, nu) du −
∫ t
0
( ˙x, ˙z)(u) du
⏐
⏐
⏐
⏐
⩽
⌊ nt
2π ⌋−1∑
k=0
⏐
⏐
⏐
⏐
∫
Ik
γ(u, nu) du −
∫
Ik
1
2π
∫ 2π
0
γ(u, v) dv du
⏐
⏐
⏐
⏐ +
+
∫
J
(
|γ(u, nu)|+ ∥γ∥C0([0,2π]×S1)
)
du
⩽
⌊ nt
2π ⌋−1∑
k=0
⏐
⏐
⏐
⏐
1
n
∫ 2π
0
γ
( v + 2kπ
n , v
)
dv −
∫
Ik
1
2π
∫ 2π
0
γ(u, v) dv du
⏐
⏐
⏐
⏐ +
+ 4π
n ∥γ∥C0([0,2π]×S1)
⩽
⌊ nt
2π ⌋−1∑
k=0
⏐
⏐
⏐
⏐
1
2π
∫
Ik
∫ 2π
0
(
γ
( v + 2kπ
n , v
)
− γ(u, v)
)
dv du
⏐
⏐
⏐
⏐ + 4π
n ∥γ∥C0([0,2π]×S1)
⩽
⌊ nt
2π
⌋ 4π2
n2 ∥∂tγ∥C0([0,2π]×S1) + 4π
n ∥γ∥C0([0,2π]×S1)
⩽ 4π
n
(
π∥∂tγ∥C0([0,2π]×S1) + ∥γ∥C0([0,2π]×S1)
)
.
By construction, the curve ( a, b, c) is Legendrian and an approximation of ( x, y, z),
provided n is large enough. □
Next we show that we can approximate closed curves by closed L egendrian
curves.
Proposition 2.4. Let υ ∈ C 0(S1, R3). For every ε > 0 there exists a Legendrian
curve η ∈ C ∞(S1, R3) such that ∥υ − η∥C0(S1) ⩽ ε.
Proof. Using standard regularization, let the curve we wish to appr oximate be
given by ( x, y, z) ∈ C ∞([0, 2π], R3), where the values of ( x, y, z) in 0 and 2 π agree
to all orders. Deﬁne g(t) := γ2
1 (t, nt). Since ∥g∥L1([0,2π]) = O(r2) as r → ∞ , we
can choose r > 0 large enough such that f := g/∥g∥L1([0,2π]) is well-deﬁned. With
the notation
I2 :=
∫ 2π
0
γ2(u, nu) du,
4

<!-- page 5 -->
we deﬁne η = (a, b, c) as follows:
(a(t), c(t)) := (x(0), z(0)) +
∫ t
0
[
γ(u, nu) − (0, I2f (u))
]
du,(3)
b(t) := ˙c(t)
˙a(t) = y(t) + γ1(t, nt)
(
2( ˙z(t) − y(t) ˙x(t))
r2 + 2 ˙x(t)2 − I2
∥g∥L1([0,2π])
)
.(4)
A straightforward computation shows that the values of ( a, b, c) in 0 and 2 π agree
to all orders, hence η ∈ C ∞(S1, R3) and it is Legendre by construction. Using (
2)
we obtain |I2|⩽ 4π2
n ∥γ2∥C1([0,2π]×S1), hence we ﬁnd using ( 4) as r → ∞ :
∥b − y∥C0([0,2π]) ⩽ ∥γ1∥C0([0,2π]×S1)
(
1 + 1
n ∥γ∥C1([0,2π]×S1)
)
O(r−2).
For the remaining components we ﬁnd ﬁnd using ( 2) and ( 3) the uniform bound
|(a(t), c(t)) − (x(t), z(t))|⩽ 4π2
n ∥γ∥C1([0,2π]×S1) + |I2|
∥g∥L1([0,2π])
∫ t
0
g(u) du
⩽ 8π2
n ∥γ∥C1([0,2π]×S1).
Choosing r large enough and n ∼ r2 concludes the proof. □
We show now how to glue together two local approximations of a curve Γ in M
on two intersecting coordinate neighborhoods. Let therefo re Uσ and Uτ in M be
coordinate patches such that U = Uσ ∩ Uτ ̸= ∅. Let Iσ and Iτ be compact intervals
such that I = Iσ ∩ Iτ contains an open neighborhood of t = 0 (after shifting the
variable t if necessary) and such that Γ( Iσ) ⊂ Uσ, Γ( Iτ ) ⊂ Uτ . Assume without
restriction that Γ is smooth and let ( x, y, z) represent Γ on U . Suppose that
(x, y, z) is approximated by Legendrian curves σ : Iσ → R3 and τ : Iτ → R3 such
that
(5) ∥σ − (x, y, z)∥C0(I) < ε 2, ∥τ − (x, y, z)∥C0(I) < ε2
for some ﬁxed 0 < ε < 1
2 . For r > 0, deﬁne R(r) to be the smallest number such
that ¯Br(0) ⊂ conv
(
R0,ε ∩ ¯BR(0)
)
. Note that R depends continuously on r and if
r > r 0 := ε/
√
1 + y(0)2, then
(6) R(r) = r
ε
√
(1 + y(0)2) (1 + (|y(0)|+ ε)2) =: r
ε w(y(0), ε).
5

<!-- page 6 -->
Choose 0 < δ < ε 2 such that [ − δ, δ] ⊂ I and such that δ∥(x, y, z)∥C1 (I) ⩽ ε2 and
deﬁne
p1 := (σ1(− δ), σ3(− δ)),
˙p1 := ( ˙σ1(− δ), ˙σ3(− δ)),
p2 := (τ1(δ), τ3(δ)),
˙p2 := ( ˙τ1(δ), ˙τ3(δ)).
From (5) and the choice of δ we obtain ˙p1, ˙p2 ∈ C ε :=
{
(u, v) ∈ R2, |v − y(0)u| ⩽
ε|u|
}
and
p2 − p1
2δ =: p ∈ B¯r(0), where ¯r = 2ε2
δ .
Since 3¯r > r 0, we can express R(3¯r) by means of formula ( 6). This will be used
in computation ( 9). We construct a path γ = ( γ1, γ2) : [ − δ, δ] → C ε as follows:
For ρ < δ/ 2, let γ|[−δ,−δ+ρ] be a continuous path from ˙ p1 to 0 and let γ|[δ−ρ,δ] be
a continuous path from 0 to ˙ p2. We construct γ such that the quotient γ2/γ1 is
well-deﬁned on [ − δ, − δ + ρ] ∪ [δ − ρ, δ] and equals y(0) in t = − δ + ρ and t = δ − ρ.
Moreover, we require that
(7)
∫ −δ+ρ
−δ
|γ(t)|dt < δε
2 and
∫ δ
δ−ρ
|γ(t)|dt < δε
2 .
On [ − δ, − δ + ρ], such a path is for example given by
t ↦→
(
1 − δ + t
ρ
) k

 ˙σ1(− δ)
y(0) ˙σ1(− δ) + ( ˙σ3(− δ) − y(0) ˙σ1(− δ))
(
1 − δ+t
ρ
) k


provided k ∈ N is suﬃciently large. We obtain
1
2(δ − ρ)
(
2δp −
∫ −δ+ρ
−δ
γ(t) dt −
∫ δ
δ−ρ
γ(t) dt
)
=: ¯p ∈ B3¯r(0)
and hence ¯p ∈ int conv(BR(3¯r)(0) ∩ R 0,ε). Using the fundamental lemma of convex
integration we let γ|[−δ+ρ,δ−ρ] be a continuous closed loop in BR(3¯r)(0)∩R 0,ε based
at 0 such that
1
2(δ − ρ)
∫ δ−ρ
−δ+ρ
γ(t) dt = ¯p.
With these deﬁnitions we obtain
1
2δ
∫ δ
−δ
γ(t) = p.
6

<!-- page 7 -->
Now we deﬁne η = (a, b, c) : [ − δ, δ] → R3 by letting b(t) := ˙c(t)/ ˙a(t), where
(a, c)(t) := p1 +
∫ t
−δ
γ(u)du.
The curve η is well-deﬁned and Legendrian by construction. It satisﬁes η(− δ) =
σ(− δ) and η(δ) = τ (δ). Moreover, ( a, c) and ( σ1, σ3) agree to ﬁrst order in t = − δ
and so do ( a, c) and ( τ1, τ3) in t = δ. From γ([− δ, δ]) ∈ C ε and the choice of δ we
ﬁnd
(8) |b(t) − y(t)|⩽ |b(t) − y(0)|+ |y(t) − y(0)|⩽ ε + δ∥y∥C1(I) < 2ε.
Using (
5), ( 6), ( 7) and the choice of δ we obtain for the remaining components
the uniform bound
(9)
|(a, c)(t) − (x, z)(t)|⩽ |p1 − (x, z)(− δ)|+
∫ t
−δ
(|γ(u)|+ |( ˙x, ˙z)(u)|) d u
⩽ ε2 + δε +
∫ δ−ρ
−δ+ρ
|γ(u)|du + 2δ∥(x, z)∥C1 (I)
⩽ 2ε + 2δR(3¯r)
⩽ ε
(
14 + 12
(
|y(0)|+ 1
2
) 2)
.
Finally, suppose υ is a continuous curve from a compact 1-manifold N (that is,
N is a compact interval or S1) into a contact 3-manifold ( M, ξ). We ﬁx some
Riemannian metric g on M . Then it follows with the bounds ( 8,9) and the com-
pactness of the domain of υ that for every ε > 0 there exists a ξ-Legendrian curve
η such that
sup
t∈N
dg(υ(t), η(t)) < ε,
where dg denotes the metric on M induced by the Riemannian metric g. In
particular, every open neighborhood of υ ∈ C 0(N, M ) – equipped with the uniform
topology – contains a Legendrian curve N → M . Since N is assumed to be
compact the uniform topology is the same as the Whitney C 0-topology, thus
proving Theorem 1.1.
3. Examples
Example 3.1 (Parallel Parking) . The trajectory of a car moving in the plane can
be thought of as a curve [0 , 2π] → S1 × R2. Denoting by ( ϕ, a, c) the natural
coordinates on S1 × R2, the angle coordinate ϕ denotes the orientation of the car
7

<!-- page 8 -->
Figure 1. The front (top) and the Lagrangian projection (bot-
tom) of the Legendrian approximation of υ.
with respect to the a-axis and the coordinates ( a, c) the position of the car in the
plane. Admissible motions of the car are curves satisfying
˙a sin ϕ = ˙c cos ϕ.
The manifold S1 × R2 together with the contact structure deﬁned by the kernel
of the 1-form θ := sin ϕ da − cos ϕ dc is a contact 3-manifold. Indeed, we have
θ ∧ dθ = − cos2ϕ dϕ ∧ da ∧ dc − sin2ϕ dϕ ∧ da ∧ dc = − dϕ ∧ da ∧ dc ̸= 0.
Applying Theorem 1.1 with b = tan ϕ gives an explicit approximation of the curve
t ↦→ (x(t), y(t), z(t)) = (0 , 0, t).
Lemma 2.3 gives the loop
γ(t, s) = 2( r cos s, cos2 s),
and hence the desired Legendrian curve
(arccot(r sec(nt)), 2rt sinc(nt), t + t sinc(2nt)) ,
provided r is large enough and n ∼ r2.
Example 3.2 (Legendrian Helix) . The Legendrian approximation of the helix
υ : [0, 2π] → R3, t ↦→ (t, cos(5t), sin(5t)),
8

<!-- page 9 -->
with n = 2
9 r2 and r = 30 is given by
a(t) = t + 3
20 sin(200t)
b(t) = 455
451 cos(5t) + 120
451 cos(5t) cos(200t)
c(t) = sin(5 t) + 459
5863 sin(195t) + 1377
18491 sin(205t) + 180
35629 sin(395t)+
+ 20
4059 sin(405t).
and produces the zig-zags and the small loops in its front and Lagrangian projec-
tions (see Figure 1).
References
[Etn05] John B. Etnyre, Legendrian and transversal knots , Handbook of knot theory, Elsevier B.
V., Amsterdam, 2005, pp. 105–185. MR 2179261 (2006j:57050)
[FT97] Dmitry Fuchs and Serge Tabachnikov, Invariants of Legendrian and transverse knots
in the standard contact space , Topology 36 (1997), no. 5, 1025–1053. MR 1445553
(99a:57006)
[Gei06] Hansj¨ org Geiges, Contact geometry , Handbook of diﬀerential geometry. Vol. II,
Elsevier/North-Holland, Amsterdam, 2006, pp. 315–382. MR 2194671 (2007c:53123)
[Gei08] , An introduction to contact topology, Cambridge Studies in Advanced Mathemat-
ics, vol. 109, Cambridge University Press, Cambridge, 2008 . MR 2397738 (2008m:57064)
[Gro86] Mikhael Gromov, Partial diﬀerential relations , Ergebnisse der Mathematik und ihrer
Grenzgebiete (3), vol. 9, Springer-Verlag, Berlin, 1986. M R 864505 (90a:58201)
[Spr10] David Spring, Convex integration theory , Modern Birkh¨ auser Classics,
Birkh¨ auser/Springer Basel AG, Basel, 2010, Solutions to t he h-principle in geom-
etry and topology, Reprint of the 1998 edition [MR1488424]. MR 3024860
Department of Mathematics, ETH Z ¨urich, R ¨amistrasse 101, 8092 Z ¨urich, Switzer-
land
E-mail address : norbert.hungerbuehler@math.ethz.ch
E-mail address : micha.wasem@math.ethz.ch
Institute for Mathematics, Goethe University Frankfurt, R obert-Mayer-Str. 10,
60325 Frankfurt am Main, Germany
E-mail address : mettler@math.uni-frankfurt.de
9