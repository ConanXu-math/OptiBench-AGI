# Optimisation of the lowest Robin eigenvalue in the exterior of a compact set

**arXiv ID:** 1608.04896v2

**Authors:** David Krejcirik, Vladimir Lotoreichik

**Abstract:** We consider the problem of geometric optimisation of the lowest eigenvalue of the Laplacian in the exterior of a compact planar set, subject to attractive Robin boundary conditions. Under either a constraint of fixed perimeter or area, we show that the maximiser within the class of exteriors of convex sets is always the exterior of a disk. We also argue why the results fail without the convexity constraint and in higher dimensions.

---

> **Note:** This text was extracted with pypdf (plain-text fallback). LaTeX formulas may be garbled. Install `marker-pdf` for better results.

<!-- page 1 -->
arXiv:1608.04896v2  [math.SP]  25 Aug 2017
OPTIMISATION OF THE LOWEST ROBIN EIGENV ALUE IN THE EXTERIOR OF A
COMPACT SET
DA VID KREJˇCI ˇRÍK AND VLADIMIR LOTOREICHIK
ABSTRACT . We consider the problem of geometric optimisation of the lo west eigenvalue of the
Laplacian in the exterior of a compact planar set, subject to attractive Robin boundary conditions.
Under either a constraint of ﬁxed perimeter or area, we show t hat the maximiser within the class
of exteriors of convex sets is always the exterior of a disk. W e also argue why the results fail
without the convexity constraint and in higher dimensions.
1. I NTRODUCTION
The isoperimetric inequality states that among all planar sets of a given perimeter, the di sk has
the largest area. This is equivalent to the isochoric inequality stating that among all planar sets
of a given area, the disk has the smallest perimeter. These tw o classical geometric optimisation
problems were known to ancient Greeks, but a ﬁrst rigorous pro of appeared only in the 19th
century (see [5] for an overview).
Going from geometric to spectral quantities, the spectral isochoric inequality states that among
all planar membranes of a given area and with ﬁxed edges, the ci rcular membrane produces
the lowest fundamental tone. It was conjectured by Lord Rayl eigh in 1877 in his famous book
The theory of sound [22], but proved only by Faber [ 14] and Krahn [ 21] almost half a century
later. Using scaling, it is easily seen that this result impl ies the spectral isoperimetric inequality
as well: among all planar membranes of a given perimeter and w ith ﬁxed edges, the circular
membrane produces the lowest fundamental tone.
The same spectral inequalities under the area or perimeter c onstraints extend to elastically
supported membranes. To be more precise, let Ω ⊂ R2 be any smooth bounded open set.
Given a real number α, consider the spectral problem for the Laplacian, subject t o Robin
boundary conditions,
(1.1)



−∆u = λu in Ω ,
∂u
∂n + α u = 0 on ∂Ω ,
where n is the outer unit normal to Ω. It is well known that ( 1.1) induces an inﬁnite number
of eigenvalues λα
1 (Ω) ≤ λα
2 (Ω) ≤ λα
3 (Ω) ≤ . . . diverging to inﬁnity . The lowest eigenvalue
admits the variational characterisation
(1.2) λα
1 (Ω) = inf
u∈ W1,2(Ω)
u̸=0
∫
Ω
|∇ u|2 + α
∫
∂Ω
|u|2
∫
Ω
|u|2
,
from which it is clear that λα
1 (Ω) is non-negative if, and only if, α is non-negative. In this
elastic regime, the spectral isoperimetric and isochoric i nequalities for the Robin Laplacian can
Date: 26 September 2016.
2010 Mathematics Subject Classiﬁcation. 35P15 (primary); 58J50 (secondary).
Key words and phrases. Robin Laplacian, negative boundary parameter, exterior of a convex set, lowest eigen-
value, spectral isoperimetric inequality , spectral isochoric inequality , parallel coordinates.
1

<!-- page 2 -->
be respectively stated as follows
α ≥ 0 min
|∂Ω|=c1
λα
1 (Ω) = λα
1 (BR1) and min
|Ω|=c2
λα
1 (Ω) = λα
1 (BR2) .(1.3)
Here BR1 and BR2 are disks of perimeter |∂BR1| = c1 and area |BR2| = c2, respectively , with c1
and c2 being two arbitrary positive constants. With an abuse of not ation, we denote by |Ω|
and |∂Ω| the 2-dimensional Lebesgue measure of Ω and the 1-dimensional Hausdorff mea-
sure of its boundary ∂Ω, respectively . In the Dirichlet case (formally correspond ing to setting
α = ∞ ), the second optimisation problem in (
1.3) is just the Rayleigh-Faber-Krahn inequality
mentioned above. Both the statements in ( 1.3) are trivial for Neumann boundary conditions
(α = 0). If α > 0 , the second result in ( 1.3) is due to Bossel [ 6] from 1986 (and Daners [ 12] from
2006, also for higher dimensions), while the ﬁrst identity ag ain follows by scaling.
Summing up, going from the ancient isoperimetric inequalit y to the most recent spectral
results on the Robin problem with positive boundary couplin g parameter, the disk turns out
to be always the extremal set for the optimisation problems u nder the area or perimeter con-
straints. Moreover, the isoperimetric and isochoric optimi sation problems are closely related.
In the last two years, however, it has been noticed that the op timisation of the Robin eigen-
value in the case of negative α is quite different. Here the story goes back to 1977 when Bare ket
conjectured that a reverse spectral isochoric inequality should hold:
Conjecture 0 (Bareket’s reverse spectral isochoric inequality [3]). For all negative α, we have
α < 0 max
|Ω|=c2
λα
1 (Ω) = λα
1 (BR2) ,(1.4)
where the maximum is taken over all bounded open connected sets Ω of a given area c2 > 0 and BR2 is
the disk of the same area as Ω, i.e. |BR2| = c2.
Notice that, contrary to (
1.3), it is natural to maximise the eigenvalue if α is negative. The
conjecture was supported by its validity for sets close to th e disk [ 3, 15] and revived both
in [ 8] and [ 13]. However, Freitas and one of the present authors have recen tly disproved the
conjecture [16]: While ( 1.4) holds for small negative values of α, it cannot hold for all values
of the boundary parameter. In fact, the annulus provides a la rger value of the lowest Robin
eigenvalue as α → −∞ . This provides the ﬁrst known example where the extremal doma in
for the lowest eigenvalue of the Robin Laplacian is not a disk . On the other hand, in the most
recent publication [2], it was shown that the reverse spectral isoperimetric inequality does hold:
Theorem 0 (Reverse spectral isoperimetric inequality [2]). For all negative α, we have
α < 0 max
|∂Ω|=c1
λα
1 (Ω) = λα
1 (BR1) ,(1.5)
where the maximum is taken over all smooth bounded open connect ed sets Ω of a given perimeter c1 > 0
and BR1 is the disk of the same perimeter as Ω, i.e. |∂BR1| = c1.
Summing up, when α is negative, the disk is still the optimiser of the reverse spectral isoperi-
metric problem, while it stops to be the optimiser of the isoc horic problem for large negative
values of α. Whether the optimiser becomes the annulus for the larger ne gative values of α
and whether Conjecture
0 holds under some geometric restrictions on Ω represent just a few
hot open problems in the recent study (see [ 2, Sec. 5.3] for conjectures supported by numerical
experiments).
In this paper, we show that both the reverse spectral isochor ic and isoperimetric inequalities
hold in the dual setting of the Robin problem in the exterior of a convex set . To this purpose, for
any open set Ω ⊂ R2, we deﬁne Ωext := R2 \ Ω.
Theorem 1. For all negative α, we have
α < 0 max
|∂Ω|=c1
Ω convex
λα
1 (Ωext) = λα
1 (Bext
R1 ) and max
|Ω|=c2
Ω convex
λα
1 (Ωext) = λα
1 (Bext
R2 ) ,(1.6)
2

<!-- page 3 -->
where the maxima are taken over all convex smooth bounded open sets Ω of a given perimeter c1 > 0 or
area c2 > 0, respectively, and BR1 and BR2 are disks of perimeter |∂BR1| = c1 and area |BR2| = c2.
It is important to mention that λα
1 (Ωext) when deﬁned by (
1.2) indeed represents a (nega-
tive) discrete eigenvalue of a self-adjoint realisation of the Robin Lapla cian in L2(Ωext). It is not
obvious because there is also the essential spectrum [0, ∞ ), but it can be shown by using the
criticality of the Laplacian in R2 and the fact that α is negative ( cf. Section 2). In fact, λα
1 (Ωext)
equals zero (the lowest point in the essential spectrum) for any domain Ω if α is non-negative,
so the optimisation problems (
1.6) are trivial in this case. At the same time, because of the exi s-
tence of the Hardy inequality in higher dimensions, λα
1 (Ωext) is also zero for all small (negative)
values of α if the dimension is greater than or equal to three, so the iden tities (
1.6) become triv-
ial in this regime, too. This is the main reason why we mostly ( but not exclusively) restrict to
planar domains in this paper. As a matter of fact, in Section 5.3 we argue that an analogue of
Theorem 1 can not hold in space dimensions greater than or equal to thre e.
We point out that the disks are the only optimisers in Theorem 1 (see Section 5.2 for the
respective argument). It is worth mentioning that the ident ities (1.6) are no longer valid if the
condition of connectedness of Ω is dropped (see Section 5.1 for a counterexample). However,
it is still unclear at the moment whether the convexity of Ω in (1.6) can be replaced by a weaker
assumption.
The organisation of this paper is as follows. In Section 2 we provide an operator-theoretic
framework for the eigenvalue problem of Robin type ( 1.1) and establish basic spectral proper-
ties in the exterior of a compact set. Section 3 is devoted to more speciﬁc results on the lowest
eigenvalue in the case of the compact set being a disk. Theore m 1 is proved in Section 4: The
isoperimetric part of the theorem follows quite straightfo rwardly by using test functions with
level lines parallel to the boundary ∂Ω, while the isochoric result is established with help of
scaling and a monotonicity result of Section 3. The method of parallel coordinates was em-
ployed also in [ 16] to establish Conjecture 0 for small values of α, however, the reader will
notice that the present implementation of the technique is q uite different and in principal does
not restrict to planar sets ( cf. Section 5.4). The paper is concluded by Section 5 with more
comments on our results and methods.
2. T HE SPECTRAL PROBLEM IN THE EXTERIOR OF A COMPACT SET
Throughout this section, Ω is an arbitrary bounded open set in R2, not necessarily connected
or convex. However, a standing assumption is that the exterior Ωext is connected. Occasionally ,
we adopt the shorthand notation Σ := ∂Ω. For simplicity , we assume thatΩ is smooth (i.e. C∞-
smooth), but less regularity is evidently needed for the maj ority of the results to hold. At the
same time, α is an arbitrary real number, not necessarily negative (unle ss otherwise stated).
We are interested in the eigenvalue problem
(2.1)



−∆u = λu in Ωext ,
− ∂u
∂n + α u = 0 on ∂Ωext .
We recall that n is the outer unit normal to Ω, that is why we have the ﬂip of sign with respect
to ( 1.1). As usual, we understand ( 2.1) as the spectral problem for the self-adjoint operator
−∆Ωext
α in L2(Ωext) associated with the closed quadratic form
(2.2) QΩext
α [u] := ∥∇ u∥2
L2(Ωext) + α ∥u∥2
L2(Σ) , D(QΩext
α ) := W1,2(Ωext) .
The boundary term is understood in the sense of traces W1,2(Ωext)֒→L2(Σ) and represents a
relatively bounded perturbation of the Neumann form QΩext
0 with the relative bound equal to
zero. Since Ω is smooth, the operator domain of −∆Ωext
α consists of functions u ∈ W2,2(Ωext)
which satisfy the Robin boundary conditions of (
2.1) in the sense of traces and the operator
acts as the distributional Laplacian ( cf. [4, Thm. 3.5] for the W2,2-regularity). We call −∆Ωext
α the
Robin Laplacian in Ωext.
3

<!-- page 4 -->
Since Ω is bounded, the embedding W1,2(Ωext)֒→L2(Ωext) is not compact. In fact, the
Robin Laplacian possesses a non-empty essential spectrum w hich equals [0, ∞ ). This property
is expected because the (essential) spectrum of the Laplaci an in the whole space R2 (i.e. Ω =
∅) equals [0, ∞ ) and removing Ω can be understood as a compact perturbation. In order to
keep the paper self-contained, we provide a proof of this cla im which relies on an explicit
construction of singular sequences and a Neumann bracketin g argument.
Proposition 1. We have σess(−∆Ωext
α ) = [ 0, ∞ ) .
Proof. First, we show the inclusion σess(−∆Ωext
α ) ⊃ [0, ∞ ) by constructing a suitable singular
sequence. For any positive integer n, let us set un(x) := ϕn(x) eik·x with an arbitrary vector
k ∈ R2 and ϕn(x) := n−1ϕ((x − nx0)/n), where ϕ is a function from C∞
0 (R2) normalised to 1
in L2(R2) and x0 := (1, 0). The prefactor in the deﬁnition of ϕn is chosen in such a way that ϕn
is normalised to 1 in L2(R2) for each n. In fact, we have
(2.3) ∥ϕn∥L2(R2) = 1 , ∥∇ ϕn∥L2(R2) = n−1∥∇ ϕ∥L2(R2) , ∥∆ϕn∥L2(R2) = n−2∥∆ϕ∥L2(R2) .
At the same time, the support of ϕn leaves any bounded set for all sufﬁciently large n. Conse-
quently , for all sufﬁciently large n, we have un ∈ C∞
0 (Ωext) ⊂ D(−∆Ωext
α ) and ∥un∥L2(Ωext) = 1.
A direct computation yields
⏐
⏐
⏐−∆Ωext
α un − |k|2un
⏐
⏐
⏐ =
⏐
⏐
⏐(−∆ϕn + 2ik · ∇ ϕn) eik·x
⏐
⏐
⏐ ≤ |∆ϕn| + 2 |k||∇ ϕn| .
Using (
2.3), we therefore have

 − ∆Ωext
α un − |k|2un

2
L2(Ωext) ≤ 2∥∆ϕn∥2
L2(R2) + 8 |k|2∥∇ ϕn∥2
L2(R2) − − −→
n→∞
0 .
Since k is arbitrary , we conclude that [0, ∞ ) ⊂ σ(−∆Ωext
α ) by [
24, Thm. 7.22]. It is clear that
[0, ∞ ) actually belongs to the essential spectrum, because the interval has no isolated points.
Second, to show the opposite inclusion σess(−∆Ωext
α ) ⊂ [0, ∞ ), we use a Neumann bracketing
argument. Let Hn be the operator that acts as −∆Ωext
α but satisﬁes an extra Neumann condition
on the circle Cn := {x ∈ R2 : |x| = n} of radius n > 0 . More speciﬁcally , Hn is the operator
associated with the form
hn[u] := ∥∇ u∥2
L2(Ωext) + α ∥u∥2
L2(Σ) , D(hn) := W1,2(Ωext \ Cn) .
Because of the domain inclusion D(hn) ⊃ D(QΩext
α ), we have −∆Ωext
α ≥ Hn and, by the min-
max principle, inf σess(−∆Ωext
α ) ≥ inf σess(Hn) for all n. Assuming that n is sufﬁciently large so
that Ω is contained in the disk Bn := {x ∈ R2 : |x| < n }, Hn decouples into an orthogonal sum
of two operators, Hn = H(1)
n ⊕ H(2)
n with respect to the decomposition L2(Ωext) = L2(Ωext ∩
Bn) ⊕ L2(Ωext \
Bn). Here H(1)
n and H(2)
n are respectively the operators in L2(Ωext ∩ Bn) and
L2(Ωext \ Bn) associated with the forms
h(1)
n [u] := ∥∇ u∥2
L2(Ωext∩ Bn) + α ∥u∥2
L2(Σ) , D(h(1)
n ) := W1,2(Ωext ∩ Bn) ,
h(2)
n [u] := ∥∇ u∥2
L2(Ωext\
Bn) , D(h(2)
n ) := W1,2(Ωext \ Bn) .
Since Ωext ∩ Bn is a smooth bounded open set, the spectrum of H(1)
n is purely discrete. Conse-
quently , infσess(−∆Ωext
α ) ≥ inf σess(H(2)
n ) ≥ inf σ(H(2)
n ) ≥ 0, where the last inequality follows by
the fact that H(2)
n is non-negative. □
Despite of the presence of essential spectrum, it still make s sense to deﬁne the lowest point
in the spectrum of −∆Ωext
α by the variational formula ( cf. (1.2))
(2.4) λα
1 (Ωext) := inf
u∈ W1,2(Ωext)
u̸=0
QΩext
α [u]
∥u∥2
L2(Ωext)
.
4

<!-- page 5 -->
However, it is not evident that it represents a discrete eige nvalue of −∆Ωext
α . Obviously , it is not
the case if α is non-negative, in which case −∆Ωext
α is non-negative and therefore its spectrum
is purely essential. The following result shows that the sit uation of negative α is different.
Proposition 2. If α < 0 and Ω is not empty, then σdisc(−∆Ωext
α ) ̸= ∅ . More speciﬁcally, λα
1 (Ωext) is
a negative discrete eigenvalue.
Proof. By Proposition
1 and (2.4), it is enough to ﬁnd a test function u ∈ W1,2(Ωext) such that
QΩext
α [u] is negative. For any positive number n, we introduce a function un : R2 → [0, 1] by
setting un(x) := ϕn(|x|) with
ϕn(r) :=









1 if r < n ,
log n2 − log r
log n2 − log n if n < r < n 2,
0 otherwise .
It is not difﬁcult to check that the restriction of un to Ωext (that we shall denote by the same
symbol) belongs to W1,2(Ωext) for every n. By employing polar coordinates, we have
∥∇ un∥2
L2(Ωext) ≤ ∥∇ un∥2
L2(R2) = 2π
∫ ∞
0
|ϕ ′
n(r)|2 r dr = 2π
∫ n2
n
1
(log n)2 r dr = 2π
log n − − −→
n→∞
0 .
On the other hand, ∥un∥2
L2(Σ) = |Σ| > 0 for all sufﬁciently large n. Since α is assumed to be
negative, it follows that QΩext
α [un] can be made negative for n large enough. □
Summing up, if α is negative, the essential spectrum of −∆Ωext
α equals the interval [0, ∞ )
and there is at least one discrete eigenvalue below 0. In particular, the lowest point λα
1 (Ωext)
in the spectrum is always a negative discrete eigenvalue. By standard methods (see, e.g., [
17,
Thm. 8.38]), it is possible to show that λα
1 (Ωext) is simple and that the corresponding eigen-
function uα
1 can be chosen to be positive in Ωext (recall that we always assume that Ωext is
connected).
It is straightforward to verify that {QΩext
α }α∈ R is a holomorphic family of forms of type (a) in
the sense of Kato [
18, Sec. VII.4]. In fact, recalling that the boundary term in ( 2.2) is relatively
bounded with respect to the Neumann form QΩext
0 with the relative bound equal to zero, one
can use [ 18, Thm. 4.8] to get the claim. Consequently , −∆Ωext
α forms a self-adjoint holomorphic
family of operators of type (B). Because of the simplicity , i t follows that α ↦→ λα
1 (Ωext) and
α ↦→ uα
1 with ∥uα
1 ∥ = 1 are real-analytic functions on (−∞ , 0).
Proposition 3. Let α < 0 and Ω ̸= ∅. Then α ↦→ λα
1 (Ωext) is a strictly concave increasing function.
Proof. For simplicity , let us set λα
1 := λα
1 (Ωext). The eigenvalue equation −∆Ωext
α uα
1 = λα
1 uα
1
means that
(2.5) ∀ϕ ∈ W1,2(Ωext) , Q Ωext
α (ϕ, uα
1 ) = λα
1 (ϕ, uα
1 )L2(Ωext) .
Differentiating the identity (
2.5) with respect to α, we easily arrive at the formula
(2.6) (∇ ϕ, ∇ ˙uα
1 )L2(Ωext) + (ϕ, uα
1 )L2(Σ) + α(ϕ, ˙uα
1 )L2(Σ) = ˙λα
1 (ϕ, uα
1 )L2(Ωext) + λα
1 (ϕ, ˙uα
1 )L2(Ωext) ,
where the dot denotes the derivative with respect to α. Notice that the differentiation below
the integral signs is permitted because ˙ uα
1 ∈ W1,2(Ωext) by standard elliptic regularity the-
ory . Moreover, differentiating the normalisation condition ∥uα
1 ∥ = 1, we get the orthogonality
property
(2.7) (uα
1 , ˙uα
1 )L2(Ωext) = 0 .
Now, substituting ϕ = uα
1 into (
2.6) and ϕ = ˙uα
1 into (
2.5) and taking the difference of the
resulting equations, we get a formula for the eigenvalue der ivative
(2.8) ˙λα
1 = ∥uα
1 ∥2
L2(Σ) > 0.
5

<!-- page 6 -->
The above inequality is strict because otherwise uα
1 would be an eigenfunction of the Dirichlet
Laplacian on Ωext corresponding to a negative eigenvalue λα
1 , which is a contradiction to the
non-negativity of the latter operator. This proves that α ↦→ λα
1 is strictly increasing.
Next, we differentiate equation (
2.8) with respect to α,
(2.9) ¨λα
1 = d
dα
(
∥uα
1 ∥2
L2(Σ)
)
= 2(uα
1 , ˙uα
1 )L2(Σ) = 2λα
1 ∥ ˙uα
1 ∥2
L2(Ωext) − 2QΩext
α [ ˙uα
1 ] < 0 .
Here the last equality employs (
2.6) with the choice ϕ = ˙uα
1 and (
2.7). The inequality follows
from the fact that λα
1 is the lowest eigenvalue of −∆Ωext
α . The above inequality is indeed strict
since otherwise ˙uα
1 would be either another eigenfunction of −∆Ωext
α corresponding to λα
1 , which
is impossible because of the simplicity , or a constant multiple of uα
1 , which would imply ˙uα
1 = 0
due to (
2.7). In the latter case ( 2.6) gives
∀ϕ ∈ C∞
0 (Ωext) , (ϕ, uα
1 )L2(Ωext) = 0 ,
and therefore uα
1 = 0, which is also a contradiction. From (
2.9) we therefore conclude that
α ↦→ λα
1 is strictly concave. □
As a consequence of Proposition
3, we get
(2.10) lim
α→−∞
λα
1 (Ωext) = − ∞ .
3. T HE LOWEST EIGENVALUE IN THE EXTERIOR OF A DISK
In this section, we establish some properties of λα
1 (Bext
R ), where BR is an open disk of radius
R > 0 . Without loss of generality , we can assume that BR is centred at the origin of R2. We
always assume that α is negative.
Using the rotational symmetry , it is easily seen that λα
1 (Bext
R ) = − k2 < 0 is the smallest
solution of the ordinary differential spectral problem
(3.1)







−r−1[rψ ′(r)] ′ = λψ(r) , r ∈ (R, ∞ ) ,
−ψ ′(R) + α ψ(R) = 0 ,
lim
r→∞
ψ(r) = 0 .
The general solution of the differential equation in (
3.1) is given by
(3.2) ψ(r) = C1K0(kr) + C2I0(kr) , C 1, C2 ∈ C ,
where K0, I0 are modiﬁed Bessel functions [ 1, Sec. 9.6]. The solution I0(kr) is excluded because
it diverges as r → ∞ , whence C2 = 0. Requiring ψ to satisfy the Robin boundary condition
at R leads us to the implicit equation
(3.3) kK1(kR) + αK0(kR) = 0
that k must satisfy as a function of α and R.
First of all, we state the following upper and lower bounds.
Proposition 4. We have
−α2 < λ α
1 (Bext
R ) < −α2 − α
R
for all negative α.
Proof. For simplicity , let us set λR := λα
1 (Bext
R ) and recall the notation λR = −k2. Using (
3.3) and
[23, Eq. 74] (with ν = 0), we have
λR = −k2 = −α2
( K0(kR)
K1(kR)
) 2
> −α2.
6

<!-- page 7 -->
This establishes the lower bound of the proposition. In the c ase α ∈ (−R−1, 0), we obtain the
upper bound from the elementary estimate
λR + α2 + α
R < α 2 + α
R = α
(
α + 1
R
)
< 0 ,
where we have used the fact that λR is negative. In the other case α ∈ (−∞ , −R−1], we get by
[23, Thm. 1] (with ν = 1/2) that
k2 = α2
( K0(kR)
K1(kR)
) 2
> α2(kR)2
1/2 + (kR)2 +
√
1/4 + (kR)2 .
The latter inequality implies
(3.4)
√
1
4R4 + k2
R2 > α 2 − 1
2R2 − k2 .
If the right-hand side of ( 3.4) is negative, then
−k2 < −α2 + 1
2R2 = −α2 + 1
R
1
2R ≤ −α2 − α
2R < −α2 − α
R ,
which is the desired inequality . If the right-hand side of ( 3.4) is non-negative, we take the
squares of both the right- and left-hand sides of ( 3.4) and obtain
1
4R4 + k2
R2 >
(
α2 − 1
2R2 − k2
) 2
=
(
α2 − 1
2R2
) 2
− 2k2
(
α2 − 1
2R2
)
+ k4 .
This inequality is equivalent to
0 > (α2 − k2)2 − α2
R2 .
Consequently ,
α2 − k2 < − α
R ,
which again yields the desired upper bound. □
We notice that analogous upper and lower bounds for λα
1 (BR) have been recently established
in [
2, Thm. 3]. Moreover, it has been shown in [ 2, Thm. 5] that R ↦→ λα
1 (BR) is strictly increasing.
Now we have a reversed monotonicity result.
Proposition 5. If α is negative, then R ↦→ λα
1 (Bext
R ) is strictly decreasing.
Proof. We follow the strategy of the proof of [
2, Thm. 5]. Computing the derivative of λR :=
λα
1 (Bext
R ) using the differential equation that λR satisﬁes, one ﬁnds ( cf. [
2, Lem. 2])
(3.5) ∂λR
∂R = − 2
R λR + α ψR(R)2
∫ ∞
R
ψR(r)2 r dr
,
where ψR(r) := K0(kr) is the eigenfunction corresponding to λR = −k2. Employing the formula
(3.6)
∫ ∞
R
K0(kr)2 r dr = r2
2
[
K0(kr)2 − K1(kr)2
] ⏐
⏐
⏐
⏐
r=∞
r=R
= R2
2
[
K1(kR)2 − K0(kR)2
]
and (3.3), we eventually arrive at the equivalent identity for the ei genvalue derivative
∂λR
∂R = − 2
R λR
λR + α2 + α
R
λR + α2 .(3.7)
The proof is concluded by recalling Proposition
4, which implies that the right-hand side is
negative. □
7

<!-- page 8 -->
4. P ROOF OF THEOREM 1
Now we are in a position to establish Theorem 1. Throughout this section, Ω ⊂ R2 is a
convex bounded open set with smooth boundary Σ := ∂Ω; then Ωext is necessarily connected.
We also assume that α is negative.
The main idea of the proof is to parameterise Ωext by means of the parallel coordinates
(4.1) L : Σ × (0, ∞ ) → Ωext : {(s, t) ↦→ s + n(s) t} ,
where n is the outer unit normal to Ω as above. Notice that L is indeed a diffeomorphism
because of the convexity and smoothness assumptions. To be more speciﬁc, the metric induced
by (4.1) acquires the diagonal form
(4.2) d L2 = (1 − κ(s) t)2 ds2 + dt2 ,
where κ := − dn is the curvature of ∂Ω. By our choice of n, the function κ is non-positive be-
cause Ω is convex (cf. [20, Thm. 2.3.2]). Consequently , the Jacobian of (4.1) given by 1 − κ(s) t is
greater than or equal to 1. In particular, it is positive and therefore L is a local diffeomorphism
by the inverse function theorem (see, e.g., [20, Thm. 0.5.1]). To see that it is a global diffeomor-
phism, notice that L is injective, because of the convexity assumption, and that L is surjective
thanks to the smoothness of Ω.
Summing up, Ωext can be identiﬁed with the product manifold Σ × (0, ∞ ) equipped with
the metric ( 4.2). Consequently , the Hilbert space L2(Ωext) can be identiﬁed with
H := L2(
Σ × (0, ∞ ), (1 − κ(s) t) ds dt
)
.
The identiﬁcation is provided by the unitary transform
(4.3) U : L2(Ωext) → H : {u ↦→ u ◦ L} .
It is thus natural to introduce the unitarily equivalent ope rator Hα := U(−∆Ωext
α )U−1, which is
the operator associated with the transformed form hα[ψ] := QΩext
α [U−1ψ], D(hα) := UD(QΩext
α ).
Of course, we have the equivalent characterisation of the lo west eigenvalue
(4.4) λα
1 (Ωext) = inf
ψ∈ D(hα)
ψ̸=0
hα[ψ]
∥ψ∥2
H
.
The set of restrictions of functions from C∞
0 (R2) to Ωext is a core of QΩext
α . Taking u from this
core, it is easily seen that ψ := Uu is a restriction of a function C∞
0 (∂Ω × R) to ∂Ω × (0, ∞ ) and
that
hα[ψ] =
∫
Σ× (0,∞)
(
|∂sψ(s, t)|2
1 − κ(s) t + |∂tψ(s, t)|2 (1 − κ(s) t)
)
ds dt + α
∫
Σ
|ψ(s, 0)|2 ds ,
∥ψ∥2
H =
∫
Σ× (0,∞)
|ψ(s, t)|2 (1 − κ(s) t) ds dt .
Restricting in (
4.4) to test functions with level lines parallel to Σ, i.e. taking ψ independent of s,
we obtain
(4.5) λα
1 (Ωext) ≤ inf
ψ∈ C∞
0 ([0,∞))
ψ̸=0
∫ ∞
0
|ψ ′(t)|2 (|Σ| + 2πt) dt + α |Σ| |ψ(0)|2
∫ ∞
0
|ψ(t)|2 (|Σ| + 2πt) dt
.
Here we have used the geometric identity
∫
Σ κ = −2π (see, e.g., [20, Cor. 2.2.2]).
Now, assume that the perimeter is ﬁxed, i.e. |Σ| = c1. Since the perimeter is the only geo-
metric quantity on which the right-hand side of ( 4.5) depends and since the eigenfunction
corresponding to λα
1 (Bext
R1
) is radially symmetric (therefore independent of s in the parallel co-
ordinates), we immediately obtain
(4.6) λα
1 (Ωext) ≤ λα
1 (Bext
R1 )
8

<!-- page 9 -->
for any Ω of the ﬁxed perimeter c1. This proves the isoperimetric optimisation result of Theo -
rem 1.
To establish the isochoric optimisation result of Theorem 1, we notice that since |Σ| = |∂BR1|,
the classical (geometric) isoperimetric inequality impli es |Ω| ≤ |BR1|, with equality if and only
if Ω is the disk. Hence, there exists BR2 ⊂ BR1 such that |BR2| = |Ω|. By ( 4.6), it is then enough to
show that λα
1 (Bext
R1
) ≤ λα
1 (Bext
R2
). This inequality (and therefore the second result of Theore m
1)
is a consequence of the more general monotonicity result of P roposition 5. □
5. C ONCLUSIONS
Let us conclude the paper by several comments on our results.
5.1. Necessity of convexity. The assumption on Ω to be convex is necessary in view of the
following simple counterexample. Let Ω ⊂ R2 be the union of two disks BR3(x1) and BR3(x2)
of the same radius R3 > 0 whose centres x1 and x2 are chosen in such a way that the closures
of the disks in R2 are disjoint. According to [ 19, Thm. 1.1] (see also [ 16, Thm. 3]), we have
λα
1 (Ωext) = − α2 − α
R3
+ o(α) , α → −∞ ,
λα
1 (Bext
R ) = − α2 − α
R + o(α) , α → −∞ .
The constraints |∂Ω| = |∂BR1| and |Ω| = |BR2| yield that R1 = 2R3 and R2 =
√
2R3, respectively .
Taking into account the above large coupling asymptotics an d the relations between the radii,
we observe that for α < 0 with sufﬁciently large |α| the reverse inequalities λα
1 (Bext
R1
) < λ α
1 (Ωext)
and λα
1 (Bext
R2
) < λ α
1 (Ωext) are satisﬁed.
We point out that, while the domain Ω of the above counterexample is disconnected, its
exterior Ωext is still connected. We leave it as an open question whether th ere exists a coun-
terexample in the class of connected non-convex domains Ω.
5.2. Uniqueness of the optimiser. In this subsection, we demonstrate that the exterior of the
disk is the unique maximiser for both the isochoric and isoperimetric spectra l optimisation
problems of Theorem
1.
Theorem 2. Let α be negative. For all convex smooth bounded open sets Ω ⊂ R2 of a ﬁxed perimeter
(respectively, ﬁxed area) different from a disk BR1 of the same perimeter (respectively, from a disk BR2 of
the same area), we have a strict inequality
α < 0 λα
1 (Ωext) < λ α
1 (Bext
R1 ) (respectively, λα
1 (Ωext) < λ α
1 (Bext
R2 )) .(5.1)
Proof. As usual, λα
1 := λα
1 (Ωext) and uα
1 denote respectively the lowest eigenvalue and the cor-
responding eigenfunction of −∆Ωext
α with α < 0 . Without loss of generality , we assume that
uα
1 is positive everywhere in Ωext. Furthermore, we introduce the auxiliary function ψ := Uuα
1
where the unitary transform U is as in (
4.3). In view of Theorem 1 and inequality ( 4.5) in its
proof, non-uniqueness of the optimiser for the spectral iso perimetric problem would necessar-
ily imply the existence of a non-circular domain Ω for which the function ψ is independent
of s; i.e. its level lines are parallel to Σ := ∂Ω. In the sequel, with a slight abuse of notation, we
use the same symbol ψ to denote the function t ↦→ ψ(t) of a single variable.
Restricting to test functions with support lying inside Ωext, the variational characterisation
of ψ implies
∀ϕ ∈ C∞
0 (Σ × (0, ∞ )) , h α(ϕ, ψ) = λα
1 (ϕ, ψ)H .
It is enough to consider real-valued test functions ϕ only . Taking into account that ψ is inde-
pendent of s, we end up with the identity
∫
Σ× (0,∞)
∂tϕ(s, t) ψ ′(t) (1 − κ(s)t) ds dt = λα
1
∫
Σ× (0,∞)
ϕ(s, t) ψ(t) (1 − κ(s)t) ds dt
9

<!-- page 10 -->
valid for all real-valued ϕ ∈ C∞
0 (Σ × (0, ∞ )). Now we restrict our attention to test functions
of the type ϕ(s, t) = ϕ1(s)ϕ2(t) ∈ C∞
0 (Ωext) with ϕ1 ∈ C∞(Σ) and ϕ2 ∈ C∞
0 ((0, ∞ )). Then the
above displayed equation reduces to
∫
Σ
ϕ1(s)
∫ ∞
0
ϕ ′
2(t) ψ ′(t) (1 − κ(s)t) dt ds = λα
1
∫
Σ
ϕ1(s)
∫ ∞
0
ϕ2(t) ψ(t) (1 − κ(s)t) dt ds.
Density of C∞(Σ) in L1(Σ) gives us
(5.2)
∀s ∈ Σ , ϕ 2 ∈ C∞
0 ((0, ∞ )) ,
∫ ∞
0
ϕ′
2(t) ψ′(t) (1 − κ(s)t) dt = λα
1
∫ ∞
0
ϕ2(t) ψ(t) (1 − κ(s)t) dt.
Since Ω is not a disk, there exist s1, s2 ∈ Σ such that κ(s1) ̸= κ(s2) (see, e.g., [
20, Prop. 1.4.3]).
Taking the difference of (5.2) with s = s1 and with s = s2, we eventually get
(5.3) ∀ϕ2 ∈ C∞
0 ((0, ∞ )) ,
∫ ∞
0
ϕ′
2(t) ψ′(t) t dt = λα
1
∫ ∞
0
ϕ2(t) ψ(t) t dt .
Let us ﬁx a function η ∈ C∞
0 ((0, ∞ )) which satisﬁes the following properties:
(i) 0 ≤ η ≤ 1,
(ii) η(t) = 1 for all t ∈ [1, 2],
(iii) supp η ⊂ [0, 3].
Furthermore, for every positive integer n, we deﬁne a function ηn ∈ C∞
0 ((0, ∞ )) by
(5.4) ηn(t) :=





η(nt), t ∈
(
0, 2
n
)
,
1, t ∈
( 2
n , n + 1
)
,
η(t − n), t ∈ (n + 1, ∞ ) .
Now we plug ϕ2 = ηnψ ∈ C∞
0 ((0, ∞ )) into (
5.3). By the dominated convergence theorem
(using that t ↦→ |ψ(t)|2 t is integrable), we obtain
(5.5)
∫ ∞
0
|ψ(t)|2 ηn(t) t dt − − −→
n→∞
∫ ∞
0
|ψ(t)|2 t dt .
The left-hand side in ( 5.3) with ϕ2 = ηnψ can be rewritten as
(5.6) In :=
∫ ∞
0
(ηnψ)′(t) ψ′(t) t dt =
∫ ∞
0
|ψ′(t)|2 ηn(t) t dt +
∫ ∞
0
ψ′(t) ψ(t) η′
n(t) t dt .
For the ﬁrst term on the right-hand side in (
5.6) we get
(5.7) I(1)
n :=
∫ ∞
0
|ψ′(t)|2 ηn(t) t dt − − −→
n→∞
∫ ∞
0
|ψ′(t)|2 t dt ,
by the dominated convergence (using that t ↦→ |ψ′(t)|2 t is integrable). The second term on the
right-hand side in ( 5.6) can be further transformed as
∫ ∞
0
ψ′(t) ψ(t) η′
n(t) t dt = n
∫ 2/n
0
ψ′(t) ψ(t) η′(nt) t dt +
∫ ∞
n+1
ψ′(t) ψ(t) η′(t − n) t dt
=
∫ 2
0
ψ′
( r
n
)
ψ
( r
n
) η′(r)r
n dr +
∫ 3
1
ψ′(r+ n) ψ(r+ n) η′(r) (r+ n) dr.
(5.8)
Again making use of the dominated convergence theorem, we ob tain
(5.9) I(2)
n :=
∫ 2
0
ψ′
( r
n
)
ψ
( r
n
) η′(r) r
n dr − − −→
n→∞
0 ;
here we implicitly employed that the integrand is uniformly bounded in n ∈ N. Observe that
⏐
⏐
⏐
⏐
⏐
∞∑
n=1
∫ 3
1
ψ′(r + n) ψ(r + n) η′(r) (r + n) dr
⏐
⏐
⏐
⏐
⏐ ≤ 2 ∥η′∥∞
∫ ∞
0
|ψ(t) ψ′(t)| t dt < ∞ ,
10

<!-- page 11 -->
where ﬁniteness of the latter integral follows from the fact t hat ψ ∈ D(hα). Therefore, we infer
(5.10) I(3)
n :=
∫ 3
1
ψ′(r + n) ψ(r + n) η′(r) (r + n) dr − − −→
n→∞
0 .
Combining the decompositions ( 5.6), (5.8) with the limits ( 5.7), (5.9), (5.10), we arrive at
(5.11) In = I(1)
n + I(2)
n + I(3)
n − − −→
n→∞
∫ ∞
0
|ψ′(t)|2 t dt .
The limits (5.5), (5.11) and the condition ( 5.3) imply
∫ ∞
0
|ψ′(t)|2 t dt = λα
1
∫ ∞
0
|ψ(t)|2 t dt .
Finally , taking into account that λα
1 is negative, we get a contradiction. This completes the
proof of the ﬁrst strict inequality in (
5.1).
To show that disk is the unique optimiser for the spectral iso choric inequality is much sim-
pler than in the isoperimetric case. Suppose that there exists a non-circular domain Ω for which
λα
1 (Ωext) = λα
1 (Bext
R2
) with |Ω| = |BR2|. Note that for BR1 with |∂BR1| = |∂Ω| we get R2 < R 1 using
the standard geometric isoperimetric inequality . Thus, Theorem
1 and Proposition 5 imply
λα
1 (Ωext) = λα
1 (Bext
R2 ) > λ α
1 (Bext
R1 ) ≥ λα
1 (Ωext) ,
which is obviously a contradiction. □
Remark 5.1. As a consequence of the ﬁrst claim in Theorem
1, we obtain the following quanti-
tative improvement upon the second inequality of ( 5.1)
(5.12) λα
1 (Ωext) ≤ λα
1 (Bext
R2 ) −
[
λα
1 (Bext
R2 ) − λα
1 (Bext
R1 )
]
,
where the radii R1 and R2 can be easily expressed through |∂Ω| and |Ω|, respectively , by virtue
of the relations |BR1| = |∂Ω| and |BR2| = |Ω|. In view of the inequality R2 < R 1, the difference
λα
1 (Bext
R2
) −λα
1 (Bext
R1
) is positive by Proposition
5, so (5.12) indeed represents a quantiﬁed version
of the reverse spectral isochoric inequality in the spirit o f [ 7, 9]. More careful analysis of the
derivative in ( 3.7) can be used to get a positive lower bound on this difference i n terms of R1
and R2.
5.3. Higher dimensions. We have already noticed that λα
1 (Ωext) is not necessarily a discrete
eigenvalue in higher dimensions. For any dimension d ≥ 3, however, there exists a critical
value α0 < 0 depending on Ω such that λα
1 (Ωext) is a discrete eigenvalue if, and only if, α < α 0,
so the optimisation problem in the exterior of a compact set b ecomes non-trivial in this regime.
In this subsection, we argue that no analogue of Theorem
1 can be expected if d ≥ 3.
To this aim we construct a simple counterexample which relies on the large coupling asymp-
totics for the lowest eigenvalue. First, we ﬁx a ball BR ⊂ Rd of arbitrary radius R > 0 . Further,
let Ω0 ⊂ Rd be the union of two disjoint balls Br(x1) and Br(x2) of the same sufﬁciently small
radius r > 0 whose centers x1 and x2 are located at a distance L > 0 . Finally , we deﬁne the
domain Ω as the convex hull of Ω0. By choosing L > 0 large enough, we can satisfy either of
the constraints |∂Ω| = |∂BR| or |Ω| = |BR|. It can be easily checked that the domain Ω has a
C1,1 boundary and that the mean curvature of ∂Ω is piecewise constant, being equal to −1/r
on the hemispheric cups and to − d−2
(d−1)r on the cylindrical face (in agreement with the rest of
this paper, we compute the mean curvature with respect to the outer normal to the bounded
set Ω). Applying [ 19, Thm. 1.1], we arrive at
λα
1 (Ωext) = − α2 − α (d − 2)
r + o(α) , α → −∞ ,
λα
1 (Bext
R ) = − α2 − α(d − 1)
R + o(α) , α → −∞ .
In view of the above asymptotics, we infer that for r < d−2
d−1 R and for α < 0 with sufﬁciently
large |α| the reverse inequality λα
1 (Bext
R ) < λ α
1 (Ωext) holds.
11

<!-- page 12 -->
We expect that a counterexample based on a ( C∞-)smooth domain can also be constructed
with additional technical efforts.
5.4. More on dimension three. The previous subsection demonstrates that, contrary to the
two-dimensional situation, the exterior of the ball can be a global maximiser neither for the
isoperimetric nor isochoric problems. Let us look at where t he technical approach of the
present paper fails in dimension three.
Let Ω be a convex smooth bounded open set in R3. In this case, the usage of parallel coordi-
nates based on Σ := ∂Ω and restricting to test functions depending only on the dist ance to the
boundary yield
(5.13) λα
1 (Ωext) ≤ inf
ψ∈ C∞
0 ([0,∞))
ψ̸=0
∫
Σ× (0,∞)
|ψ ′(t)|2 (1 − 2M(s) t + K(s) t2) dΣ dt + α
∫
Σ
|ψ(0)|2 dΣ
∫
Σ× (0,∞)
|ψ(t)|2 (1 − 2M(s) t + K(s) t2) dΣ dt
.
Here dΣ := |g|1/2(s) ds is the surface measure of Σ, with g being the Riemannian metric of Σ
induced by the embedding of Σ in R3, and K and M denote respectively the Gauss curvature
and the mean curvature of Σ (see [ 11] for more geometric details). K is an intrinsic quantity ,
while M is non-positive when computed with respect to our choice (ou ter to Ω) of the normal
vector ﬁeld n.
By deﬁnition,
∫
Σ 1 dΣ equals the total area |Σ| of Σ, while
∫
Σ K(s) dΣ = 4π by the Gauss-
Bonnet theorem for closed surfaces diffeomorphic to the sph ere (see [ 20, Thm. 6.3.5]). The
quantity MΣ :=
∫
Σ |M(s)| dΣ is known as the half of the total mean curvature of Σ (see [ 10,
§ 28.1.3]). Moreover, we have [ 10, § 19] MΣ = 2π b(Ω), where b(Ω) is the mean width of Ω.
Consequently ,
(5.14) λα
1 (Ωext) ≤ inf
ψ∈ C∞
0 ([0,∞))
ψ̸=0
∫ ∞
0
|ψ ′(t)|2 (|Σ| + 2 MΣ t + 4π t2) dt + α |Σ| |ψ(0)|2
∫ ∞
0
|ψ(t)|2 (|Σ| + 2 MΣ t + 4π t2) dt
.
To get a reverse spectral isoperimetric inequality as in the planar case above, we would need
in addition to the constraint |∂Ω| = c1 also require that the mean width b(Ω) is ﬁxed. However,
the Minkowski quadratic inequality for cross-sectional mea sures (cf. [10, § 20.2]), M2
Σ ≥ 4π |Σ|,
with equality only if Ω is a ball, implies that the two simultaneous constraints are possible only
if either the class of admissible domains excludes the ball o r the class of admissible domains
consists of the ball only . In the ﬁrst case our method is not app licable, while in the second case
the method can be applied but it yields a trivial statement.
Acknowledgments. D.K. was partially supported by FCT (Portugal) through proje ct PTDC/-
MAT-CAL/4334/2014. V .L. acknowledges the support by the grant No. 17-01706S of the Czech
Science Foundation (GA ˇCR).
REFERENCES
1. M. S. Abramowitz and I. A. Stegun, eds., Handbook of mathematical functions , Dover, New York, 1965.
2. P . R. S. Antunes, P . Freitas, and D. Krejˇ ciˇ rík, Bounds and extremal domains for Robin eigenvalues with nega tive
boundary parameter, Adv . Calc. Var., to appear; preprint on
arXiv:1605.08161 [math.SP]; doi:10.1515/acv-2015-
0045.
3. M. Bareket, On an isoperimetric inequality for the ﬁrst eigenvalue of a b oundary value problem, SIAM J. Math. Anal. 8
(1977), 280–287.
4. J. Behrndt, M. Langer, V . Lotoreichik, and J. Rohleder, Quasi boundary triples and semibounded self-
adjoint extensions , Proc. Roy . Soc. Edinburgh Sect. A., to appear; preprint on arXiv:1504.03885 [math.SP];
doi:10.1017/S0308210516000421.
5. V . Blåsjö, The isoperimetric problem, Am. Math. Mon. 112 (2005), 526–566.
6. M.-H. Bossel, Membranes élastiquement liées: Extension du théoréme de Ra yleigh-Faber-Krahn et de l’inégalité de
Cheeger, C. R. Acad. Sci. Paris Sér. I Math. 302 (1986), 47–50.
12

<!-- page 13 -->
7. L. Brasco and A. Pratelli, Sharp stability of some spectral inequalities , Geom. Funct. Anal. 22 (2012), 107–135.
8. F. Brock and D. Daners, Conjecture concerning a Faber-Krahn inequality for Robin p roblems, Oberwolfach Rep. 4
(2007), 1022–1023, Open Problem in Mini-Workshop: Shape An alysis for Eigenvalues (Organized by: D. Bucur,
G. Buttazzo, A. Henrot).
9. D. Bucur, V . Ferone, C. Nitsch, and C. Trombetti, The quant itative Faber-Krahn inequality for the Robin Lapla-
cian, preprint on arXiv:1611.06704 [math.AP].
10. Yu. D. Burago and V . A. Zalgaller, Geometric inequalities, Springer-V erlag, Berlin Heidelberg, 1988.
11. G. Carron, P . Exner, and D. Krejˇ ciˇ rík,T opologically nontrivial quantum layers, J. Math. Phys. 45 (2004), 774–784.
12. D. Daners, A Faber-Krahn inequality for Robin problems in any space dim ension, Math. Ann. 335 (2006), 767–785.
13. , Principal eigenvalues for generalised indeﬁnite Robin pro blems, Pot. Anal. 38 (2013), 1047–1069.
14. G. Faber, Beweis dass unter allen homogenen Membranen von gleicher Fl äche und gleicher Spannung die kreisförmige
den tiefsten Grundton gibt , Sitz. bayer. Akad. Wiss. (1923), 169–172.
15. V . Ferone, C. Nitsch, and C. Trombetti, On a conjectured reversed Faber-Krahn inequality for a Stek lov-type Laplacian
eigenvalue, Commun. Pure Appl. Anal. 14 (2015), 63–81.
16. P . Freitas and D. Krejˇ ciˇ rík,The ﬁrst Robin eigenvalue with negative boundary parameter , Adv . Math. 280 (2015),
322–339.
17. D. Gilbarg and N. S. Trudinger, Elliptic partial differential equations of second order , Springer-V erlag, Berlin, 1983.
18. T. Kato, Perturbation theory for linear operators , Springer-V erlag, Berlin, 1966.
19. H. Kovaˇ rík and K. Pankrashkin, On the p-Laplacian with Robin boundary conditions and boundary tra ce theorems,
Calc. Var. PDE 56 (2017), 49.
20. W. Klingenberg, A course in differential geometry , Springer-V erlag, New York, 1978.
21. E. Krahn, Über eine von Rayleigh formulierte Minimaleigenschaft des Kreises, Math. Ann. 94 (1924), 97–100.
22. J. W. S. Rayleigh, The theory of sound, Macmillan, London, 1877, 1st edition (reprinted: Dover, N ew York (1945)).
23. J. Segura, Bounds for ratios of modiﬁed Bessel functions and associate d Turán-type inequalities , J. Math. Anal. Appl.
374 (2011), 516–528.
24. J. Weidmann, Linear operators in Hilbert spaces , Springer-V erlag, New York Inc., 1980.
DEPARTMENT OF MATHEMATICS , F ACULTY OF NUCLEAR SCIENCES AND PHYSICAL ENGINEERING , C ZECH
TECHNICAL UNIVERSITY IN PRAGUE , T ROJANOVA 13, 12000 P RAGUE 2, C ZECH REPUBLIC
E-mail address: david.krejcirik@fjfi.cvut.cz
DEPARTMENT OF THEORETICAL PHYSICS , N UCLEAR PHYSICS INSTITUTE , C ZECH ACADEMY OF SCIENCES ,
25068 ˇREŽ , C ZECH REPUBLIC
E-mail address: lotoreichik@ujf.cas.cz
13